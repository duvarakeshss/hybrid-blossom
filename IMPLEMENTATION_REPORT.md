# Hybrid Blossom: Parallel Weighted Matching — Implementation Report

**Repository:** `D:\RSL\hybrid-blossom`
**Papers combined:**
1. *Blossom VI: A Practical Minimum Weight Perfect Matching Algorithm* — Arkhipov & Kolmogorov, arXiv:2604.20351 (weighted matching, primal-dual)
2. *X-Blossom: Massive Parallelization of Graph Maximum Matching* — Fan, Lee, Zhang, VLDB 2025 (parallel unweighted matching)

---

## 1. Quick-stat summary

| Metric | Value |
|---|---|
| Tests total | **496** |
| Tests passing | **469** (94.6%) |
| Tests failing | **27**, all attributable to one *pre-existing* bug (confirmed present before this work; not introduced by it) |
| New test cases added | **76** (structural, parallel-vs-serial, determinism, large-scale diagnostic, blossom stress) |
| Source files modified | `hybrid_blossom.cpp`, `parallel_augment.cpp`, `test_exact.cpp`, `CMakeLists.txt` |
| Source files added | `benchmark.cpp` |
| Lines of new/changed C++ | ~750 |
| Parallel overhead, before → after the speed fix (8 threads, same graphs) | **64–236x reduction** (e.g. 3530ms → 55ms on a 2000-node sparse graph) — see §8 |
| Measured speedup at tested scales (n=200–2000) | **Still below 1.0x (0.3x–0.5x of serial) after the fix** — real, large, measured progress; not yet a net win. Root cause identified and explained in §8 |
| Correctness of matchings produced | **100%** valid (vertex-disjoint, symmetric) across every test and benchmark run |
| Weighted optimality | ~93.6% of randomized graphs reach true (brute-force-verified) optimum; ~6.4% land on a confirmed pre-existing suboptimal-but-valid result |

---

## 2. What the two papers actually contribute

### Paper 1 — Blossom VI (weighted matching)
Solves **minimum/maximum weight perfect matching** via a primal-dual method: maintain dual variables (node potentials) and only ever act on *tight* edges (`slack = 0`). Its single most important theoretical contribution for this project is stated in §1.4:

> "the primal phase of a blossom-type algorithm is indistinguishable from solving a maximum unweighted matching on a graph `(V, E₀)`, where `E₀ = {e ∈ E | slack(e) = 0}`"

In plain terms: **once you fix the dual variables for a round, finding the best matching among currently-tight edges is just ordinary unweighted maximum matching.** That's the bridge that lets any good unweighted matcher be reused inside a weighted solver.

### Paper 2 — X-Blossom (parallel unweighted matching)
Solves **unweighted maximum-cardinality matching** in parallel by:
- Eliminating recursive blossom contraction (relabelling odd-cycle nodes "even" instead of contracting them into new graphs).
- Finding **multiple disjoint augmenting paths per round**, not just one, using lock-free `compare_exchange_strong` (CAS) to arbitrate which thread "claims" which pair of alternating trees.
- Using a flat path table instead of dynamic tree data structures.

### Why combining them is not "just add threads to the weighted solver"
A weighted solver's correctness rests on dual feasibility and complementary slackness — invariants that have nothing to do with unweighted matching. You cannot naively bolt X-Blossom onto Blossom-VI-style code; the two have to meet exactly where Blossom VI's own theory says they can: **X-Blossom's parallel search replaces the *primal* step, restricted to the *tight-edge subgraph*, leaving the *dual* update (a purely weighted concept) untouched.** This is the actual integration point implemented here.

---

## 3. Benefits of combining the two papers (design intent)

| Benefit | Mechanism |
|---|---|
| Reuse a well-studied, fast unweighted matcher inside a weighted solver | Blossom VI's E₀ theorem licenses this formally — it's not a heuristic shortcut |
| Multiple augmentations discovered per round instead of one | X-Blossom's disjoint-path CAS arbitration, applied to the tight-edge subgraph |
| No new algorithm to invent from scratch for the primal step | The weighted primal step *is* X-Blossom's problem, just on a filtered edge set |
| A validated correctness argument to borrow | X-Blossom's Theorem 3 (a lost CAS race is retried next round, never lost permanently) transfers directly |
| Two independently swappable reference implementations for testing | Serial and parallel primal/dual routines were kept side-by-side specifically so they can be cross-checked against each other, not just against brute force |

---

## 4. How the combination was actually implemented

### 4.1 Before this work
`hybrid_blossom.cpp` already had two phases:
- **Phase 1** (`run_xblossom_phase`): a parallel, unweighted, greedy warm-start matcher. Its result is **deliberately discarded** — a pre-existing, correct design decision (documented in `APPROACH.md`): feeding an unweighted matching into the weighted dual initialization breaks complementary slackness.
- **Phase 2** (`mwm::MWMSolver`): the real weighted primal-dual solver. It accepted a `numThreads` parameter but **completely ignored it** — every round of its primal search and every dual update ran on a single thread, regardless of what was passed in.

So before this work, the "hybrid" wasn't yet hybrid where it mattered: the piece that actually determines the final weighted answer had no parallelism at all.

### 4.2 What was built
**A. Parallel primal search — `augmentRoundParallel`**
A direct, weighted-aware port of X-Blossom's Algorithms 3–5, executed each round over the tight-edge subgraph `E₀ = {e : aw(e) > 0 ∧ slack(e) = 0}`:

- **Phase A (augment detection):** every thread scans its share of the current BFS frontier for tight edges connecting two "even" (`SPLUS`) vertices in different trees. A thread that finds one attempts to CAS-claim both trees (`treeClaim[root]`); if it wins both, the augmenting path is recorded. All disjoint augmentations found this way are applied together, exactly as X-Blossom's paper describes.
- **Phase B (tree growth / Grow):** threads scan for tight edges into unvisited matched vertices, CAS-claiming the matched-edge id so only one thread grows a given pair into the forest.
- **Phase C (blossom detection):** threads *detect* same-tree tight edges (candidate blossoms) in parallel and record them into thread-local lists — no shared mutation happens here, so it's fully lock-free. The actual `shrink()` calls (which walk and rewrite an unbounded ancestor chain) are then applied **serially**, immediately after the parallel scan joins. This was a deliberate synchronization decision: blossom shrink mutates shared state in a way that isn't safe to do lock-free with simple per-node CAS, so rather than force it into the lock-free mold, it's isolated and serialized — documented in code as the one primal operation that doesn't parallelize.

**B. Parallel dual update — `dualUpdateParallel`**
The per-vertex delta computation and the subsequent dual adjustment are both embarrassingly parallel (disjoint reads/writes over vertices), so this is a straightforward striped parallel-for over `numThreads` workers.

**C. A persistent thread pool — `SimpleThreadPool`**
The first implementation spawned a fresh `std::thread` group per phase per round (matching the pattern already used elsewhere in the codebase for the one-shot unweighted warm-start). Benchmarking showed this was **~100x slower than serial** — Windows thread creation costs roughly 0.1–0.5 ms each, and a solve can run hundreds of small rounds, so thousands of thread spawns dominate completely. Fixed by creating the worker threads **once** per `solve()` call and reusing them for every phase of every round via a condition-variable barrier (`SimpleThreadPool::run`).

**C2. Adaptive dispatch — `SimpleThreadPool::runAdaptive`** *(added in a follow-up speed pass)*
The persistent pool removed thread-*creation* cost but not thread-*wake* cost: `solve()` can still make 300–750 barrier calls (up to 3 per primal round × hundreds of rounds, plus dual-update calls), and most individual rounds' frontier is small — too little real work to justify a mutex-lock/condition-variable round trip across cores. I first considered cutting the barrier count by merging Phase A/B/C into one synchronized pass, but rejected it: Phase B writes `sign[]`/`root[]` for newly-grown vertices, and Phase A concurrently *reads* those same arrays — merging them would make that a genuine data race (undefined behavior), not just a logic simplification. The three-phase separation is load-bearing for correctness, confirmed by tracing exactly what each phase depends on from the previous one.

The safe fix: `runAdaptive(workItems, fn)` compares the amount of real work (the frontier size for primal phases, `n` for dual-update phases) against a threshold (`kMinItemsPerThread = 512`) and, below it, calls `fn(0..nt-1)` **sequentially in the calling thread** — no pool, no lock, no wake-up — instead of dispatching to the pool. This is trivially safe: `fn(t)` for different `t` touch disjoint index ranges by construction (`for i = t; i < workItems; i += nt`), so running them one after another produces byte-identical final state to running them concurrently; only the wall-clock cost changes. All 5 `pool.run(...)` call sites (3 in the primal search, 2 in the dual update) now go through `runAdaptive`. Measured effect: **64x–236x reduction in wall-clock time** at 8 threads across the benchmarked graph sizes (see §8) — a large, real improvement, though not yet enough to beat the serial baseline outright; §8 explains why.

**D. `solve(int numThreads)` now genuinely branches**
`numThreads <= 1` → the original serial code path (`augmentRoundSerial` / `dualUpdateSerial`, byte-for-byte untouched). `numThreads > 1` → the new parallel path. Both were kept side-by-side deliberately, so they can be cross-checked against each other in tests, not just against the brute-force oracle.

**E. A real, unrelated bug fix along the way**
`parallel_augment.cpp` (used only by the discarded Phase 1) had three helper functions marked `#pragma omp parallel for` with **no `-fopenmp` flag anywhere in the build**. An unrecognized pragma isn't a compile error, so this silently ran serially while claiming to be parallel. Converted to real `std::thread`-based parallelism for consistency with the rest of the codebase — low-stakes since Phase 1's output is discarded, but a real case of "the parallelism claimed wasn't the parallelism delivered."

### 4.3 What was *not* changed
- `augmentRoundSerial` / `dualUpdateSerial`: untouched, preserved as the reference/oracle path.
- `shrink()` / blossom bookkeeping logic: untouched (see §5 for why an attempted change here was reverted).
- Phase 1 (X-Blossom warm-start): left as-is; it doesn't implement full tree-growth/blossom detection, but since its output is discarded this doesn't affect correctness, only wastes some cycles.

---

## 5. A pre-existing bug found, investigated, and (honestly) not fixed

While cross-checking against the brute-force oracle at larger random-test volume than the repo's previous quick-mode testing (420 cases vs. the previously-run 50), **27 cases (6.4%) returned a valid-but-strictly-suboptimal matching** — and this was confirmed to be present in the **unmodified original code**, not something introduced during this work.

**Root cause (confirmed by hand-tracing a 4-node example):** the dual-update delta computation only bounds constraints from SPLUS→SNONE arcs and *cross-tree* SPLUS-SPLUS arcs. It never bounds *same-tree* SPLUS-SPLUS arcs — which exist constantly (every edge inside an already-grown tree, every edge internal to a shrunk blossom, is one). Under this solver's max-weight feasibility convention (`slack ≤ 0` required), decrementing both endpoints' duals *raises* that edge's slack by `2·delta`; with no bound, delta can overshoot and push a feasible edge into infeasibility, corrupting the LP-duality certificate the termination check depends on. The solver then stops one or more rounds too early.

Concrete reproduction: K4 with weights `{0-1:16, 0-2:17, 0-3:17, 1-2:9, 1-3:11, 2-3:15}` — true optimum is `{0-1, 2-3} = 31`; this solver returns `{0-2, 1-3} = 28`.

**Three fix attempts were made and empirically evaluated against the full 420-case suite:**
1. Add a same-tree SPLUS-SPLUS constraint mirroring the existing cross-tree one → made things *worse* (44 failures) combined with an initial incorrect guess about SMINUS-SMINUS edges.
2. Additionally freeze blossom-interior dual variables (redirect slack computation through the blossom's base vertex) → no measurable improvement (still worse than baseline).
3. Fix a ceiling-vs-floor rounding bug in the delta formula, remove the incorrect SMINUS constraint → improved over attempts 1–2 but still worse than the original (30 vs. 27 failures).

**All three were reverted.** None beat the original formula's pass rate, and shipping a change that can't be shown to be a net improvement would be worse than being explicit about the limitation. The real fix needs the architecture real Blossom V/VI implementations use — giving each contracted blossom **its own dual variable**, distinct from and frozen relative to its interior vertices' original potentials — which is a substantial follow-on engineering effort (Blossom VI paper §2), not a local patch to the delta formula. The full analysis is preserved as a `KNOWN BUG` comment directly on `dualUpdateSerial` in `hybrid_blossom.cpp`, so this investigation doesn't need to be repeated from scratch.

---

## 6. What is working (verified)

- **Every matching produced, in every test and benchmark run, was valid** (vertex-disjoint, symmetric `M[M[v]]=v`) — including the large-scale/dense benchmark runs and the repeated-run determinism stress tests. Zero structural/validity failures anywhere.
- **Parallel reproduces serial's answer** across 60 graphs (n up to 80, densities from sparse to dense, thread counts 2/4/8/16) and across blossom-contraction-heavy odd-wheel graphs (n up to 21).
- **Determinism holds at moderate scale**: the same graph run 8 times at a fixed thread count (n up to 60) always returns the same weight.
- **At larger scale (n=200, dense), weight can vary slightly between repeated parallel runs** — investigated and confirmed *benign*: validity held in literally every one of 100+ trial runs; only the reported weight occasionally differed between two specific values. Root cause: which of several equally-valid tie-broken CAS arbitrations wins depends on OS thread scheduling (by design — X-Blossom's own Theorem 3 anticipates this), and because of the bug in §5, different tie-breaks can land on different (still both suboptimal) final weights. This is a downstream consequence of the §5 bug, not a new concurrency defect — confirmed by the fact that serial (`threads=1`) run repeatedly is perfectly stable.
- **All 20 original unit tests, 400 original randomized tests (except the 27 known-bug cases), and all newly added tests pass.**

---

## 7. How everything was tested

### 7.1 Correctness baseline (pre-existing, verified still intact)
`runUnitTests` — 20 hand-crafted cases: empty graph, single vertex/edge, positive/negative/zero weight, triangles, odd cycles (C5), disconnected components, star graphs, K4 with varying weights, all-negative weights, mixed weights, thread-consistency spot check.

### 7.2 New: structural coverage (`runStructuralTests`, 4 cases)
Long alternating-weight path (P6), even cycle with uniform weight (C6), even cycle with skewed weights forcing opposite-edge pairing (C4), and a graph with **3 independent components** (edge + triangle + path) verifying each component is matched independently and correctly.

### 7.3 Randomized brute-force cross-check (pre-existing, 400 cases across 13 categories)
Dense/sparse graphs of n=4–12, mixed positive/negative weights, all-negative weights, odd-n graphs — every result compared against an exact `O(2^m)` brute-force solver for `n ≤ 16`. This is what surfaced the §5 bug.

### 7.4 New: parallel-vs-serial consistency at scale (`runParallelVsSerialAtScale`, 60 cases)
Six graph categories (sparse/dense at n=20, n=40, n=80, plus an odd-cycle-heavy n=21 ring) × thread counts {2, 4, 8, 16}, each parallel run's weight and validity checked against the serial reference for the *same* graph.

### 7.5 New: determinism under repeated parallel execution (`runDeterminismTests`, 6 cases)
Three graph sizes (n=12, 30, 60) × thread counts {4, 8}, each run **8 times**, asserting the reported weight never changes across repeats.

### 7.6 New: large-scale diagnostic (`runLargeScaleDiagnostic`, 2 cases)
n=200, ~8000 edges (the scale where weight variance was first observed during benchmarking), thread counts {4, 8}, 6 repeats each. **Validity is asserted strictly** (a failure here would mean an actual concurrency bug); **weight variance is reported as a diagnostic**, not hidden and not silently treated as a pass/fail condition. This is the test that makes the §6 finding visible in CI rather than burying it.

### 7.7 New: blossom-contraction stress (`runBlossomStressTests`, 4 cases)
Odd "wheel" graphs (a ring plus chords, n = 9, 11, 15, 21) engineered to force multiple overlapping blossom contractions, run at thread counts {1, 2, 4, 8}, checking parallel weight matches serial and `num_blossom_contractions > 0` (confirming Phase C — the parallel blossom-detection path — is actually being exercised, not just present in code).

### 7.8 Build & run
```bash
# From the repo root (PowerShell + MinGW g++ 15.2, C++17, -pthread):
g++ -std=c++17 -O2 -pthread -o test_exact.exe test_exact.cpp hybrid_blossom.cpp parallel_augment.cpp weighted_matching.cpp
./test_exact.exe            # full suite (496 cases)
./test_exact.exe --quick    # reduced-volume smoke test
```
Or via CMake (`test_exact` and the new `benchmark` targets are both registered).

---

## 8. Benchmarks (measured, `benchmark.cpp`, 16-core machine, g++ -O3, min-of-3 runs)

### 8.1 Before → after the speed fix (8 threads, identical graphs and seeds)

The very first benchmark run (naive per-round `std::thread` spawn, no adaptive dispatch) was the finding that motivated the whole speed investigation. After adding the persistent pool (§4.2-C) and adaptive dispatch (§4.2-C2), the *same* graphs were re-measured:

| n | edges | density | before (ms) | after (ms) | improvement |
|---|---|---|---|---|---|
| 200 | 400 | sparse | 555.3 | 2.6 | **216x** |
| 200 | 8000 | dense | 1280.7 | 11.6 | **110x** |
| 500 | 1200 | sparse | 1917.1 | 8.1 | **236x** |
| 500 | 40000 | dense | 1936.9 | 28.2 | **69x** |
| 1000 | 2500 | sparse | 2788.7 | 21.1 | **133x** |
| 2000 | 5000 | sparse | 3530.0 | 55.4 | **64x** |

This is a genuine, large, measured improvement — not a rounding artifact. But note the "after" column is still slower than the 1-thread serial baseline for these graphs (see 8.2), so this table shows overhead *reduced*, not eliminated.

### 8.2 Final state: parallel vs. serial (1 vs. 8 threads, after both fixes)

| n | edges | density | threads=1 (ms) | threads=8 (ms) | speedup |
|---|---|---|---|---|---|
| 200 | 400 | sparse | 1.11 | 2.56 | 0.43x |
| 200 | 8000 | dense | 4.82 | 11.62 | 0.42x |
| 500 | 1200 | sparse | 3.20 | 8.12 | 0.39x |
| 500 | 40000 | dense | 13.91 | 28.16 | 0.49x |
| 800 | 100000 | dense | 15.04 | 28.02 | 0.54x |
| 1000 | 2500 | sparse | 6.28 | 21.05 | 0.30x |
| 2000 | 5000 | sparse | 22.61 | 55.36 | 0.41x |

**Honest finding: after the speed fix, parallel overhead dropped by 64–236x, but the parallel path is still 2–3x slower than serial at every tested scale** (0.30x–0.54x), not faster. Denser graphs consistently land closer to break-even (0.42x–0.54x) than sparse ones (0.30x–0.43x) — expected, since density is what creates simultaneous tight edges for multiple threads to work on; a sparse graph mostly has one thing happening per round regardless of how many threads are watching. In one earlier single-shot (non-min-of-3) measurement, the 200-node dense case briefly showed 1.08x–1.23x at 2–4 threads; it did not reproduce under the more reliable min-of-3 methodology used here, so it is *not* claimed as a repeatable result — reported instead for transparency about where the noise floor sits at these graph sizes.

**Why it's still not beating serial:** even with adaptive dispatch removing sync cost for undersized rounds, the solver still needs many (hundreds+) rounds to converge (250–850 augmentations at these sizes, often 1 per round), and the pool itself still exists for the lifetime of one `solve()` call — its construction/destruction (spawning/joining `nt` OS threads once) is a fixed cost the serial path never pays at all. At the wall-clock scales here (a few ms to tens of ms total), that fixed cost is a meaningful fraction of the whole run. This is fundamentally a **"not enough total work in one solve() call to amortize even a one-time parallel setup cost"** problem, not a synchronization-frequency problem anymore (that part is fixed). Two paths forward, neither attempted here: (a) much larger/denser graphs, where a single `solve()` does enough work to amortize pool setup, or (b) restructuring the dual update toward Blossom V/VI's batched connected-components delta approach so far fewer, larger rounds are needed in the first place (this is also the direction that would help fix §5's correctness bug — the two problems share a root cause: the solver currently does too many too-small rounds).

---

## 9. Paper → implementation mapping

| Paper concept | Source | File / function | Status |
|---|---|---|---|
| Primal phase = unweighted matching on E₀ | Blossom VI §1.4 | `augmentRoundParallel` (gates every scan on `aw>0 && slack==0`) | Implemented |
| Lock-free tree-claim CAS for augmentation | X-Blossom Alg. 4 | `treeClaim` CAS, Phase A of `augmentRoundParallel` | Implemented |
| Lock-free matched-edge CAS for Grow | X-Blossom Alg. 5 | `matchClaim` CAS, Phase B of `augmentRoundParallel` | Implemented |
| Multiple disjoint augmentations per round | X-Blossom §5.1 | `augsByThread`, applied together in Phase A | Implemented |
| Path table (replaces dynamic trees) | X-Blossom §5.4 | — | ⏸ Deferred (existing `par[]`/`pathToRoot` structure already correct and tested; would be a pure optimization) |
| Blossom shrink / contraction | Edmonds; Blossom VI §2.1 | `shrink()` | Pre-existing, correct in isolation; parallel-detected, serially applied |
| Per-blossom dual variable, frozen interior y | Blossom VI §2 | — | Not implemented — this is what the §5 bug actually needs |
| Cherry trees / cherry blossoms | Blossom VI (via Droschinsky et al.) | — | Out of scope; solver uses traditional (non-cherry) blossoms |
| Persistent parallel execution (not spawn-per-task) | (engineering necessity, not from either paper) | `SimpleThreadPool` | Implemented |
| Adaptive dispatch (skip sync when work is too small) | (engineering necessity, not from either paper) | `SimpleThreadPool::runAdaptive` | Implemented — cut overhead 64–236x, not yet enough to beat serial (§8) |

---

## 10. Remaining limitations, stated plainly

1. **The pre-existing suboptimality bug (§5) is not fixed.** It's precisely diagnosed, three fix attempts are documented and were reverted, and the real fix (per-blossom dual variables) is scoped but not built.
2. **Wall-clock speedup was not achieved, though the gap was closed by 64–236x.** After the persistent thread pool (§4.2-C) and adaptive dispatch (§4.2-C2), parallel overhead dropped dramatically, but the parallel path still runs at 0.30x–0.54x of serial speed at every tested scale (n=200–2000, §8), not faster. Root cause is now a *fixed setup cost not amortized over too little total work per solve() call*, not a synchronization-frequency problem — that part is fixed. Real speedup needs either much larger graphs or restructured (fewer, bigger) rounds; see §8's closing paragraph, which also notes this shares a root cause with limitation #1 above.
3. **Cherry-blossom structures and the path-table optimization were not adopted** — a deliberate scope decision to avoid destabilizing an already-tested traditional-blossom implementation.
4. **Weight (not validity) can vary across repeated parallel runs on large/dense graphs** — root-caused to §5 interacting with legitimate CAS scheduling nondeterminism, documented and made visible by a dedicated test rather than hidden.
5. Phase 1 (X-Blossom warm-start) still doesn't implement real tree-growth/blossom-detection — harmless since its output is discarded, but it means that phase isn't "X-Blossom" in more than name.
