# Hybrid Blossom: A Parallel Primal-Dual Approach to Maximum-Weight Matching

**RS Lab**
**Status: complete draft — one open correctness limitation, precisely characterized (§6)**

---

## Abstract

We present Hybrid Blossom, an implementation of maximum-weight matching on general
undirected graphs that combines two lines of prior work: Blossom VI's primal-dual
framework for weighted matching (Arkhipov & Kolmogorov, arXiv:2604.20351) and
X-Blossom's lock-free parallel search for unweighted matching (Fan, Lee, Zhang, VLDB
2025). The integration point is Blossom VI's own theorem that, once dual variables
are fixed for a round, the primal step of a weighted blossom algorithm is exactly an
unweighted maximum-matching problem restricted to the tight-edge subgraph — which
licenses reusing X-Blossom's parallel search as that primal step. We built both a
serial reference implementation and a parallel one (persistent thread pool, adaptive
dispatch), verified correctness against a brute-force oracle on hundreds of randomized
graphs, and benchmarked both across graph sizes from n=200 to n=20,000. We report
three honest findings rather than three successes: (1) matchings are always
structurally valid, but a specific, precisely diagnosed dual-update bug produces
strictly suboptimal (not just imperfect) results on roughly 6% of random graphs; (2)
parallel execution does not beat serial execution at the graph sizes this project has
historically tested (n≤2000), but efficiency climbs steadily as graph size grows,
reaching approximate parity by n≈20,000; (3) a deeper investigation into the
correctness bug — including a new, reverted fix attempt — identifies a structural
reason a single-scalar-dual-per-vertex representation cannot host the necessary fix,
narrowing exactly what a correct fix requires.

---

## 1. Problem Statement

**Maximum-Weight Matching (MWM).** Given an undirected graph G=(V,E) with integer
edge weights w(u,v), find a matching M ⊆ E — a set of edges sharing no vertex — that
maximizes Σ_{e∈M} w(e). The matching need not be perfect; only vertices for which
some incident edge would improve the objective need be covered, and edges with
non-positive weight are never included in an optimal solution.

This is polynomial-time solvable (unlike the NP-hard weighted set-packing problems it
superficially resembles) via Edmonds' 1965 primal-dual blossom method, but a
correct, efficient implementation is notoriously difficult to get right — the
"blossom" (odd-cycle contraction) machinery is where most implementations
accumulate subtle bugs, ours included (§6).

## 2. Related Work

- **Edmonds (1965)**, "Paths, Trees, and Flowers" — the original blossom algorithm
  for maximum matching, introducing odd-cycle contraction as the mechanism that
  makes greedy augmenting-path search correct on general (non-bipartite) graphs.
- **Kolmogorov (2009)**, "Blossom V" — a widely used, carefully engineered
  implementation of minimum-weight perfect matching; the reference architecture for
  how real implementations give each contracted blossom its own dual variable
  (relevant directly to §6 below).
- **Arkhipov & Kolmogorov, arXiv:2604.20351**, "Blossom VI" — the primal-dual
  weighted-matching formulation this project follows directly, in particular its
  §1.4 theorem that the primal step of a weighted blossom algorithm restricted to
  tight edges is an ordinary unweighted matching problem.
- **Fan, Lee, Zhang, VLDB 2025**, "X-Blossom" — a lock-free parallel algorithm for
  unweighted maximum-cardinality matching using CAS-based tree claiming and a flat
  path-table representation instead of dynamic tree data structures; the source of
  this project's parallel primal search.

## 3. Algorithm

### 3.1 Dual variables and slack

Every vertex v carries a dual variable `dual[v] = 2·y[v]` (kept as `2y` so all
arithmetic stays integral). For an arc from u to v with weight w, the slack is

```
slack(u,v) = 4·w − dual[u] − dual[v]
```

Weights and duals are scaled by 4 (rather than the more familiar 2) purely to keep
every intermediate value an integer through the `⌈−slack/2⌉`-style delta
computations in dual updates (§3.4). An edge is *tight* (`slack = 0`) when it is
eligible to be added to the matching. Feasibility, for this max-weight convention,
requires `slack ≤ 0` for every edge and `slack = 0` for every currently matched
edge (complementary slackness).

Dual initialization sets `dual[v] = 2 · max{w(u,v) : (u,v) incident to v, w>0}`,
which guarantees feasibility from the start.

### 3.2 One augmentation round

Every exposed (unmatched) vertex with positive dual becomes the root of an
alternating-forest search (Edmonds' `SPLUS`/`SMINUS` labeling). The search follows
only tight, positive-weight arcs:

- **SPLUS → unlabeled, unmatched vertex**: an augmenting path exists; flip the
  matching along it.
- **SPLUS → unlabeled, matched vertex**: grow the tree (label the vertex `SMINUS`,
  its mate `SPLUS`, enqueue the mate).
- **SPLUS → SPLUS, different tree**: cross-tree augmenting path.
- **SPLUS → SPLUS, same tree**: an odd cycle has closed — shrink it into a single
  blossom (relabel every `SMINUS` vertex on the cycle as `SPLUS`, so the search can
  continue through it).

Multiple non-conflicting augmentations can be applied in a single round.

### 3.3 Blossom shrinking

The odd cycle discovered by a same-tree `SPLUS`-`SPLUS` collision is contracted by
finding the lowest common ancestor of the two colliding vertices in the alternating
forest and relabeling every vertex on the cycle to point at that ancestor
(`blossomBase[v] = LCA`). Crucially, in this implementation blossoms exist **only
for the remainder of the current round** — `blossomBase` is reset to the identity
at the start of every new `augmentRound()` call, and no blossom-*expansion*
operation is implemented (§6 explains why this matters).

### 3.4 Dual update

When a round's search exhausts without finding an augmenting path, dual variables
are adjusted so that new edges become tight: `SPLUS` duals *decrease* by δ,
`SMINUS` duals *increase* by δ. δ is the minimum of three constraints — (a) the
smallest `−slack` needed to make some `SPLUS→unlabeled` arc exactly tight, (b) the
smallest `⌈−slack/2⌉` needed to keep a cross-tree `SPLUS→SPLUS` arc from
overshooting into infeasibility, and (c) `dual[u]` itself, so no dual is driven
negative. δ = 0 signals optimality.

### 3.5 Why the X-Blossom warm-start is discarded

An earlier version of this project fed X-Blossom's fast unweighted matching in as a
warm start for the weighted solver. This was removed: an unweighted matching's
edges generally have nonzero slack under the weighted dual initialization, which
breaks complementary slackness and causes the dual-update loop to terminate
prematurely (δ=0 on the first call). The weighted solver now always starts from an
empty matching. X-Blossom's actual contribution to this project is **not** a warm
start — it is the parallel search algorithm reused, unmodified in its logic, as the
primal step described in §3.2, restricted to the tight-edge subgraph. This
narrower, more accurate framing of "hybrid" (a shared search algorithm across two
different admissibility filters, not a two-stage pipeline) is the one this paper
uses; earlier project documentation described a warm-start pipeline that no longer
reflects what the code does, and should be treated as superseded by this paper.

## 4. Implementation

### 4.1 Parallel primal search

`augmentRoundParallel` is a direct, weighted-aware port of X-Blossom's tree-claiming
algorithm, executed each round over the tight-edge subgraph
`E₀ = {e : w(e)>0 ∧ slack(e)=0}` — Blossom VI's own §1.4 theorem is what licenses
this substitution. Three bulk-synchronous phases per round:

- **Phase A (augment)**: each thread scans its share of the frontier for tight
  `SPLUS`-`SPLUS` collisions across trees; a thread that finds one CAS-claims both
  trees' root slots, and only the winner records the augmenting path.
- **Phase B (grow)**: threads scan for tight arcs into unvisited matched vertices,
  CAS-claiming the matched-edge id so exactly one thread grows a given pair.
- **Phase C (blossom detection)**: threads *detect* same-tree collisions in
  parallel (thread-local, no shared mutation) and the actual `shrink()` calls are
  applied serially afterward — the one primal operation deliberately not
  parallelized, since walking an unbounded ancestor chain isn't safe under simple
  per-node CAS.

`dualUpdateParallel` stripes the embarrassingly-parallel per-vertex delta
computation and adjustment across threads.

### 4.2 Making the parallelism real

A first implementation spawned a fresh `std::thread` group per phase per round —
matching the (single-shot) pattern already used elsewhere in the codebase — and
measured roughly 100x slower than serial, since a `solve()` call can run hundreds
of small rounds and OS thread creation costs ~0.1–0.5ms each on Windows. This was
fixed with a persistent `SimpleThreadPool` (workers created once per `solve()`
call, reused via a condition-variable barrier) plus adaptive dispatch
(`runAdaptive`): when a round's frontier is smaller than `512 × numThreads`, the
striped loop runs sequentially in the calling thread instead of paying a
mutex/condition-variable round trip. This cut wall-clock parallel overhead by
64–236x across benchmarked graph sizes (§7.1).

### 4.3 What is deliberately unimplemented

Cherry-blossom structures and a path-table replacement for the current
`par[]`/`pathToRoot()` representation were both scoped out as pure optimizations
that would add risk without being necessary for correctness. More importantly,
**blossom expansion is not implemented at all** — blossoms are contracted within a
round and discarded (not un-contracted) at the start of the next. This is a
legitimate simplification given blossoms never need to persist across rounds in
this design, but it is also directly implicated in the correctness limitation
below.

## 5. Correctness Methodology

Every produced matching is checked for structural validity (vertex-disjoint,
symmetric `M[M[v]]=v`) on every test and benchmark run — this has never failed,
including under 100+ repeated parallel runs on n=200 dense graphs. Weight
*optimality* is checked against a brute-force `O(2^m)` reference solver for n≤16,
across 13 categories (dense/sparse/mixed-sign/odd-cycle-heavy/all-negative graphs).
Parallel-vs-serial consistency is additionally checked at n up to 80 (thread counts
2/4/8/16) and determinism is checked via 8 repeated runs per configuration at fixed
thread counts.

Current results: 469 of 496 total test cases pass. All 27 failures are **valid but
strictly suboptimal** matchings — never invalid ones. This is the correctness
limitation investigated in depth below.

## 6. The Correctness Limitation

### 6.1 Symptom and root cause

Cross-checking against the brute-force oracle at higher volume than earlier smoke
testing (420 randomized cases) surfaced 27 valid-but-suboptimal results (6.4%),
confirmed present in the codebase prior to any of the parallelization work in this
project — not a regression introduced by it. Concrete counterexample: K4 with
weights `{0-1:16, 0-2:17, 0-3:17, 1-2:9, 1-3:11, 2-3:15}`. True optimum is
`{0-1, 2-3} = 31`; the solver returns `{0-2, 1-3} = 28`.

Root cause (§3.4's constraint (b)): the delta computation bounds cross-tree
`SPLUS`-`SPLUS` arcs but never bounds a **same-tree, different-blossom**
`SPLUS`-`SPLUS` arc — which exists whenever two separate branches of the same
alternating tree carry a non-tight edge between them that hasn't yet triggered a
shrink (shrinks only fire on *tight* same-tree collisions during the search; a
negative-slack same-tree edge is invisible to the search entirely). Left
unconstrained, δ can overshoot such an edge's slack from feasible-negative to
infeasible-positive, corrupting the LP-duality certificate the termination check
relies on and causing the solver to stop — and therefore accept a matching — one or
more rounds too early.

### 6.2 What has been tried

Four fix attempts have now been made and evaluated against the full randomized
suite. **All four were reverted; none improved on the 27-failure baseline.**

1. Add a same-tree `SPLUS`-`SPLUS` bound mirroring the existing cross-tree one,
   combined with an (independently incorrect) guess about `SMINUS`-`SMINUS` edges
   → 44 failures.
2. Additionally freeze blossom-interior duals by redirecting slack computation
   through the blossom's base vertex → no measurable improvement over (1).
3. Fix a ceiling/floor rounding error in the delta formula and remove the
   incorrect `SMINUS` constraint from (1) → improved to 30 failures, still worse
   than baseline.
4. **(This work.)** Key the same-tree bound on `blossomBase` identity rather than
   tree/root identity — i.e. bound δ for *any* `SPLUS`-`SPLUS` arc connecting two
   distinct blossoms, whether same-tree or cross-tree, while correctly excluding
   edges strictly interior to one already-shrunk blossom. This is a closer match
   to the textbook LP-duality granularity (distinct *pseudonodes*, not distinct
   *trees*) and does fix the K4 counterexample's first divergence point when
   traced by hand. Full-suite result: **30 failures** — worse than baseline.

### 6.3 Why attempt 4 fails: a structural finding, not just an empirical one

Hand-tracing attempt 4 through to termination (rather than only to its first
divergence point) exposes a second, independent problem that the earlier three
attempts' framing did not isolate. Once a blossom forms during a round, every one
of its member vertices is labeled `SPLUS` (by construction — that is what
shrinking does). The dual update applies δ **directly and uniformly to every
`SPLUS` vertex's own `dual[v]`**, blossom-interior or not. This means an edge
*internal* to a freshly-formed blossom — including, critically, a blossom's own
internal **matched** edge — has both endpoints shift by δ simultaneously, moving
its slack away from 0 by exactly 2δ and breaking complementary slackness on it.

This is not a missing bound that a fifth attempt could add: with a single scalar
`dual[v]` per vertex, there is no assignment of δ that simultaneously (a) correctly
shifts a blossom member's *external*-facing edges toward feasibility and (b) leaves
that same member's *internal* (intra-blossom) edges unchanged, because both
depend on the identical shift to the identical variable. The real fix needs two
separate quantities per vertex: a **frozen** per-vertex potential `y[v]` that stops
moving once v is absorbed into a blossom, and a separate **per-blossom** dual
`z[B]` that moves instead — with slack computed as
`w − (y[u] + Σ ancestor z's) − (y[v] + Σ ancestor z's)`, and edges strictly
interior to an already-contracted blossom excluded from ever being re-examined at
all (their tightness is guaranteed structurally at contraction time via Edmonds'
theorem, not re-verified via a slack formula). This is exactly the per-blossom
dual-variable architecture that Blossom V and this project's own Blossom VI
reference (arXiv:2604.20351, §2) use, and it requires genuine data-structure
changes (a `blossomDual` array, freeze/fold semantics on contraction, and — since
this codebase never expands a blossom once formed — at minimum verifying that
never expanding is still sound under the new representation) rather than a
one-line change to the existing delta formula. We regard this as the accurate,
final scope statement of the fix, sharper than the "give each blossom its own dual
variable" note in earlier project documentation, which did not identify *why* a
flat representation is structurally incapable of hosting even a well-targeted
partial fix.

### 6.4 Disposition

Given four independent, empirically-tested attempts — the most recent
additionally supported by a structural (not just empirical) argument for why
delta-formula patches cannot close this gap — we are treating the per-blossom
dual-variable rewrite as required, scoped, follow-on work rather than attempting a
fifth local patch. The current formula (pre-attempt-4, restored) remains in the
codebase as the best-performing of the five variants tested, with the full
investigation preserved as inline documentation at `dualUpdateSerial` in
`hybrid_blossom.cpp` so this analysis does not need to be repeated from scratch by
whoever undertakes the rewrite.

## 7. Performance Evaluation

### 7.1 Parallel overhead: before and after the threading fix (§4.2)

Same graphs and seeds, 8 threads, before vs. after the persistent thread pool and
adaptive dispatch:

| n | edges | density | before (ms) | after (ms) | improvement |
|---|---|---|---|---|---|
| 200 | 400 | sparse | 555.3 | 2.6 | 216x |
| 200 | 8000 | dense | 1280.7 | 11.6 | 110x |
| 500 | 1200 | sparse | 1917.1 | 8.1 | 236x |
| 500 | 40000 | dense | 1936.9 | 28.2 | 69x |
| 1000 | 2500 | sparse | 2788.7 | 21.1 | 133x |
| 2000 | 5000 | sparse | 3530.0 | 55.4 | 64x |

### 7.2 Parallel vs. serial at n=200–2000 (after the fix)

| n | edges | density | threads=1 (ms) | threads=8 (ms) | speedup |
|---|---|---|---|---|---|
| 200 | 400 | sparse | 1.11 | 2.56 | 0.43x |
| 200 | 8000 | dense | 4.82 | 11.62 | 0.42x |
| 500 | 1200 | sparse | 3.20 | 8.12 | 0.39x |
| 500 | 40000 | dense | 13.91 | 28.16 | 0.49x |
| 800 | 100000 | dense | 15.04 | 28.02 | 0.54x |
| 1000 | 2500 | sparse | 6.28 | 21.05 | 0.30x |
| 2000 | 5000 | sparse | 22.61 | 55.36 | 0.41x |

At these sizes, closing the *overhead* gap by 64–236x is not the same as closing
the *speedup* gap: the parallel path is still 2–3x slower than serial everywhere
tested. The identified reason is a fixed cost, not a synchronization-frequency
problem: a `solve()` call still needs hundreds of small rounds to converge (often
one augmentation per round), and the thread pool's own construction/destruction
(spawning/joining `nt` OS threads once per `solve()` call) is a cost the serial
path never pays at all — at wall-clock scales of single-digit-to-tens of
milliseconds, that fixed cost is a meaningful fraction of the whole run.

### 7.3 New data: does the gap close as graphs grow?

If the bottleneck is fixed per-`solve()`-call overhead not amortized over enough
work, efficiency should improve as graphs grow, since more real work exists to
absorb the same fixed cost. We tested this directly on synthetic sparse random
graphs (average degree 6) at n=5,000/10,000/20,000, best-of-3 wall-clock timing of
the weighted phase only:

| n | threads=1 | threads=8 | speedup@8 | best speedup |
|---|-----------|-----------|-----------|--------------|
| 5,000 | 105ms | 235ms | 0.45x | 0.54x (@4) |
| 10,000 | 327ms | 611ms | 0.54x | 0.64x (@4) |
| 20,000 | 761ms | 778ms | 0.98x | **1.03x (@4)** |

This confirms the diagnosis: efficiency climbs steadily with n and reaches
approximate parity (4 threads slightly *beating* serial for the first time) by
n≈20,000. A separate n=50,000 trial was excluded from this table: for that
particular random graph, the initial tight-edge matching alone happened to reach a
stable fixed point with zero dual-update rounds needed (an artifact of that
specific random instance, not a general property of n=50,000 graphs), so it
finished in single-digit milliseconds regardless of thread count and is not a
valid comparison point. Whether a graph requires many dual-update rounds appears
to depend on structural properties of the specific instance more than on n alone
— both n=5,000 and n=20,000 needed ~500–520 dual-update rounds despite differing
by 4x in vertex count.

### 7.4 Implication for further acceleration (e.g. GPU/CUDA)

Given CPU-thread parallelism only reaches break-even around n≈20,000, and a CUDA
kernel launch has a materially higher fixed-latency floor than the
`std::condition_variable` wake this thread pool already pays per round, a GPU port
of the same round structure would be expected to need substantially larger graphs
than n≈20,000 to break even — plausibly two or more orders of magnitude larger —
and would still be bounded by the same "too many too-small rounds" structural
issue as CPU threading is. We did not pursue a CUDA implementation for this
reason: restructuring toward fewer, larger rounds (§8) is the higher-leverage next
step regardless of which hardware executes the parallel primitives.

## 8. Limitations and Future Work

1. **The correctness bug (§6) is not fixed.** It is precisely diagnosed down to a
   structural argument for why the existing delta-formula representation cannot
   host a fix; the real fix (frozen per-vertex potentials + separate per-blossom
   duals, §6.3) is scoped but not built.
2. **Wall-clock parallel speedup is not demonstrated at the sizes historically
   tested (n≤2000)**, though the overhead gap was closed 64–236x, and new data
   (§7.3) shows efficiency reaching parity by n≈20,000 and continuing to improve
   with n.
3. **The round-granularity problem underlies both limitations above.** §7.2's
   closing analysis and §6.1's root cause share a common origin: the solver
   currently does many small rounds where Blossom V/VI's batched,
   connected-components-style delta update would do fewer, larger ones. Fixing
   this is likely to simultaneously improve both the correctness picture (fewer
   opportunities for the same-tree bound to be exercised incorrectly) and the
   performance picture (more work per round to amortize parallel dispatch cost).
   This is the highest-leverage next step identified by this work.
4. **Cherry-blossom structures and the path-table optimization are out of
   scope**, a deliberate decision to avoid destabilizing an already-tested
   traditional-blossom implementation.
5. **The X-Blossom phase's role in the "hybrid" architecture has been narrowed**
   from an earlier warm-start design (removed; see §3.5) to a shared search
   algorithm reused under a different admissibility filter. Earlier project
   documentation describing a warm-start pipeline should be treated as superseded
   by §3.5 of this paper.

## 9. Conclusion

Hybrid Blossom demonstrates that X-Blossom's parallel search algorithm can be
reused, without modification to its core logic, as the primal step of a weighted
primal-dual matching solver — exactly as Blossom VI's own theory predicts. The
engineering required to make that parallelism *net beneficial* (a persistent
thread pool, adaptive dispatch) is substantial and is now built and measured. Two
honest results remain open rather than resolved: a precisely characterized
correctness bug whose fix is now scoped down to a specific, necessary
architectural change, and a performance profile that improves with scale but does
not yet win at the scales this project has tested to date. Both point at the same
underlying fix — restructuring the primal-dual loop toward fewer, larger rounds —
identified here as the single highest-leverage next step for this project.

## References

- Edmonds, J. (1965). "Paths, Trees, and Flowers." *Canadian Journal of
  Mathematics*.
- Kolmogorov, V. (2009). "Blossom V: A New Implementation of a Minimum Cost
  Perfect Matching Algorithm." *Mathematical Programming Computation*.
- Arkhipov, D. & Kolmogorov, V. "Blossom VI: A Practical Minimum Weight Perfect
  Matching Algorithm." arXiv:2604.20351.
- Fan, W., Lee, K., Zhang, R. (2025). "X-Blossom: Massive Parallelization of Graph
  Maximum Matching." *VLDB 2025*.
