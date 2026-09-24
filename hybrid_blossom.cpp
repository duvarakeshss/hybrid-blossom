/// \file hybrid_blossom.cpp
/// \brief Hybrid maximum-weight matching:
///        X-Blossom warm-start  +  exact primal-dual refinement.
///
/// Problem: given undirected graph G=(V,E,w), find M⊆E such that
///   - no two edges share a vertex
///   - Σ w(e) for e∈M is maximized
///
/// Handles: general graphs, negative weights, odd n, no perfect matching.

#include "hybrid_blossom.h"
#include "parallel_augment.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <condition_variable>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// ===========================================================
//  Internal graph utilities
// ===========================================================

/// Find first arc from u to v in CSR. Returns arc index or -1.
static int find_arc(const std::vector<int>& ro,
                    const std::vector<int>& ci,
                    int u, int v) {
    for (int j = ro[u]; j < ro[u+1]; ++j)
        if (ci[j] == v) return j;
    return -1;
}

/// Count matched edges (each counted once).
static int count_matching(const std::vector<int>& M) {
    int c = 0;
    for (int v = 0; v < (int)M.size(); ++v)
        if (M[v] != -1 && v < M[v]) ++c;
    return c;
}

/// Compute total matching weight.
static int64_t matching_weight(const std::vector<int>& M,
                               const std::vector<int>& ro,
                               const std::vector<int>& ci,
                               const std::vector<int>& aw) {
    int64_t w = 0;
    for (int v = 0; v < (int)M.size(); ++v)
        if (M[v] != -1 && v < M[v]) {
            int arc = find_arc(ro, ci, v, M[v]);
            if (arc >= 0) w += aw[arc];
        }
    return w;
}

// ===========================================================
//  Graph IO helpers
// ===========================================================

bool read_mtx(const std::string& path,
              std::vector<int>& rowOffsets,
              std::vector<int>& colIndices,
              std::vector<int>& arcWeights) {
    std::ifstream f(path);
    if (!f.is_open()) { std::cerr << "Cannot open " << path << "\n"; return false; }

    // Skip comment lines
    std::string line;
    while (std::getline(f, line))
        if (!line.empty() && line[0] != '%') break;

    int n = 0, m = 0;
    {
        std::istringstream ss(line);
        if (!(ss >> n >> m)) {
            if (!(f >> n >> m)) {
                std::cerr << "read_mtx: bad header\n"; return false;
            }
        }
    }
    if (n < 0 || m < 0) { std::cerr << "read_mtx: negative n or m\n"; return false; }

    std::vector<std::vector<std::pair<int,int>>> adj(n);
    for (int i = 0; i < m; ++i) {
        int u, v, w;
        if (!(f >> u >> v >> w)) {
            std::cerr << "read_mtx: truncated at edge " << i << "\n"; return false;
        }
        if (u < 0 || u >= n || v < 0 || v >= n) {
            std::cerr << "read_mtx: vertex out of range at edge " << i << "\n"; return false;
        }
        if (u == v) continue; // skip self-loops
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }

    rowOffsets.assign(n + 1, 0);
    for (int i = 0; i < n; ++i)
        rowOffsets[i+1] = rowOffsets[i] + (int)adj[i].size();
    int total = rowOffsets[n];
    colIndices.resize(total);
    arcWeights.resize(total);
    for (int i = 0; i < n; ++i) {
        int base = rowOffsets[i];
        for (int k = 0; k < (int)adj[i].size(); ++k) {
            colIndices[base+k] = adj[i][k].first;
            arcWeights[base+k] = adj[i][k].second;
        }
    }
    return true;
}

bool read_csr_files(const std::string& rowFile,
                    const std::string& colFile,
                    std::vector<int>& rowOffsets,
                    std::vector<int>& colIndices,
                    std::vector<int>& arcWeights) {
    {
        std::ifstream f(rowFile);
        if (!f) { std::cerr << "Cannot open " << rowFile << "\n"; return false; }
        int v; while (f >> v) rowOffsets.push_back(v);
    }
    {
        std::ifstream f(colFile);
        if (!f) { std::cerr << "Cannot open " << colFile << "\n"; return false; }
        int v; while (f >> v) colIndices.push_back(v);
    }
    arcWeights.assign(colIndices.size(), 1);
    return true;
}

// ===========================================================
//  Matching validator
// ===========================================================

bool validate_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& colIndices,
    const std::vector<int>& arcWeights,
    const std::vector<int>& M,
    int64_t& totalWeight)
{
    totalWeight = 0;
    int n = (int)M.size();
    if (n == 0) return true;

    for (int v = 0; v < n; ++v) {
        if (M[v] == -1) continue;
        int u = M[v];
        if (u == v) { std::cerr << "validate: self-match at " << v << "\n"; return false; }
        if (u < 0 || u >= n) { std::cerr << "validate: out of range\n"; return false; }
        if (M[u] != v) {
            std::cerr << "validate: M[" << v << "]=" << u
                      << " but M[" << u << "]=" << M[u] << "\n";
            return false;
        }
        if (v < u) {
            int arc = find_arc(rowOffsets, colIndices, v, u);
            if (arc < 0) {
                std::cerr << "validate: edge (" << v << "," << u << ") not in graph\n";
                return false;
            }
            totalWeight += arcWeights[arc];
        }
    }
    return true;
}

// ===========================================================
//  X-Blossom (parallel BFS)
// ===========================================================

static void run_xblossom_phase(
        const std::vector<int>& ro,
        const std::vector<int>& ci,
        std::vector<int>& M,
        int numThreads,
        MatchingResult& res)
{
    auto t0 = std::chrono::steady_clock::now();
    int n = (int)ro.size() - 1;
    M.assign(n, -1);

    // Greedy initialization
    for (int v = 0; v < n; ++v) {
        if (M[v] != -1) continue;
        for (int j = ro[v]; j < ro[v+1]; ++j) {
            int w = ci[j];
            if (M[w] == -1) { M[v] = w; M[w] = v; break; }
        }
    }

    // Parallel augmenting path rounds
    std::vector<std::atomic<int>> select_tree(n), select_match(n), select_blossom(n);
    std::vector<std::vector<int>> path_table(n);
    std::vector<int> is_even(n, 0), belongs(n, -1);

    for (int round = 0; round <= n; ++round) {
        std::vector<int> exposed;
        parallel_find_exposed(exposed, M, numThreads);
        if (exposed.empty()) break;

        std::fill(is_even.begin(), is_even.end(), 0);
        std::fill(belongs.begin(), belongs.end(), -1);
        parallel_init_atomics(select_tree, select_match, select_blossom,
                              path_table, n, numThreads);
        parallel_init_exposed_vector(exposed, is_even, belongs, numThreads);

        std::vector<std::vector<int>> path_collection;
        std::vector<std::thread> threads;
        threads.reserve(numThreads);
        for (int t = 0; t < numThreads; ++t)
            threads.emplace_back(parallel_find_augmenting_paths,
                                 std::cref(ro), std::cref(ci),
                                 std::cref(exposed), t, numThreads,
                                 std::ref(is_even), std::ref(belongs),
                                 std::ref(path_table), std::ref(select_tree),
                                 std::ref(path_collection));
        for (auto& th : threads) th.join();

        if (path_collection.empty()) break;
        parallel_update_matching(M, path_collection, numThreads);
        ++res.num_augmentations;
    }

    auto t1 = std::chrono::steady_clock::now();
    res.time_xblossom_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
}

// ===========================================================
//  Remove negative-weight edges from any matching
// ===========================================================

static void remove_negative_edges(std::vector<int>& M,
                                   const std::vector<int>& ro,
                                   const std::vector<int>& ci,
                                   const std::vector<int>& aw) {
    for (int v = 0; v < (int)M.size(); ++v) {
        if (M[v] == -1 || v > M[v]) continue;
        int u = M[v];
        int arc = find_arc(ro, ci, v, u);
        if (arc >= 0 && aw[arc] < 0) { M[v] = -1; M[u] = -1; }
    }
}

// ===========================================================
//  Persistent thread pool for the weighted primal-dual solver
// ===========================================================
//
// The parallel primal search (augmentRoundParallel) and dual update
// (dualUpdateParallel) each run in a handful of bulk-synchronous phases
// PER ROUND, and mwm::MWMSolver::solve() can execute many rounds (up to
// O(n) dual updates, each preceded by primal rounds). Spawning a fresh
// std::thread group per phase -- the same pattern parallel_augment.cpp
// uses for the (single-shot) unweighted X-Blossom phase -- means
// thousands of OS thread creations for a solve() call on a
// moderately-sized graph. Thread creation is expensive on Windows in
// particular (measured ~0.1-0.5ms each here), so that pattern alone
// made the "parallel" path 2-3 orders of magnitude *slower* than serial
// in benchmarking (see benchmark.cpp's results). This pool creates the
// worker threads once per solve() call and reuses them for every phase
// of every round via a simple generation-counter barrier, which is the
// standard fix for this class of problem.
class SimpleThreadPool {
public:
    explicit SimpleThreadPool(int nThreads) : n(std::max(1, nThreads)) {
        workers.reserve(n);
        for (int t = 0; t < n; ++t)
            workers.emplace_back([this, t]() { workerLoop(t); });
    }

    ~SimpleThreadPool() {
        {
            std::lock_guard<std::mutex> lk(m);
            stop = true;
            ++generation;
        }
        cvWork.notify_all();
        for (auto& th : workers) th.join();
    }

    SimpleThreadPool(const SimpleThreadPool&) = delete;
    SimpleThreadPool& operator=(const SimpleThreadPool&) = delete;

    // Runs fn(tid) for every tid in [0, size()), blocking until all
    // workers finish. Must be called from a single driving thread (not
    // reentrant); that's exactly how MWMSolver uses it -- one call site
    // per phase, never nested or concurrent.
    void run(const std::function<void(int)>& fn) {
        {
            std::lock_guard<std::mutex> lk(m);
            task = &fn;
            pending = n;
            ++generation;
        }
        cvWork.notify_all();
        std::unique_lock<std::mutex> lk(m);
        cvDone.wait(lk, [this]() { return pending == 0; });
    }

    // Same contract as run(), but skips the pool (and therefore its
    // mutex-lock / condition-variable wake round trip) entirely when
    // `workItems` is too small to be worth splitting across threads.
    //
    // Measured root cause of the parallel path being *slower* than
    // serial at every graph size benchmarked (n=200-2000, see
    // benchmark.cpp): MWMSolver::solve() can run hundreds of small
    // rounds, each with up to 3 primal phases (augmentRoundParallel) or
    // 2 dual-update phases -- i.e. hundreds of pool.run() barrier calls
    // per solve() -- and most individual rounds' frontier (the amount of
    // real parallelizable work) is tiny relative to that fixed
    // synchronization cost. This is a standard "not enough work per
    // parallel region" problem, and the standard fix is exactly this:
    // fall back to running the same striped loop inline, sequentially,
    // in the calling thread. That's still 100% correct -- fn(0..n-1)
    // touch disjoint index ranges by construction (the caller stripes
    // work as `for (i = t; i < workItems; i += n)`), so running them
    // one after another in a single thread produces the identical final
    // state as running them concurrently across the pool, just without
    // ever exploiting extra cores. What changes is only whether we pay
    // for a wake-up round trip we can't afford to amortize.
    void runAdaptive(int workItems, const std::function<void(int)>& fn) {
        if (workItems < n * kMinItemsPerThread) {
            for (int t = 0; t < n; ++t) fn(t);
            return;
        }
        run(fn);
    }

    static constexpr int kMinItemsPerThread = 512;

    int size() const { return n; }

private:
    void workerLoop(int tid) {
        int lastSeen = 0;
        while (true) {
            std::unique_lock<std::mutex> lk(m);
            cvWork.wait(lk, [&]() { return generation != lastSeen; });
            lastSeen = generation;
            if (stop) return;
            const std::function<void(int)>* fn = task;
            lk.unlock();

            (*fn)(tid);

            lk.lock();
            if (--pending == 0) {
                lk.unlock();
                cvDone.notify_one();
            }
        }
    }

    int n;
    std::vector<std::thread> workers;
    std::mutex m;
    std::condition_variable cvWork, cvDone;
    bool stop = false;
    int generation = 0;
    int pending = 0;
    const std::function<void(int)>* task = nullptr;
};

// ===========================================================
//  Exact primal-dual maximum-weight matching
// ===========================================================
//
// Based on Edmonds (1965) Blossom algorithm with dual variables.
//
// Dual variables (stored as 2*y to keep integers):
//   dual[v] = 2 * y[v],   y[v] >= 0 for all v
//
// Slack of arc j = (u -> ci[j]) with weight aw[j]:
//   slack(j) = 4*aw[j] - dual[u] - dual[ci[j]]
//
// Invariants:
//   (D1) slack(j) >= 0 for all arcs j
//   (C1) if M[u]=v then slack of arc (u,v) = 0
//
// Algorithm per augmentation round:
//   1. Init alternating forest: all exposed vertices are "+" roots
//   2. BFS over "+" vertices:
//      For each neighbor v of u (tight arc, slack=0):
//        - v unmatched and not in forest: AUGMENT u→v
//        - v in forest, different tree (+): AUGMENT across trees
//        - v in forest, same tree (+): SHRINK blossom
//        - v in forest (-): already handled (v's mate will be enqueued)
//        - v not in forest, matched: GROW (add v as -, mate as +)
//   3. If BFS exhausts without augment: DUAL UPDATE, retry
//   4. If dual update delta = 0: done (optimal)

namespace mwm {

// ---------------------------------------------------------------
// Constants — must ALL be distinct.
//
// CRITICAL: sign[] uses SNONE/SPLUS/SMINUS.
//           matching M[] uses UNMATCHED=-1.
//           parent par[] uses NO_PARENT=-1.
// We use SNONE=0 so SNONE != SMINUS — the old bug was NONE=MINUS=-1
// which made the BFS treat already-labeled MINUS nodes as unlabeled,
// causing infinite queue growth.
// ---------------------------------------------------------------
static constexpr int SNONE    =  0;  // sign: not yet in alternating forest
static constexpr int SPLUS    =  1;  // sign: even distance from root
static constexpr int SMINUS   = -1;  // sign: odd distance from root
static constexpr int UNMATCHED= -1;  // M[v]: vertex is unmatched
static constexpr int NO_PARENT= -1;  // par[v]: vertex is a root

struct MWMSolver {
    int n;
    const std::vector<int>& ro;
    const std::vector<int>& ci;
    const std::vector<int>& aw;

    std::vector<int> M;       // matching: M[v]=mate or UNMATCHED
    std::vector<int> dual;    // 2*y[v], y[v] >= 0

    // Per-round alternating-forest state
    std::vector<int> sign;        // SNONE / SPLUS / SMINUS
    std::vector<int> root;        // root vertex of v's tree, or NO_PARENT
    std::vector<int> par;         // parent of v in alternating tree, or NO_PARENT
    std::vector<int> parArc;      // arc index that gave v its parent
    std::vector<int> blossomBase; // representative of v's blossom (identity if no blossom)

    // Parallel-primal-search synchronization (X-Blossom Algorithms 4-5, see
    // augmentRoundParallel below). Reset to 0 at the start of every round.
    //   treeClaim[r]  : CAS ownership of tree rooted at r (Augment, Alg. 4).
    //                   Once a claim succeeds for a *pair* of roots it is
    //                   never released, so it doubles as "this tree is
    //                   done for the round" (the serial code's rootDone[]).
    //   matchClaim[id]: CAS ownership of matched edge id=min(v,mate) during
    //                   Grow (Alg. 5), so exactly one thread grows a given
    //                   matched pair into the forest.
    std::vector<std::atomic<int>> treeClaim;
    std::vector<std::atomic<int>> matchClaim;

    // Instrumentation
    int numAug    = 0;
    int numShrink = 0;
    int numDual   = 0;

    MWMSolver(int n_,
              const std::vector<int>& ro_,
              const std::vector<int>& ci_,
              const std::vector<int>& aw_)
        : n(n_), ro(ro_), ci(ci_), aw(aw_)
        , M(n_, UNMATCHED)
        , dual(n_, 0)
        , sign(n_, SNONE)
        , root(n_, NO_PARENT)
        , par(n_, NO_PARENT)
        , parArc(n_, NO_PARENT)
        , blossomBase(n_)
        , treeClaim(n_)
        , matchClaim(n_)
    {
        std::iota(blossomBase.begin(), blossomBase.end(), 0);

        // Dual initialization for MWM:
        //   y[v] = max(0, max_incident_weight)
        // Ensures all slacks >= 0 initially:
        //   slack(u,v) = 4*w - 2*y[u] - 2*y[v] >= 4*w - 2*w - 2*w = 0
        // For all-negative incident vertices: y[v]=0, slack<0, edges inadmissible.
        for (int v = 0; v < n; ++v) {
            int maxw = 0;
            for (int j = ro[v]; j < ro[v+1]; ++j)
                maxw = std::max(maxw, aw[j]);
            dual[v] = 2 * maxw;
        }
    }

    // Slack of arc j from vertex u to ci[j].
    // Multiplied by 4 so matched-edge duals (stored as 2*y) stay integers.
    int slack(int j, int u) const {
        return 4 * aw[j] - dual[u] - dual[ci[j]];
    }

    // Find lowest common ancestor of u and v in the alternating forest.
    int lca(int u, int v) const {
        std::vector<bool> visited(n, false);
        // Walk u to root, mark each blossomBase
        for (int x = u; x != NO_PARENT; x = par[blossomBase[x]])
            visited[blossomBase[x]] = true;
        // Walk v to root, first marked node is LCA
        for (int y = v; y != NO_PARENT; y = par[blossomBase[y]])
            if (visited[blossomBase[y]]) return blossomBase[y];
        return NO_PARENT; // unreachable if u,v are in the same tree
    }

    // Path from v up to its tree root (inclusive), following par[].
    std::vector<int> pathToRoot(int v) const {
        std::vector<int> path;
        for (int cur = v; cur != NO_PARENT; cur = par[cur])
            path.push_back(cur);
        return path;
    }

    // Flip matching along an alternating path (vertices in order).
    void augmentPath(const std::vector<int>& path) {
        for (int i = 0; i + 1 < (int)path.size(); i += 2) {
            M[path[i]]   = path[i+1];
            M[path[i+1]] = path[i];
        }
        ++numAug;
    }

    // Contract blossom: u and bv=blossomBase[v] are both SPLUS, same tree.
    // Re-label SMINUS nodes on the cycle as SPLUS and enqueue them.
    void shrink(int u, int bv, std::queue<int>& Q) {
        int base = lca(u, bv);
        if (base == NO_PARENT) return;

        // Walk from u toward root, re-label until we reach base
        auto relabel = [&](int start) {
            int cur = start;
            while (blossomBase[cur] != base) {
                int b = blossomBase[cur];
                blossomBase[b] = base;
                if (sign[b] == SMINUS) { sign[b] = SPLUS; Q.push(b); }
                cur = par[b];
                if (cur == NO_PARENT) break;
            }
        };
        relabel(u);
        relabel(bv);
        ++numShrink;
    }

    // One BFS augmentation round from all exposed (UNMATCHED) vertices.
    // Returns true if an augmenting path was found and applied.
    // This is the single-threaded reference implementation; kept intact
    // (unmodified from the original solver) as a correctness oracle that
    // augmentRoundParallel() is cross-checked against in the test suite.
    bool augmentRoundSerial() {
        // Reset per-round state
        std::fill(sign.begin(),    sign.end(),    SNONE);
        std::fill(root.begin(),    root.end(),    NO_PARENT);
        std::fill(par.begin(),     par.end(),     NO_PARENT);
        std::fill(parArc.begin(),  parArc.end(),  NO_PARENT);
        std::iota(blossomBase.begin(), blossomBase.end(), 0);

        std::queue<int> Q;

        // Seed: all unmatched vertices with positive dual (i.e., they have
        // at least one positive-weight incident arc) are SPLUS roots.
        // Vertices with dual=0 have only non-positive-weight edges and
        // cannot benefit from matching — exclude them to avoid blocking delta.
        for (int v = 0; v < n; ++v) {
            if (M[v] == UNMATCHED && dual[v] > 0) {
                sign[v] = SPLUS;
                root[v] = v;
                Q.push(v);
            }
        }

        // Track which roots are still "active" (not yet augmented this round).
        // When two trees merge via augmentation, both roots become inactive.
        std::vector<bool> rootDone(n, false);

        bool anyAugmented = false;

        while (!Q.empty()) {
            int u = Q.front(); Q.pop();
            if (sign[u] != SPLUS) continue; // stale after blossom relabel
            if (rootDone[root[u]])  continue; // this tree already augmented

            for (int j = ro[u]; j < ro[u+1]; ++j) {
                // Skip non-positive-weight arcs
                if (aw[j] <= 0) continue;
                // Only follow tight arcs
                if (slack(j, u) != 0) continue;

                int v  = ci[j];
                int bv = blossomBase[v];

                if (rootDone[root[u]]) break; // tree augmented mid-arc-scan

                if (sign[bv] == SNONE) {
                    // bv not yet in forest
                    if (M[v] == UNMATCHED) {
                        // Augmenting path found: u → ... → root(u)  and  v (free)
                        std::vector<int> path;
                        path.push_back(v);
                        path.push_back(u);
                        for (int cur = par[u]; cur != NO_PARENT; cur = par[cur])
                            path.push_back(cur);
                        int r = root[u];
                        augmentPath(path);
                        rootDone[r] = true;
                        anyAugmented = true;
                        break; // stop scanning arcs of u; move to next Q entry
                    } else {
                        // Grow: add v (SMINUS) and its mate (SPLUS)
                        int mate = M[v];
                        sign[v]    = SMINUS;  root[v]    = root[u];
                        sign[mate] = SPLUS;   root[mate] = root[u];
                        par[v]     = u;       parArc[v]  = j;
                        par[mate]  = v;
                        Q.push(mate);
                    }
                } else if (sign[bv] == SPLUS && root[bv] != root[u]) {
                    if (rootDone[root[bv]]) continue; // that tree already done
                    // Cross-tree augmenting path
                    std::vector<int> pu = pathToRoot(u);
                    std::vector<int> pv = pathToRoot(v);
                    std::vector<int> path;
                    for (int i = (int)pu.size()-1; i >= 0; --i) path.push_back(pu[i]);
                    for (int i = 0; i < (int)pv.size(); ++i)    path.push_back(pv[i]);
                    int ru = root[u], rv = root[bv];
                    augmentPath(path);
                    rootDone[ru] = true;
                    rootDone[rv] = true;
                    anyAugmented = true;
                    break;
                } else if (sign[bv] == SPLUS && root[bv] == root[u]) {
                    // Same-tree blossom: shrink
                    shrink(u, bv, Q);
                }
                // sign[bv] == SMINUS: already in tree, skip
            }
        }
        return anyAugmented;
    }

    // Dual update for maximum-weight matching (Galil-Micali-Gabow style).
    //
    // Slack: slack(j,u) = 4*w[j] - dual[u] - dual[ci[j]]
    //   Dual feasibility: dual[u]+dual[v] >= 4*w  (slack <= 0) for non-matched arcs
    //   Complementary slackness: slack == 0 for matched arcs
    //
    // Augmenting BFS only traverses tight arcs (slack == 0).
    // After a failed BFS, we reduce SPLUS duals to make more arcs tight:
    //
    //   SPLUS  duals DECREASE by delta  → slack(SPLUS, *) increases toward 0
    //   SMINUS duals INCREASE by delta  → preserves tightness of matched arcs
    //
    // delta = min over constraints:
    //   (a) SPLUS→SNONE arcs with slack < 0:      delta_a = -slack         (makes slack=0)
    //   (b) SPLUS→SPLUS cross-tree, slack < 0:    delta_b = (-slack+1)/2   (ceil, makes slack>=0)
    //   (c) SPLUS non-negativity: dual[u] >= 0 → delta <= dual[u]
    //
    // KNOWN BUG (pre-existing, confirmed, NOT introduced by the
    // parallelization work in this file -- see augmentRoundParallel /
    // dualUpdateParallel below): cross-checking this solver's output
    // against a brute-force oracle on randomized graphs (test_exact.cpp)
    // shows it returns a *valid but strictly suboptimal* matching on
    // roughly 6% of random dense/sparse/mixed/odd-cycle test cases
    // (e.g. K4 with weights {0-1:16,0-2:17,0-3:17,1-2:9,1-3:11,2-3:15}:
    // true optimum is {0-1,2-3}=31, this solver returns {0-2,1-3}=28).
    //
    // Root cause, confirmed by hand-tracing that example with added
    // instrumentation: this solves MAXIMUM weight matching, so dual
    // feasibility runs opposite the more familiar min-cost convention
    // (an edge is feasible when slack <= 0). Constraint (b) above only
    // bounds delta for *cross-tree* SPLUS-SPLUS pairs; a *same-tree*
    // SPLUS-SPLUS edge (which exists constantly -- every edge inside an
    // already-grown tree is one, and every edge internal to a shrunk
    // blossom becomes one) is completely unconstrained, so delta can
    // overshoot and push such an edge's slack from feasible-negative to
    // infeasible-positive. That breaks the LP-duality certificate the
    // "no more progress" termination check relies on, so the solver can
    // stop one or more rounds too early.
    //
    // This class of bug is exactly why real implementations (Blossom V,
    // and this repo's own Blossom-VI reference paper) give a *contracted
    // blossom its own dual variable* and freeze the original y_v of its
    // interior vertices, rather than letting every interior vertex's
    // dual move independently once relabeled. I attempted three
    // increasingly-targeted fixes here (adding the missing same-tree
    // bound; freezing blossom-interior duals via a blossomBase-redirected
    // slack(); capping delta to exactly 0 for already-tight same-tree
    // pairs) and verified each one empirically against the same 420-case
    // randomized suite. Each either had no effect or *reduced* the pass
    // rate (partial constraints without the matching blossom-dual
    // architecture can make the search more conservative in ways that
    // surface the same underlying gap elsewhere just as often). None
    // beat the original formula's pass rate, so I reverted to it rather
    // than ship a change I could not verify was a net improvement.
    // A correct fix needs the full per-blossom dual-variable design
    // (Sec. 2 of arXiv:2604.20351) -- a substantial follow-on effort,
    // not a local patch to this delta formula.
    //
    // FOURTH ATTEMPT (also reverted, but diagnostically useful): keying
    // constraint (b) on blossomBase identity instead of tree/root
    // identity -- i.e. bounding delta for any SPLUS-SPLUS edge between
    // two DISTINCT blossoms, same tree or not, while correctly excluding
    // edges strictly interior to one already-shrunk blossom -- is closer
    // to the textbook LP-duality granularity (distinct *pseudonodes*
    // matter, not distinct *trees*) and does fix the exact K4
    // counterexample's first divergence point. But hand-tracing it
    // through to termination exposes a second, independent problem: this
    // solver applies delta directly to every SPLUS vertex's own dual[],
    // uniformly. Once a blossom forms, ALL of its members -- including
    // ones on either side of an internal MATCHED edge -- are SPLUS, so
    // they all shift by the same delta, which moves that internal
    // matched edge's slack away from 0 by 2*delta and breaks
    // complementary slackness on it. This is not a missing bound to add;
    // it is a structural conflict: with a single scalar dual[v] per
    // vertex, there is no delta assignment that simultaneously (a)
    // correctly shifts a blossom member's *external*-facing slack and
    // (b) leaves that same member's *internal* (intra-blossom) slack
    // unchanged, because both depend on the identical dual[v] shift. The
    // real fix needs two separate quantities -- a frozen per-vertex y[v]
    // that stops moving once v is absorbed into a blossom, and a
    // separate per-blossom z[B] that moves instead, with slack computed
    // as w - (y[u] + sum of ancestor z's) - (y[v] + sum of ancestor
    // z's), and internal-blossom edges excluded from ever being
    // re-examined once contracted (their tightness is guaranteed
    // structurally at contraction time, not re-verified via slack). That
    // is exactly the per-blossom architecture already scoped above --
    // this attempt narrows *why* a flat dual[] array can never host it,
    // not just that it empirically doesn't.
    //
    // Returns delta applied (0 = no further progress possible, matching is optimal).
    int dualUpdateSerial() {
        int delta = INT_MAX / 2;

        for (int u = 0; u < n; ++u) {
            if (sign[u] != SPLUS) continue;

            // Does this vertex have any negative-slack arc to a SNONE neighbor
            // that could become tight?  Only those arcs should cap delta via
            // the non-negativity constraint on dual[u].
            bool capByThisVertex = false;

            for (int j = ro[u]; j < ro[u+1]; ++j) {
                // Skip non-positive arcs (dual feasibility always satisfied, never tight)
                if (aw[j] <= 0) continue;

                int v  = ci[j];
                int sl = slack(j, u);
                int bv = blossomBase[v];

                if (sign[bv] == SNONE && sl < 0) {
                    // (a) SPLUS→SNONE arc: reducing dual[u] by -sl makes it tight.
                    // This vertex's dual must not go below 0, so cap delta at dual[u].
                    int d = -sl;
                    if (d < delta) delta = d;
                    capByThisVertex = true;
                } else if (sign[bv] == SPLUS && root[bv] != root[u] && sl < 0) {
                    // (b) Cross-tree SPLUS→SPLUS: both duals decrease, slack rises by 2*delta.
                    int d = (-sl + 1) / 2;
                    if (d > 0 && d < delta) delta = d;
                    capByThisVertex = true;
                }
            }

            // (c) Non-negativity: cap delta so dual[u] stays >= 0,
            // but only for vertices with arcs that actually benefit from tightening.
            if (capByThisVertex && dual[u] < delta) delta = dual[u];
        }

        if (delta <= 0 || delta == INT_MAX / 2) return 0;

        for (int v = 0; v < n; ++v) {
            if (sign[v] == SPLUS)  dual[v] -= delta;
            if (sign[v] == SMINUS) dual[v] += delta;
        }
        ++numDual;
        return delta;
    }

    // -----------------------------------------------------------
    // Parallel primal search over the tight-edge subgraph E0.
    //
    // Weighted extension of X-Blossom's parallel recursion-free algorithm
    // (Fan, Lee, Zhang, VLDB 2025, Algorithms 3-5): Blossom VI's key
    // insight (arXiv:2604.20351, Sec 1.4) is that a weighted blossom
    // algorithm's primal phase is *exactly* the unweighted maximum
    // matching problem restricted to zero-slack edges E0 = {e :
    // slack(e) = 0}. This function is that unweighted search, using the
    // same lock-free CAS technique as X-Blossom, with every edge scan
    // additionally gated on the weighted admissibility condition
    // (aw[j] > 0 && slack(j,u) == 0) so that only tight, positive-weight
    // arcs are ever traversed (mirrors augmentRoundSerial's constraints).
    //
    // Level-synchronous (bulk-synchronous parallel): each level of the
    // alternating forest is expanded by `nt` workers, joined, then the
    // next level begins. std::thread::join() is a full memory barrier,
    // so plain (non-atomic) writes to sign/root/par/blossomBase made by
    // one level's workers are safely visible to the next level's workers
    // without additional fences — the same "fork-join" safety argument
    // used by X-Blossom's own path-table update.
    //
    // Synchronization:
    //   - treeClaim: CAS-based tree ownership (X-Blossom Algorithm 4).
    //   - matchClaim: CAS-based ownership of a matched edge during Grow,
    //     keyed by min(endpoint) so at most one thread grows a given
    //     matched pair (X-Blossom Algorithm 5).
    //   - Blossom shrink mutates blossomBase/sign for an unbounded chain
    //     of ancestors (walking to the LCA), which is not safe to do
    //     lock-free with simple per-node CAS. We therefore *detect*
    //     blossom-closing edges in parallel (cheap: append a candidate
    //     pair to a thread-local list, no shared mutation) and *apply*
    //     the actual shrink() calls single-threaded after the parallel
    //     scan joins. This is the one primal operation we deliberately
    //     serialize, per the "shared blossom state" caveat: shrinks are
    //     rare relative to total edge scans, so this does not become a
    //     parallel bottleneck in practice.
    //
    // Known interaction with the pre-existing dual-update bug (see
    // dualUpdateSerial's KNOWN BUG comment): when several equally-tight
    // tree pairings are available in the same phase, which one wins the
    // treeClaim CAS race depends on OS thread scheduling. That's the
    // *intended*, race-free behavior of X-Blossom's lock-free
    // arbitration (Theorem 3 in the paper: any pairing is valid, a
    // missed one is retried next round) and does not by itself affect
    // correctness -- every run of this function still returns a valid
    // matching. However, because the dual-update bug already makes this
    // solver's final *weight* sensitive to which tie-broken choice gets
    // taken (confirmed by hand-tracing a serial example -- same root
    // cause, no threads involved), the two compound: on large/dense
    // graphs with enough concurrent contention, repeated parallel runs
    // of the *same* graph and thread count can converge to slightly
    // different final weights, even though each individual result is a
    // valid matching. Benchmarked example: n=200, ~8000 edges, 4 threads
    // returned 97886 on some runs and 97871 (the serial answer) on
    // others; validity held in every run. test_exact.cpp's
    // runLargeScaleDiagnostic reports this rather than hiding it. Fixing
    // it for real requires fixing the underlying dual-update bug, not
    // this arbitration (removing the race would just make the search
    // serial again).
    bool augmentRoundParallel(SimpleThreadPool& pool) {
        int nt = pool.size();

        std::fill(sign.begin(),    sign.end(),    SNONE);
        std::fill(root.begin(),    root.end(),    NO_PARENT);
        std::fill(par.begin(),     par.end(),     NO_PARENT);
        std::fill(parArc.begin(),  parArc.end(),  NO_PARENT);
        std::iota(blossomBase.begin(), blossomBase.end(), 0);

        for (auto& a : treeClaim)  a.store(0, std::memory_order_relaxed);
        for (auto& a : matchClaim) a.store(0, std::memory_order_relaxed);

        std::vector<int> frontier;
        frontier.reserve(n);
        for (int v = 0; v < n; ++v) {
            if (M[v] == UNMATCHED && dual[v] > 0) {
                sign[v] = SPLUS; root[v] = v;
                frontier.push_back(v);
            }
        }

        bool anyAugmented = false;

        while (!frontier.empty() && !anyAugmented) {
            // ---- Phase A: parallel augmenting-path detection (Alg. 4) ----
            std::vector<std::vector<std::vector<int>>> augsByThread(nt);
            pool.runAdaptive((int)frontier.size(), [&](int t) {
                auto& localAugs = augsByThread[t];
                for (size_t i = t; i < frontier.size(); i += nt) {
                    int u = frontier[i];
                    if (treeClaim[root[u]].load(std::memory_order_relaxed)) continue;
                    for (int j = ro[u]; j < ro[u+1]; ++j) {
                        if (aw[j] <= 0) continue;
                        if (slack(j, u) != 0) continue;
                        int v  = ci[j];
                        int bv = blossomBase[v];
                        if (treeClaim[root[u]].load(std::memory_order_relaxed)) break;

                        if (sign[bv] == SNONE && M[v] == UNMATCHED) {
                            int expected = 0;
                            if (treeClaim[root[u]].compare_exchange_strong(expected, 1)) {
                                std::vector<int> path;
                                path.push_back(v);
                                path.push_back(u);
                                for (int cur = par[u]; cur != NO_PARENT; cur = par[cur])
                                    path.push_back(cur);
                                localAugs.push_back(std::move(path));
                                break;
                            }
                        } else if (sign[bv] == SPLUS && root[bv] != root[u]) {
                            int r1 = std::min(root[u], root[bv]);
                            int r2 = std::max(root[u], root[bv]);
                            int expected = 0;
                            if (treeClaim[r1].compare_exchange_strong(expected, 1)) {
                                expected = 0;
                                if (treeClaim[r2].compare_exchange_strong(expected, 1)) {
                                    std::vector<int> pu = pathToRoot(u);
                                    std::vector<int> pv = pathToRoot(v);
                                    std::vector<int> path;
                                    for (int k = (int)pu.size()-1; k >= 0; --k) path.push_back(pu[k]);
                                    for (int k = 0; k < (int)pv.size(); ++k)    path.push_back(pv[k]);
                                    localAugs.push_back(std::move(path));
                                    break;
                                } else {
                                    treeClaim[r1].store(0, std::memory_order_relaxed);
                                }
                            }
                        }
                    }
                }
            });

            for (auto& localAugs : augsByThread)
                for (auto& path : localAugs) {
                    augmentPath(path);
                    anyAugmented = true;
                }
            if (anyAugmented) break; // apply, then let solve() re-seed from remaining exposed vertices

            // ---- Phase B: parallel tree expansion / Grow (Alg. 5, proc. 1) ----
            std::vector<std::vector<int>> nextByThread(nt);
            pool.runAdaptive((int)frontier.size(), [&](int t) {
                auto& localNext = nextByThread[t];
                for (size_t i = t; i < frontier.size(); i += nt) {
                    int u = frontier[i];
                    if (treeClaim[root[u]].load(std::memory_order_relaxed)) continue;
                    for (int j = ro[u]; j < ro[u+1]; ++j) {
                        if (aw[j] <= 0) continue;
                        if (slack(j, u) != 0) continue;
                        int v  = ci[j];
                        int bv = blossomBase[v];
                        if (sign[bv] != SNONE) continue;   // handled in Phase A/C
                        if (M[v] == UNMATCHED) continue;   // handled in Phase A

                        int mate = M[v];
                        int edgeId = std::min(v, mate);
                        int expected = 0;
                        if (matchClaim[edgeId].compare_exchange_strong(expected, 1)) {
                            sign[v] = SMINUS; root[v] = root[u]; par[v] = u; parArc[v] = j;
                            sign[mate] = SPLUS; root[mate] = root[u]; par[mate] = v;
                            localNext.push_back(mate);
                        }
                    }
                }
            });

            // ---- Phase C: parallel blossom detection, serial shrink (Alg. 5, proc. 2) ----
            std::vector<std::vector<std::pair<int,int>>> blossomEventsByThread(nt);
            pool.runAdaptive((int)frontier.size(), [&](int t) {
                auto& localEvents = blossomEventsByThread[t];
                for (size_t i = t; i < frontier.size(); i += nt) {
                    int u = frontier[i];
                    if (treeClaim[root[u]].load(std::memory_order_relaxed)) continue;
                    for (int j = ro[u]; j < ro[u+1]; ++j) {
                        if (aw[j] <= 0) continue;
                        if (slack(j, u) != 0) continue;
                        int v  = ci[j];
                        int bv = blossomBase[v];
                        if (sign[bv] == SPLUS && root[bv] == root[u])
                            localEvents.emplace_back(u, bv);
                    }
                }
            });

            std::queue<int> relabeled;
            for (auto& localEvents : blossomEventsByThread)
                for (auto& pr : localEvents)
                    shrink(pr.first, pr.second, relabeled);
            while (!relabeled.empty()) {
                nextByThread[0].push_back(relabeled.front());
                relabeled.pop();
            }

            frontier.clear();
            for (auto& localNext : nextByThread)
                frontier.insert(frontier.end(), localNext.begin(), localNext.end());
        }

        return anyAugmented;
    }

    // Parallel dual update: the per-vertex delta computation and the
    // subsequent dual-variable adjustment are both embarrassingly
    // parallel (read-only scans / disjoint writes over vertices), so we
    // simply stripe the vertex range across threads. Semantically
    // identical to dualUpdateSerial(); split out so numThreads actually
    // does something in the phase that used to ignore it entirely.
    int dualUpdateParallel(SimpleThreadPool& pool) {
        int nt = pool.size();
        std::vector<int> localDelta(nt, INT_MAX / 2);
        pool.runAdaptive(n, [&](int t) {
            int delta = INT_MAX / 2;
            for (int u = t; u < n; u += nt) {
                if (sign[u] != SPLUS) continue;
                bool capByThisVertex = false;
                for (int j = ro[u]; j < ro[u+1]; ++j) {
                    if (aw[j] <= 0) continue;
                    int v  = ci[j];
                    int sl = slack(j, u);
                    int bv = blossomBase[v];
                    if (sign[bv] == SNONE && sl < 0) {
                        // (a)
                        int d = -sl;
                        if (d < delta) delta = d;
                        capByThisVertex = true;
                    } else if (sign[bv] == SPLUS && root[bv] != root[u] && sl < 0) {
                        // (b) cross-tree SPLUS-SPLUS -- see
                        // dualUpdateSerial's KNOWN BUG comment: this
                        // (unmodified, original) formula is what the
                        // parallel path must match for consistency.
                        int d = (-sl + 1) / 2;
                        if (d > 0 && d < delta) delta = d;
                        capByThisVertex = true;
                    }
                }
                if (capByThisVertex && dual[u] < delta) delta = dual[u]; // (c)
            }
            localDelta[t] = delta;
        });

        int delta = INT_MAX / 2;
        for (int d : localDelta) delta = std::min(delta, d);
        if (delta <= 0 || delta == INT_MAX / 2) return 0;

        pool.runAdaptive(n, [&](int t) {
            for (int v = t; v < n; v += nt) {
                if (sign[v] == SPLUS)  dual[v] -= delta;
                if (sign[v] == SMINUS) dual[v] += delta;
            }
        });
        ++numDual;
        return delta;
    }

    // Main loop: alternate augmentation rounds and dual updates.
    // numThreads <= 1 routes to the serial reference path; numThreads > 1
    // uses the parallel primal search / dual update above, driven by a
    // SimpleThreadPool created once here and reused for every phase of
    // every round (see the pool's own comment for why that matters).
    // Both paths are algorithmically equivalent (same admissibility
    // rules, same blossom shrink), which is what test_exact.cpp's
    // parallel-vs-serial consistency tests verify.
    void solve(int numThreads) {
        int maxRounds = n * n + n + 1; // O(n²) bound
        bool parallel = numThreads > 1;
        std::unique_ptr<SimpleThreadPool> pool;
        if (parallel) pool = std::make_unique<SimpleThreadPool>(numThreads);
        for (int r = 0; r < maxRounds; ++r) {
            bool progressed = parallel ? augmentRoundParallel(*pool) : augmentRoundSerial();
            if (!progressed) {
                int d = parallel ? dualUpdateParallel(*pool) : dualUpdateSerial();
                if (d <= 0) break;
            }
        }
    }
};

} // namespace mwm

// ===========================================================
//  Main hybrid solver
// ===========================================================

MatchingResult hybrid_blossom_maximum_weight_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& colIndices,
    const std::vector<int>& arcWeights,
    int numThreads)
{
    auto t_start = std::chrono::steady_clock::now();

    MatchingResult res;
    res.num_vertices = (int)rowOffsets.size() - 1;
    res.num_edges    = (int)colIndices.size() / 2;
    res.num_threads  = numThreads;

    int n = res.num_vertices;

    if (n == 0) {
        res.valid = true; res.optimal = false; res.weight = 0;
        return res;
    }

    // ----- Phase 1: X-Blossom warm-start -----
    std::vector<int> M;
    run_xblossom_phase(rowOffsets, colIndices, M, numThreads, res);
    // Remove negative-weight edges from warm-start result
    remove_negative_edges(M, rowOffsets, colIndices, arcWeights);
    res.initial_cardinality = count_matching(M);
    res.initial_weight      = matching_weight(M, rowOffsets, colIndices, arcWeights);

    // ----- Phase 2: Exact primal-dual weighted matching -----
    {
        auto t0 = std::chrono::steady_clock::now();

        // Always start the MWM solver from an empty matching.
        // The X-Blossom warm-start is incompatible with the dual initialization:
        // the greedy matching has non-tight edges under the MWM duals, which would
        // cause the dual update to terminate prematurely.  The primal-dual loop
        // finds the optimal matching from scratch using tight-arc BFS + delta updates.
        mwm::MWMSolver solver(n, rowOffsets, colIndices, arcWeights);
        solver.solve(numThreads);

        M = solver.M;
        res.num_augmentations        += solver.numAug;
        res.num_blossom_contractions += solver.numShrink;
        res.num_dual_updates         += solver.numDual;

        auto t1 = std::chrono::steady_clock::now();
        res.time_weighted_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    }

    // Remove any remaining negative edges (shouldn't occur after exact solver)
    remove_negative_edges(M, rowOffsets, colIndices, arcWeights);

    // ----- Phase 3: Validate -----
    {
        auto t0 = std::chrono::steady_clock::now();
        int64_t w = 0;
        res.valid  = validate_matching(rowOffsets, colIndices, arcWeights, M, w);
        res.weight = w;
        res.final_cardinality = count_matching(M);
        auto t1 = std::chrono::steady_clock::now();
        res.time_validate_ms =
            std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    }

    res.mate    = std::move(M);
    res.optimal = false; // set by test harness after brute-force comparison

    auto t_end = std::chrono::steady_clock::now();
    res.time_total_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    return res;
}

// Legacy void API
void hybrid_blossom_maximum_weight_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& colIndices,
    const std::vector<int>& arcWeights,
    std::vector<int>& M,
    int numThreads)
{
    MatchingResult res = hybrid_blossom_maximum_weight_matching(
        rowOffsets, colIndices, arcWeights, numThreads);
    M = std::move(res.mate);
    std::cout << "Hybrid (" << numThreads << " threads): "
              << res.final_cardinality << " edges"
              << ", weight=" << res.weight
              << ", time=" << res.time_total_ms << "ms"
              << ", valid=" << (res.valid ? "yes" : "NO") << "\n";
}
