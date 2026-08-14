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
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iostream>
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
//  X-Blossom warm-start (parallel BFS)
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
    bool augmentRound() {
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
    // Returns delta applied (0 = no further progress possible, matching is optimal).
    int dualUpdate() {
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

    // Main loop: alternate augmentation rounds and dual updates.
    void solve() {
        int maxRounds = n * n + n + 1; // O(n²) bound
        for (int r = 0; r < maxRounds; ++r) {
            if (!augmentRound()) {
                if (dualUpdate() <= 0) break;
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
        solver.solve();

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
