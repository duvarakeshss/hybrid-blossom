/// \file test_exact.cpp
/// \brief Exact correctness testing for hybrid maximum-weight matching.
///
/// Contains:
///   1. Brute-force exact MWM solver (for n <= 16)
///   2. Deterministic unit tests
///   3. Randomized test generator
///   4. Reference comparison (Hybrid vs. Brute-force)
///
/// Build (from hybrid-x-blossom/ directory):
///   g++ -std=c++17 -O2 -pthread -o test_exact \
///       test_exact.cpp hybrid_blossom.cpp parallel_augment.cpp \
///       weighted_matching.cpp
///
/// Run:   ./test_exact
///        ./test_exact --quick

#include "hybrid_blossom.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

// ============================================================
//  Edge type alias
// ============================================================
using Edge3 = std::tuple<int,int,int>; // (u, v, weight)

inline int  eu(const Edge3& e) { return std::get<0>(e); }
inline int  ev(const Edge3& e) { return std::get<1>(e); }
inline int  ew(const Edge3& e) { return std::get<2>(e); }

// ============================================================
//  Brute-force exact MWM (exponential, n <= 16)
// ============================================================

struct SmallGraph {
    int n;
    std::vector<Edge3> edges; // undirected

    int64_t maxWeightMatching() const {
        int m = (int)edges.size();
        int64_t best = 0; // empty matching is always feasible

        for (int mask = 1; mask < (1 << m); ++mask) {
            std::vector<bool> used(n, false);
            bool valid = true;
            int64_t w = 0;

            for (int i = 0; i < m && valid; ++i) {
                if (!(mask & (1 << i))) continue;
                int u = eu(edges[i]), v = ev(edges[i]), wt = ew(edges[i]);
                if (used[u] || used[v]) { valid = false; break; }
                used[u] = used[v] = true;
                w += wt;
            }
            if (valid) best = std::max(best, w);
        }
        return best;
    }
};

/// Convert CSR to SmallGraph (keep only edges u < v to avoid duplicates)
static SmallGraph toSmallGraph(int n,
                               const std::vector<int>& ro,
                               const std::vector<int>& ci,
                               const std::vector<int>& aw) {
    SmallGraph g;
    g.n = n;
    for (int u = 0; u < n; ++u)
        for (int j = ro[u]; j < ro[u+1]; ++j) {
            int v = ci[j];
            if (u < v) g.edges.push_back(std::make_tuple(u, v, aw[j]));
        }
    return g;
}

// ============================================================
//  Graph builder: edge list → CSR
// ============================================================

static void buildCSR(int n,
                     const std::vector<Edge3>& edges,
                     std::vector<int>& ro,
                     std::vector<int>& ci,
                     std::vector<int>& aw)
{
    std::vector<std::vector<std::pair<int,int>>> adj(n);
    for (int i = 0; i < (int)edges.size(); ++i) {
        int u = eu(edges[i]), v = ev(edges[i]), w = ew(edges[i]);
        if (u == v) continue;
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }
    ro.assign(n + 1, 0);
    for (int i = 0; i < n; ++i)
        ro[i+1] = ro[i] + (int)adj[i].size();
    ci.resize(ro[n]);
    aw.resize(ro[n]);
    for (int i = 0; i < n; ++i)
        for (int k = 0; k < (int)adj[i].size(); ++k) {
            ci[ro[i]+k] = adj[i][k].first;
            aw[ro[i]+k] = adj[i][k].second;
        }
}

// ============================================================
//  Test statistics
// ============================================================

struct TestStats {
    int total = 0, passed = 0, failed = 0;

    void record(bool ok) { ++total; if (ok) ++passed; else ++failed; }

    void print() const {
        std::cout << "\n=== Test Results ===\n"
                  << "Total:  " << total  << "\n"
                  << "Passed: " << passed << "\n"
                  << "Failed: " << failed << "\n";
    }
};

// ============================================================
//  Core test runner
// ============================================================

static bool runTest(const std::string& name,
                    int n,
                    const std::vector<Edge3>& edges,
                    int64_t expectedWeight, // LLONG_MIN → compute from brute-force
                    TestStats& stats,
                    int numThreads = 1)
{
    std::vector<int> ro, ci, aw;
    buildCSR(n, edges, ro, ci, aw);

    // Brute-force reference
    int64_t bruteForce;
    if (expectedWeight == LLONG_MIN) {
        if (n > 16) {
            std::cout << "[SKIP] " << name << " (n=" << n << " > 16)\n";
            return true;
        }
        SmallGraph g = toSmallGraph(n, ro, ci, aw);
        bruteForce = g.maxWeightMatching();
    } else {
        bruteForce = expectedWeight;
    }

    // Hybrid solver
    MatchingResult res = hybrid_blossom_maximum_weight_matching(ro, ci, aw, numThreads);

    bool weightOk   = (res.weight == bruteForce);
    bool matchingOk = res.valid;
    bool ok         = matchingOk && weightOk;

    if (!ok) {
        std::cout << "[FAIL] " << name
                  << "  expected=" << bruteForce
                  << "  got=" << res.weight
                  << "  valid=" << (matchingOk ? "yes" : "NO") << "\n";
        if (n <= 8) {
            std::cout << "       Edges:";
            for (int i = 0; i < (int)edges.size(); ++i)
                std::cout << " (" << eu(edges[i]) << "-" << ev(edges[i])
                          << ",w=" << ew(edges[i]) << ")";
            std::cout << "\n";
            std::cout << "       Mate:";
            for (int v = 0; v < n; ++v)
                std::cout << " " << v << "->" << res.mate[v];
            std::cout << "\n";
        }
    } else {
        std::cout << "[PASS] " << name
                  << "  weight=" << res.weight << "\n";
    }
    stats.record(ok);
    return ok;
}

// ============================================================
//  Deterministic unit tests
// ============================================================

static void runUnitTests(TestStats& stats) {
    std::cout << "\n--- Unit Tests ---\n";

    // 1. Empty graph
    {
        std::vector<int> ro = {0}, ci, aw;
        MatchingResult res = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 1);
        bool ok = res.valid && res.weight == 0 && res.mate.empty();
        std::cout << (ok ? "[PASS]" : "[FAIL]") << " Empty graph (0 vertices)\n";
        stats.record(ok);
    }

    // 2. Single vertex
    runTest("Single vertex, no edges", 1, {}, 0, stats);

    // 3. Two vertices, positive weight
    runTest("2v positive w=10",   2, {std::make_tuple(0,1,10)},  10, stats);

    // 4. Two vertices, negative weight — should NOT match
    runTest("2v negative w=-10",  2, {std::make_tuple(0,1,-10)},  0, stats);

    // 5. Two vertices, zero weight — empty matching is equally optimal
    runTest("2v zero weight",     2, {std::make_tuple(0,1,0)},    0, stats);

    // 6. Triangle, all weight 5 — best is one edge = 5
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,5),
                                std::make_tuple(1,2,5),
                                std::make_tuple(0,2,5)};
        runTest("Triangle all w=5", 3, e, 5, stats);
    }

    // 7. Triangle: max-cardinality != max-weight
    //    Edge(0,1)=100 >> others; only one match possible in triangle
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,100),
                                std::make_tuple(1,2,1),
                                std::make_tuple(0,2,1)};
        runTest("Triangle max-card != max-weight", 3, e, 100, stats);
    }

    // 8. Path 4 nodes, heaviest edge in middle
    //    0-1(w=1), 1-2(w=100), 2-3(w=1) → match {1,2} only, weight=100
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,1),
                                std::make_tuple(1,2,100),
                                std::make_tuple(2,3,1)};
        runTest("Path-4 heavy middle edge", 4, e, 100, stats);
    }

    // 9. Two disjoint edges
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,5), std::make_tuple(2,3,7)};
        runTest("Two disjoint edges", 4, e, 12, stats);
    }

    // 10. C5 (odd cycle), all w=1 — can match at most 2 edges
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,1), std::make_tuple(1,2,1),
                                std::make_tuple(2,3,1), std::make_tuple(3,4,1),
                                std::make_tuple(4,0,1)};
        runTest("C5 odd cycle all w=1", 5, e, 2, stats);
    }

    // 11. Disconnected graph — two components matched independently
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,10), std::make_tuple(2,3,20)};
        runTest("Disconnected 2 components", 4, e, 30, stats);
    }

    // 12. Star K_{1,4}: center=0, only 1 edge can be matched
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,1), std::make_tuple(0,2,1),
                                std::make_tuple(0,3,1), std::make_tuple(0,4,1)};
        runTest("Star K_1,4 all w=1", 5, e, 1, stats);
    }

    // 13. K4 with varying weights — best matching: (0,3)=10 + (1,2)=9 = 19
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,5), std::make_tuple(0,2,6),
                                std::make_tuple(0,3,10), std::make_tuple(1,2,9),
                                std::make_tuple(1,3,4), std::make_tuple(2,3,3)};
        runTest("K4 varying weights", 4, e, 19, stats);
    }

    // 14. All-negative weights → empty matching
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,-5), std::make_tuple(1,2,-3),
                                std::make_tuple(2,3,-1), std::make_tuple(0,3,-4)};
        runTest("All negative weights", 4, e, 0, stats);
    }

    // 15. Mixed weights — best is single edge (0,1)=5
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,5),
                                std::make_tuple(0,2,-3),
                                std::make_tuple(1,2,2)};
        runTest("Mixed weights triangle", 3, e, 5, stats);
    }

    // 16. Odd n=3, one heavy edge
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,100), std::make_tuple(1,2,1)};
        runTest("Odd n=3 one heavy edge", 3, e, 100, stats);
    }

    // 17. Isolated vertex: n=3, only edge (0,1)
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,7)};
        runTest("Isolated vertex n=3", 3, e, 7, stats);
    }

    // 18. Weighted conflict: must choose between (0,2)=9 or (0,1)=5+(1,2) impossible
    //     Graph: 0-1(w=5), 0-2(w=9). Best = match (0,2)=9
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,5), std::make_tuple(0,2,9)};
        runTest("Weighted conflict", 3, e, 9, stats);
    }

    // 19. Perfect matching exists and is optimal
    //     4-cycle: 0-1(w=1), 1-2(w=10), 2-3(w=1), 3-0(w=10) → {1-2, 3-0} = 20
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,1), std::make_tuple(1,2,10),
                                std::make_tuple(2,3,1), std::make_tuple(3,0,10)};
        runTest("4-cycle perfect matching", 4, e, 20, stats);
    }

    // 20. Thread consistency: 1 vs 4 threads on K4
    {
        std::vector<Edge3> e = {std::make_tuple(0,1,5), std::make_tuple(0,2,6),
                                std::make_tuple(0,3,10), std::make_tuple(1,2,9),
                                std::make_tuple(1,3,4), std::make_tuple(2,3,3)};
        std::vector<int> ro, ci, aw;
        buildCSR(4, e, ro, ci, aw);
        MatchingResult r1 = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 1);
        MatchingResult r4 = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 4);
        bool ok = r1.valid && r4.valid && r1.weight == r4.weight;
        std::cout << (ok ? "[PASS]" : "[FAIL]")
                  << " Thread consistency (1 vs 4 threads)"
                  << "  w1=" << r1.weight << " w4=" << r4.weight << "\n";
        stats.record(ok);
    }
}

// ============================================================
//  Randomized tests
// ============================================================

static std::vector<Edge3> randomGraph(int n, int seed, int maxEdges, int minW, int maxW) {
    std::mt19937 rng(static_cast<unsigned>(seed));
    std::uniform_int_distribution<int> wDist(minW, maxW);
    std::uniform_int_distribution<int> vDist(0, n > 1 ? n-1 : 0);

    std::vector<std::vector<bool>> has(n, std::vector<bool>(n, false));
    std::vector<Edge3> edges;

    int attempts = maxEdges * 4;
    for (int i = 0; i < attempts && (int)edges.size() < maxEdges; ++i) {
        int u = vDist(rng), v = vDist(rng);
        if (u == v || has[u][v]) continue;
        has[u][v] = has[v][u] = true;
        edges.push_back(std::make_tuple(u, v, wDist(rng)));
    }
    return edges;
}

static void runRandomTests(TestStats& stats, int numTests) {
    std::cout << "\n--- Randomized Tests (" << numTests << ") ---\n";

    struct Category {
        std::string name;
        int n, maxE, minW, maxW;
    };

    std::vector<Category> cats = {
        {"Dense+4",    4,  6,   1,  20},
        {"Dense+6",    6,  10,  1,  50},
        {"Dense+8",    8,  15,  1, 100},
        {"Sparse+10",  10,  8,  1,  50},
        {"Sparse+12",  12, 10,  1, 100},
        {"Mixed4",     4,   5, -10, 10},
        {"Mixed6",     6,   8, -20, 20},
        {"Mixed8",     8,  12,  -5, 15},
        {"AllNeg4",    4,   6, -20, -1},
        {"AllNeg6",    6,   8, -10, -1},
        {"OddN5",      5,   8,   1, 30},
        {"OddN7",      7,  10,   1, 30},
        {"Sparse10",   10,  4,   1, 50},
    };

    int perCat = std::max(1, numTests / (int)cats.size());

    for (int ci = 0; ci < (int)cats.size(); ++ci) {
        const Category& cat = cats[ci];
        for (int i = 0; i < perCat; ++i) {
            int seed = i * 7919 + ci * 1337;
            auto edges = randomGraph(cat.n, seed, cat.maxE, cat.minW, cat.maxW);
            std::string name = cat.name + "/s" + std::to_string(seed);
            runTest(name, cat.n, edges, LLONG_MIN, stats);
        }
    }
}

// ============================================================
//  Thread consistency tests
// ============================================================

static void runThreadTests(TestStats& stats) {
    std::cout << "\n--- Thread Consistency Tests ---\n";
    for (int seed = 0; seed < 10; ++seed) {
        auto edges = randomGraph(6, seed * 999, 10, 1, 50);
        std::vector<int> ro, ci, aw;
        buildCSR(6, edges, ro, ci, aw);

        MatchingResult r1 = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 1);
        MatchingResult r2 = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 2);
        MatchingResult r4 = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 4);
        MatchingResult r8 = hybrid_blossom_maximum_weight_matching(ro, ci, aw, 8);

        bool ok = r1.valid && r2.valid && r4.valid && r8.valid &&
                  r1.weight == r2.weight &&
                  r1.weight == r4.weight &&
                  r1.weight == r8.weight;

        std::cout << (ok ? "[PASS]" : "[FAIL]")
                  << " seed=" << seed
                  << "  w1=" << r1.weight
                  << " w2=" << r2.weight
                  << " w4=" << r4.weight
                  << " w8=" << r8.weight << "\n";
        stats.record(ok);
    }
}

// ============================================================
//  main
// ============================================================

int main(int argc, char* argv[]) {
    bool quick = (argc >= 2 && std::string(argv[1]) == "--quick");
    int numRandom = quick ? 50 : 400;

    std::cout << "=== Hybrid Blossom Exact Correctness Tests ===\n";
    std::cout << "Mode: " << (quick ? "quick" : "full") << "\n\n";

    TestStats stats;
    runUnitTests(stats);
    runRandomTests(stats, numRandom);
    runThreadTests(stats);
    stats.print();

    return (stats.failed == 0) ? 0 : 1;
}
