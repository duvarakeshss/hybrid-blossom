/// \file benchmark.cpp
/// \brief Serial-vs-parallel benchmark for the weighted hybrid Blossom solver.
///
/// Measures, for a matrix of (graph size x density x thread count):
///   - wall-clock runtime
///   - matching weight and cardinality (to confirm parallel doesn't
///     trade correctness for speed -- weight must match the serial run)
///   - speedup and parallel efficiency relative to the 1-thread run
///   - number of augmentations and blossom contractions performed
///
/// This exercises exactly the code path this project's "hybrid" claim
/// rests on: mwm::MWMSolver's primal search (augmentRoundParallel) and
/// dual update (dualUpdateParallel), reached through the public
/// hybrid_blossom_maximum_weight_matching() API with numThreads > 1.
///
/// Build:
///   g++ -std=c++17 -O3 -pthread -o benchmark \
///       benchmark.cpp hybrid_blossom.cpp parallel_augment.cpp weighted_matching.cpp
///
/// Run:
///   ./benchmark            (prints a Markdown table to stdout)

#include "hybrid_blossom.h"
#include <chrono>
#include <cstdio>
#include <random>
#include <tuple>
#include <vector>

using Edge3 = std::tuple<int,int,int>;

static std::vector<Edge3> randomGraph(int n, int seed, int maxEdges, int minW, int maxW) {
    std::mt19937 rng(static_cast<unsigned>(seed));
    std::uniform_int_distribution<int> wDist(minW, maxW);
    std::uniform_int_distribution<int> vDist(0, n > 1 ? n - 1 : 0);
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

static void buildCSR(int n, const std::vector<Edge3>& edges,
                      std::vector<int>& ro, std::vector<int>& ci, std::vector<int>& aw) {
    std::vector<std::vector<std::pair<int,int>>> adj(n);
    for (auto& e : edges) {
        int u = std::get<0>(e), v = std::get<1>(e), w = std::get<2>(e);
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }
    ro.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) ro[i + 1] = ro[i] + (int)adj[i].size();
    ci.resize(ro[n]);
    aw.resize(ro[n]);
    for (int i = 0; i < n; ++i)
        for (size_t k = 0; k < adj[i].size(); ++k) {
            ci[ro[i] + k] = adj[i][k].first;
            aw[ro[i] + k] = adj[i][k].second;
        }
}

struct BenchRow {
    int n, m, threads;
    double ms;
    int64_t weight;
    int cardinality;
    int augmentations;
    int shrinks;
};

// Runs `repeats` times and keeps the fastest -- at these graph sizes a
// single measurement is dominated by OS scheduling jitter (background
// processes, first-touch memory effects, thread wake latency variance),
// so min-of-N is the standard way to get a trustworthy wall-clock number
// without needing a much bigger (and much slower to benchmark) graph.
static BenchRow runOnce(const std::vector<int>& ro, const std::vector<int>& ci,
                         const std::vector<int>& aw, int n, int threads, int repeats = 3) {
    BenchRow best{};
    double bestMs = 1e18;
    for (int r = 0; r < repeats; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        MatchingResult res = hybrid_blossom_maximum_weight_matching(ro, ci, aw, threads);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        if (ms < bestMs) {
            bestMs = ms;
            best = BenchRow{n, (int)ci.size() / 2, threads, ms, res.weight,
                             res.final_cardinality, res.num_augmentations,
                             res.num_blossom_contractions};
        }
    }
    return best;
}

int main() {
    struct SizeCfg { int n, maxE, minW, maxW; const char* densityLabel; };
    std::vector<SizeCfg> sizes = {
        {200,   400,  1, 1000, "sparse"},
        {200,  8000,  1, 1000, "dense"},
        {500,  1200,  1, 1000, "sparse"},
        {500, 40000,  1, 1000, "dense"},
        {1000, 2500,  1, 1000, "sparse"},
        {2000, 5000,  1, 1000, "sparse"},
        {800, 100000, 1, 1000, "dense"},
    };
    std::vector<int> threadCounts = {1, 2, 4, 8};

    std::printf("| n | edges | density | threads | time (ms) | speedup | efficiency | weight | cardinality | augmentations | shrinks |\n");
    std::printf("|---|-------|---------|---------|-----------|---------|------------|--------|-------------|---------------|---------|\n");

    for (auto& cfg : sizes) {
        auto edges = randomGraph(cfg.n, cfg.n * 7919 + cfg.maxE, cfg.maxE, cfg.minW, cfg.maxW);
        std::vector<int> ro, ci, aw;
        buildCSR(cfg.n, edges, ro, ci, aw);

        double baselineMs = 0.0;
        int64_t baselineWeight = 0;
        for (int t : threadCounts) {
            BenchRow row = runOnce(ro, ci, aw, cfg.n, t);
            if (t == 1) { baselineMs = row.ms; baselineWeight = row.weight; }
            double speedup = baselineMs > 0 ? baselineMs / row.ms : 1.0;
            double efficiency = speedup / t;
            bool weightOk = row.weight == baselineWeight;
            std::printf("| %d | %d | %s | %d | %.2f | %.2fx | %.1f%% | %lld%s | %d | %d | %d |\n",
                        row.n, row.m, cfg.densityLabel, t, row.ms, speedup, efficiency * 100.0,
                        (long long)row.weight, weightOk ? "" : " (MISMATCH!)",
                        row.cardinality, row.augmentations, row.shrinks);
        }
    }
    return 0;
}
