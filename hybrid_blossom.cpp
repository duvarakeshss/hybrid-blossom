#include "hybrid_blossom.h"
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <climits>

/// \brief Hybrid solver state
struct State {
    int n;
    const std::vector<int>& ro;
    const std::vector<int>& ci;
    const std::vector<int>& ew;

    std::vector<int> M;           // matching
    std::vector<int> dual;        // dual vars (amortized)
    std::vector<int> slack;      // slack per edge

    std::mutex mtx;

    State(const std::vector<int>& rowOffsets,
          const std::vector<int>& columnIndices,
          const std::vector<int>& edgeWeights)
        : ro(rowOffsets), ci(columnIndices), ew(edgeWeights)
    {
        n = static_cast<int>(ro.size()) - 1;
        int m = static_cast<int>(ci.size()) / 2;
        M.assign(n, -1);
        dual.assign(n, 0);
        slack.assign(m, 0);
    }
};

/// \brief Find exposed nodes (M[v] == -1)
static void find_exposed(const std::vector<int>& M, std::vector<int>& out) {
    out.clear();
    for (int i = 0; i < static_cast<int>(M.size()); ++i) {
        if (M[i] == -1) out.push_back(i);
    }
}

/// \brief Init: make slacks non-negative (Blossom VI Init)
///
/// For each vertex v, dual[v] = 2 * min(edge weight)
static void init_slacks(State& s, int threads) {
    int chunk = (s.n + threads - 1) / threads;
    std::vector<std::thread> th;

    for (int t = 0; t < threads; ++t) {
        th.emplace_back([&, t, chunk]() {
            int start = t * chunk;
            int end = std::min(start + chunk, s.n);
            for (int v = start; v < end; ++v) {
                int min_w = INT32_MAX;
                for (int j = s.ro[v]; j < s.ro[v + 1]; ++j) {
                    int edge = j / 2;
                    if (edge < static_cast<int>(s.ew.size())) {
                        min_w = std::min(min_w, s.ew[edge]);
                    }
                }
                s.dual[v] = 4 * min_w;
            }
        });
    }
    for (auto& t : th) t.join();
}

/// \brief Try to match two exposed nodes via tight edge
///
/// slack = dual[v] + dual[w] - 2*weight
/// tight when slack == 0
static void match_exposed(State& s) {
    for (int v = 0; v < s.n; ++v) {
        if (s.M[v] != -1) continue;
        for (int j = s.ro[v]; j < s.ro[v + 1]; ++j) {
            int w = s.ci[j];
            if (s.M[w] != -1) continue;
            int edge = j / 2;
            if (edge < static_cast<int>(s.ew.size())) {
                int slack = 4 * s.ew[edge] - s.dual[v] - s.dual[w];
                if (slack == 0) {
                    s.M[v] = w;
                    s.M[w] = v;
                    break; // only break inner loop — match ALL pairs
                }
            }
        }
    }
}

/// \brief Increase dual of exposed nodes by min slack
///
/// Blossom VI: for each tree, find min slack → increase dual by that
/// \brief Increase dual of exposed nodes by min slack
///
/// Blossom VI: for each tree, dual[v] += min_slack where slack = 4*weight - dual
static void update_duals(State& s) {
    for (int v = 0; v < s.n; ++v) {
        if (s.M[v] != -1) continue;
        int min_s = INT32_MAX;
        for (int j = s.ro[v]; j < s.ro[v + 1]; ++j) {
            int edge = j / 2;
            if (edge < static_cast<int>(s.ew.size())) {
                // Same slack formula as match_exposed
                int sl = 4 * s.ew[edge] - s.dual[v] - s.dual[s.ci[j]];
                min_s = std::min(min_s, sl);
            }
        }
        s.dual[v] += min_s;
    }
}

/// \brief Validate matching
static bool valid(const std::vector<int>& M) {
    int n = static_cast<int>(M.size());
    for (int v = 0; v < n; ++v) {
        if (M[v] != -1) {
            if (M[v] < 0 || M[v] >= n) return false;
            if (M[M[v]] != v) return false;
        }
    }
    return true;
}

/// \brief Count matching edges
static int count(const std::vector<int>& M) {
    int c = 0;
    int n = static_cast<int>(M.size());
    for (int v = 0; v < n; ++v) {
        if (M[v] != -1 && v < M[v]) ++c;
    }
    return c;
}

/// \brief Main hybrid solver
void hybrid_blossom_maximum_weight_matching(
        const std::vector<int>& rowOffsets,
        const std::vector<int>& columnIndices,
        const std::vector<int>& edgeWeights,
        std::vector<int>& M,
        int num_threads) {

    State s(rowOffsets, columnIndices, edgeWeights);
    auto t0 = std::chrono::steady_clock::now();

    // 1. Init slacks
    init_slacks(s, num_threads);

    // 2. Match exposed with tight edges
    match_exposed(s);

    // 3. Find remaining exposed
    std::vector<int> exposed;
    find_exposed(s.M, exposed);

    // 4. While exposed remain: increase duals until new tight edge found
    int iter = 0;
    while (!exposed.empty() && iter < s.n) {
        ++iter;

        // Increase dual of each exposed node until an incident edge is tight
        // Blossom VI: for each tree, dual[v] += min_slack where slack = 4*weight - dual
        for (int v : exposed) {
            int min_s = INT32_MAX;
            for (int j = s.ro[v]; j < s.ro[v + 1]; ++j) {
                int w = s.ci[j];
                // Consider all edges — min slack to any neighbor
                // (Blossom VI: each vertex pays min slack to all incident edges)
                int edge = j / 2;
                if (edge < static_cast<int>(s.ew.size())) {
                    int sl = 4 * s.ew[edge] - s.dual[v] - s.dual[w];
                    min_s = std::min(min_s, sl);
                }
            }
            if (min_s < INT32_MAX) {
                s.dual[v] += min_s;
            }
        }

        // Try to match exposed nodes via now-tight edges
        for (int v : exposed) {
            if (s.M[v] != -1) continue;
            for (int j = s.ro[v]; j < s.ro[v + 1]; ++j) {
                int w = s.ci[j];
                if (s.M[w] != -1) continue;
                int edge = j / 2;
                if (edge < static_cast<int>(s.ew.size())) {
                    int slack = 4 * s.ew[edge] - s.dual[v] - s.dual[w];
                    if (slack == 0) {
                        s.M[v] = w;
                        s.M[w] = v;
                        break;
                    }
                }
            }
        }

        find_exposed(s.M, exposed);
    }

    auto t1 = std::chrono::steady_clock::now();
    int ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    M = s.M;
    bool ok = valid(M);

    std::cout << "Hybrid (" << num_threads << " threads): "
              << count(M) << " edges in " << ms << "ms"
              << " | valid: " << (ok ? "yes" : "no") << "\n";
}