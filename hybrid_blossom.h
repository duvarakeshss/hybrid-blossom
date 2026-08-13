#ifndef HYBRID_BLOSSOM_H
#define HYBRID_BLOSSOM_H

#include <vector>
#include <mutex>
#include <atomic>
#include <thread>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <cstdint>

/// \brief Edge in CSR format
struct HybridEdge {
    int head;
    int tail;
    int weight;
    int slack_quadrupled;

    HybridEdge(int h, int t, int w) : head(h), tail(t), weight(w),
                                       slack_quadrupled(4 * w) {}
};

/// \brief Node state
struct HybridNode {
    int tree;
    int matched_edge;
    int minus_parent;
    int old_tree;
    int old_matched;
    bool is_even;
    bool is_alive;
    bool plus;
    bool old_plus;
    int receptacle;

    HybridNode() : tree(-1), matched_edge(-1), minus_parent(-1),
                   old_tree(-1), old_matched(-1), is_even(false),
                   is_alive(true), plus(false), old_plus(false),
                   receptacle(-1) {}
};

/// \brief Timer
struct StopWatch {
    std::chrono::milliseconds total{0};
    std::chrono::milliseconds augment{0};
    std::chrono::milliseconds expand{0};
    std::chrono::milliseconds dual{0};

    void reset() {
        total = std::chrono::milliseconds(0);
        augment = std::chrono::milliseconds(0);
        expand = std::chrono::milliseconds(0);
        dual = std::chrono::milliseconds(0);
    }
};

/// \brief Main API — hybrid maximum-weight matching
///
/// Uses X-Blossom's parallel path-table search within Blossom VI's
/// dual-variable framework. This is a weighted (not unweighted) solver.
void hybrid_blossom_maximum_weight_matching(
        const std::vector<int>& rowOffsets,
        const std::vector<int>& columnIndices,
        const std::vector<int>& edgeWeights,
        std::vector<int>& M,
        int numThreads);

#endif //HYBRID_BLOSSOM_H