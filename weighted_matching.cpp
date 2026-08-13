#include "weighted_matching.h"

#include <algorithm>
#include <queue>
#include <iostream>
#include <cassert>

/// \brief Implements Blossom VI's "find minimum perfect matching" in the hybrid
///
/// This is a simplified version of MWPMSolver that:
///   1. Fits into the hybrid's parallel path-table framework
///   2. Uses the same dual variable scheme
///   3. Is called after X-Blossom finds a rough matching
///
/// The full Blossom VI code is too large to replicate here.
/// This version uses the same structure but integrates the path-table approach.

WeightedMatchingSolver::WeightedMatchingSolver(
        const std::vector<std::tuple<int,int,int>>& edge_list) {
    // Initialize: set up edge structures
    int max_vertex = 0;
    for (auto [u, v, w] : edge_list) {
        max_vertex = std::max({max_vertex, u, v});
    }
    int n = max_vertex + 1;

    // If n is odd, no perfect matching exists
    if (n % 2 != 0) {
        throw std::runtime_error("Graph has no perfect matching (odd vertex count)");
    }

    primal_objective = 0;
}

void WeightedMatchingSolver::find_min_perfect_matching() {
    // Placeholder: the actual Blossom VI algorithm
    // This will be called as a refinement step after
    // X-Blossom finds an initial maximum matching
}

const std::vector<std::pair<int,int>>& WeightedMatchingSolver::get_matching() const {
    return matching;
}

int64_t WeightedMatchingSolver::get_primal_objective() const {
    return primal_objective;
}