/// \file weighted_matching.cpp
/// \brief Blossom VI-style weighted matching solver stub.
///
/// The actual weighted matching is now implemented in hybrid_blossom.cpp
/// (mwm::MWMSolver).  This file keeps the WeightedMatchingSolver class
/// compilable for any external code that still references it.

#include "weighted_matching.h"
#include <stdexcept>
#include <algorithm>

WeightedMatchingSolver::WeightedMatchingSolver(
        const std::vector<std::tuple<int,int,int>>& edge_list)
    : primal_objective(0)
{
    int max_vertex = 0;
    for (const auto& e : edge_list)
        max_vertex = std::max({max_vertex,
                               std::get<0>(e),
                               std::get<1>(e)});
    int n = max_vertex + 1;
    (void)n; // unused in stub
}

void WeightedMatchingSolver::find_min_perfect_matching() {
    // Stub: exact solver is in mwm::MWMSolver (hybrid_blossom.cpp)
}

const std::vector<std::pair<int,int>>& WeightedMatchingSolver::get_matching() const {
    return matching;
}

int64_t WeightedMatchingSolver::get_primal_objective() const {
    return primal_objective;
}
