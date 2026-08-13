#ifndef WEIGHTED_MATCHING_H
#define WEIGHTED_MATCHING_H

#include <vector>
#include <cstdint>
#include <utility>
#include <tuple>

/// \brief Blossom VI-style weighted minimum perfect matching solver
///
/// This is the "Blossom VI" side of the hybrid:
///   - maintains dual variables (amortized on edges)
///   - uses edge heaps for O(E log V) complexity
///   - supports negative weights
///
/// Key difference from X-Blossom: this works with weights.
/// X-Blossom finds ANY maximum matching (unweighted).
/// This finds the MINIMUM-WEIGHT perfect matching.
struct WeightedMatchingSolver {
    int64_t primal_objective;
    int64_t dual_objective;

    std::vector<std::pair<int,int>> matching;

    /// \brief Construct from a list of (head, tail, weight) tuples
    explicit WeightedMatchingSolver(
            const std::vector<std::tuple<int,int,int>>& edge_list);

    /// \brief Find the minimum-weight perfect matching
    void find_min_perfect_matching();

    /// \brief Get the matching result
    const std::vector<std::pair<int,int>>& get_matching() const;

    /// \brief Get the primal objective (total weight)
    int64_t get_primal_objective() const;
};

#endif // WEIGHTED_MATCHING_H