#ifndef PARALLEL_AUGMENT_H
#define PARALLEL_AUGMENT_H

#include <vector>
#include <atomic>
#include <thread>
#include <mutex>

/// \brief Parallel augmenting path finder, using X-Blossom's technique
///
/// This implements the core X-Blossom parallelism:
///   - Multiple threads divide the adjacency list
///   - Each thread uses std::atomic compare_exchange to claim trees
///   - Path tables are flat std::vector<int> arrays (no recursion)
///
/// \param rowOffsets     CSR row offsets
/// \param columnIndices  CSR column indices
/// \param nodes_vector   Vector of nodes to search
/// \param index          Thread partition index
/// \param num_threads    Total threads
/// \param is_even        Even/odd parity per vertex
/// \param belongs        Tree membership per vertex
/// \param path_table     Path table (array of flat vectors)
/// \param select_tree    Atomic flags for tree selection
/// \param path_collection Output: collected augmenting paths
void parallel_find_augmenting_paths(
        const std::vector<int>& rowOffsets,
        const std::vector<int>& columnIndices,
        const std::vector<int>& nodes_vector,
        int index,
        int num_threads,
        std::vector<int>& is_even,
        std::vector<int>& belongs,
        std::vector<std::vector<int>>& path_table,
        std::vector<std::atomic<int>>& select_tree,
        std::vector<std::vector<int>>& path_collection);

/// \brief Initialize atomic flags in parallel
void parallel_init_atomics(
        std::vector<std::atomic<int>>& select_tree,
        std::vector<std::atomic<int>>& select_match,
        std::vector<std::atomic<int>>& select_blossom,
        std::vector<std::vector<int>>& path_table,
        int nodes,
        int num_threads);

/// \brief Find exposed (unmatched) nodes in parallel
void parallel_find_exposed(
        std::vector<int>& exposed,
        const std::vector<int>& M,
        int num_threads);

/// \brief Initialize is_even and belongs for exposed nodes
void parallel_init_exposed_vector(
        const std::vector<int>& exposed,
        std::vector<int>& is_even,
        std::vector<int>& belongs,
        int num_threads);

/// \brief Apply matching update: flip matched/unmatched along a path
void parallel_update_matching(
        std::vector<int>& M,
        const std::vector<std::vector<int>>& path_collection,
        int num_threads);

/// \brief Find the full path from a path table
std::vector<int> find_path_from_table(
        const std::vector<std::vector<int>>& path_table,
        int v);

#endif // PARALLEL_AUGMENT_H