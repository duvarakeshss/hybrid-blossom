#include "parallel_augment.h"
#include <set>

/// \brief X-Blossom's core parallel augmenting path search
///
/// This replicates parFindAugmentingPathNoRecursionUpdatePathTable from x_blossom.cpp
///
/// Each thread scans its partition of the adjacency list.
/// When it finds an even vertex connected to an even vertex in a different tree,
/// it uses compare_exchange_strong (atomic) to claim both trees.
/// If both are claimed, it builds the augmenting path from the path table.
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
        std::vector<std::vector<int>>& path_collection) {

    int n = static_cast<int>(nodes_vector.size());
    int estimated_size = static_cast<int>(static_cast<double>(n) / num_threads * 1.25);
    std::vector<int> local_path;
    local_path.reserve(estimated_size);

    // Thread-local path construction
    for (int i = index; i < n; i += num_threads) {
        int v = nodes_vector[i];
        int start = rowOffsets[v];
        int end = rowOffsets[v + 1];

        for (int j = start; j < end; ++j) {
            int w = columnIndices[j];
            int tree_v = belongs[v];
            int tree_w = belongs[w];

            if (is_even[w] && tree_v != tree_w && tree_v != -1 && tree_w != -1) {
                // Found two vertices in different trees
                int min_tree = std::min(tree_v, tree_w);
                int max_tree = std::max(tree_v, tree_w);
                int expected = 0;

                // Atomic CAS: try to claim both trees
                if (select_tree[min_tree].compare_exchange_strong(expected, 1)) {
                    expected = 0;
                    if (select_tree[max_tree].compare_exchange_strong(expected, 1)) {
                        // Both trees claimed → build path
                        std::vector<int> pv = find_path_from_table(path_table, v);
                        std::vector<int> pw = find_path_from_table(path_table, w);

                        // Copy paths in X-Blossom order
                        for (int s = static_cast<int>(pv.size()) - 1; s >= 0; --s) {
                            local_path.push_back(pv[s]);
                        }
                        for (int t = 0; t < static_cast<int>(pw.size()); ++t) {
                            local_path.push_back(pw[t]);
                        }

                        // Reset the atomic locks
                        select_tree[max_tree] = 0;
                    } else {
                        // Failed to claim second tree → back out
                        select_tree[min_tree] = 0;
                    }
                }
            }
        }
    }

    // Push collected path to shared collection
    if (!local_path.empty()) {
        static std::mutex path_mutex;
        std::lock_guard<std::mutex> guard(path_mutex);
        path_collection.push_back(local_path);
    }
}

std::vector<int> find_path_from_table(
        const std::vector<std::vector<int>>& path_table,
        int v) {
    std::vector<int> path;
    path.push_back(v);
    if (path_table[v].empty()) return path;

    int cur = v;
    while (!path_table[cur].empty()) {
        path.insert(path.end(), path_table[cur].begin(), path_table[cur].end());
        cur = path.back();
    }
    return path;
}

void parallel_init_atomics(
        std::vector<std::atomic<int>>& select_tree,
        std::vector<std::atomic<int>>& select_match,
        std::vector<std::atomic<int>>& select_blossom,
        std::vector<std::vector<int>>& path_table,
        int nodes,
        int num_threads) {
    int chunk = (nodes + num_threads - 1) / num_threads;
#pragma omp parallel for
    for (int t = 0; t < num_threads; ++t) {
        int start = t * chunk;
        int end = std::min(start + chunk, nodes);
        for (int i = start; i < end; ++i) {
            select_tree[i] = 0;
            select_match[i] = 0;
            select_blossom[i] = 0;
            path_table[i].clear();
        }
    }
}

void parallel_find_exposed(
        std::vector<int>& exposed,
        const std::vector<int>& M,
        int num_threads) {
    exposed.clear();
    int n = static_cast<int>(M.size());
    int chunk = (n + num_threads - 1) / num_threads;
    std::mutex exp_mutex;
    std::vector<std::vector<int>> local(n);

#pragma omp parallel for
    for (int t = 0; t < num_threads; ++t) {
        int start = t * chunk;
        int end = std::min(start + chunk, n);
        for (int i = start; i < end; ++i) {
            if (M[i] == -1) local[t].push_back(i);
        }
    }
    for (auto& l : local)
        exposed.insert(exposed.end(), l.begin(), l.end());
}

void parallel_init_exposed_vector(
        const std::vector<int>& exposed,
        std::vector<int>& is_even,
        std::vector<int>& belongs,
        int num_threads) {
    int n = static_cast<int>(exposed.size());
    int chunk = (n + num_threads - 1) / num_threads;
#pragma omp parallel for
    for (int t = 0; t < num_threads; ++t) {
        int start = t * chunk;
        int end = std::min(start + chunk, n);
        for (int i = start; i < end; ++i) {
            int v = exposed[i];
            is_even[v] = 1;
            belongs[v] = v;
        }
    }
}

void parallel_update_matching(
        std::vector<int>& M,
        const std::vector<std::vector<int>>& path_collection,
        int num_threads) {
    for (const auto& path : path_collection) {
        for (int i = 0; i < static_cast<int>(path.size()) - 1; i += 2) {
            int v = path[i];
            int w = path[i + 1];
            M[v] = w;
            M[w] = v;
        }
    }
}

/// \brief Check if a path-collection-based matching is valid
static bool validate_matching_impl(const std::vector<int>& M, int n) {
    std::set<std::pair<int,int>> edges;
    for (int v = 0; v < n; ++v) {
        if (M[v] != -1) {
            int w = M[v];
            if (w < 0 || w >= n) return false;
            if (M[w] != v) return false;
            if (v < w) edges.insert({v, w});
        }
    }
    return true;
}