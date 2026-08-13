#include "hybrid_blossom.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <cstdlib>
#include <sstream>

/// \brief Read a file of integers into a vector
void read_ints(const std::string& path, std::vector<int>& out) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }
    int v;
    while (f >> v) out.push_back(v);
}

/// \brief Build edge weights from CSR arrays
///
/// The CSR format: rowOffsets array, columnIndices array.
/// Each edge appears twice (once per direction).
/// We assign weight 1 to each edge.
void build_unweighted(const std::vector<int>& rowOffsets,
                      const std::vector<int>& columnIndices,
                      std::vector<int>& edgeWeights) {
    edgeWeights.assign(columnIndices.size() / 2, 1);
}

/// \brief Build CSR from .mtx file (with edge weights)
///
/// Format: n m
///         u v w
///         ...
void read_mtx(const std::string& path,
              std::vector<int>& rowOffsets,
              std::vector<int>& columnIndices,
              std::vector<int>& edgeWeights) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "Cannot open " << path << "\n";
        std::exit(1);
    }

    int n, m;
    f >> n >> m;

    // Build adjacency list
    std::vector<std::vector<int>> adj(n);
    std::vector<int> weights;

    for (int i = 0; i < m; ++i) {
        int u, v, w;
        f >> u >> v >> w;
        adj[u].push_back(v);
        adj[v].push_back(u);
        weights.push_back(w);
    }

    // Build CSR
    rowOffsets.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        rowOffsets[i + 1] = rowOffsets[i] + static_cast<int>(adj[i].size());
        for (int v : adj[i]) {
            columnIndices.push_back(v);
        }
    }

    edgeWeights = weights;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " <rowOffsets.txt> <columnIndices.txt> <num_threads>\n"
                  << "  " << argv[0] << " --mtx <graph.mtx> <num_threads>\n";
        return 1;
    }

    int num_threads = 4;
    std::vector<int> rowOffsets, columnIndices, edgeWeights;

    if (std::string(argv[1]) == "--mtx") {
        // Read weighted graph from .mtx
        read_mtx(argv[2], rowOffsets, columnIndices, edgeWeights);
        num_threads = (argc > 3) ? std::stoi(argv[3]) : 4;
    } else {
        // Read CSR from rowOffsets + columnIndices files
        read_ints(argv[1], rowOffsets);
        read_ints(argv[2], columnIndices);
        num_threads = (argc > 3) ? std::stoi(argv[3]) : 4;
        build_unweighted(rowOffsets, columnIndices, edgeWeights);
    }

    int n = static_cast<int>(rowOffsets.size()) - 1;
    std::cout << "Graph: " << n << " vertices, "
              << edgeWeights.size() << " weighted edges\n";

    // Run
    std::vector<int> M;
    auto start = std::chrono::steady_clock::now();
    hybrid_blossom_maximum_weight_matching(rowOffsets, columnIndices, edgeWeights, M, num_threads);
    auto end = std::chrono::steady_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Total: " << ms << "ms\n";

    return 0;
}