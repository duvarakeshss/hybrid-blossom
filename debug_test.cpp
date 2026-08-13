#include "hybrid_blossom.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

int main() {
    // Read mtx
    std::ifstream f("weighted_test.mtx");
    int n, m;
    f >> n >> m;
    std::vector<std::vector<int>> adj(n);
    std::vector<int> weights;

    for (int i = 0; i < m; ++i) {
        int u, v, w;
        f >> u >> v >> w;
        adj[u].push_back(v);
        adj[v].push_back(u);
        weights.push_back(w);
    }

    std::cout << "=== ADJ LIST ===\n";
    for (int i = 0; i < n; ++i) {
        std::cout << i << ": ";
        for (int v : adj[i]) std::cout << v << " ";
        std::cout << "\n";
    }

    // Build CSR
    std::vector<int> rowOffsets(n + 1);
    std::vector<int> columnIndices;
    rowOffsets[0] = 0;
    for (int i = 0; i < n; ++i) {
        rowOffsets[i + 1] = rowOffsets[i] + adj[i].size();
        for (int v : adj[i]) {
            columnIndices.push_back(v);
        }
    }

    std::cout << "\n=== CSR ===\n";
    std::cout << "rowOffsets: ";
    for (int r : rowOffsets) std::cout << r << " ";
    std::cout << "\n";

    // Run hybrid
    std::vector<int> M;
    hybrid_blossom_maximum_weight_matching(rowOffsets, columnIndices, weights, M, 4);

    std::cout << "\n=== MATCHING ===\n";
    for (int v = 0; v < n; ++v) {
        std::cout << v << " -> " << M[v] << "\n";
    }
    return 0;
}
