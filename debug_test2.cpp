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

    // Build CSR from mtx
    std::vector<std::vector<int>> adj(n);
    std::vector<int> weights;
    std::vector<int> rowOffsets(n + 1);
    std::vector<int> columnIndices;

    for (int i = 0; i < m; ++i) {
        int u, v, w;
        f >> u >> v >> w;
        adj[u].push_back(v);
        adj[v].push_back(u);
        weights.push_back(w);
    }

    rowOffsets[0] = 0;
    for (int i = 0; i < n; ++i) {
        rowOffsets[i + 1] = rowOffsets[i] + adj[i].size();
        for (int v : adj[i]) {
            columnIndices.push_back(v);
        }
    }

    // Run hybrid
    std::vector<int> M;
    hybrid_blossom_maximum_weight_matching(rowOffsets, columnIndices, weights, M, 4);

    // Print
    for (int v = 0; v < n; ++v) {
        std::cout << "M[" << v << "] = " << M[v] << "\n";
    }

    bool ok = true;
    for (int v = 0; v < n; ++v) {
        if (M[v] != -1 && (M[M[v]] != v)) {
            std::cerr << "Invalid: " << v << " -> " << M[v] << " but M[" << M[v] << "] = " << M[M[v]] << "\n";
            ok = false;
        }
    }
    std::cout << "Valid: " << (ok ? "yes" : "no") << "\n";
    return 0;
}
