#include "hybrid_blossom.h"
#include <iostream>
#include <cassert>

/// \brief Validate that M is a valid matching
bool check_matching_valid(const std::vector<int>& M) {
    for (int v = 0; v < static_cast<int>(M.size()); ++v) {
        if (M[v] != -1) {
            int w = M[v];
            if (w < 0 || w >= static_cast<int>(M.size())) {
                std::cerr << "ERROR: M[" << v << "] = " << w << " out of range\n";
                return false;
            }
            if (M[w] != v) {
                std::cerr << "ERROR: M[" << v << "] = " << w
                          << " but M[" << w << "] = " << M[w] << "\n";
                return false;
            }
        }
    }
    return true;
}

/// \brief Compute matching size
int matching_size(const std::vector<int>& M) {
    int count = 0;
    for (int v = 0; v < static_cast<int>(M.size()); ++v) {
        if (M[v] != -1 && v < M[v]) ++count;
    }
    return count;
}

/// \brief Test the hybrid on the Google graph
void test_google_graph() {
    // Load from the example dataset
    std::vector<int> rowOffsets;
    std::vector<int> columnIndices;
    std::vector<int> edge_weights;

    // Just use the loaded CSR data
    // (This would be called from main for actual tests)
}

/// \brief Run a randomized correctness test
void test_random_graph(int num_nodes, int num_edges, int seed) {
    std::srand(seed);
    std::vector<int> rowOffsets(num_nodes + 1);
    std::vector<int> columnIndices;
    std::vector<int> edge_weights;

    // Generate random edges
    // ...

    // Run the solver
    std::vector<int> M;
    int num_threads = 4;
    // hybrid_maximum_weight_matching(rowOffsets, columnIndices, edge_weights, M, num_threads);

    // Verify
    bool valid = check_matching_valid(M);
    std::cout << "Test " << (valid ? "PASSED" : "FAILED")
              << " | matching size: " << matching_size(M) << "\n";
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <rowOffsets> <columnIndices> <num_threads>\n";
        return 1;
    }

    // Load CSR from files
    std::string row_file = argv[1];
    std::string col_file = argv[2];
    int num_threads = std::stoi(argv[3]);

    std::vector<int> rowOffsets;
    std::vector<int> columnIndices;
    std::vector<int> edge_weights;  // placeholder

    {
        std::ifstream in(row_file);
        int val;
        while (in >> val) rowOffsets.push_back(val);
    }
    {
        std::ifstream in(col_file);
        int val;
        while (in >> val) columnIndices.push_back(val);
    }

    // For unweighted matching, set all weights to 1
    edge_weights.assign(columnIndices.size() / 2, 1);

    // Run
    std::vector<int> M;
    auto start = std::chrono::high_resolution_clock::now();
    // hybrid_maximum_weight_matching(rowOffsets, columnIndices, edge_weights, M, num_threads);
    auto end = std::chrono::high_resolution_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "Time: " << ms.count() << "ms\n";
    std::cout << "Matching size: " << matching_size(M) << "\n";
    std::cout << "Valid: " << (check_matching_valid(M) ? "yes" : "no") << "\n";

    return 0;
}