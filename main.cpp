/// \file main.cpp
/// \brief CLI entry point for hybrid maximum-weight matching.
///
/// Usage:
///   ./hybrid_blossom --mtx <graph.mtx> <num_threads>
///   ./hybrid_blossom <rowOffsets.txt> <columnIndices.txt> <num_threads>

#include "hybrid_blossom.h"
#include <iostream>
#include <string>
#include <vector>
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage:\n"
                  << "  " << argv[0] << " --mtx <graph.mtx> <num_threads>\n"
                  << "  " << argv[0] << " <rowOffsets.txt> <columnIndices.txt> <num_threads>\n";
        return 1;
    }

    std::vector<int> rowOffsets, colIndices, arcWeights;
    int numThreads = 4;
    bool ok = false;

    if (std::string(argv[1]) == "--mtx") {
        if (argc < 4) {
            std::cerr << "Missing num_threads\n";
            return 1;
        }
        auto t0 = std::chrono::steady_clock::now();
        ok = read_mtx(argv[2], rowOffsets, colIndices, arcWeights);
        auto t1 = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
        std::cout << "Graph loaded in " << ms << "ms\n";
        numThreads = std::stoi(argv[3]);
    } else {
        if (argc < 4) {
            std::cerr << "Missing num_threads\n";
            return 1;
        }
        auto t0 = std::chrono::steady_clock::now();
        ok = read_csr_files(argv[1], argv[2], rowOffsets, colIndices, arcWeights);
        auto t1 = std::chrono::steady_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
        std::cout << "Graph loaded in " << ms << "ms\n";
        numThreads = std::stoi(argv[3]);
    }

    if (!ok) return 1;

    int n = (int)rowOffsets.size() - 1;
    int m = (int)colIndices.size() / 2;
    std::cout << "Graph: " << n << " vertices, " << m << " edges, "
              << numThreads << " threads\n";

    // Run hybrid solver
    MatchingResult res = hybrid_blossom_maximum_weight_matching(
        rowOffsets, colIndices, arcWeights, numThreads);

    // Print instrumentation
    std::cout << "\n--- Results ---\n";
    std::cout << "Final cardinality:    " << res.final_cardinality << "\n";
    std::cout << "Final weight:         " << res.weight << "\n";
    std::cout << "Valid matching:       " << (res.valid ? "YES" : "NO") << "\n";
    std::cout << "\n--- Timing ---\n";
    std::cout << "X-Blossom phase:      " << res.time_xblossom_ms  << "ms\n";
    std::cout << "Weighted phase:       " << res.time_weighted_ms  << "ms\n";
    std::cout << "Validation:           " << res.time_validate_ms  << "ms\n";
    std::cout << "Total:                " << res.time_total_ms     << "ms\n";
    std::cout << "\n--- Algorithm stats ---\n";
    std::cout << "Initial cardinality:  " << res.initial_cardinality << "\n";
    std::cout << "Initial weight:       " << res.initial_weight << "\n";
    std::cout << "Augmentations:        " << res.num_augmentations << "\n";
    std::cout << "Blossom contractions: " << res.num_blossom_contractions << "\n";
    std::cout << "Blossom expansions:   " << res.num_blossom_expansions << "\n";
    std::cout << "Dual updates:         " << res.num_dual_updates << "\n";

    if (!res.valid) {
        std::cerr << "\nERROR: matching is INVALID\n";
        return 2;
    }

    return 0;
}
