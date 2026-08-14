#ifndef HYBRID_BLOSSOM_H
#define HYBRID_BLOSSOM_H

#include <vector>
#include <cstdint>
#include <string>

// ============================================================
//  Public result type
// ============================================================

/// \brief Result of hybrid maximum-weight matching
struct MatchingResult {
    std::vector<int> mate;   ///< mate[v] = matched vertex, or -1
    int64_t          weight; ///< total weight of matched edges
    bool             valid;  ///< mate is a valid matching
    bool             optimal;///< weight equals brute-force optimum (set by tests)

    // --- Instrumentation fields ---
    long long time_load_ms       = 0;
    long long time_xblossom_ms   = 0;
    long long time_weighted_ms   = 0;
    long long time_validate_ms   = 0;
    long long time_total_ms      = 0;

    int  num_vertices            = 0;
    int  num_edges               = 0;
    int  num_threads             = 0;
    int  initial_cardinality     = 0;
    int  final_cardinality       = 0;
    int64_t initial_weight       = 0;

    int  num_augmentations       = 0;
    int  num_blossom_contractions= 0;
    int  num_blossom_expansions  = 0;
    int  num_dual_updates        = 0;
};

// ============================================================
//  Main API
// ============================================================

/// \brief Compute maximum-weight matching on an undirected graph.
///
/// Uses a parallel X-Blossom warm-start followed by an exact
/// primal-dual weighted matching refinement.
///
/// Works for:
///   - general (non-bipartite) graphs
///   - bipartite graphs
///   - positive, zero, and negative edge weights
///   - graphs without perfect matchings
///   - odd numbers of vertices
///   - isolated vertices
///
/// \param rowOffsets    CSR row offsets (size n+1)
/// \param colIndices    CSR column indices (each undirected edge
///                      appears twice: u->v and v->u)
/// \param arcWeights    Weight for each arc: arcWeights[j] = weight
///                      of arc j (same weight for both directions
///                      of the same undirected edge)
/// \param numThreads    Parallel threads for X-Blossom phase (>=1)
/// \return              MatchingResult with mate[], weight, valid, instrumentation
MatchingResult hybrid_blossom_maximum_weight_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& colIndices,
    const std::vector<int>& arcWeights,
    int numThreads);

/// \brief Legacy void API — fills M in-place.
///        Prints a summary line. Calls the result-returning overload.
void hybrid_blossom_maximum_weight_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& colIndices,
    const std::vector<int>& arcWeights,
    std::vector<int>& M,
    int numThreads);

// ============================================================
//  Graph-building helpers (used by main.cpp and tests)
// ============================================================

/// \brief Build CSR + per-arc weights from a .mtx edge list.
///
/// .mtx format:  first line "n m", then m lines "u v w" (0-based).
/// Self-loops are silently ignored.
/// Each undirected edge (u,v,w) produces arcs u->v and v->u,
/// both with weight w.
///
/// Returns false on I/O or format error.
bool read_mtx(const std::string& path,
              std::vector<int>& rowOffsets,
              std::vector<int>& colIndices,
              std::vector<int>& arcWeights);

/// \brief Build CSR from plain rowOffsets / columnIndices text files.
///        Arc weights are set to 1 (unweighted).
bool read_csr_files(const std::string& rowFile,
                    const std::string& colFile,
                    std::vector<int>& rowOffsets,
                    std::vector<int>& colIndices,
                    std::vector<int>& arcWeights);

// ============================================================
//  Matching validator (standalone, usable from tests)
// ============================================================

/// \brief Full matching validator.
///
/// Checks:
///   1. Symmetry:  M[v]==-1  OR  M[M[v]]==v
///   2. Valid arc: (v, M[v]) must exist in the CSR graph
///   3. No self-match:  M[v] != v
///   4. Edge uniqueness
///   5. Computes total weight (int64_t, each edge counted once)
///
/// \param rowOffsets   CSR row offsets
/// \param colIndices   CSR column indices
/// \param arcWeights   CSR arc weights
/// \param M            Matching vector
/// \param totalWeight  [out] total matched weight
/// \return true if matching is valid
bool validate_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& colIndices,
    const std::vector<int>& arcWeights,
    const std::vector<int>& M,
    int64_t& totalWeight);

#endif // HYBRID_BLOSSOM_H
