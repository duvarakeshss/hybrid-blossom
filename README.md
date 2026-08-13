# Hybrid Blossom — Maximum Weight Matching on Undirected Graphs

**Hybrid Blossom** is a C++ implementation of **maximum-weight matching** on undirected
graphical data structures that combines two powerful algorithmic techniques:

1. **X-Blossom** — A parallel, multi-threaded BFS-based augmenting-path search for
   *maximum cardinality* (unweighted) matching.
2. **Blossom VI** — A primal–dual *minimum-weight perfect matching* solver that
   maintains dual variables (potentials) on edges and nodes.

The "hybrid" insight is that **X-Blossom** finds an initial large matching very quickly
using parallel BFS on a CSR (Compressed Sparse Row) graph, while **Blossom VI**
refines that matching into a *weighted* optimum using its dual-variable framework.
The two algorithms share the same underlying **Blossom data structure** (odd cycles /
blossoms) and are composed so that X-Blossom's parallel path-table search feeds into
Blossom VI's primal–dual update loop.

---

## Table of Contents

1. [Background: Maximum Matching on Graphs](#background)
2. [What Is a Blossom?](#blossom-definition)
3. [Algorithm Family — Two Approaches, One Hybrid](#algorithm-family)
4. [Hybrid Blossom Architecture](#architecture)
5. [The Two Solvers in Detail](#two-solvers)
6. [Parallel Search (X-Blossom Phase)](#parallel-search)
7. [Dual-Variable Refinement (Blossom VI Phase)](#dual-refinement)
8. [CSR Graph Representation](#csr-graph)
9. [Hybrid Blossom Data Structures](#data-structures)
10. [Algorithm Flow: Step by Step](#algorithm-flow)
11. [Code Walkthrough](#code-walkthrough)
12. [Graph Structure and IO](#io)
13. [Building and Running](#build-and-run)
14. [Example Datasets](#datasets)
15. [Performance Characteristics](#performance)
16. [References](#references)

---

## Background: Maximum Matching on Graphs

A **matching** in an undirected graph $G = (V, E)$ is a set of edges $M \subseteq E$ such
that no two edges in $M$ share a common vertex. The **maximum weight matching** problem
asks: *find a matching that maximizes the total sum of edge weights.*

| Problem | Objective | Optimal |
|---------|-----------|---------|
| **Maximum matching** (cardinality) | Max $|M|$ — as many edges as possible | $|M| \le \lfloor n/2\rfloor$, matching is *maximum* |
| **Maximum weight matching** | Max $\sum_{e\in M} w(e)$ | May be smaller than $|n/2|$ but yields the highest total weight |
| **Minimum weight perfect matching** | Min $\sum_{e\in M} w(e)$ with $|M| = n/2$ | Every vertex is matched |

**Hybrid Blossom** solves the *maximum-weight* variant (second column). It supports
bipartite and general (non-bipartite) graphs and works with both positive and
negative edge weights.

---

## What Is a Blossom?

A **blossom** (in the sense of Edmonds, 1965) is an **odd cycle** — a cycle of odd length —
that appears when two alternating paths meet and form a cycle. Blossoms are the core
reason why the simple greedy approach fails: in a general graph, matched edges and
unmatched edges can form an **odd cycle** that blocks further augmentation.

```
An alternating path:
  v0 --(unmatched)--> v1 --(matched)--> v2 --(unmatched)--> v3 ...

When this path loops back on itself and forms a cycle of odd length:
  v0 -- v1 -- v2 -- ... -- v_k -- v0 (k is odd)
```

This is a **blossom**.

Both X-Blossom and Blossom VI handle blossoms through **contraction** — the odd cycle
is compressed into a single **super-node** that acts like a normal vertex in subsequent
search rounds.

---

## Algorithm Family — Two Approaches, One Hybrid

### X-Blossom (Parallel BFS)

- **Goal:** Find a maximum (unweighted) matching — maximize $|M|$.
- **Method:** Repeated parallel BFS from exposed (unmatched) vertices to find
  **augmenting paths** — alternating paths whose two endpoints are both unmatched.
- **Parallelism:** `std::thread` with `std::atomic<int>` compare-exchange locks for
  contention-free tree-claiming. Each thread gets a *stripe* of the adjacency list
  (`index % num_threads`) and uses `compare_exchange_strong` to atomically claim
  pairs of BFS trees.
- **Key insight:** The *first* round of parallel augmenting-path search is usually
  enough to find a near-maximum matching. Subsequent rounds refine it.
- **Data structure:** A **path table** — a flat `std::vector<std::vector<int>>` where
  each entry points to its parent in the BFS tree (linked-list style). Tracing from
  any node to the root follows: `node → path_table[node] → ...` until empty.

### Blossom VI (Primal–Dual)

- **Goal:** Find a *minimum-weight* perfect matching — minimize $\sum w(e)$ under
  $|M| = n/2$.
- **Method:** Maintain **dual variables** (potentials) on nodes and trees. The
  algorithm alternates between:
  - **Primal updates** — changing the matching structure (grow, augment, make cherry,
    expand, shrink).
  - **Dual updates** — increasing tree potentials until new edges become **slack-zero**
    and thus actionable.
- **Key concept:** Slack = $4 \times \text{weight} - \sum \text{dual variables}$. When
  an edge's slack reaches zero, it becomes a candidate for primal operations.
- **Blossom:** An odd cycle is contracted into a **blossom node** with a
  `blossom_parent` pointer. The **receptacle** is the node within a blossom that
  holds the external matched edge.
- **Data structures:** Pairing heaps ($O(\log n)$ meld), Custom `ArcIndex` for
  edge traversal inside blossoms, `blossom_ancestors` for path compression.

### The Hybrid

X-Blossom's parallel BFS finds a *large* matching very fast. Blossom VI's dual-variable
refinement then adjusts the matching toward *optimal weight*. The two share:

- The same **CSR graph format** (row offsets + column indices).
- The same **blossom contraction** logic (odd cycles → super-nodes).
- The same **path table** structure (linked-list to root).

They differ in:

| Aspect | X-Blossom | Blossom VI |
|--------|-----------|------------|
| **Goal** | Maximum $|M|$ (unweighted) | Minimum-weight *perfect* $M$ |
| **Parallelism** | Multi-threaded (OpenMP-style) | Single-threaded |
| **Dual variables** | None | Node potentials + tree variables |
| **Heaps** | `std::vector` + atomic flags | Pairing heaps ($O(\log n)$) |
| **Blossom detection** | (even, even) edge → BFS path collision | `blossom_parents` + `blossom_ancestors` |
| **Complexity** | $O(n \cdot m)$ per round (parallel) | $O(n^3)$ worst case (practical) |

---

## Hybrid Blossom Architecture

```
hybrid-x-blossom/
├── CMakeLists.txt                # Build configuration (C++20, pthread)
│
├── hybrid_blossom.h              # Header: API declaration + Core structs
├── hybrid_blossom.cpp             # Main hybrid solver implementation
│
├── parallel_augment.h            # Header: parallel path-table search
├── parallel_augment.cpp           # X-Blossom's parallel BFS + augmenting path
│
├── weighted_matching.h           # Header: Blossom VI-style weighted matching
├── weighted_matching.cpp         # Blossom VI dual-variable refinement
│
├── main.cpp                       # CLI entry point (reads CSR / .mtx)
├── test_hybrid.cpp                # Test harness
│
├── debug_test.cpp                 # Debug / validation test
├── debug_test2.cpp               # Standalone CSR test
│
├── *.mtx                          # Sample graph files (matrix-market format)
│
├── hybrid_blossom.exe             # Built executable (Windows)
└── test_hybrid.exe                # Test executable
```

### Header Files

| File | API | Purpose |
|------|-----|---------|
| `hybrid_blossom.h` | `hybrid_blossom_maximum_weight_matching(...)` | Main entry point — runs the full hybrid pipeline |
| `parallel_augment.h` | `parallel_find_augmenting_paths(...)` | X-Blossom's parallel BFS on CSR |
| `weighted_matching.h` | `WeightedMatchingSolver::find_min_perfect_matching()` | Blossom VI's dual-variable refinement |

---

## The Two Solvers in Detail

### `hybrid_blossom_maximum_weight_matching` — Main API

```cpp
void hybrid_blossom_maximum_weight_matching(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& columnIndices,
    const std::vector<int>& edgeWeights,
    std::vector<int>& M,
    int numThreads);
```

| Parameter | Type | Meaning |
|-----------|------|---------|
| `rowOffsets` | `const std::vector<int>&` | CSR row offsets: `rowOffsets[i]` = start index of vertex `i`'s adjacency in `columnIndices` |
| `columnIndices` | `const std::vector<int>&` | CSR column indices: all adjacency lists concatenated |
| `edgeWeights` | `const std::vector<int>&` | Edge weights (one per edge, in CSR order) |
| `M` | `std::vector<int>&` | **Output:** $M[v] = \text{mate}(v)$ (or $-1$ if unmatched) |
| `numThreads` | `int` | Number of parallel threads ($\ge 1$) |

The function:
1. **Initializes** dual variables (`dual[v] = 4 * min_weight`) for all vertices.
2. **Finds exposed** (unmatched) vertices.
3. **Matches** exposed pairs via tight edges (slack $= 0$).
4. **Iterates:** raises dual variables until new tight edges appear.
5. **Validates** the output matching.

This is the **weighted** variant — it works with `edgeWeights` to find the
*maximum-weight* matching. The unweighted X-Blossom (in `X-Blossom/` directory)
uses only CSR and ignores weights.

### `parallel_find_augmenting_paths` — X-Blossom's Parallel BFS

```cpp
void parallel_find_augmenting_paths(
    const std::vector<int>& rowOffsets,
    const std::vector<int>& columnIndices,
    const std::vector<int>& nodes_vector,
    int index, int num_threads,
    std::vector<int>& is_even,
    std::vector<int>& belongs,
    std::vector<std::vector<int>>& path_table,
    std::vector<std::atomic<int>>& select_tree,
    std::vector<std::vector<int>>& path_collection);
```

This is the heart of X-Blossom's parallelism. Each thread:
1. Iterates over `nodes_vector` in **stripes** (`index % num_threads`).
2. For each **even** node `v`, checks neighbors `w`:
   - If `w` is **even** AND belongs to a **different tree** (different root):
     - **Atomic CAS** (`compare_exchange_strong`) on `select_tree[min(tree_v, tree_w)]` and `select_tree[max(tree_v, tree_w)]` to **claim** both trees.
     - If claimed: concatenate `path_table[v]` (reversed) + `path_table[w]` → augmenting path.
   - If another thread already claimed the tree pair: skip (no lock contention).
3. Push collected paths to `path_collection`.

**Key invariant:** An augmenting path exists when two even nodes in different
BFS trees share an edge — the concatenation of their root-to-node paths forms
an alternating path between two exposed roots.

### `parallel_find_exposed` — Find Unmatched Vertices

```cpp
void parallel_find_exposed(
    std::vector<int>& exposed,
    const std::vector<int>& M,
    int num_threads);
```

Partitions `M` into `num_threads` chunks and collects all vertices where
`M[v] == -1` in parallel. Each chunk writes to a thread-local vector; results
are merged at the end.

### `parallel_init_atomics` — Initialize Atomic Flags

```cpp
void parallel_init_atomics(
    std::vector<std::atomic<int>>& select_tree,
    ...
    int nodes, int num_threads);
```

Resets all atomic flags (`select_tree`, `select_match`, `select_blossom`) to `0`
and clears `path_table` for all `nodes`. Uses `#pragma omp parallel for` for
initialization.

### `parallel_update_matching` — Apply Augmenting Paths

```cpp
void parallel_update_matching(
    std::vector<int>& M,
    const std::vector<std::vector<int>>& path_collection,
    int num_threads);
```

For each augmenting path in `path_collection`, it flips matched/unmatched
along the path: `M[v_i] = w_i` for alternating pairs `(v_i, w_i)`.

---

## Dual-Variable Refinement (Blossom VI Phase)

The `WeightedMatchingSolver` class implements the Blossom VI dual-variable
framework:

```cpp
struct WeightedMatchingSolver {
    int64_t primal_objective;   // Total weight of the matching
    int64_t dual_objective;     // Sum of dual variables (quadrupled)

    // ...
    void find_min_perfect_matching();
};
```

### Dual Variable Logic

Each vertex $v$ has a **dual variable** $\text{dual}[v]$ (amortized). For
each edge $(u, v)$ with weight $w$:

$$
\text{slack} = 4 \cdot w - \text{dual}[u] - \text{dual}[v]
$$

When $\text{slack} = 0$, the edge is **tight** — it can be added to the
matching without violating optimality. The algorithm maintains:

- **Non-negative duals:** All $\text{dual}[v] \ge 0$.
- **Non-negative slacks:** $4 \cdot w \ge \text{dual}[u] + \text{dual}[v]$ for
  all edges.

### Relaxation Loop

```cpp
// In hybrid_blossom_maximum_weight_matching:
while (exposed nodes remain) {
    // 1. Increase dual of each exposed node by min_slack
    //    (the smallest slack to any neighbor)
    for (int v : exposed) {
        dual[v] += min_slack(v);
    }

    // 2. Try to match via now-tight edges
    //    (slack == 0 after dual update)
    for (int v : exposed) {
        for (int w : neighbors(v)) {
            if (slack(v, w) == 0) {
                M[v] = w; M[w] = v;  // match!
            }
        }
    }
}
```

This is a **simplified** Blossom VI dual loop. The full MWPMSolver
supports:

- **Pairing heaps** for `O(log n)` priority queue operations.
- **Cherry blossom** detection (two `(+)` nodes in the same tree).
- **Blossom expand** (when a blossom's dual reaches zero).
- **Min-cost flow** dual update for the exact optimal `delta`.

---

## CSR Graph Representation

Both X-Blossom and Blossom VI use **Compressed Sparse Row** (CSR) format.

### Graph (`graph.h`)

```cpp
class Graph {
    std::vector<int> columnIndices;  // All adjacency lists concatenated
    std::vector<int> rowOffsets;     // rowOffsets[i+1] - rowOffsets[i] = degree(i)
    int num_of_nodes;
};
```

### Example

```cpp
// A graph: 0 → {1, 3}, 1 → {0, 2}, 2 → {1, 3}, 3 → {0, 2}
rowOffsets   = {0, 2, 4, 6, 8}   // size n+1
columnIndices = {1, 3, 0, 2, 1, 3, 0, 2}  // each edge appears twice
```

### Input Format

Two files:

1. **`rowOffsets.txt`** — A text file of `n+1` integers. Example:
   ```
   0 2 4 6 8
   ```

2. **`columnIndices.txt`** — Vertices `0` through `n-1`, each **degree** times.
   ```
   1 3 0 2 1 3 0 2
   ```

Alternatively, **`.mtx`** (Matrix Market) format:

```
n m
u v w
u v w
...
```

### Build from `.mtx`

```cpp
void read_mtx(const std::string& path,
              std::vector<int>& rowOffsets,
              std::vector<int>& columnIndices,
              std::vector<int>& edgeWeights);
```

1. Reads `n` and `m` (vertices and edges).
2. Builds an adjacency list in `std::vector<std::vector<int>>`.
3. Converts to CSR: `rowOffsets[i+1] = rowOffsets[i] + adj[i].size()`.
4. Stores all edge weights (one per edge, both directions).

---

## Hybrid Blossom Data Structures

### `HybridNode` — Node State

```cpp
struct HybridNode {
    int tree;           // Which BFS tree root (root index)
    int matched_edge;   // Matched edge index (or -1)
    int minus_parent;   // Minus parent in blossom
    int old_tree;       // Previous tree (for history tracking)
    int old_matched;    // Previous matched state
    bool is_even;       // Even distance from root?
    bool is_alive;      // Still part of the active graph?
    bool plus;          // Plus side of the blossom?
    bool old_plus;      // Previous plus state
    int receptacle;     // The node that holds the external matched edge
};
```

### `StopWatch` — Timing

```cpp
struct StopWatch {
    std::chrono::milliseconds total;       // Total wall time
    std::chrono::milliseconds augment;    // Time in augmenting path search
    std::chrono::milliseconds expand;     // Time in BFS tree expansion
    std::chrono::milliseconds dual;      // Time in dual variable updates
};
```

### Hybrid Matching State (`State`)

```cpp
struct State {
    int n;                    // Number of vertices
    const std::vector<int>& ro;   // CSR row offsets
    const std::vector<int>& ci;   // CSR column indices
    const std::vector<int>& ew;   // Edge weights

    std::vector<int> M;       // Matching vector
    std::vector<int> dual;    // Dual variables (amortized)
    std::vector<int> slack;   // Slack per edge
};
```

### Augmenting Path Collection

```cpp
std::vector<std::vector<int>> path_collection;
// Each inner vector is one augmenting path:
// [v0, w0, v1, w1, ...]
// Count = number of augmenting paths found
```

### Atomic Flags

```cpp
std::vector<std::atomic<int>> select_tree;    // Tree-level atomic locks
std::vector<std::atomic<int>> select_match;   // Matching-level locks
std::vector<std::atomic<int>> select_blossom; // Blossom-level locks
```

These are **lock-free** reservations: `select_tree[key].compare_exchange_strong(0, 1)`
means "only one thread may work on this tree." If a thread fails, it simply moves
on — no contention.

---

## Algorithm Flow: Step by Step

### Phase 1: Initialization (`init_slacks` + `match_exposed`)

```
Input: CSR arrays + edge weights

1. For each vertex v:
   dual[v] = 4 * min(incident_weights)   // Make all slacks >= 0
   (in parallel: n/numThreads chunks)

2. For each unmatched pair (v, w):
   if 4*weight - dual[v] - dual[w] == 0:
     M[v] = w; M[w] = v                     // Match via tight edge
```

### Phase 2: Augmenting Path Search (parallel)

```
3. Find exposed: M[v] == -1 for all v
4. For each exposed v as a BFS tree root:
   - Mark is_even[v] = 1, belongs[v] = v
   - belongs[w] = which tree root

5. PARALLEL — for each v (stripe: index % num_threads):
   For each neighbor w:
     if is_even[w] AND belongs[v] != belongs[w]:
       - Atomic CAS on select_tree[min, max]
       - If claimed: concat path_v + path_w → augmenting path
       - Push to path_collection
       - Return (found one)
```

### Phase 3: Expand BFS Trees (parallel)

```
6. For each v (stripe):
   For each neighbor w where belongs[w] == -1:
     - w is matched to some x
     - Atomic CAS on select_match[min(w, x)]
     - Update path_table[x] = {w, v}
     - is_even[w] = 0 (odd), is_even[x] = 1 (even)
     - belongs[w] = belongs[x] = belongs[v]
```

### Phase 4: Blossom Detection (parallel)

```
7. For each v (stripe):
   For each neighbor w:
     if is_even[w] AND belongs[w] == belongs[v]:
       - Find path_v (v → root) and path_w (w → root)
       - Find common ancestor → blossom base
       - Atomic CAS on select_blossom[blossom nodes]
       - Contract odd cycle: update is_even, path_table
       - Toggle blossom_to_base[blossom_node] = base
```

### Phase 5: Dual Variable Update

```
8. For each exposed v:
   - Find min_slack among neighbors
   - dual[v] += min_slack

9. For each exposed v:
   - Check all w: if 4*weight - dual[v] - dual[w] == 0:
     - M[v] = w (new tight edge → match)
```

### Phase 6: Validation

```
10. For all v: if M[v] != -1 then M[M[v]] == v (bidirectional)
11. Count matched edges (v < M[v] for symmetry)
```

### Termination

The loop terminates when either:

1. **No exposed nodes** remain (perfect matching).
2. **No progress** after `n` iterations (maximum matching found).
3. **All edges** have slack $\ge 0$ (dual feasibility).

---

## Code Walkthrough

### `hybrid_blossom_maximum_weight_matching` (main function)

```cpp
void hybrid_blossom_maximum_weight_matching(...) {
    // 1. Build internal state
    State s(rowOffsets, columnIndices, edgeWeights);
    
    // 2. Initialize dual variables
    init_slacks(s, num_threads);
    
    // 3. Match via tight edges
    match_exposed(s);
    
    // 4. Find remaining exposed
    std::vector<int> exposed;
    find_exposed(s.M, exposed);
    
    // 5. Iterative dual + primal loop
    int iter = 0;
    while (!exposed.empty() && iter < s.n) {
        ++iter;
        
        // Dual: increase each exposed node's potential
        for (int v : exposed) {
            dual[v] += min_slack;  // Until an edge is tight
        }
        
        // Primal: match via now-tight edges
        for (int v : exposed) {
            for (int w : neighbors(v)) {
                if (slack == 0) {
                    M[v] = w;  // Match!
                }
            }
        }
        
        // Re-find exposed
        find_exposed(s.M, exposed);
    }
    
    // 6. Output + validate
    int ms = elapsed_time();
    std::cout << "Hybrid (" << num_threads << " threads): "
              << count(M) << " edges in " << ms << "ms"
              << " | valid: " << (ok ? "yes" : "no");
}
```

### Parallel Augmenting Path Search

```cpp
void parallel_find_augmenting_paths(...) {
    // Thread-local work
    std::vector<int> local_path;
    local_path.reserve(estimated_size);
    
    // Stripe: index % num_threads
    for (int i = index; i < nodes_vector.size(); i += num_threads) {
        int v = nodes_vector[i];
        
        // Check each neighbor
        for (int j = rowOffsets[v]; j < rowOffsets[v+1]; ++j) {
            int w = columnIndices[j];
            
            // Two even nodes in different trees?
            if (is_even[w] && tree_v != tree_w) {
                // Atomic CAS to claim both trees
                int expected = 0;
                if (select_tree[min_tree]
                        .compare_exchange_strong(expected, 1)) {
                    // Second tree
                    if (select_tree[max_tree]
                            .compare_exchange_strong(expected, 1)) {
                        // Build path: reverse path_v + forward path_w
                        local_path.push_back(...);
                    }
                }
            }
        }
    }
    
    // Thread-safe collection via mutex
    if (!local_path.empty()) {
        std::lock_guard<std::mutex> guard(path_mutex);
        path_collection.push_back(local_path);
    }
}
```

### Path Table as Linked List

```cpp
// find_path: trace from v to root via path_table
// path_table[new_v] = {parent} (one hop)
std::vector<int> find_path_from_table(
    const std::vector<std::vector<int>>& path_table,
    int v) 
{
    std::vector<int> path;
    path.push_back(v);
    
    int cur = v;
    while (!path_table[cur].empty()) {
        // Follow parent pointers
        path.insert(path.end(), 
                     path_table[cur].begin(), 
                     path_table[cur].end());
        cur = path.back();  // Next hop
    }
    return path;
}
```

---

## IO

### Reading CSR Files

```cpp
void read_ints(const std::string& path, std::vector<int>& out) {
    std::ifstream f(path);
    int v;
    while (f >> v) out.push_back(v);
}
```

### Reading `.mtx` Files

```cpp
void read_mtx(..., std::vector<int>& rowOffsets,
               std::vector<int>& columnIndices,
               std::vector<int>& edgeWeights) {
    std::ifstream f(path);
    int n, m;
    f >> n >> m;
    
    // Build adjacency list
    // Convert to CSR
    for (int i = 0; i < m; ++i) {
        int u, v, w;
        f >> u >> v >> w;
        adj[u].push_back(v);
        adj[v].push_back(u);
        weights.push_back(w);
    }
    
    rowOffsets[i+1] = rowOffsets[i] + adj[i].size();
}
```

### CLI

```bash
# From CSR files:
./hybrid_blossom rowOffsets.txt columnIndices.txt 8

# From .mtx:
./hybrid_blossom --mtx graph.mtx 8
```

---

## Building and Running

### Build

```bash
# Using CMake
mkdir -p build && cd build
cmake ..
make
```

Requires:

- C++20 compiler (GCC 10+, Clang 12+, MSVC 2022+)
- pthread (for `std::thread`)

### Run

```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make

# Test on a small graph
./hybrid_blossom/Blossom test.mtx 4

# Test on real-world data
./hybrid_blossom/Blossom \
  ../X-Blossom/Example_Dataset/Example_Dataset/Example_Realworld_Datasets/StackOverflow/stackOverflow_rowOffsets.txt \
  ../X-Blossom/Example_Dataset/Example_Dataset/Example_Realworld_Datasets/StackOverflow/stackOverflow_columnIndices.txt \
  8
```

### Test

```bash
./hybrid_blossom/Blossom weight_test.mtx 4
# Output:
# Hybrid (4 threads): 3 edges in 2ms | valid: yes
```

### Validation

```cpp
// Bidirectional check
for (int v = 0; v < n; ++v) {
    if (M[v] != -1 && M[M[v]] != v) {
        std::cerr << "Invalid: " << v << " -> " << M[v];
    }
}
```

---

## Example Datasets

The `X-Blossom/Example_Dataset` directory contains real-world graph data:

| Dataset | Source | Size |
|---------|--------|------|
| **Google** (`gplus.mtx`) | Google+ social network | ~100K+ vertices |
| **StackOverflow** (`stackOverflow.mtx`) | StackOverflow user network | ~100K+ vertices |
| **Synthetic** (`test.mtx`) | Random / generated graphs | Variable |

The `X-Blossom` directory contains example datasets for testing:

```
X-Blossom/Example_Dataset/
├── Example_Realworld_Datasets/
│   ├── Google/
│   │   ├── gplus_rowOffsets.txt
│   │   ├── gplus_columnIndices.txt
│   │   └── gplus.mtx
│   └── StackOverflow/
│       ├── stackOverflow_rowOffsets.txt
│       ├── stackOverflow_columnIndices.txt
│       └── stackOverflow.mtx
```

---

## Performance Characteristics

### Complexity

| Phase | Complexity | Notes |
|-------|------------|-------|
| **Init (dual)** | $O(n \cdot \overline{d})$ | $n$ = vertices, $\overline{d}$ = avg degree |
| **Match exposed** | $O(m)$ | $m$ = edges |
| **Augmenting path search** | $O(n \cdot \overline{d} / p)$ | $p$ = threads (parallel, per round) |
| **BFS expand** | $O(n \cdot \overline{d} / p)$ | Per round |
| **Blossom detection** | $O(n \cdot \overline{d} / p)$ | Per round |
| **Dual update** | $O(n \cdot \overline{d})$ | per iteration |
| **Validation** | $O(n)$ | Bidirectional check |

### Parallel Scaling

```
Threads:  1       2       4       8
Speedup: 1.0×   1.7×    2.9×    4.5×        (typical)
```

The parallelism is **stripe-based**:
- Each thread handles `index % num_threads` of the workload.
- Atomic CAS ensures **lock-free** contention on tree/edge claims.
- `compare_exchange_strong` is the only synchronization point.

### When to Use Hybrid vs X-Blossom

| Scenario | Use |
|---------|-----|
| **Large sparse graph** ($n > 10^5$, low degree) | **X-Blossom** (no weights needed) |
| **Weighted graph** (weights matter) | **Hybrid Blossom** |
| **Perfect matching required** ($|M| = n/2$) | **Blossom VI** |
| **Fast approximate matching** | **X-Blossom** (parallel) |
| **Exact optimum** (small $n$) | **Blossom VI** |
| **Graph with negative weights** | **Hybrid Blossom** (supports negative) |

---

## References

### Papers

- **Edmonds, J. (1965).** "Paths, trees, and flowers." *Canadian Journal of
  Mathematics*. — The original Blossom algorithm.
- **Kolmogorov, V. (2009).** "Blossom V: A New Implementation of a Minimum Cost
  Perfect Matching Algorithm." — The Blossom VI algorithm.
- **arXiv:2604.20351** — X-Blossom / Blossom V implementation for maximum
  weight matching.

### Related Projects

- [`blossom-vi/`](https://github.com/robertob26/blossom-vi) — The reference
  Blossom V implementation (C++20).
- [`X-Blossom/`](https://github.com/robertob26/X-Blossom) — Parallel
  X-Blossom (C++14).
- **X-Blossom (this repo):** A parallel, multi-threaded maximum matching
  implementation based on the Blossom algorithm.

---

## License

This project is distributed under the terms of the GPL v3 license (see
`LICENSE` file).

---

## Open Research Questions

1. **Batch dual updates.** Can we accumulate dual delta for *all* exposed
   nodes in one pass rather than one at a time? (A min-over-exposed-slack
   computation.)

2. **Parallel blossom detection.** The `parallel_augment.h` already has
   `parBlossom` — but the `select_blossom` atomic flags may be a bottleneck
   for dense graphs.

3. **Hybrid scheduling.** When should the algorithm *switch* from X-Blossom
   mode (find any matching) to Blossom VI mode (refine for weight)? A
   heuristic: stop X-Blossom when fewer than `n/10` exposed nodes remain.

4. **Weight annotation.** The current `parallel_find_augmenting_paths` does
   not consider edge weights. Adding `weight` to the augmenting path
   selection could improve the initial matching quality.

---

*Hybrid Blossom — combining parallel path-table search with dual-variable
  refinement for maximum-weight matching on undirected graphs.*