cmake_minimum_required(VERSION 3.10)
project(hybrid_x_blossom)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build flags: -pthread for std::thread, -O3 for performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -O3 -march=native -Wall -Wextra")

# ---------- Shared library of core algorithm sources ----------
set(CORE_SOURCES
    hybrid_blossom.cpp
    hybrid_blossom.h
    parallel_augment.cpp
    parallel_augment.h
    weighted_matching.cpp
    weighted_matching.h
)

# ---------- Main CLI binary ----------
add_executable(hybrid_blossom
    main.cpp
    ${CORE_SOURCES}
)

# ---------- Exact correctness test binary ----------
add_executable(test_exact
    test_exact.cpp
    ${CORE_SOURCES}
)

# ---------- Enable CTest ----------
enable_testing()
add_test(NAME ExactTests COMMAND test_exact)
add_test(NAME QuickTests COMMAND test_exact --quick)
