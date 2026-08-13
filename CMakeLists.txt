cmake_minimum_required(VERSION 3.10)
project(hybrid_x_blossom)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -O3 -march=native")

add_executable(hybrid_blossom
        main.cpp
        hybrid_blossom.cpp
        hybrid_blossom.h
        weighted_matching.cpp
        weighted_matching.h
        parallel_augment.cpp
        parallel_augment.h
        test_hybrid.cpp
)