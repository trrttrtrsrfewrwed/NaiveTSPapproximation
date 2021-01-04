cmake_minimum_required(VERSION 3.17)
project(untitled1)

set(CMAKE_CXX_STANDARD 14)

add_executable(untitled1 main.cpp EuclideanGraph.h EuclideanGraph.cpp Vector.h Vector.cpp Matching/BinaryHeap.h Matching/BinaryHeap.cpp Matching/Globals.h Matching/Graph.h Matching/Graph.cpp Matching/Matching.h Matching/Matching.cpp)

set_target_properties( untitled1
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )