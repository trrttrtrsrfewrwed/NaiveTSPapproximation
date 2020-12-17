cmake_minimum_required(VERSION 3.17)
project(untitled1)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(lib)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(untitled1 main.cpp EuclideanGraph.h EuclideanGraph.cpp Vector.h Vector.cpp Matching/BinaryHeap.h Matching/BinaryHeap.cpp Matching/Globals.h Matching/Graph.h Matching/Graph.cpp Matching/Matching.h Matching/Matching.cpp)

# 'Google_Tests_run' is the target name
# 'test1.cpp tests2.cpp' are source files with tests
add_executable(Google_Tests_run test1.cpp)
target_link_libraries(Google_Tests_run gtest gtest_main)