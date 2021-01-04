#include <iostream>
#include <vector>
#include "EuclideanGraph.h"
#include <chrono>

int main() {
    size_t N;
    std::cin >> N;

    std::vector<Vector> points(N);
    for (size_t i = 0; i < N; ++i) {
        size_t unused;
        std::cin >> unused;
        std::cin >> points[i];
    }

    EuclideanGraph graph(points);
    std::string method;

    std::cin >> method;

    graph.FindTSPApproximation(method);
    return 0;
}
