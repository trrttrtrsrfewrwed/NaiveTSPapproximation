#include <iostream>
#include <vector>
#include "Vector.cpp"
#include "Graph.h"

int main() {
    size_t N;
    std::cin >> N;

    std::vector<Vector> points(N);
    for (size_t i = 0; i < N; ++i) {
        size_t unused;
        std::cin >> unused;
        std::cin >> points[i];
    }

    Graph graph(points);
    graph.FindTSPApproximation();
    return 0;
}
