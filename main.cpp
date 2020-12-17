#include <iostream>
#include <vector>
#include "EuclideanGraph.h"

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
    std::cout << "Running ... \n";
    graph.FindTSPApproximation("christofides naive");
    graph.FindTSPApproximation("christofides");
    graph.calcAntigreedy();
    return 0;
}
