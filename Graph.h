#pragma once

#include <vector>
#include <set>
#include "Vector.cpp"

class Graph {
public:
    explicit Graph(std::vector<Vector> points);

    /**
     * Находит путь, приближающий оптимальный путь в задаче коммивояжёра
     */
    void FindTSPApproximation();
private:
    /**
     * Представление ребра в графе
     */
    struct Edge {
        Edge() = default;
        Edge(size_t first_idx, size_t second_idx, double weight);
        size_t first_idx_;
        size_t second_idx_;
        double weight_;
    };

    /**
     * Находит минимальное остовное дерево
     * https://clck.ru/L6AY7 - описание используемого алгоритма Прима
     *
     * @return набор ребер, образующий минимальное остовное дерево
     */
    std::multiset<Graph::Edge> BuildMst();

    /**
     * @return минимальное остовное дерево + жадное приближение совершенного сочетания минимального веса
     */
    std::multiset<Graph::Edge> BuildMSTWithPerfectMatching();

    /**
     * Находит степени вершин в дереве
     * @param tree - дерево
     * @return степени каждой из вершин дерева
     */
    static std::vector<size_t> FindDegrees(const std::multiset<Graph::Edge>& tree);

    /**
     * Находит эйлеров цикл
     * https://cutt.ly/XhkERxR - описание алгоритма
     * @param graph - (мульти)граф, в котором существует эйлеров цикл
     * @return - номера вершин в порядке обхода эйлерова цикла
     */
    static std::vector<size_t> FindEulerCycle(const std::multiset<Graph::Edge>& graph);

    /**
     * Преобразует эйлеров цикл в гамильтонов, путём пропуска повторяющихся вершин
     * @param euler_cycle
     * @return
     */
    static std::vector<size_t> CovertEulerToHamilton(const std::vector<size_t>& euler_cycle);

    std::vector<Vector> points_;
};
