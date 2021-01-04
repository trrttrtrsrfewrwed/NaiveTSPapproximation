#ifndef UNTITLED1_EUCLIDEANGRAPH_H
#define UNTITLED1_EUCLIDEANGRAPH_H

#include <vector>
#include <set>
#include "Vector.h"

class EuclideanGraph {
public:
    explicit EuclideanGraph(std::vector<Vector> points);

    /**
     * Находит путь, приближающий оптимальный путь в задаче коммивояжёра
     * @param method_name - название используемого алгоритма
     */
    void FindTSPApproximation(const std::string& method_name = "christofides");
private:
    /**
     * Представление ребра в графе
     */
    struct Edge {
        Edge() = default;
        Edge(size_t first_idx, size_t second_idx, double weight);
        bool operator<(const Edge &another) const;
        bool operator>(const Edge &another) const;
        size_t first_idx_{};
        size_t second_idx_{};
        double weight_{};
    };

    /**
     * Находит минимальное остовное дерево
     * https://clck.ru/L6AY7 - описание используемого алгоритма Прима
     *
     * @return набор ребер, образующий минимальное остовное дерево
     */
    std::multiset<EuclideanGraph::Edge> BuildMst();

    /**
     *
     * @param odds номера вершин нечётных степеней в мин остове
     * @return совершенное сочетание минимального веса
     */
    std::multiset<EuclideanGraph::Edge> findPerfectMatching(const std::vector<size_t>& odds);

    /**
     *
     * @param odds номера вершин нечётных степеней в мин остове
     * @param degrees степени вершин графа
     * @return жадное приближение совершенного сочетания минимального веса
     */
    std::multiset<EuclideanGraph::Edge> findNaiveMatching(const std::vector<size_t>& odds, std::vector<size_t>& degrees);

    /**
     *
     * @param naive - флаг, указывающий, какой алгоритм поиска паросочетания используется, наивный или нет
     * @return минимальное остовное дерево + паросочетание
     */
    std::multiset<EuclideanGraph::Edge> BuildMSTWithMatching(bool naive);


    /**
     * Находит степени вершин в дереве
     * @param tree - дерево
     * @return степени каждой из вершин дерева
     */
    static std::vector<size_t> FindDegrees(const std::multiset<EuclideanGraph::Edge>& tree);

    /**
     * Находит эйлеров цикл
     * https://cutt.ly/XhkERxR - описание алгоритма
     * @param graph - (мульти)граф, в котором существует эйлеров цикл
     * @return - номера вершин в порядке обхода эйлерова цикла
     */
    static std::vector<size_t> FindEulerCycle(const std::multiset<EuclideanGraph::Edge>& graph);

    /**
     * Преобразует эйлеров цикл в гамильтонов, путём пропуска повторяющихся вершин
     * @param euler_cycle
     * @return
     */
    static std::vector<size_t> ConvertEulerToHamilton(const std::vector<size_t>& euler_cycle);

    std::vector<Vector> points_;

    /**
     *  Находит путь, приближающий оптимальный путь в задаче коммивояжёра, с помощью антижадного алгоритма
     * @return вершины пути
     */
    std::vector<size_t> CalcAntigreedy();

    /**
     *
     * @param currState - текущее состояние рёбер (добавлено/не зафиксировано/удалено)
     * @param edges - ребра графа в порядке убывания веса
     * @param adj_matrix - матрица смежности, в которой для пары индексов хранится номер ребра в массиве edges
     * @return выполнены ли 2 условия:
     * 1 - нет циклов, содержащих не все вершины.
     * 2 - каждой вершине инцидентно не более 2 добавленных рёбер и не более (n - 3) удаленных рёбер
     */
    bool CheckGraph(const std::vector<size_t>& currState, const std::vector<EuclideanGraph::Edge>& edges, const std::vector<std::vector<size_t>>& adj_matrix);

    /**
     * Изменяет текущее состояние рёбер:
     *  1. Удаление возможности появления мелких циклов.
     *  2. Удаление ненужных ребер в тех случаях, когда среди инцидентных вершине ребер уже
     *  есть два фиксированных ребра.
     *  3. Фиксирование ребер, если у конкретной вершины осталось всего два инцидентных
     *  ей не удаленных ребра
     * @param currState - текущее состояние рёбер (добавлено/не зафиксировано/удалено)
     * @param edges - ребра графа в порядке убывания веса
     * @param adj_matrix - матрица смежности, в которой для пары индексов хранится номер ребра в массиве edges
     */
    void Processing(std::vector<size_t>& currState, const std::vector<EuclideanGraph::Edge>& edges, const std::vector<std::vector<size_t>>& adj_matrix);

    /**
     * Удаляет возможности появления мелких циклов
     * @param n_start - начальная вершинап добавленного ребра
     * @param n_end - конечная вершина добавленного ребра     * @param edges - ребра графа в порядке убывания веса

     * @param currState - текущее состояние рёбер (добавлено/не зафиксировано/удалено)
     * @param adj_matrix - матрица смежности, в которой для пары индексов хранится номер ребра в массиве edges
     */
    void RemovingSmallCycles(size_t n_start, size_t n_end, std::vector<size_t> &currState,
                             const std::vector<std::vector<size_t>> &adj_matrix);
};
#endif //UNTITLED1_EUCLIDEANGRAPH_H