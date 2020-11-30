#include "Graph.h"

#include <utility>
#include <unordered_set>
#include <algorithm>
#include <map>
#include <stack>
#include <iostream>

Graph::Graph(std::vector<Vector> points) : points_(std::move(points)) {

}

void Graph::FindTSPApproximation() {
    std::vector<size_t> result = CovertEulerToHamilton(FindEulerCycle(BuildMSTWithPerfectMatching()));

    double length = 0;
    std::cout << "Путь: ";
    for (size_t i = 0; i < points_.size(); ++i) {
        size_t from = result[i];
        std::cout << result[i] << " ";
        size_t to = result[(i + 1) % points_.size()];
        length += (points_[from] - points_[to]).Norm();
    }
    std::cout << std::endl;

    std::cout << "Длина: " << length << std::endl;

}

std::multiset<Graph::Edge> Graph::BuildMst() {
    std::multiset<Graph::Edge> result;

    // Вершины, находящиеся в мин остове на данном шаге алгоритма
    std::unordered_set<size_t> mst_vertices;

    // Ребра, соединяющие вершины простроенной части мин остова с оставшимися вершинами в графе
    std::set<Edge> connecting;

    // Добавляем первую вершину в мин остов
    mst_vertices.insert(0);
    for (size_t i = 1; i < points_.size(); ++i) {
        connecting.insert({0, i, (points_[i] - points_[0]).Norm()});
    }

    // Добавляем ребро минимального веса, соединяющее вершины,
    // одна из которых лежит в построенной части мин остова, а другая - нет,
    // пока не будет построен мин остов
    for (size_t i = 1; i < points_.size(); ++i) {
        while (true) {
            Edge min = *connecting.begin();

            connecting.erase(min);

            size_t from = min.first_idx_;
            size_t to = min.second_idx_;

            if (mst_vertices.find(from) == mst_vertices.end() || mst_vertices.find(to) == mst_vertices.end()) {
                result.insert({from, to, (points_[from] - points_[to]).Norm()});
                if (mst_vertices.find(from) == mst_vertices.end() ) {
                    mst_vertices.insert(from);
                    for (size_t j = 0; j < points_.size(); ++j) {
                        if (mst_vertices.find(j) != mst_vertices.end()) {
                            continue;
                        }
                        connecting.insert({from, j, (points_[from] - points_[j]).Norm()});
                    }
                } else {
                    mst_vertices.insert(to);
                    for (size_t j = 0; j < points_.size(); ++j) {
                        if (mst_vertices.find(j) != mst_vertices.end()) {
                            continue;
                        }
                        connecting.insert({to, j, (points_[to] - points_[j]).Norm()});
                    }
                }
                break;
            }
        }
    }
    return result;
}

std::multiset<Graph::Edge> Graph::BuildMSTWithPerfectMatching() {
    std::multiset<Graph::Edge> mst = BuildMst();
    std::vector<size_t> degrees = FindDegrees(mst);

    // Находим вершины нечётной степени в мин остове
    std::vector<size_t> odds;
    for (size_t i = 0; i < degrees.size(); ++i) {
        if (degrees[i] % 2 == 1) {
            odds.push_back(i);
        }
    }
    size_t odds_cnt = odds.size();

    // Находим ребра, соединяющие вершины нечётных степеней
    std::vector<Edge> odd_edges;

    for (size_t i = 0; i < odds.size(); ++i) {
        for (size_t j = i + 1; j < odds.size(); ++j) {
            odd_edges.emplace_back(odds[i], odds[j], (points_[odds[i]] - points_[odds[j]]).Norm());
        }
    }

    // Упорядочим по возрастанию ребра, соединяющие вершины нечётных степеней
    std::sort(odd_edges.begin(), odd_edges.end());

    // Находим жадно приближение совершенного паросочетания минимального веса
    for (size_t i = 0; i < odd_edges.size() && odds_cnt > 0; ++i) {
        size_t from = odd_edges[i].first_idx_;
        size_t to = odd_edges[i].second_idx_;

        if (degrees[from] % 2 == 1 && degrees[to] % 2 == 1) {
            ++degrees[from];
            ++degrees[to];
            odds_cnt -= 2;
            mst.insert({from, to, (points_[from] - points_[to]).Norm()});
        }
    }
    return mst;
}

std::vector<size_t> Graph::FindDegrees(const std::multiset<Graph::Edge>& tree) {
    // степени вершин
    std::vector<size_t> degrees(tree.size() + 1, 0);
    for (auto edge: tree) {
        ++degrees[edge.first_idx_];
        ++degrees[edge.second_idx_];
    }
    return degrees;
}

std::vector<size_t> Graph::FindEulerCycle(const std::multiset<Graph::Edge> &graph) {
    std::vector<size_t> result;

    // Заполним матрицу смежности для графа
    std::vector<std::map<size_t, size_t>> adj_matrix(graph.size() + 1);
    for (auto e: graph) {
        ++adj_matrix[e.first_idx_][e.second_idx_];
        ++adj_matrix[e.second_idx_][e.first_idx_];
    }

    std::stack<size_t> stack;
    stack.push(0);
    while (!stack.empty()) {
        size_t curr = stack.top();
        // Если из данной вершины есть непосещенное ребро, добавляем её в стек и проходим дальше по ребру
        if (!adj_matrix[curr].empty()) {
            size_t next = adj_matrix[curr].begin()->first;
            // "Посещаем" ребро - удаляем его из матрицы смежности
            if (adj_matrix[curr][next] == 1) {
                adj_matrix[curr].erase(next);
            } else {
                --adj_matrix[curr][next];
            }
            if (adj_matrix[next][curr] == 1) {
                adj_matrix[next].erase(curr);
            } else {
                --adj_matrix[next][curr];
            }
            stack.push(next);
        } else { // Иначе записываем вершины в ответ, пока снова не встретим вершину, у которой не все ребра посещены
            result.push_back(curr);
            stack.pop();
        }
    }
    return result;
}

std::vector<size_t> Graph::CovertEulerToHamilton(const std::vector<size_t> &euler_cycle) {
    std::unordered_set<size_t> hamilton_vertices;
    std::vector<size_t> result;
    for (size_t vertex: euler_cycle) {
        if (hamilton_vertices.find(vertex) == hamilton_vertices.end()) {
            hamilton_vertices.insert(vertex);
            result.push_back(vertex);
        }
    }
    return result;
}


Graph::Edge::Edge(size_t first_idx, size_t second_idx, double weight) : first_idx_(first_idx), second_idx_(second_idx),
                                                                        weight_(weight) {
}
