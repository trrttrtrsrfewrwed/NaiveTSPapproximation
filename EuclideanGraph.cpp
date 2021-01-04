
#include "EuclideanGraph.h"
#include "./Matching/Matching.h"
#include <unordered_set>
#include <algorithm>
#include <map>
#include <stack>
#include <iostream>
#include <chrono>
#include <queue>

EuclideanGraph::EuclideanGraph(std::vector<Vector> points) : points_(std::move(points)) {

}

void EuclideanGraph::FindTSPApproximation(const std::string& method_name) {
    std::vector<size_t> result;

    auto start = std::chrono::steady_clock::now();
    if (method_name == "christofides") {
        result = ConvertEulerToHamilton(FindEulerCycle(BuildMSTWithMatching(false)));
    } else if (method_name == "christofides naive") {
        result = ConvertEulerToHamilton(FindEulerCycle(BuildMSTWithMatching(true)));
    } else if (method_name == "antigreedy") {
        result = CalcAntigreedy();
    }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;

    std::cout << "Time: " << chrono::duration <double, milli> (diff).count() << " ms" << endl;

    double length = 0;
    std::cout << "Path: ";
    for (size_t i = 0; i < points_.size(); ++i) {
        size_t from = result[i];
        std::cout << result[i] << " ";
        size_t to = result[(i + 1) % points_.size()];
        length += (points_[from] - points_[to]).Norm();
    }
    std::cout << std::endl;

    std::cout << "Length: " << length << std::endl;

}

// Методы для алгоритма Кристофидеса

std::multiset<EuclideanGraph::Edge> EuclideanGraph::BuildMst() {
    std::multiset<EuclideanGraph::Edge> result;

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
                if (mst_vertices.find(from) == mst_vertices.end()) {
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

std::multiset<EuclideanGraph::Edge> EuclideanGraph::findPerfectMatching(const std::vector<size_t>& odds) {
    std::multiset<EuclideanGraph::Edge> result;

    // Используем готовый код для поиска совершенного паросочетания минимального веса
    // Источник: https://github.com/dilsonpereira/Minimum-Cost-Perfect-Matching/tree/a65c16ae27b398fe4b97ed3471eb067abb2d7949
    Graph O(odds.size());
    vector<double> costO;

    for (size_t i = 0; i < odds.size(); ++i) {
        for (size_t j = i + 1; j < odds.size(); ++j) {
            O.AddEdge(i, j);
            costO.push_back((points_[odds[i]] - points_[odds[j]]).Norm());
        }
    }

    Matching M(O);
    pair< list<int>, double > p = M.SolveMinimumCostPerfectMatching(costO);
    list<int> matching = p.first;

    for (int & it : matching) {
        pair<int, int> edge = O.GetEdge(it);
        size_t from = odds[edge.first];
        size_t to = odds[edge.second];

        result.insert({from, to, (points_[from] - points_[to]).Norm()});
    }

    return result;
}

std::multiset<EuclideanGraph::Edge> EuclideanGraph::findNaiveMatching(const std::vector<size_t>& odds, std::vector<size_t>& degrees) {
    std::multiset<EuclideanGraph::Edge> result;
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
            result.insert({from, to, (points_[from] - points_[to]).Norm()});
        }
    }
    return result;
}

std::multiset<EuclideanGraph::Edge> EuclideanGraph::BuildMSTWithMatching(bool naive) {
    std::multiset<EuclideanGraph::Edge> mst = BuildMst();
    std::vector<size_t> degrees = FindDegrees(mst);

    // Находим вершины нечётной степени в мин остове
    std::vector<size_t> odds;
    for (size_t i = 0; i < degrees.size(); ++i) {
        if (degrees[i] % 2 == 1) {
            odds.push_back(i);
        }
    }
    // Добавляем паросочетание к мин остову
    std::multiset<EuclideanGraph::Edge> matching = (naive)? findNaiveMatching(odds, degrees): findPerfectMatching(odds);
    mst.insert(matching.begin(), matching.end());
    return mst;
}

std::vector<size_t> EuclideanGraph::FindDegrees(const std::multiset<EuclideanGraph::Edge> &tree) {
    // степени вершин
    std::vector<size_t> degrees(tree.size() + 1, 0);
    for (auto edge: tree) {
        ++degrees[edge.first_idx_];
        ++degrees[edge.second_idx_];
    }
    return degrees;
}

std::vector<size_t> EuclideanGraph::FindEulerCycle(const std::multiset<EuclideanGraph::Edge> &graph) {
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

std::vector<size_t> EuclideanGraph::ConvertEulerToHamilton(const std::vector<size_t> &euler_cycle) {
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

// Методы для антижадного алгоритма

std::vector<size_t> EuclideanGraph::CalcAntigreedy() {
    // лог изменений
    std::stack<std::vector<size_t>> log;

    // ребра графа в порядке убывания
    std::vector<EuclideanGraph::Edge> edges;
    for (size_t from = 0; from < points_.size(); ++from) {
        for (size_t to = from + 1; to < points_.size(); ++to) {
            edges.emplace_back(from, to, (points_[from] - points_[to]).Norm());
        }
    }
    std::sort(edges.begin(), edges.end(), greater<>());

    // заполним матрицу смежности для графа номерами рёбер в порядке убывания
    std::vector<std::vector<size_t>> adj_matrix(points_.size(), std::vector<size_t>(points_.size()));
    for (size_t i = 0; i < edges.size(); ++i) {
        size_t from = edges[i].first_idx_;
        size_t to = edges[i].second_idx_;

        adj_matrix[from][to] = i;
        adj_matrix[to][from] = i;
    }

    // текущее состояние графа (все ребра не зафиксированы)
    std::vector<size_t> currState(edges.size(), 1);

    // Добавляем начальное состояние в лог
    log.push(currState);

    // Пока есть незафиксированные ребра, находим незафиксированное ребро максимального веса
    int max_edge_idx = std::find(currState.begin(), currState.end(), 1) - currState.begin();
    while (max_edge_idx != currState.size()) {
        currState = log.top();
        log.pop();

        // создаем состояние фиксации и состояние удаления ребра
        std::vector<size_t> deleteEdge = currState;
        std::vector<size_t> fixEdge = currState;

        deleteEdge[max_edge_idx] = 0;
        fixEdge[max_edge_idx] = 2;

        // Обрабатываем новые состояния
        Processing(deleteEdge, edges, adj_matrix);
        Processing(fixEdge, edges, adj_matrix);

        if (CheckGraph(fixEdge, edges, adj_matrix)) {
            log.push(fixEdge);
        }

        if (CheckGraph(deleteEdge, edges, adj_matrix)) {
            log.push(deleteEdge);
        }

        max_edge_idx = std::find(log.top().begin(), log.top().end(), 1) - log.top().begin();
    }

    std::vector<std::vector<size_t>> neighbours(points_.size());

    for (size_t i = 0; i < log.top().size(); ++i) {
        if (log.top()[i] == 2) {
            std::cout << i << " " << edges[i].first_idx_ << " " << edges[i].second_idx_ << "\n";
            size_t from = edges[i].first_idx_;
            size_t to = edges[i].second_idx_;
            neighbours[from].push_back(to);
            neighbours[to].push_back(from);
        }
    }

    std::vector<size_t> result;
    result.push_back(0);

    size_t prev = 0;
    size_t curr = neighbours[prev][0];
    for (size_t i = 1; i < points_.size(); ++i) {
        result.push_back(curr);
        if (neighbours[curr][0] == prev) {
            prev = curr;
            curr = neighbours[curr][1];
        } else {
            prev = curr;
            curr = neighbours[curr][0];
        }
    }

    return result;
}

bool EuclideanGraph::CheckGraph(const std::vector<size_t>& currState, const std::vector<EuclideanGraph::Edge>& edges, const std::vector<std::vector<size_t>>& adj_matrix) {
    for (size_t i = 0; i < points_.size(); ++i) {
        size_t fixed_cnt = 0;
        size_t deleted_cnt = 0;
        // Проходимся по смежным рёбрам
        for (size_t j = 0; j < points_.size(); ++j) {
            if (j == i) {
                continue;
            }
            size_t edge_number = adj_matrix[i][j];
            // Нашли добавленное ребро
            if (currState[edge_number] == 2) {
                ++fixed_cnt;
            }
            // Нашли удаленное ребро
            if (currState[edge_number] == 0) {
                ++deleted_cnt;
            }
        }

        if (fixed_cnt > 2 || deleted_cnt > points_.size() - 1 - 2) {
            return false;
        }
    }
    // В случае, если рёбер хотя бы число вершин - 1, нужно проверить связность
    if (std::count(currState.begin(), currState.end(), 2) >= points_.size() - 1) {
        std::vector<std::vector<size_t>> neighbours(points_.size());
        std::vector<size_t> visited(points_.size(), 0);

        // запомним соседей каждой вершины
        for (size_t i = 0; i < currState.size(); ++i) {
            if (currState[i] == 2) {
                size_t from = edges[i].first_idx_;
                size_t to = edges[i].second_idx_;
                neighbours[from].push_back(to);
                neighbours[to].push_back(from);
            }
        }

        std::queue<size_t> q;
        q.push(0);
        visited[0] = 1;

        while (!q.empty()) {
            size_t curr = q.front();
            q.pop();
            for (size_t neighbour: neighbours[curr]) {
                if (!visited[neighbour]) {
                    q.push(neighbour);
                    visited[neighbour] = 1;
                }
            }
        }
        if (std::count(visited.begin(), visited.end(), 1) < points_.size()) {
            return false;
        }
    }
    return true;
}

void EuclideanGraph::Processing(std::vector<size_t>& currState, const std::vector<EuclideanGraph::Edge>& edges, const std::vector<std::vector<size_t>>& adj_matrix) {
    std::vector<size_t> copyState;
    // Будем улучшать состояние, пока оно не стабилизируется
    do {
        // копируем текущее состояние
        copyState = currState;

        // Для каждой вершины, у которой добавлены два инцидентных ребра, удаляем все остальные ребра
        // Для каждой вершины, у которой удалены все ребра, кроме двух, добавляем эти два ребра
        for (size_t i = 0; i < points_.size(); ++i) {
            size_t fixed_cnt = 0;
            size_t deleted_cnt = 0;

            // Проходимся по смежным рёбрам
            for (size_t j = 0; j < points_.size(); ++j) {
                if (j == i) {
                    continue;
                }
                size_t edge_number = adj_matrix[i][j];
                // Нашли добавленное ребро
                if (currState[edge_number] == 2) {
                    ++fixed_cnt;
                }
                // Нашли удаленное ребро
                if (currState[edge_number] == 0) {
                    ++deleted_cnt;
                }
            }
            // Если два ребра добавлены, удаляем остальные
            if (fixed_cnt == 2) {
                for (size_t edge_number: adj_matrix[i]) {
                    if (currState[edge_number] != 2) {
                        currState[edge_number] = 0;
                    }
                }
            }
            // Если все ребра, кроме двух, удалены, добавляем оставшиеся 2 ребра
            if (deleted_cnt == (points_.size() - 1) - 2) {
                for (size_t edge_number: adj_matrix[i]) {
                    if (currState[edge_number] != 0) {
                        currState[edge_number] = 2;
                    }
                }
            }
        }
        // Если фиксированных рёбер меньше, чем количество вершин - 1, необходимо удалить возможности появления мини-циклов
        if (std::count(currState.begin(), currState.end(), 2) < points_.size() - 1) {
            for (size_t i = 0; i < currState.size(); ++i) {
                // Для каждого добавленного ребра
                if (currState[i] == 2) {
                    Edge edge = edges[i];
                    RemovingSmallCycles(edge.first_idx_, edge.second_idx_, currState, adj_matrix);
                }
            }
        }
    } while (copyState != currState);
}

void EuclideanGraph::RemovingSmallCycles(size_t n_start, size_t n_end, std::vector<size_t>& currState, const std::vector<std::vector<size_t>>& adj_matrix) {
    // Посещенные вершины
    std::vector<int> visited(points_.size(), 0);
    visited[n_start] = 1;
    visited[n_end] = 1;
    int curr = 0;
    do {
        // номер ребра, соединяющего текущую вершину с конечной
        int edge_to_end_number = adj_matrix[curr][n_end];

        // Если мы уже добавили ребро из текущей вершины в конечную, не должно быть ребра из неё же в начальную
        if (!visited[curr] && currState[edge_to_end_number] == 2) {
         visited[curr] = 1;

         // номер ребра, соединяющего текущую вершину с начальной
         int edge_to_start_number = adj_matrix[curr][n_start];
         // Удалим данное ребро
         currState[edge_to_start_number] = 0;

         // обновляем конечную вершину
         n_end = curr;
         curr = -1;
        }
        ++curr;
    } while ( curr < points_.size());
}
//----------------

EuclideanGraph::Edge::Edge(size_t first_idx, size_t second_idx, double weight) : first_idx_(first_idx), second_idx_(second_idx),
                                                                                 weight_(weight) {
}

bool EuclideanGraph::Edge::operator<(const EuclideanGraph::Edge &another) const {
    if (weight_ < another.weight_) {
        return true;
    }
    if (weight_ == another.weight_) {
        if (first_idx_ < another.first_idx_) {
            return true;
        }
        if (first_idx_ == another.first_idx_) {
            return second_idx_ < another.second_idx_;
        }
    }
    return false;

}

bool EuclideanGraph::Edge::operator>(const EuclideanGraph::Edge &another) const {
    if (weight_ > another.weight_) {
        return true;
    }
    if (weight_ == another.weight_) {
        if (first_idx_ > another.first_idx_) {
            return true;
        }
        if (first_idx_ == another.first_idx_) {
            return second_idx_ > another.second_idx_;
        }
    }
    return false;
}
