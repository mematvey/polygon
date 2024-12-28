#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <algorithm>
#include <cmath>
#include <chrono>


struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& pair) const {
        // Вычисляем хэщ для каждого элемента отдельно и комбинируем их
        auto hash1 = std::hash<T1>{}(pair.first);
        auto hash2 = std::hash<T2>{}(pair.second);
        return hash1 ^ (hash2 << 1);
    }
};


struct Node {
    std::pair<long double, long double> coord;
    std::vector<std::pair<Node*, long double>> neighbours;

    Node(std::pair<long double, long double> location) : coord(location) {}
};


struct Graph {
    std::vector<Node*> nodes;
    std::unordered_map<std::pair<long double, long double>, Node*, pair_hash> nodes_map;

    // Удаляем все, что создали через new
    ~Graph() {
        for (auto node : nodes) {
            for (auto& [neighbour, _] : node->neighbours) {
                delete neighbour;
            }
            node->neighbours.clear();
            delete node; 
        }
        nodes.clear();

        nodes_map.clear();
    }
};


Node* create_or_get_node(Graph& graph, const std::pair<long double, long double>& coord) {
    auto existNode = graph.nodes_map.find(coord);
    if (existNode != graph.nodes_map.end()) {
        return existNode->second;
    }

    Node* node = new Node(coord);
    graph.nodes.push_back(node);
    graph.nodes_map[coord] = node;
    return node;
}


void create_node_and_edges(Graph& graph, const std::string& str) {
    std::stringstream sstr(str);
    std::string segment;

    std::vector<long double> result;
    while (std::getline(sstr, segment, ',')) {
        size_t dotCommaPos = segment.find(';');

        // Если есть еще ';', то делим по ним дальше
        if (dotCommaPos != std::string::npos) {
            if (!segment.substr(dotCommaPos + 1).empty()) {
                result.push_back(std::stold(segment.substr(0, dotCommaPos)));
                result.push_back(std::stold(segment.substr(dotCommaPos + 1)));
            } else {
                result.push_back(std::stold(segment.substr(0, dotCommaPos)));
            }
        } else {
            size_t colonPos = segment.find(':');

            // Проверяем, присутствует ли ':' в segment
            if (colonPos != std::string::npos) {
                result.push_back(std::stold(segment.substr(0, colonPos)));
                result.push_back(std::stold(segment.substr(colonPos + 1)));
            } else {
                result.push_back(std::stold(segment));
            }
        }
    }

    if (result.size() < 2) {
        throw std::invalid_argument("Incorrect format for Node. Given only 1 coordinate.");
    }

    Node* node = create_or_get_node(graph, {result[0], result[1]});

    size_t i = 2;
    while (i + 2 < result.size()) {
        Node* neighbour = create_or_get_node(graph, {result[i], result[i + 1]});
        long double distance = result[i + 2];
        node->neighbours.push_back({neighbour, distance});
        i += 3;
    }
}


Graph create_graph(const std::vector<std::string>& sparsed_graph) {
    Graph graph;

    for (const std::string& str : sparsed_graph) {
        create_node_and_edges(graph, str);
    }

    return graph;
}


template <typename T>
std::vector<T> read_array_from_file(const char filename[]) {
    std::vector<T> data;

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file");
    }

    std::string value;
    while (std::getline(file, value)) {
        T value_processed;

        if constexpr (std::is_same<T, int>::value) {
            value_processed = std::stoi(value);
        } else if constexpr (std::is_same<T, float>::value) {
            value_processed = std::stof(value);
        } else {
            value_processed = value;
        }

        data.push_back(value_processed);
    }

    return data;
}


Node* find_closest_node(const Graph& graph, long double lon, long double lat) {
    // Функция для поиска ближайшего узла в графе по координатам
    long double min_distance = 99999999;
    Node* node_found = nullptr;

    for (auto node : graph.nodes) {
        long double distance = std::sqrt(
            std::pow(node->coord.first - lon, 2) + 
            std::pow(node->coord.second - lat, 2)
        );

        if (distance < min_distance) {
            node_found = node;
            min_distance = distance;
        }
    }

    return node_found;
}


void print_path(const std::vector<Node*>& path) {
    if (path.empty()) {
        std::cout << "No path found!" << '\n';
        return;
    }

    for (const auto& node : path) {
        std::cout << "(" << node->coord.first << ", " << node->coord.second << ")" << '\n';
    }
}


// Перегружаем функцию для взвешенного графа
void print_path(const std::vector<std::pair<Node*, long double>>& path) {
    if (path.empty()) {
        std::cout << "No path found!" << '\n';
        return;
    }

    long double totalDistance = 0;
    for (auto& [node, distance] : path) {
        totalDistance += distance;
        std::cout << "(" << node->coord.first << ", " << node->coord.second << ")" << '\n';
    }

    std::cout << "Total distance: " << totalDistance << '\n';
}



std::vector<Node*> BFS(Graph* graph, Node* start, Node* end) {
    if (start == nullptr || end == nullptr) {
        std::cout << "Incorrect start or end vertex" << '\n';
        return {};
    } 

    if (start == end){
        return {start};
    }

    std::queue<Node*> line;
    std::unordered_set<Node*> visited;
    std::unordered_map<Node*, Node*> parents;

    line.push(start);
    parents[start] = nullptr;

    while (!line.empty()) {
        Node* node = line.front();
        line.pop();

        if (visited.find(node) != visited.end()){
            continue;
        }
        visited.insert(node);

        if (node == end) {
            // Собираем путь и потом переворачиваем его, тк он в обратном порядке
            std::vector<Node*> path;
            for (Node* current = end; current != nullptr; current = parents[current]) {
                path.push_back(current);
            }

            std::reverse(path.begin(), path.end());
            return path;
        }

        // Распаковываем соседей красиво через auto. Расстояния нам не нужны
        for (auto& [neighbour, _]: node->neighbours) {
            if (visited.find(neighbour) == visited.end()) {
                line.push(neighbour);
                parents[neighbour] = node;
            }
        }
    }

    return {};
}


std::vector<Node*> DFS(Graph* graph, Node* start, Node* end) {
    if (start == nullptr || end == nullptr) {
        std::cout << "Incorrect start or end vertex" << '\n';
        return {};
    } 

    if (start == end){
        return {start};
    }

    std::stack<Node*> line;
    std::unordered_set<Node*> visited;
    std::unordered_map<Node*, Node*> parents;

    line.push(start);
    parents[start] = nullptr;

    while (!line.empty()) {
        Node* node = line.top();
        line.pop();

        if (visited.find(node) != visited.end()){
            continue;
        }
        visited.insert(node);

        if (node == end) {
            // Собираем путь и потом переворачиваем его, тк он в обратном порядке
            std::vector<Node*> path;
            for (Node* current = end; current != nullptr; current = parents[current]) {
                path.push_back(current);
            }

            std::reverse(path.begin(), path.end());
            return path;
        }

        // Распаковываем соседей красиво через auto. Расстояния нам не нужны
        for (auto& [neighbour, _]: node->neighbours) {
            if (visited.find(neighbour) == visited.end()) {
                line.push(neighbour);
                parents[neighbour] = node;
            }
        }
    }

    return {};
}


std::vector<std::pair<Node*, long double>> Dijkstra(Graph* graph, Node* start, Node* end) {
    if (start == nullptr || end == nullptr) {
        std::cout << "Incorrect start or end vertex" << '\n';
        return {};
    }

    std::priority_queue<
        std::pair<long double, Node*>, 
        std::vector<std::pair<long double, Node*>>, 
        std::greater<>
    > pq;

    std::unordered_map<Node*, long double> dist;
    std::unordered_map<Node*, Node*> parents;
    std::unordered_set<Node*> visited;

    for (auto node : graph->nodes) {
        dist[node] = std::numeric_limits<long double>::infinity();
    }

    pq.push({0, start});
    dist[start] = 0;
    parents[start] = nullptr;

    while (!pq.empty()) {
        auto [current_dist, node] = pq.top();
        pq.pop();

        // Пропускаем уже посещенные вершины
        if (visited.find(node) != visited.end()) {
            continue;
        }
        visited.insert(node);

        if (node == end) {
            std::vector<std::pair<Node*, long double>> path;
            // Собираем путь
            for (Node* current = end; current != nullptr; current = parents[current]) {
                path.push_back({current, dist[current]});
            }
            std::reverse(path.begin(), path.end());
            return path;
        }

        for (auto& [neighbour, weight] : node->neighbours) {
            if (visited.find(neighbour) != visited.end()){
              continue;
            }

            // Делаем операции с вершиной только если удалось уменьшить расстояние
            long double new_dist = current_dist + weight;
            if (new_dist < dist[neighbour]) {
                dist[neighbour] = new_dist;
                parents[neighbour] = node;
                pq.push({new_dist, neighbour});
            }
        }
    }

    return {};
}


long double heuristic(const Node& a, const Node& b) {
    return std::sqrt(
        std::pow(a.coord.first - b.coord.first, 2) + 
        std::pow(a.coord.second - b.coord.second, 2)
    );
}


std::vector<std::pair<Node*, long double>> AStar(Graph* graph, Node* start, Node* end) {
    if (start == nullptr || end == nullptr) {
        std::cout << "Incorrect start or end vertex" << '\n';
        return {};
    }

    std::priority_queue<
        std::pair<long double, Node*>, 
        std::vector<std::pair<long double, Node*>>, 
        std::greater<>
    > pq;

    std::unordered_map<Node*, Node*> parents;
    std::unordered_map<Node*, long double> dist;
    std::unordered_set<Node*> visited;

    for (Node* node : graph->nodes) {
        dist[node] = std::numeric_limits<long double>::infinity();
    }

    dist[start] = 0;
    // Тк наша функция для вычисления эвристики принимает ссылки, 
    // то мы передаем объект, на который указывает указатель
    pq.push({heuristic(*start, *end), start});
    parents[start] = nullptr;

    while (!pq.empty()) {
        auto [_, node] = pq.top();
        pq.pop();

        if (visited.find(node) != visited.end()) {
            continue; 
        }
        visited.insert(node);

        if (node == end) {
            // Собираем путь и потом переворачиваем его, тк он в обратном порядке
            std::vector<std::pair<Node*, long double>> path;
            for (Node* node = end; node != nullptr; node = parents[node]) {
                path.push_back({node, dist[node]});
            }
            std::reverse(path.begin(), path.end());
            return path;
        }

        for (auto& [neighbour, weight] : node->neighbours) {
            if (visited.find(neighbour) != visited.end()) {
                continue; 
            }

            long double newDistance = dist[node] + weight;
            
            if (newDistance < dist[neighbour]) {
                parents[neighbour] = node;
                dist[neighbour] = newDistance;
                pq.push({newDistance + heuristic(*neighbour, *end), neighbour});
            }
        }
    }

    return {};
}


void print_distance(const std::vector<std::pair<Node*, long double>>& path) {
    if (path.empty()) {
        std::cout << "No path found!" << '\n';
        return;
    }

    long double totalDistance = 0;
    for (auto& [node, distance] : path) {
        totalDistance += distance;
    }

    std::cout << "Total distance: " << totalDistance << '\n';
}


int main() {
    const char fileName[] = "./spb_graph.txt";

    Graph graph = create_graph(read_array_from_file<std::string>(fileName));

    Node* startNode = find_closest_node(graph, 30.308108, 59.957238);
    Node* endNode = find_closest_node(graph, 30.49996, 59.936558);

    // BFS
    auto start_bfs = std::chrono::high_resolution_clock::now();
    auto bfsPath = BFS(&graph, startNode, endNode);
    auto end_bfs = std::chrono::high_resolution_clock::now();
    auto duration_bfs = std::chrono::duration<double>(end_bfs - start_bfs).count();
    std::cout << "BFS: " << duration_bfs << " seconds\n";

    // DFS
    auto start_dfs = std::chrono::high_resolution_clock::now();
    auto dfsPath = DFS(&graph, startNode, endNode);
    auto end_dfs = std::chrono::high_resolution_clock::now();
    auto duration_dfs = std::chrono::duration<double>(end_dfs - start_dfs).count();
    std::cout << "DFS: " << duration_dfs << " seconds\n";

    // Dijkstra
    auto start_dijkstra = std::chrono::high_resolution_clock::now();
    auto dijkstraPath = Dijkstra(&graph, startNode, endNode);
    auto end_dijkstra = std::chrono::high_resolution_clock::now();
    auto duration_dijkstra = std::chrono::duration<double>(end_dijkstra - start_dijkstra).count();
    std::cout << "Dijkstra: " << duration_dijkstra << " seconds\n";
    print_distance(dijkstraPath);
    
    // A*
    auto start_astar = std::chrono::high_resolution_clock::now();
    auto aStarPath = AStar(&graph, startNode, endNode);
    auto end_astar = std::chrono::high_resolution_clock::now();
    auto duration_astar = std::chrono::duration<double>(end_astar - start_astar).count();
    std::cout << "A*: " << duration_astar << " seconds\n";
    print_distance(aStarPath);

    return 0;
}
