#include "graph.h"

Graph::Graph(set<int> v, int INPUT_N) {
    this->V = v;
    this->adj.assign(INPUT_N, vector<int>(INPUT_N, 0));
    this->is_edge.assign(INPUT_N, vector<bool>(INPUT_N, false));
}

void Graph::add_edge(int from, int to, int weight) {
    adj[from][to] = weight;
    is_edge[from][to] = true;
}
