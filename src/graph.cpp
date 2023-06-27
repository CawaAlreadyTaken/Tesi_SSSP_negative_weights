#include "graph.h"
#include "utils.h"

Graph::Graph(set<int> v) {
    this->V = v;
    this->adj.assign(INPUT_N, vector<int>(INPUT_N, 0));
    this->edges.assign(INPUT_N, vector<int>());
}

void Graph::add_edge(int from, int to, int weight) {
    adj[from][to] = weight;
    edges[from].push_back(to);
}
