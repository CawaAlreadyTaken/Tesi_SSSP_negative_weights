#include "graph.h"
#include "utils.h"

Graph::Graph(set<int> v) {
    this->V = v;
    this->adj.assign(INPUT_N, vector<int>());
    this->edges.assign(INPUT_N, vector<int>());
}

void Graph::add_edge(int from, int to, int weight) {
    adj[from].push_back(weight);
    edges[from].push_back(to);
}
