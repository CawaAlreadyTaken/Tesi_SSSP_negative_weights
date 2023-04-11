#include "graph.h"

Graph::Graph(int n, int m) {
    this->n = n;
    this->m = m;
    this->adj.assign(n, vector<int>());
    this->weights.assign(n, vector<int>());
}

void Graph::add_edge(int from, int to, int weight) {
    assert(from < n);
    assert(to < n);
    adj[from].push_back(to);
    weights[from].push_back(weight);
}
