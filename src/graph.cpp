#include "graph.h"

Graph::Graph(int n, int m) {
    this->n = n;
    this->m = m;
    this->adj.assign(n, vector<int>());
}

void Graph::add_edge(int from, int to) {
    assert(from < n);
    assert(to < n);
    adj[from].push_back(to);
}

SSSP_Result sssp(Graph* graph) {
    // TODO implement this (facile a dirsi, e' praticamente l'intera tesi)
}
