#include "graph.h"

Graph::Graph(vector<int> v) {
    this->V = v;
    this->adj.assign(MAX_N, vector<int>(MAX_N, 0));
    this->is_edge.assign(MAX_N, vector<bool>(MAX_N, false));
}

void Graph::add_edge(int from, int to, int weight) {
    assert(from < MAX_N);
    assert(to < MAX_N);
    adj[from][to] = weight;
    is_edge[from][to] = true;
}
