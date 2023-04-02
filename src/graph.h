#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h> // TODO: use more specific libraries

using namespace std;

class Graph {
    int n, m;
    public:
    vector<vector<int>> adj;
    Graph(int n, int m);
    void add_edge(int from, int to);
};

struct SSSP_Result {
    bool has_negative_cycle;
    Graph* shortest_paths_tree;
};

SSSP_Result sssp(Graph* graph);

#endif

