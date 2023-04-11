#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h> // TODO: use more specific libraries

using namespace std;

class Graph {
    public:
    int n, m;
    vector<vector<int>> adj;
    vector<vector<int>> weights;
    Graph(int n, int m);
    void add_edge(int from, int to, int weight);
};

struct SSSP_Result {
    bool has_negative_cycle;
    Graph* shortest_paths_tree;
};

#endif

