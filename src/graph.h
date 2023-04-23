#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h> // TODO: use more specific libraries

using namespace std;

class Graph {
    public:
    vector<vector<int>> adj;
    vector<vector<bool>> is_edge;
    set<int> V;
    Graph(set<int> v, int INPUT_N);
    void add_edge(int from, int to, int weight);
};

struct SSSP_Result {
    bool has_negative_cycle;
    Graph* shortest_paths_tree;
};

#endif

