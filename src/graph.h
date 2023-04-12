#ifndef GRAPH_H
#define GRAPH_H

#include <bits/stdc++.h> // TODO: use more specific libraries

using namespace std;

int INPUT_N;
int MAX_N = 10000;

class Graph {
    public:
    vector<vector<int>> adj;
    vector<vector<bool>> is_edge;
    vector<int> V;
    Graph(vector<int> v);
    void add_edge(int from, int to, int weight);
};

struct SSSP_Result {
    bool has_negative_cycle;
    Graph* shortest_paths_tree;
};

#endif

