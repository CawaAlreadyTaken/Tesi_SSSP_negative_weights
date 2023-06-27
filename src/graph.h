#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <set>

using namespace std;

class Graph {
    public:
    Graph(set<int> v);
    vector<vector<int>> adj;
    vector<vector<int>> edges;
    set<int> V;
    void add_edge(int from, int to, int weight);
};

struct SSSP_Result {
    bool has_negative_cycle;
    Graph* shortest_paths_tree;
};

#endif

