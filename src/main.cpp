#include <bits/stdc++.h> // TODO: use more specific libraries

using namespace std;

class Graph {
    int n, m;
    public:
    vector<vector<int>> adj;
    Graph(int n, int m) {
        this->n = n;
        this->m = m;
        this->adj.assign(n, vector<int>());
    }

    void add_edge(int from, int to) {
        assert(from < n);
        assert(to < n);
        adj[from].push_back(to);
    }
};

struct SSSP_Result {
    bool has_negative_cycle;
    Graph* shortest_paths_tree;
};

SSSP_Result sssp(Graph* graph) {
    // TODO implement this (facile a dirsi, e' praticamente l'intera tesi)
}

void print_shortest_path_tree(SSSP_Result result) {
    vector<vector<int>> shortest_paths = result.shortest_paths_tree->adj;
    for (int i = 0; i < shortest_paths.size(); i++) {
        cout << i << ":\t";
        for (int nodo_adiacente : shortest_paths[i]) {
            cout << nodo_adiacente << "\t";
        }
        cout << endl;
    }
}

/* INPUT FORMAT:
* n m
* from_0 to_0
* from_1 to_1
* from_2 to_2
* ...
*/

int main() {
    // n: number of vertices
    // m: number of edges
    int n, m;
    cin >> n >> m;

    // Create graph with given constraints
    Graph* graph = new Graph(n, m);

    // Receive in input all edges
    for (int i = 0; i < m; i++) {
        int from, to; cin >> from >> to;
        graph->add_edge(from, to);
    }

    // END INPUT

    SSSP_Result sssp_result = sssp(graph);

    if (sssp_result.has_negative_cycle) {
        // TODO: Maybe print the negative cycle
        cout << "Graph constains a negative cycle" << endl;
    } else {
        print_shortest_path_tree(sssp_result);
    }
}


