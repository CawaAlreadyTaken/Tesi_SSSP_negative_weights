#include <bits/stdc++.h> // TODO: use more specific libraries
#include "graph.h"
#include "sssp.h"

using namespace std;

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
* n m s
* from_0 to_0 weight_0
* from_1 to_1 weight_1
* from_2 to_2 weight_2
* ...
*/

int main() {
    // n: number of vertices
    // m: number of edges
    // s: index of source
    int n, m, s;
    cin >> n >> m >> s;

    // Create graph with given constraints
    Graph* graph = new Graph(n, m);

    // Receive in input all edges
    for (int i = 0; i < m; i++) {
        int from, to, weight; cin >> from >> to >> weight;
        graph->add_edge(from, to, weight);
    }

    // END INPUT

    SSSP_Result sssp_result = SPmain(graph, s);

    if (sssp_result.has_negative_cycle) {
        // TODO: Maybe print the negative cycle
        cout << "Graph constains a negative cycle" << endl;
    } else {
        print_shortest_path_tree(sssp_result);
    }
}
