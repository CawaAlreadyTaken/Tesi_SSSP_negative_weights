#include "graph.h"
#include "sssp.h"
#include "utils.h"

#include <iostream>

using namespace std;

int INPUT_N;

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
    INPUT_N = n+1;

    // Create graph with given constraints
    set<int> vertices;
    for (int i = 1; i <= n; i++)
        vertices.insert(i);
    Graph* graph = new Graph(vertices, INPUT_N);

    // Receive in input all edges
    for (int i = 0; i < m; i++) {
        int from, to, weight; cin >> from >> to >> weight;
        graph->add_edge(from, to, weight);
    }

    // END INPUT

    // DEBUG PRINT
    csacademy_printGraph(graph);

    SSSP_Result sssp_result = SPmain(graph, s, INPUT_N);

    if (sssp_result.has_negative_cycle) {
        // TODO: Maybe print the negative cycle
        cout << "Graph constains a negative cycle" << endl;
    } else {
        print_shortest_path_tree(sssp_result);
    }
}
