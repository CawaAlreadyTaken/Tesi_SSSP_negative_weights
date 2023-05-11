#include "graph.h"
#include "sssp.h"
#include "utils.h"

#include <iostream>

using namespace std;

/* INPUT FORMAT:
* n m s
* from_0 to_0 weight_0
* from_1 to_1 weight_1
* from_2 to_2 weight_2
* ... (m lines)
*/

int main() {

    /* INPUT */
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
    Graph* graph = new Graph(vertices);

    // Receive in input all edges
    for (int i = 0; i < m; i++) {
        int from, to, weight; cin >> from >> to >> weight;
        graph->add_edge(from, to, weight);
    }
    /* END INPUT */

    log(false, 0, "Input graph:");
    printGraph(graph, 0, false);

    SSSP_Result sssp_result = SPmain(graph, s);

    if (sssp_result.has_negative_cycle) {
        // TODO: Maybe print the negative cycle
        log(false, 0, "[RESULT] Graph constains a negative cycle");
    } else {
        print_shortest_path_tree(sssp_result);
    }
}
