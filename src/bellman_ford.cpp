#include <iostream>
#include <vector>
#include <set>
#include <cstdint>

using namespace std;

vector<int> bellman_ford(int n, set<pair<int, pair<int, int>>> & edges, int s) {
    // Initialize distances
    vector<int> distances(n, INT32_MAX);
    distances[s] = 0;

    // Relax edges
    for (int i = 0; i < n; i++) {
        for (auto [weight, ft] : edges) {
            auto [from, to] = ft;
            if (distances[from] != INT32_MAX && distances[to] + weight < distances[to])
                distances[to] = distances[from] + weight;
        }
    }

    // Return distances
    return distances;
}

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
    n++;
    set<pair<int, pair<int, int>>> edges; // {weight, {from, to}}

    // Receive in input all edges
    for (int i = 0; i < m; i++) {
        int from, to, weight; cin >> from >> to >> weight;
        edges.insert({weight, {from, to}});
    }
    /* END INPUT */

    // Run Bellman-Ford algorithm
    vector<int> distances = bellman_ford(n, edges, s);
}
