#include "utils.h"

int roundB(int b) {
    // round up b to the nearest power of 2
    if (b==1)
        return b;
    int k = 1;
    while (true) {
        k*=2;
        if (k >= b)
            break;
    }
    return k;
}

SSSP_Result dijkstra(Graph* g, int s) {
    SSSP_Result result;
    vector<vector<int>> result_adj(g->n);
    vector<vector<int>> result_weights(g->n);
    int result_edges_number = 0;
    // If dijkstra is called, then there is no negative weight
    // If this is wrong, exit
    for (int i = 0; i < g->n; i++) {
        for (int w : g->weights[i])
            assert(w>=0);
    }
    result.has_negative_cycle = false;

    vector<bool> confirmed(g->n, false);
    confirmed[s] = true;

    vector<int> distances(g->n);
    distances[s] = 0;

    priority_queue<pair<int, pair<int, int>>> pq; // <weight*-1, <indexFrom, indexTo>>
    for (int i = 0; i < g->adj[s].size(); i++) {
        pq.push({g->weights[s][i]*-1, {s, g->adj[s][i]}});
    }

    while (!pq.empty()) {
        auto top = pq.top();
        pq.pop();
        int weight = top.first*-1;
        int indexFrom = top.second.first;
        int indexTo = top.second.second;

        // If already confirmed, skip
        if (confirmed[indexTo])
            continue;

        // Otherwise, we found a new distance to confirm
        distances[indexTo] = distances[indexFrom]+weight;
        confirmed[indexTo] = true;
        result_adj[indexFrom].push_back(indexTo);
        result_weights[indexFrom].push_back(weight);
        result_edges_number++;

        // Add all his adj not yet confirmed to the priority queue
        for (int i = 0; i < g->adj[indexTo].size(); i++) {
            if (!confirmed[g->adj[indexTo][i]]) {
                pair<int, pair<int, int>> newNode = {g->weights[indexTo][i]*-1, {indexTo, g->adj[indexTo][i]}};
                pq.push(newNode);
            }
        }
    }
    
    // Compose result tree
    result.shortest_paths_tree = new Graph(g->n, result_edges_number);
    result.shortest_paths_tree->adj = result_adj;
    result.shortest_paths_tree->weights = result_weights;

    return result;
}
