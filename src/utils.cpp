#include "utils.h"
#include "graph.h"

double d_min(double a, double b) {
    if (a < b)
        return a;
    return b;
}

set<int> getRandomVertices(Graph* g, int k) {
    set<int> vertices;
    while (k--) {
        int randomIndex = rand() % g->V.size();
        vertices.insert(g->V[randomIndex]);
    }
    return vertices;
}

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
    vector<vector<int>> result_adj(MAX_N, vector<int>(MAX_N));
    vector<vector<bool>> result_is_edge(MAX_N, vector<bool>(MAX_N, false));
    // If dijkstra is called, then there is no negative weight
    // If this is wrong, exit
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j])
                assert(g->adj[i][j]>=0);
        }
    }
    result.has_negative_cycle = false;

    vector<bool> confirmed(MAX_N, false);
    confirmed[s] = true;

    // TODO maybe remove this?
    vector<int> distances(MAX_N);
    distances[s] = 0;

    priority_queue<pair<int, pair<int, int>>> pq; // <weight*-1, <indexFrom, indexTo>>
    for (int nodo:g->V) {
        if (g->is_edge[s][nodo])
            pq.push({g->adj[s][nodo]*-1, {s, nodo}});
    }

    while (!pq.empty()) {
        auto top = pq.top();
        pq.pop();
        int weight = top.first*-1;
        int vertexFrom = top.second.first;
        int vertexTo = top.second.second;

        // If already confirmed, skip
        if (confirmed[vertexTo])
            continue;

        // Otherwise, we found a new distance to confirm
        distances[vertexTo] = distances[vertexFrom]+weight;
        confirmed[vertexTo] = true;
        result_adj[vertexFrom][vertexTo] = weight;
        result_is_edge[vertexFrom][vertexTo] = true;

        // Add all his adj not yet confirmed to the priority queue
        for (int nodo:g->V) {
            if (g->is_edge[vertexTo][nodo] && !confirmed[nodo]) {
                pair<int, pair<int, int>> newNode = {(distances[vertexTo]+g->adj[vertexTo][nodo])*-1, {vertexTo, nodo}};
                pq.push(newNode);
            }
        }
    }
    
    // Compose result tree
    result.shortest_paths_tree = new Graph(g->V);
    result.shortest_paths_tree->adj = result_adj;
    result.shortest_paths_tree->is_edge = result_is_edge;

    return result;
}
