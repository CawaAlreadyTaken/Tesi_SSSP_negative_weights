#include "utils.h"
#include "graph.h"
#include "sssp.h"

double d_min(double a, double b) {
    if (a < b)
        return a;
    return b;
}

set<int> getRandomVertices(Graph* g, int k, int INPUT_N) {
    set<int> vertices;
    while (k--) {
        int randomVert = rand() % INPUT_N;
        if (g->V.find(randomVert) != g->V.end()) {
            vertices.insert(randomVert);
        } else {
            g->V.insert(randomVert);
            if (g->V.find(randomVert++) != g->V.end())
                vertices.insert(*(g->V.find(randomVert++)));
            else
                vertices.insert(*(g->V.find(randomVert--)));
            g->V.erase(randomVert);
        }
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

SSSP_Result dijkstra(Graph* g, int s, int INPUT_N) {
    SSSP_Result result;
    vector<vector<int>> result_adj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> result_is_edge(INPUT_N, vector<bool>(INPUT_N, false));
    // If dijkstra is called, then there is no negative weight
    // If this is wrong, exit
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j])
                assert(g->adj[i][j]>=0);
        }
    }
    result.has_negative_cycle = false;

    vector<bool> confirmed(INPUT_N, false);
    confirmed[s] = true;

    // TODO maybe remove this?
    vector<int> distances(INPUT_N);
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
    result.shortest_paths_tree = new Graph(g->V, INPUT_N);
    result.shortest_paths_tree->adj = result_adj;
    result.shortest_paths_tree->is_edge = result_is_edge;

    return result;
}

Graph* induced_graph(Graph* g, set<int> vertices, int INPUT_N) {
    Graph* result = new Graph(vertices, INPUT_N);
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j] && vertices.find(i)!=vertices.end() && vertices.find(j)!=vertices.end()) {
                result->adj[i][j] = g->adj[i][j];
                result->is_edge[i][j] = true;
            }
        }
    }
    return result;
}

set<pair<int, int>> intersect(set<pair<int, int>> a, set<pair<int, int>> b) {
    set<pair<int, int>> result;
    for (auto i:a) {
        if (b.find(i) != b.end())
            result.insert(i);
    }
    return result;
}

Graph* subtractVertices(Graph* g, set<int> vertices, int INPUT_N) {
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    set<int> resultV;
    for (int i = 0; i < g->V.size(); i++) {
        if (vertices.find(i) != vertices.end())
            continue;
        resultV.insert(i);
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j] && vertices.find(j)==vertices.end()) {
                resultAdj[i][j] = g->adj[i][j];
                resultIsEdge[i][j] = true;
            }
        }
    }
    Graph* result = new Graph(resultV, INPUT_N);
    result->adj = resultAdj;
    result->is_edge = resultIsEdge;
    return result;
}

Graph* subtractEdges(Graph* g, set<pair<int, int>> edges, int INPUT_N) {
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    set<int> resultV;
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j] && edges.find({i, j})==edges.end()) {
                resultV.insert(i);
                resultV.insert(j);
                resultAdj[i][j] = g->adj[i][j];
                resultIsEdge[i][j] = true;
            }
        }
    }
    Graph* result = new Graph(resultV, INPUT_N);
    result->adj = resultAdj;
    result->is_edge = resultIsEdge;
    return result;
}

set<pair<int, int>> fromMatrixToSet(vector<vector<bool>> isEdge) {
    set<pair<int, int>> result;
    for (int i = 0; i < isEdge.size(); i++) {
        for (int j = 0; j < isEdge.size(); j++) {
            if (isEdge[i][j])
                result.insert({i, j});
        }
    }
    return result;
}

bool isSubset(set<int> a, set<int> b) {
    for (auto i:a) {
        if (b.find(i) == b.end())
            return false;
    }
    return true;
}

Graph* applyPriceFunction(Graph* g, PriceFunction p, int INPUT_N) {
    Graph* result = new Graph(g->V, INPUT_N);
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j]) {
                result->adj[i][j] = g->adj[i][j]+p.prices[i]-p.prices[j];
                result->is_edge[i][j] = true;
            }
        }
    }
    return result;
}

Graph* addIntegerToEdges(Graph* g, int e, int INPUT_N) {
    Graph* result = new Graph(g->V, INPUT_N);
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j]) {
                result->adj[i][j] = g->adj[i][j]+e;
                result->is_edge[i][j] = true;
            }
        }
    }
    return result;
}

Graph* mergeGraphs(Graph* g1, Graph* g2, int INPUT_N) {
    set<int> resultVertices;
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    for (int i = 0; i < g1->V.size(); i++) {
        resultVertices.insert(i);
        for (int j = 0; j < g1->V.size(); j++) {
            if (g1->is_edge[i][j]) {
                resultAdj[i][j] = g1->adj[i][j];
                resultIsEdge[i][j] = true;
            }
        }
    }
    for (int i = 0; i < g2->V.size(); i++) {
        resultVertices.insert(i);
        for (int j = 0; j < g2->V.size(); j++) {
            if (g2->is_edge[i][j]) {
                resultAdj[i][j] = g1->adj[i][j];
                resultIsEdge[i][j] = true;
            }
        }
    }
    Graph* result = new Graph(resultVertices, INPUT_N);
    result->adj = resultAdj;
    result->is_edge = resultIsEdge;
    return result;
}

vector<set<int>> computeSCCs(Graph* g, int INPUT_N) {
    vector<set<int>> result;
    vector<bool> visited(INPUT_N, false);
    stack<int> s;
    for (int v : g->V) {
        DFS(g, v, visited, s);
    }
    visited.assign(INPUT_N, false);
    Graph* gT = transpose(g, INPUT_N);
    while (!s.empty()) {
        int v = s.top();
        s.pop();
        if (!visited[v]) {
            set<int> component;
            DFSaddComp(gT, v, visited, component);
            result.push_back(component);
        }
    }
    return result;
}

Graph* transpose(Graph* g, int INPUT_N) {
    Graph* result = new Graph(g->V, INPUT_N);
    for (int i = 0; i < g->V.size(); i++) {
        for (int j = 0; j < g->V.size(); j++) {
            if (g->is_edge[i][j]) {
                result->adj[j][i] = g->adj[i][j];
                result->is_edge[j][i] = true;
            }
        }
    }

    return result;
}

void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s) {
    visited[v] = true;
    for (int i = 0; i < g->V.size(); i++) {
        if (g->is_edge[v][i] && !visited[i]) {
            DFS(g, i, visited, s);
        }
    }
    s.push(v);
}

void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component) {
    visited[v] = true;
    component.insert(v);
    for (int i = 0; i < g->V.size(); i++) {
        if (g->is_edge[v][i] && !visited[i]) {
            DFSaddComp(g, i, visited, component);
        }
    }
}
