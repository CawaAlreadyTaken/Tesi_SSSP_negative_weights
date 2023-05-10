#include "utils.h"
#include "graph.h"
#include "sssp.h"

bool terminateLDD;

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
            auto it = g->V.find(randomVert);
            ++it;
            if (it != g->V.end())
                vertices.insert(*it);
            else
                vertices.insert(*(--(--it)));
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
    //cout << "[PARTIAL RESULT] Dijkstra on graph:" << endl;
    //csacademy_printGraph(g);

    /* DEBUG INPUT REQUIREMENTS */
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j])
                assert(g->adj[i][j]>=0);
        }
    }
    /* END DEBUG INPUT REQUIREMENTS */

    SSSP_Result result;
    vector<vector<int>> result_adj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> result_is_edge(INPUT_N, vector<bool>(INPUT_N, false));
    result.has_negative_cycle = false;

    vector<bool> confirmed(INPUT_N, false);
    confirmed[s] = true;

    // DEBUG
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
        distances[vertexTo] = weight;
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
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j] && vertices.find(i)!=vertices.end() && vertices.find(j)!=vertices.end()) {
                result->adj[i][j] = g->adj[i][j];
                result->is_edge[i][j] = true;
            }
        }
    }
    return result;
}

set<pair<int, int>> edgesUnion(set<pair<int, int>> a, set<pair<int, int>> b) {
    for (auto i:b) {
        a.insert(i);
    }
    return a;
}

Graph* subtractVertices(Graph* g, set<int> vertices, int INPUT_N) {
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    set<int> resultV;
    for (int i : g->V) {
        if (vertices.find(i) != vertices.end())
            continue;
        resultV.insert(i);
        for (int j : g->V) {
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
    Graph* result = new Graph(g->V, INPUT_N);
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j] && edges.find({i, j})==edges.end()) {
                result->adj[i][j] = g->adj[i][j];
                result->is_edge[i][j] = true;
            }
        }
    }
    return result;
}

set<pair<int, int>> fromMatrixToSet(vector<vector<bool>> isEdge) {
    set<pair<int, int>> result;
    for (int i = 0; i < isEdge.size(); i++) {
        for (int j = 0; j < isEdge.size(); j++) {
            if (isEdge[i][j]) {
                result.insert({i, j});
            }
        }
    }
    return result;
}

bool isSubset(set<int> a, set<int> b) {
    for (int i:a) {
        if (b.find(i) == b.end())
            return false;
    }
    return true;
}

Graph* applyPriceFunction(Graph* g, PriceFunction p, int INPUT_N) {
    Graph* result = new Graph(g->V, INPUT_N);
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j]) {
                result->adj[i][j] = g->adj[i][j]+p.prices[i]-p.prices[j];
                result->is_edge[i][j] = true;
            }
        }
    }
    return result;
}

Graph* addIntegerToNegativeEdges(Graph* g, int e, int INPUT_N) {
    Graph* result = new Graph(g->V, INPUT_N);
    result->adj = g->adj;
    result->is_edge = g->is_edge;
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j] && g->adj[i][j] < 0) {
                result->adj[i][j] = g->adj[i][j]+e;
            }
        }
    }
    return result;
}

Graph* mergeGraphs(Graph* g1, Graph* g2, int INPUT_N) {
    set<int> resultVertices;
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    for (int i : g1->V) {
        resultVertices.insert(i);
        for (int j : g1->V) {
            if (g1->is_edge[i][j]) {
                resultAdj[i][j] = g1->adj[i][j];
                resultIsEdge[i][j] = true;
            }
        }
    }
    for (int i : g2->V) {
        resultVertices.insert(i);
        for (int j : g2->V) {
            if (g2->is_edge[i][j]) {
                resultAdj[i][j] = g2->adj[i][j];
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
    cerr << "[DEBUG] Computing SCCs..." << endl;
    cerr << "[DEBUG] Graph has " << g->V.size() << " vertices" << endl;
    vector<set<int>> result;
    vector<bool> visited(INPUT_N, false);
    stack<int> s;
    for (int v : g->V) {
        if (!visited[v])
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
    cerr << "[DEBUG] SCCs computed. SCCs number: " << result.size() << endl;
    return result;
}

Graph* transpose(Graph* g, int INPUT_N) {
    Graph* result = new Graph(g->V, INPUT_N);
    for (int i : g->V) {
        for (int j : g->V) {
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
    for (int i : g->V) {
        if (g->is_edge[v][i] && !visited[i]) {
            DFS(g, i, visited, s);
        }
    }
    s.push(v);
}

void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component) {
    visited[v] = true;
    component.insert(v);
    for (int i : g->V) {
        if (g->is_edge[v][i] && !visited[i]) {
            DFSaddComp(g, i, visited, component);
        }
    }
}

set<pair<int, int>> edgesMinusEdges(set<pair<int, int>> a, set<pair<int, int>> b) {
    set<pair<int, int>> result;
    for (auto i:a) {
        if (b.find(i) == b.end())
            result.insert(i);
    }
    return result;
}

vector<set<int>> topologicalSort(vector<set<int>> SCCs, Graph* graph, vector<int> fromVertixToSCC) {

    vector<vector<int>> sccAdj(SCCs.size(), vector<int>());
    for (int v : graph->V) {
        for (int w : graph->V) {
            if (graph->is_edge[v][w] && fromVertixToSCC[v] != fromVertixToSCC[w]) {
                sccAdj[fromVertixToSCC[v]].push_back(fromVertixToSCC[w]);
            }
        }
    }

    stack<int> s;
    vector<set<int>> result;
    vector<bool> visited(SCCs.size(), false);

    for (int i = 0; i < SCCs.size(); i++) {
        if (!visited[i])
            topoDFS(i, visited, s, sccAdj);
    }

    while (!s.empty()) {
        result.push_back(SCCs[s.top()]);
        s.pop();
    }
    return result;
}

void topoDFS(int index, vector<bool>& visited, stack<int>& s, vector<vector<int>>& sccAdj) {
    visited[index] = true;
    
    for (int k : sccAdj[index]) {
        if (!visited[k])
            topoDFS(k, visited, s, sccAdj);
    }

    s.push(index);
}

void printGraph(Graph* g) {
    cerr << "[DEBUG] PrintGraph:" << endl;
    cerr << "[DEBUG] Graph has " << g->V.size() << " vertices:" << endl;
    cerr << "[DEBUG] ";
    for (int i : g->V) {
        cerr << i << ' ';
    }
    cerr << endl;
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j]) {
                cerr << "[DEBUG] " << i << " -> " << j << " : " << g->adj[i][j] << endl;
            }
        }
    }
    cerr << endl << endl;
}

void csacademy_printGraph(Graph* g) {
    for (int i : g->V) {
        cout << i << endl;
    }
    for (int i : g->V) {
        for (int j : g->V) {
            if (g->is_edge[i][j])
                cout << i << ' ' << j << ' ' << g->adj[i][j] << endl;
        }
    }
}

bool checkConstantOutDegree(Graph* graph) {
    int shared_out_degree = -1;
    for (int v : graph->V) {
        int out_degree = 0;
        for (int u : graph->V) {
            if (graph->is_edge[v][u]) {
                out_degree++;
            }
        }
        if (shared_out_degree == -1)
            shared_out_degree = out_degree;
        else if (shared_out_degree != out_degree)
            return false;
    }
    return true;
}

Graph * addDummySource(Graph* g, int INPUT_N) {
    Graph* graph = new Graph(g->V, INPUT_N);
    graph->adj = g->adj;
    graph->is_edge = g->is_edge;

    graph->V.insert(0);
    for (int i : graph->V) {
        graph->adj[0][i] = 0;
        graph->is_edge[0][i] = true;
    }

    return graph;
}

void print_shortest_path_tree(SSSP_Result result) {
    vector<vector<int>> shortest_paths = result.shortest_paths_tree->adj;
    for (int i = 0; i < shortest_paths.size(); i++) {
        cout << "[RESULT] Dal nodo " << i << ":" << endl;
        for (int j = 0; j < shortest_paths[i].size(); j++) {
            if (result.shortest_paths_tree->is_edge[i][j])
                cout << "[RESULT] al nodo " << j << ": " << result.shortest_paths_tree->adj[i][j] << endl;
        }
        cout << endl;
    }
}
