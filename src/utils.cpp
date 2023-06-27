#include "utils.h"
#include "graph.h"
#include "sssp.h"

bool terminateLDD;
int INPUT_N;

void log(bool _cerr, int depth, string message) {
    if (_cerr) {
        cerr << "[DEBUG] ";
        for (int i=0; i<depth; ++i)
            cerr << '\t';
        cerr << message << endl;
    } else {
        cout << "[SSSP] ";
        for (int i=0; i<depth; ++i)
            cout << '\t';
        cout << message << endl;
    }
}

double d_min(double a, double b) {
    if (a < b)
        return a;
    return b;
}

set<int> getRandomVertices(Graph* g, int k) {
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

SSSP_Result dijkstra(Graph* g, int s) {
    log(false, 0, "Executing Dijkstra on graph:");
    printGraph(g, 0, false);

    /* DEBUG INPUT REQUIREMENTS */
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            assert(g->adj[i][j]>=0);
        }
    }
    /* END DEBUG INPUT REQUIREMENTS */

    SSSP_Result result;
    vector<vector<int>> result_adj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> result_is_edge(INPUT_N, vector<bool>(INPUT_N, false));
    vector<vector<int>> result_edges(INPUT_N, vector<int>());
    result.has_negative_cycle = false;

    vector<bool> confirmed(INPUT_N, false);
    confirmed[s] = true;

    // DEBUG
    vector<int> distances(INPUT_N);
    distances[s] = 0;

    priority_queue<pair<int, pair<int, int>>> pq; // <weight*-1, <indexFrom, indexTo>>
    for (int nodo:g->edges[s]) {
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
        result_edges[vertexFrom].push_back(vertexTo);

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
    result.shortest_paths_tree->edges = result_edges;

    return result;
}

Graph* induced_graph(Graph* g, set<int> vertices) {
    Graph* result = new Graph(vertices);
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            if (vertices.find(i)!=vertices.end() && vertices.find(j)!=vertices.end()) {
                result->adj[i][j] = g->adj[i][j];
                result->is_edge[i][j] = true;
                result->edges[i].push_back(j);
            }
        }
    }
    return result;
}

set<pair<int, int>> edgesUnion(set<pair<int, int>> a, set<pair<int, int>>& b) {
    for (auto i:b) {
        a.insert(i);
    }
    return a;
}

Graph* subtractVertices(Graph* g, set<int> vertices) {
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    vector<vector<int>> resultEdges(INPUT_N, vector<int>());
    set<int> resultV;
    for (int i : g->V) {
        if (vertices.find(i) != vertices.end())
            continue;
        resultV.insert(i);
        for (int j : g->edges[i]) {
            if (vertices.find(j)==vertices.end()) {
                resultAdj[i][j] = g->adj[i][j];
                resultIsEdge[i][j] = true;
                resultEdges[i].push_back(j);
            }
        }
    }
    Graph* result = new Graph(resultV);
    result->adj = resultAdj;
    result->is_edge = resultIsEdge;
    result->edges = resultEdges;
    delete g;
    return result;
}

Graph* subtractEdges(Graph* g, set<pair<int, int>> edges) {
    Graph* result = new Graph(g->V);
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            if (edges.find({i, j})==edges.end()) {
                result->adj[i][j] = g->adj[i][j];
                result->is_edge[i][j] = true;
                result->edges[i].push_back(j);
            }
        }
    }
    return result;
}

set<pair<int, int>> fromGraphToSetOfEdges(Graph* g) {
    set<pair<int, int>> result;
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            result.insert({i, j});
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

Graph* applyPriceFunction(Graph* g, PriceFunction p) {
    Graph* result = new Graph(g->V);
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            result->adj[i][j] = g->adj[i][j]+p.prices[i]-p.prices[j];
            result->is_edge[i][j] = true;
            result->edges[i].push_back(j);
        }
    }
    delete g;
    return result;
}

Graph* addIntegerToNegativeEdges(Graph* g, int e) {
    Graph* result = new Graph(g->V);
    result->adj = g->adj;
    result->is_edge = g->is_edge;
    result->edges = g->edges;
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            if (g->adj[i][j] < 0) {
                result->adj[i][j] = g->adj[i][j]+e;
            }
        }
    }
    return result;
}

Graph* mergeGraphs(Graph* g1, Graph* g2) {
    set<int> resultVertices;
    vector<vector<int>> resultAdj(INPUT_N, vector<int>(INPUT_N));
    vector<vector<bool>> resultIsEdge(INPUT_N, vector<bool>(INPUT_N, false));
    vector<vector<int>> resultEdges(INPUT_N, vector<int>());
    for (int i : g1->V) {
        resultVertices.insert(i);
        for (int j : g1->edges[i]) {
            resultAdj[i][j] = g1->adj[i][j];
            resultIsEdge[i][j] = true;
            resultEdges[i].push_back(j);
        }
    }
    for (int i : g2->V) {
        resultVertices.insert(i);
        for (int j : g2->edges[i]) {
            resultAdj[i][j] = g2->adj[i][j];
            resultIsEdge[i][j] = true;
            // This works because mergegraphs will always be called on graphs with disjoint vertex sets
            resultEdges[i].push_back(j);
        }
    }
    Graph* result = new Graph(resultVertices);
    result->adj = resultAdj;
    result->is_edge = resultIsEdge;
    result->edges = resultEdges;
    delete g1;
    delete g2;
    return result;
}

vector<set<int>> computeSCCs(Graph* g, int depth) {
    log(true, depth, "Computing SCCs...");
    string s_log = "Graph has " + to_string(g->V.size()) + " vertices";
    log(true, depth, s_log);
    vector<set<int>> result;
    vector<bool> visited(INPUT_N, false);
    stack<int> s;
    for (int v : g->V) {
        if (!visited[v])
            DFS(g, v, visited, s);
    }
    visited.assign(INPUT_N, false);
    Graph* gT = transpose(g);
    while (!s.empty()) {
        int v = s.top();
        s.pop();
        if (!visited[v]) {
            set<int> component;
            DFSaddComp(gT, v, visited, component);
            result.push_back(component);
        }
    }
    s_log = "SCCs computed. SCCs number: " + to_string(result.size());
    log(true, depth, s_log);
    delete gT;
    delete g;
    return result;
}

Graph* transpose(Graph* g) {
    Graph* result = new Graph(g->V);
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            result->adj[j][i] = g->adj[i][j];
            result->is_edge[j][i] = true;
            result->edges[j].push_back(i);
        }
    }

    return result;
}

void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s) {
    visited[v] = true;
    for (int i : g->edges[v]) {
        if (!visited[i]) {
            DFS(g, i, visited, s);
        }
    }
    s.push(v);
}

void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component) {
    visited[v] = true;
    component.insert(v);
    for (int i : g->edges[v]) {
        if (!visited[i]) {
            DFSaddComp(g, i, visited, component);
        }
    }
}

set<pair<int, int>> edgesMinusEdges(set<pair<int, int>>& a, set<pair<int, int>>& b) {
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

void printGraph(Graph* g, int depth, bool cerr) {
    log(cerr, depth, "PrintGraph:");
    log(cerr, depth, "Graph has " + to_string(g->V.size()) + " vertices:");
    string s = "";
    for (int i : g->V) {
        s += to_string(i) + " ";
    }
    log(cerr, depth, s);
    for (int i : g->V) {
        for (int j : g->edges[i]) {
            log(cerr, depth, to_string(i) + " -> " + to_string(j) + " : " + to_string(g->adj[i][j]));
        }
    }
}

bool checkConstantOutDegree(Graph* graph) {
    for (int v : graph->V) {
        if (graph->edges[v].size() > 2) {
            return false;
        }
    }
    return true;
}

Graph* addDummySource(Graph* g) {
    Graph* graph = new Graph(g->V);
    graph->adj = g->adj;
    graph->is_edge = g->is_edge;
    graph->edges = g->edges;

    graph->V.insert(0);
    for (int i : graph->V) {
        graph->adj[0][i] = 0;
        graph->is_edge[0][i] = true;
        graph->edges[0].push_back(i);
    }

    return graph;
}

Graph* onlyEdgesInsideSCCs(Graph* g, vector<set<int>> SCCs) {
    Graph* graph = new Graph(g->V);

    for (int i : graph->V) {
        for (int j : g->edges[i]) {
            for (set<int>& k : SCCs) {
                if (k.find(i) != k.end()) {
                    if (k.find(j) == k.end())
                        break;
                    else {
                        graph->adj[i][j] = g->adj[i][j];
                        graph->is_edge[i][j] = true;
                        graph->edges[i].push_back(j);
                    }
                }
            }
        }
    }

    return graph;
}

void print_shortest_path_tree(SSSP_Result result) {
    log(false, 0, "Shortest paths tree:");
    vector<vector<int>> adj = result.shortest_paths_tree->adj;
    vector<vector<int>> edges = result.shortest_paths_tree->edges;
    for (int i = 0; i < adj.size(); i++) {
        string s = "Dal nodo " + to_string(i) + ":";
        log(false, 0, s);
        for (int j : edges[i]) {
            string s = "al nodo " + to_string(j) + ": " + to_string(adj[i][j]);
            log(false, 0, s);
        }
        log(false, 0, "");
    }
}
