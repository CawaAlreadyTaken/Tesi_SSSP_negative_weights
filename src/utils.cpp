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

long long l_max(long long a, long long b) {
    return a > b;
}

double d_min(double a, double b) {
    return a < b;
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

long long roundB(long long b) {
    // round up b to the nearest power of 2
    if (b==1)
        return b;
    long long k = 1;
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
        for (int j = 0; j < g->edges[i].size(); j++) {
            assert(g->adj[i][j]>=0);
        }
    }
    /* END DEBUG INPUT REQUIREMENTS */

    SSSP_Result result;
    vector<vector<long long>> result_adj(INPUT_N, vector<long long>());
    vector<vector<int>> result_edges(INPUT_N, vector<int>());
    result.has_negative_cycle = false;

    vector<bool> confirmed(INPUT_N, false);
    confirmed[s] = true;

    priority_queue<pair<long long, pair<int, int>>> pq; // <weight*-1, <indexFrom, indexTo>>
    for (int i = 0; i < g->edges[s].size(); i++) {
        int nodo = g->edges[s][i];
        pq.push({g->adj[s][i]*-1, {s, nodo}});
    }

    while (!pq.empty()) {
        auto top = pq.top();
        pq.pop();
        long long weight = top.first*-1;
        int vertexFrom = top.second.first;
        int vertexTo = top.second.second;

        // If already confirmed, skip
        if (confirmed[vertexTo])
            continue;

        // Otherwise, we found a new distance to confirm
        confirmed[vertexTo] = true;
        result_adj[vertexFrom].push_back(weight);
        result_edges[vertexFrom].push_back(vertexTo);

        // Add all his adj not yet confirmed to the priority queue
        for (int i = 0; i < g->edges[vertexTo].size(); i++) {
            int nodo = g->edges[vertexTo][i];
            if (!confirmed[nodo]) {
                pair<long long, pair<int, int>> newNode = {(weight+g->adj[vertexTo][i])*-1, {vertexTo, nodo}};
                pq.push(newNode);
            }
        }
    }

    // Compose result tree
    result.shortest_paths_tree = new Graph(g->V);
    result.shortest_paths_tree->adj = result_adj;
    result.shortest_paths_tree->edges = result_edges;

    return result;
}

Graph* induced_graph(Graph* g, set<int>& vertices) {
    Graph* result = new Graph(vertices);
    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            if (vertices.find(i)!=vertices.end() && vertices.find(g->edges[i][j])!=vertices.end()) {
                result->adj[i].push_back(g->adj[i][j]);
                result->edges[i].push_back(g->edges[i][j]);
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

Graph* subtractVertices(Graph* g, set<int>& vertices) {
    vector<vector<long long>> resultAdj(INPUT_N, vector<long long>());
    vector<vector<int>> resultEdges(INPUT_N, vector<int>());
    set<int> resultV;
    for (int i : g->V) {
        if (vertices.find(i) != vertices.end())
            continue;
        resultV.insert(i);
        for (int j = 0; j < g->edges[i].size(); j++) {
            if (vertices.find(g->edges[i][j])==vertices.end()) {
                resultAdj[i].push_back(g->adj[i][j]);
                resultEdges[i].push_back(g->edges[i][j]);
            }
        }
    }
    Graph* result = new Graph(resultV);
    result->adj = resultAdj;
    result->edges = resultEdges;
    delete g;
    return result;
}

Graph* subtractEdges(Graph* g, set<pair<int, int>>& edges) {
    Graph* result = new Graph(g->V);
    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            if (edges.find({i, g->edges[i][j]})==edges.end()) {
                result->adj[i].push_back(g->adj[i][j]);
                result->edges[i].push_back(g->edges[i][j]);
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

bool isSubset(set<int>& a, set<int>& b) {
    for (int i:a) {
        if (b.find(i) == b.end())
            return false;
    }
    return true;
}

Graph* applyPriceFunction(Graph* g, PriceFunction p) {
    Graph* result = new Graph(g->V);
    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            result->adj[i].push_back(g->adj[i][j]+p.prices[i]-p.prices[g->edges[i][j]]);
            result->edges[i].push_back(g->edges[i][j]);
        }
    }
    delete g;
    return result;
}

Graph* addIntegerToNegativeEdges(Graph* g, long long e) {
    Graph* result = new Graph(g->V);
    result->adj = g->adj;
    result->edges = g->edges;
    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            if (g->adj[i][j] < 0) {
                result->adj[i].push_back(g->adj[i][j]+e);
                result->edges[i].push_back(g->edges[i][j]);
            }
        }
    }
    return result;
}

vector<set<int>> computeSCCs(Graph* g, int depth) {
    log(true, depth, "Computing SCCs...");
    string s_log = "Graph has " + to_string(g->V.size()) + " vertices";
    log(true, depth, s_log);
    vector<set<int>> result;
    vector<bool> visited(INPUT_N, false);
    stack<int> s;
    /*
    stack<int> s1;
    for (int v : g->V) {
        if (!visited[v]) {
            s1.push(v);
            while (!s1.empty()) {
                int current = s1.top();
                s1.pop();
                if (current < 0) {
                    s.push(-current);
                    continue;
                }
                s1.push(-v);
                if (!visited[current]) {
                    visited[current] = true;
                    for (int i : g->edges[current]) {
                        if (!visited[i]) {
                            s1.push(i);
                        }
                    }
                }
            }
        }
    }
    */
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
    cout << "f" << endl;
    s_log = "SCCs computed. SCCs number: " + to_string(result.size());
    log(true, depth, s_log);
    delete gT;
    delete g;
    return result;
}

Graph* transpose(Graph* g) {
    Graph* result = new Graph(g->V);
    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            result->edges[g->edges[i][j]].push_back(i);
        }
    }

    return result;
}

void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s) {
    visited[v] = true;
    for (int i : g->edges[v]) {
        cout << i << endl;
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

vector<set<int>> topologicalSort(vector<set<int>>& SCCs, Graph* graph, vector<int>& fromVertixToSCC) {

    vector<vector<int>> sccAdj(SCCs.size(), vector<int>());
    for (int v : graph->V) {
        for (int w : graph->edges[v]) {
            if (fromVertixToSCC[v] != fromVertixToSCC[w]) {
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
        for (int j = 0; j < g->edges[i].size(); j++) {
            log(cerr, depth, to_string(i) + " -> " + to_string(g->edges[i][j]) + " : " + to_string(g->adj[i][j]));
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
    graph->edges = g->edges;

    graph->V.insert(0);
    for (int i : graph->V) {
        if (i != 0) {
            graph->adj[0].push_back(0);
            graph->edges[0].push_back(i);
        }
    }

    return graph;
}

Graph* onlyEdgesInsideSCCs(Graph* g, vector<set<int>>& SCCs) {
    Graph* graph = new Graph(g->V);

    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            for (set<int>& k : SCCs) {
                if (k.find(i) != k.end()) {
                    if (k.find(g->edges[i][j]) == k.end())
                        break;
                    else {
                        graph->adj[i].push_back(g->adj[i][j]);
                        graph->edges[i].push_back(g->edges[i][j]);
                    }
                }
            }
        }
    }

    return graph;
}

map<pair<int, int>, long long> createSupportMap(Graph* g) {
    map<pair<int, int>, long long> supportMap;
    for (int i : g->V) {
        for (int j = 0; j < g->edges[i].size(); j++) {
            supportMap[make_pair(i, g->edges[i][j])] = g->adj[i][j];
        }
    }
    return supportMap;
}

void print_shortest_path_tree(SSSP_Result result) {
    log(false, 0, "Shortest paths tree:");
    vector<vector<long long>> adj = result.shortest_paths_tree->adj;
    vector<vector<int>> edges = result.shortest_paths_tree->edges;
    for (int i = 0; i < adj.size(); i++) {
        if (adj[i].size() == 0)
            continue;
        string s = "Dal nodo " + to_string(i) + ":";
        log(false, 0, s);
        for (int j = 0; j < edges[i].size(); j++) {
            string s = "al nodo " + to_string(edges[i][j]) + ": " + to_string(adj[i][j]);
            log(false, 0, s);
        }
        log(false, 0, "");
    }
}
