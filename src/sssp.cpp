#include "sssp.h"
#include "graph.h"
#include "utils.h"
#include <cstdint>
#include <random>

pair<set<int>, set<pair<int, int>>> ballIn(Graph* graph, int v, int R, int INPUT_N) {
    // In order to do ballIn, we can reverse the edges and call ballOut
    Graph* graph_copy = new Graph(graph->V, INPUT_N);
    for (int i : graph->V) {
        for (int j : graph->V) {
            if (graph->is_edge[i][j])
                graph_copy->add_edge(j, i, graph->adj[i][j]);
        }
    }
    return ballOut(graph_copy, v, R, INPUT_N, true);
}

pair<set<int>, set<pair<int, int>>> ballOut(Graph* graph, int v, int R, int INPUT_N, bool fromBallIn) {
    // basically dijkstra in order to find vertices closer than R
    set<int> ris;
    set<pair<int, int>> boundary;

    vector<bool> confirmed(INPUT_N, false);
    confirmed[v] = true;
    ris.insert(v);

    priority_queue<pair<int, pair<int, int>>> pq; // <weight*-1, <indexFrom, indexTo>>
    for (int nodo : graph->V) {
        if (graph->is_edge[v][nodo])
            pq.push({graph->adj[v][nodo]*-1, {v, nodo}});
    }

    while (!pq.empty()) {
        auto top = pq.top();
        pq.pop();
        int weight = top.first*-1;
        int vertexFrom = top.second.first;
        int vertexTo = top.second.second;

        if (weight > R) { // Difference from dijstra
            if (fromBallIn) {
                boundary.insert({vertexTo, vertexFrom});
            } else {
                boundary.insert({vertexFrom, vertexTo});
            }
            continue;
        }

        if (confirmed[vertexTo])
            continue;

        confirmed[vertexTo] = true;
        ris.insert(vertexTo);

        for (int nodo : graph->V) {
            if (graph->is_edge[vertexTo][nodo] && !confirmed[nodo]) {
                pair<int, pair<int, int>> newNode = {weight*-1, {vertexTo, nodo}};
                pq.push(newNode);
            }
        }
    }

    return {ris, boundary};
}

set<pair<int, int>> LDD(Graph* graph, int D, int INPUT_N) {
    // Save a copy of the original graph, since we are going to modify it
    Graph* g0 = new Graph(graph->V, INPUT_N);
    g0->adj = graph->adj;
    g0->is_edge = graph->is_edge;

    set<pair<int, int>> Erem;
    // Phase 1: mark vertices as light or heavy
    float k = 27*log(INPUT_N);  // TODO: change this
    set<int> S = getRandomVertices(graph, k, INPUT_N);

    map<int, set<int>> ballInIntersec;
    map<int, set<int>> ballOutIntersec;

    // Da ricontrollare
    for (int x : S) {
        set<int> bIn = ballIn(graph, x, D/4, INPUT_N).first;
        set<int> bOut = ballOut(graph, x, D/4, INPUT_N).first;
        for (int z : bIn)
            ballOutIntersec[z].insert(x);
        for (int z : bOut)
            ballInIntersec[z].insert(x);
    }

    vector<int> in_light;
    vector<int> out_light;
    vector<int> heavy;

    for (int v : graph->V) {
        if (ballInIntersec[v].size() <= 0.6*k)
            in_light.push_back(v);
        else if (ballOutIntersec[v].size() <= 0.6*k)
            out_light.push_back(v);
        else
            heavy.push_back(v);
    }
    // Phase 2: Carve out balls until no light vertices remain
    default_random_engine generator;
    while (in_light.size()) {
        int v = in_light.back(); // Check unmarkation
        double p = d_min(1, 80.0*log2(INPUT_N)/D);
        geometric_distribution<int> distribution(p);
        int Rv = distribution(generator);
        auto result = ballIn(graph, v, Rv, INPUT_N);
        set<int> newBallIn = result.first;
        set<pair<int, int>> Ebound = result.second;
        if (Rv > D/4 || newBallIn.size() > 0.7*graph->V.size()) {
            //return Erem = E(G) and terminate
            return fromMatrixToSet(graph->is_edge);
        }
        set<pair<int, int>> Erecurs = LDD(induced_graph(graph, newBallIn, INPUT_N), D, INPUT_N);
        Erem = intersect(intersect(Erem, Ebound), Erecurs);
        graph = subtractVertices(graph, newBallIn, INPUT_N);
    }

    while (out_light.size()) {
        int v = out_light.back(); // Check unmarkation
        double p = d_min(1, 80.0*log2(INPUT_N)/D);
        geometric_distribution<int> distribution(p);
        int Rv = distribution(generator);
        auto result = ballOut(graph, v, Rv, INPUT_N);
        set<int> newBallOut = result.first;
        set<pair<int, int>> Ebound = result.second;
        if (Rv > D/4 || newBallOut.size() > 0.7*graph->V.size()) {
            //return Erem = E(G) and terminate TODO maybe terminate means really terminate
            return fromMatrixToSet(graph->is_edge);
        }
        set<pair<int, int>> Erecurs = LDD(induced_graph(graph, newBallOut, INPUT_N), D, INPUT_N);
        Erem = intersect(intersect(Erem, Ebound), Erecurs);
        graph = subtractVertices(graph, newBallOut, INPUT_N);
    }

    // Clean Up: check that remaining vertices  have small weak diameter in initial input graph G0
    // TODO: check if terminate means really terminate
    int v = *graph->V.begin();
    set<int> ballInTest = ballIn(g0, v, D/2, INPUT_N).first;
    if (!isSubset(ballInTest, graph->V))
        return fromMatrixToSet(graph->is_edge);
    set<int> ballOutTest = ballIn(g0, v, D/2, INPUT_N).first;
    if (!isSubset(ballOutTest, graph->V))
        return fromMatrixToSet(graph->is_edge);
    
    return Erem;
}

PriceFunction elimNeg(Graph *graph, int s, int INPUT_N) {

    set<pair<int, int>> e = fromMatrixToSet(graph->is_edge);
    set<pair<int, int>> eNeg;
    for (int v:graph->V) {
        for (int u:graph->V) {
            if (graph->is_edge[v][u]) {
                eNeg.insert({v, u});
            }
        }
    }
    set<pair<int, int>> e_minus_eNeg = edgesMinusEdges(e, eNeg);

    vector<int> prices(INPUT_N, 0);
    for (int v:graph->V) {
        prices[v] = INT32_MAX;
    }
    prices[s] = 0;
    priority_queue<pair<int, int>> Q;
    Q.push({0, s});
    set<int> marked;
    while (!Q.empty()) {
        while (!Q.empty()) {
            int v = Q.top().second;
            int d = -Q.top().first;
            Q.pop();
            if (prices[v] != d)
                continue;
            marked.insert(v);
            for (pair<int, int> edge : e_minus_eNeg) {
                if (d + graph->adj[edge.first][edge.second] < prices[edge.second]) {
                    Q.push({-prices[edge.second], edge.second});
                    prices[edge.second] = d + graph->adj[edge.first][edge.second];
                }
            }
        }

        for (int v : marked) {
            for (pair<int, int> edge : eNeg) {
                if (prices[edge.first] + graph->adj[edge.first][edge.second] < prices[edge.second]) {
                    Q.push({-prices[edge.second], edge.second});
                    prices[edge.second] = prices[edge.first] + graph->adj[edge.first][edge.second];
                }
            }
        }
    }

    PriceFunction result;
    result.prices = prices;
    // DEBUG 
    for (int i = 0; i < result.prices.size(); i++) {
        assert(result.prices[i] != INT32_MAX);
    }
    return result;
}

PriceFunction PriceFunction::sum(PriceFunction a, PriceFunction b) {
    assert (a.prices.size() == b.prices.size());
    PriceFunction c;
    for (int i = 0; i < a.prices.size(); i++) {
        c.prices.push_back(a.prices[i]+b.prices[i]);
    }
    return c;
}

PriceFunction scaleDown(Graph *graph, int delta, int B, int s, int INPUT_N) {
    PriceFunction Phi2;
    Graph * graph_B_Phi2;
    if (delta > 2) {
        int d = delta/2;
        Graph* graph_B_pos = new Graph(graph->V, INPUT_N);
        graph_B_pos->adj = graph->adj;
        graph_B_pos->is_edge = graph->is_edge;

        for (int i : graph_B_pos->V) {
            for (int j : graph_B_pos->V) {
                if (graph_B_pos->is_edge[i][j]) {
                    if (graph_B_pos->adj[i][j] < 0)
                        graph_B_pos->adj[i][j] = max(graph_B_pos->adj[i][j]+B, 0);
                }
            }
        }

        // phase 0: Decompose V to SCCs V1, V2... with weak diameter dB in G
        set<pair<int, int>> Erem = LDD(graph_B_pos, d*B, INPUT_N);
        Graph* graph_B = addIntegerToEdges(graph, B, INPUT_N);
        Graph* graph_B_rem = subtractEdges(graph_B, Erem, INPUT_N);
        vector<set<int>> SCCs = computeSCCs(graph_B_rem, INPUT_N);
        // phase 1: Make edges inside the SCCs G^B[V_i] non-negative
        Graph* H = induced_graph(graph, *SCCs.begin(), INPUT_N);
        for (set<int> SCC : SCCs) {
            if (SCC == *SCCs.begin())
                continue;
            H = mergeGraphs(H, induced_graph(graph, SCC, INPUT_N), INPUT_N);
        }
        PriceFunction Phi1 = scaleDown(H, delta/2, B, s, INPUT_N);
        // phase 2: Make all edges in G^B \ E^rem non-negative
        Graph* graph_B_Phi1 = applyPriceFunction(graph_B, Phi1, INPUT_N);
        Graph* graph_B_Phi1_rem = subtractEdges(graph_B_Phi1, Erem, INPUT_N);
        PriceFunction psi = FixDAGEdges(graph_B_Phi1_rem, SCCs, INPUT_N);
        Phi2 = PriceFunction::sum(Phi1, psi);
    } else {
        Phi2.prices.assign(INPUT_N, 0);
    }
    // phase 3: Make all edges in G^B non-negative
    graph_B_Phi2 = applyPriceFunction(addIntegerToEdges(graph, B, INPUT_N), Phi2, INPUT_N);
    PriceFunction psi_first = elimNeg(graph_B_Phi2, s, INPUT_N);
    PriceFunction Phi3 = PriceFunction::sum(Phi2, psi_first);

    return Phi3;
}

SSSP_Result SPmain(Graph* g_in, int s_in, int INPUT_N) {
    Graph* g_up = new Graph(g_in->V, INPUT_N);
    g_up->adj = g_in->adj;
    g_up->is_edge = g_in->is_edge;
    for (int i : g_up->V) {
        for (int j: g_up->V) {
            if (g_up->is_edge[i][j])
                g_up->adj[i][j]*=2*g_in->V.size();
        }
    }

    // Round B up to nearest power of 2
    int B = 2*g_in->V.size();
    B = roundB(B);

    // Identity price function
    PriceFunction Phi0;
    Phi0.prices.assign(INPUT_N, 0);

    for (int i = 1; pow(2, i)<=B; i++) {
        Graph* graph_B_phi0 = new Graph(g_up->V, INPUT_N);
        graph_B_phi0->adj = g_up->adj;
        graph_B_phi0->is_edge = g_up->is_edge;
        // Apply the price function
        for (int j : graph_B_phi0->V) {
            for (int k : graph_B_phi0->V) {
                if (graph_B_phi0->is_edge[j][k])
                    graph_B_phi0->adj[j][k] += Phi0.prices[j] - Phi0.prices[k];
            }
        }
        // Call scaleDown function
        PriceFunction Psi0 = scaleDown(graph_B_phi0, g_in->V.size(), B/pow(2, i), s_in, INPUT_N);
        Phi0 = PriceFunction::sum(Phi0, Psi0);
    }

    Graph* g_star = new Graph(g_up->V, INPUT_N);
    g_star->adj = g_up->adj;
    g_star->is_edge = g_up->is_edge;
    for (int j : g_star->V) {
        for (int k : g_star->V) {
            if (g_star->is_edge[j][k])
                g_star->adj[j][k] += Phi0.prices[j] - Phi0.prices[k] + 1;
        }
    }

    SSSP_Result result = dijkstra(g_star, s_in, INPUT_N);
    return result;
}

PriceFunction FixDAGEdges(Graph* graph, vector<set<int>> SCCs, int INPUT_N) {
    PriceFunction Phi;
    Phi.prices.assign(INPUT_N, 0);
    // Sort SCCs vector in topological order
    vector<int> old_fromVertixToSCC(INPUT_N, -1);
    for (int i = 0; i < SCCs.size(); i++) {
        for (int v : SCCs[i]) {
            old_fromVertixToSCC[v] = i;
        }
    }
    vector<set<int>> SCCs_topo = topologicalSort(SCCs, graph, old_fromVertixToSCC);

    vector<int> fromVertixToSCC(INPUT_N, -1);
    for (int i = 0; i < SCCs_topo.size(); i++) {
        for (int v : SCCs_topo[i]) {
            fromVertixToSCC[v] = i;
        }
    }

    vector<int> mu_j(SCCs_topo.size(), 0);
    for (int v : graph->V) {
        for (int w : graph->V) {
            if (graph->is_edge[v][w] && graph->adj[v][w] < 0) {
                if (fromVertixToSCC[v] != fromVertixToSCC[w]) {
                    mu_j[fromVertixToSCC[w]] = min(mu_j[fromVertixToSCC[w]], graph->adj[v][w]);
                }
            }
        }
    }

    int mLast = 0;
    for (int j = 1; j < SCCs_topo.size(); j++) {
        mLast += mu_j[j];
        for (int v : SCCs_topo[j]) {
            Phi.prices[v] = mLast;
        }
    }

    return Phi;
}
