#include "sssp.h"
#include "graph.h"
#include "utils.h"

set<int> ballIn(Graph* graph, int v, int R) {
    // In order to do ballIn, we can reverse the edges and call ballOut
    Graph* graph_copy = new Graph(graph->V);
    for (int i : graph->V) {
        for (int j : graph->V) {
            if (graph->is_edge[i][j])
                graph_copy->add_edge(j, i, graph->adj[i][j]);
        }
    }
    return ballOut(graph_copy, v, R);
}

set<int> ballOut(Graph* graph, int v, int R) {
    // basically dijkstra in order to find vertices closer than R
    set<int> ris;

    vector<bool> confirmed(MAX_N, false);
    confirmed[v] = true;
    ris.insert(v);

    // TODO: maybe can remove indexFrom
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

        if (weight > R) // Difference from dijstra
            break;

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

    return ris;
}

set<int> LDD(Graph* graph, int D) {
    // Save a copy of the original graph, since we are going to modify it
    Graph* g0 = new Graph(graph->V);
    g0->adj = graph->adj;
    g0->is_edge = graph->is_edge;

    set<int> Erem;
    // Phase 1: mark vertices as light or heavy
    float k = c*ln(INPUT_N);
    set<int> S = getRandomVertices(graph, k);
    /*
    map<int, set<int>> sBallIns;
    map<int, set<int>> sBallOuts;
    map<int, set<int>> ballInIntersec;
    map<int, set<int>> ballOutIntersec;

    for (auto x : S) {
        if (sBallIns[x].empty()) {
            sBallIns[x] = ballIn(graph, x, D/4);
            sBallOuts[x] = ballOut(graph, x, D/4);
        }
    }

    for (auto v : graph->V) {
        for (auto)
    }
    */

    map<int, set<int>> ballInIntersec;
    map<int, set<int>> ballOutIntersec;

    // Da ricontrollare
    for (int x : S) {
        set<int> bIn = ballIn(graph, x, D/4);
        set<int> bOut = ballOut(graph, x, D/4);
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
    while (in_light.size()) {
        int v = in_light.back(); // Check unmarkation
        double p = d_min(1, 80.0*log2(INPUT_N)/D);
        int Rv = sampleGeo(p);
        set<int> newBallIn = ballIn(graph, v, Rv);
        set<int> Ebound = getBoundariesIn(newBallIn);
        if (Rv > D/4 || newBallIn.size() > 0.7*graph->V.size()) {
            //return Erem = E(G) and terminate
        }
        set<int> Erecurs = LDD(induced_graph(graph, newBallIn), D);
        Erem = intersect(intersect(Erem, Ebound), Erecurs);
        graph = subtract(graph, newBallIn);
    }

    while (out_light.size()) {
        int v = out_light.back(); // Check unmarkation
        double p = d_min(1, 80.0*log2(INPUT_N)/D);
        int Rv = sampleGeo(p);
        set<int> newBallOut = ballOut(graph, v, Rv);
        set<int> Ebound = getBoundariesOut(newBallOut);
        if (Rv > D/4 || newBallOut.size() > 0.7*graph->V.size()) {
            //return Erem = E(G) and terminate
        }
        set<int> Erecurs = LDD(induced_graph(graph, newBallOut), D);
        Erem = intersect(intersect(Erem, Ebound), Erecurs);
        graph = subtract(graph, newBallOut);
    }

    // Clean Up: check that remaining vertices  have small weak diameter in initial input graph G0
    // TODO
    
    return Erem;
}

PriceFunction elimNeg(Graph *graph) {
}

PriceFunction PriceFunction::sum(PriceFunction a, PriceFunction b) {
    assert (a.prices.size() == b.prices.size());
    PriceFunction c;
    for (int i = 0; i < a.prices.size(); i++) {
        c.prices.push_back(a.prices[i]+b.prices[i]);
    }
    return c;
}

PriceFunction scaleDown(Graph *graph, int delta, int B) {
    PriceFunction Phi2;
    Graph * graph_B_Phi2 = new Graph(graph->V);
    graph_B_Phi2->adj = graph->adj;
    graph_B_Phi2->is_edge = graph_B_Phi2->is_edge;

    // Add B to negative edges TODO check that questo sia davvero da fare prima di applicare la price function
    for (int i : graph_B_Phi2->V) {
        for (int j : graph_B_Phi2->V) {
            if (graph_B_Phi2->is_edge[i][j]) {
                if (graph_B_Phi2->adj[i][j] < 0)
                    graph_B_Phi2->adj[i][j]+=B;
            }
        }
    }

    if (delta > 2) {
        int d = delta/2;  // TODO check this
        Graph* graph_B_pos = new Graph(graph->V);
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
        set<int> Erem = LDD(graph_B_pos, d*B);
        int SCCs = computeSCCs();
        // phase 1: Make edges inside the SCCs G^B[V_i] non-negative
        Graph* H = new Graph();
        PriceFunction Phi1 = scaleDown(H, delta/2, B);
        // phase 2: Make all edges in G^B \ E^rem non-negative
        PriceFunction psi = FixDAGEdges(graph_B_Phi1_rem, SCCs);
        Phi2 = PriceFunction::sum(Phi1, psi);
        // Add Phi2 to graph_B_Phi2 (this should be in phase 3 but if phi2 is all 0 we can avoid doing it)
        for (int i : graph_B_Phi2->V) {
            for (int j : graph_B_Phi2->V) {
                if (graph_B_Phi2->is_edge[i][j]) {
                    graph_B_Phi2->adj[i][j]+=Phi2.prices[i]-Phi2.prices[j];
                }
            }
        }
    } else {
        Phi2.prices.assign(MAX_N, 0);
    }
    // phase 3: Make all edges in G^B non-negative
    PriceFunction psi_first = elimNeg(graph_B_Phi2);
    PriceFunction Phi3 = PriceFunction::sum(Phi2, psi_first);

    return Phi3;
}

SSSP_Result SPmain(Graph* g_in, int s_in) {
    Graph* g_up = new Graph(g_in->V);
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
    Phi0.prices.assign(MAX_N, 0);

    for (int i = 1; pow(2, i)<=B; i++) {
        Graph* graph_B_phi0 = new Graph(g_up->V);
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
        PriceFunction Psi0 = scaleDown(graph_B_phi0, g_in->V.size(), B/pow(2, i));
        Phi0 = PriceFunction::sum(Phi0, Psi0);
    }

    Graph* g_star = new Graph(g_up->V);
    g_star->adj = g_up->adj;
    g_star->is_edge = g_up->is_edge;
    for (int j : g_star->V) {
        for (int k : g_star->V) {
            if (g_star->is_edge[j][k])
                g_star->adj[j][k] += Phi0.prices[j] - Phi0.prices[k] + 1;
        }
    }

    SSSP_Result result = dijkstra(g_star, s_in);
    return result;
}
