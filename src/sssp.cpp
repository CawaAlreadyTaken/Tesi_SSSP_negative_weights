#include "sssp.h"
#include "graph.h"
#include "utils.h"

set<int> LDD(Graph* graph, int D) {

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
    Graph * graph_B_Phi2 = new Graph(graph->n, graph->m);
    graph_B_Phi2->adj = graph->adj;
    // Add B to negative edges TODO check that questo sia davvero da fare prima di applicare la price function
    for (int i = 0; i < graph->n; i++) {
        assert(graph->adj[i].size() == graph->weights[i].size());
        for (int j = 0; j < graph->adj[i].size(); j++) {
            if (graph->weights[i][j] < 0)
                graph_B_Phi2->weights[i].push_back(graph->weights[i][j]+B);
            else
                graph_B_Phi2->weights[i].push_back(graph->weights[i][j]);
        }
    }

    if (delta > 2) {
        int d = delta/2;  // TODO check this
        Graph* graph_B_pos = new Graph(graph->n, graph->m);
        graph_B_pos->adj = graph->adj;
        graph_B_pos->weights.resize(graph->n);

        for (int i = 0; i < graph_B_pos->n; i++) {
            for (int j = 0; j < graph_B_pos->adj[i].size(); j++) {
                if (graph->weights[i][j] < 0)
                    graph_B_pos->weights[i].push_back(max(graph->weights[i][j]+B, 0));
                else
                    graph_B_pos->weights[i].push_back(graph->weights[i][j]);
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
        for (int i = 0; i < graph_B_Phi2->n; i++) {
            assert(graph_B_Phi2->adj[i].size() == graph_B_Phi2->weights[i].size());
            for (int j = 0; j < graph_B_Phi2->weights[i].size(); j++) {
                graph_B_Phi2->weights[i][j] = graph_B_Phi2->weights[i][j] + Phi2.prices[i] - Phi2.prices[graph_B_Phi2->adj[i][j]];
            }
        }
    } else {
    Phi2.prices.assign(graph->n, 0);
}
    // phase 3: Make all edges in G^B non-negative
    PriceFunction psi_first = elimNeg(graph_B_Phi2);
    PriceFunction Phi3 = PriceFunction::sum(Phi2, psi_first);

    return Phi3;
}

SSSP_Result SPmain(Graph* g_in, int s_in) {
    Graph* g_up = new Graph(g_in->n, g_in->m);
    vector<vector<int>> w_up;
    w_up.assign(g_in->n, vector<int>());
    for (int i = 0; i < g_in->n; i++) {
        for (int weight : g_in->weights[i]) {
            w_up[i].push_back(weight*2*g_in->n);
        }
    }
    g_up->adj = g_in->adj;
    g_up->weights = w_up;

    // Round B up to nearest power of 2
    int B = 2*g_in->n;
    B = roundB(B);

    // Identity price function
    PriceFunction Phi0;
    Phi0.prices.assign(g_in->n, 0);

    for (int i = 1; pow(2, i)<=B; i++) {
        Graph* graph_B_phi0 = new Graph(g_up->n, g_up->m);
        graph_B_phi0->adj = g_up->adj;
        graph_B_phi0->weights.resize(g_up->n);
        // Apply the price function
        for (int j = 0; j < g_up->n; j++) {
            for (int k = 0; k < g_up->weights[j].size(); k++) {
                graph_B_phi0->weights[j].push_back(g_up->weights[j][k]+Phi0.prices[j]-Phi0.prices[g_up->adj[j][k]]);
            }
        }
        // Call scaleDown function
        PriceFunction Psi0 = scaleDown(graph_B_phi0, g_in->n, B/pow(2, i));
        Phi0 = PriceFunction::sum(Phi0, Psi0);
    }

    Graph* g_star = new Graph(g_up->n, g_up->m);
    g_star->adj = g_up->adj;
    g_star->weights.assign(g_star->n, vector<int>());
    for (int i = 0; i < g_up->n; i++) {
        for (int j = 0; j < g_up->weights[i].size(); j++) {
            int weight = g_up->weights[i][j];
            g_star->weights[i].push_back(weight+Phi0.prices[i]-Phi0.prices[g_up->adj[i][j]]+1);
        }
    }

    SSSP_Result result = dijkstra(g_star, s_in);
    return result;
}
