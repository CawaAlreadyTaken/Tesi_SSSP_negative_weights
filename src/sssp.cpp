#include "sssp.h"
#include "graph.h"
#include "utils.h"

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
    int B = 2*g_in->n;
    B = roundB(B);
    // Identity price function
    PriceFunction Phi0;
    Phi0.prices.assign(g_in->n, 0);
    for (int i = 1; pow(2, i)<=B; i++) {
        Graph* graph_B_phi0 = new Graph(g_up->n, g_up->m);
        graph_B_phi0->adj = g_up->adj;
        PriceFunction Psi0 = scaleDown(graph_B_phi0, g_in->n, B/pow(2, i));
        Phi0 = PriceFunction::sum(Phi0, Psi0);
    }
    Graph* g_star = new Graph(g_up->n, g_up->m);
    g_star->adj = g_up->adj;
    g_star->weights.assign(g_star->n, vector<int>());
    for (int i = 0; i < g_up->n; i++) {
        for (int j = 0; j < g_up->weights[i].size(); j++) {
            int weight = g_up->weights[i][j];
            g_star->weights[i].push_back(weight+Phi0.prices[g_up->adj[i][j]]-Phi0.prices[i]+1);
        }
    }
    SSSP_Result result = dijkstra(g_star, s_in);
}

PriceFunction elimNeg(Graph *graphA) {

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
    if (delta > 2) {
        int d = delta/2;  // TODO check this
        // phase 0
        // phase 1
        // phase 2
    } else {
        Phi2.prices.assign(2, 1);
    }
    // phase 3
    PriceFunction psi = elimNeg(graph_B_Phi2);
    PriceFunction Phi3 = PriceFunction::sum(Phi2, psi);
    return Phi3;
}
