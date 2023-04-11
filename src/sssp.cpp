#include "sssp.h"
#include "graph.h"

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
