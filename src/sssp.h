#ifndef SSSP_H
#define SSSP_H

#include <bits/stdc++.h> // TODO: use more specific libraries
#include "graph.h"

using namespace std;

class PriceFunction {
    // Will be an integral price function

    public:
    vector<int> prices;
    static PriceFunction sum(PriceFunction a, PriceFunction b);
};

PriceFunction scaleDown(Graph* graph, int delta, int B);
PriceFunction elimNeg(Graph* graph);
SSSP_Result SPmain(Graph* g_in, int s_in);

#endif

