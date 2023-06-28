#ifndef SSSP_H
#define SSSP_H

#include "graph.h"

using namespace std;

class PriceFunction {
    // Will be an integral price function
    public:
    vector<long long> prices;
    static PriceFunction sum(PriceFunction a, PriceFunction b);
};

SSSP_Result SPmain(Graph* g_in, int s_in);
PriceFunction scaleDown(Graph* graph, long long delta, long long B, int depth);
PriceFunction elimNeg(Graph* graph);
PriceFunction FixDAGEdges(Graph* g, vector<set<int>>& v);
pair<set<int>, set<pair<int, int>>> ballOut(Graph* graph, int v, int R, bool fromBallIn=false);
pair<set<int>, set<pair<int, int>>> ballIn(Graph* graph, int v, int R);
set<pair<int, int>> LDD(Graph* graph, int D, int depth);

#endif

