#ifndef SSSP_H
#define SSSP_H

#include "graph.h"

using namespace std;

class PriceFunction {
    // Will be an integral price function
    public:
    vector<int> prices;
    static PriceFunction sum(PriceFunction a, PriceFunction b);
};

SSSP_Result SPmain(Graph* g_in, int s_in);
PriceFunction scaleDown(Graph* graph, int delta, int B, int depth);
PriceFunction elimNeg(Graph* graph);
set<pair<int, int>> LDD(Graph* graph, int D, int depth);
pair<set<int>, set<pair<int, int>>> ballOut(Graph* graph, int v, int R, bool fromBallIn=false);
pair<set<int>, set<pair<int, int>>> ballIn(Graph* graph, int v, int R);
PriceFunction FixDAGEdges(Graph* g, vector<set<int>> v);

#endif

