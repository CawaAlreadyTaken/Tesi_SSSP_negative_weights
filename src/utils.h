#ifndef UTILS_H
#define UTILS_H

#include "graph.h"
#include "sssp.h"

#include <iostream>
#include <stack>
#include <random>
#include <map>
#include <cassert>
#include <queue>
#include <chrono>

using namespace std;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

extern bool terminateLDD;
extern int INPUT_N;

int roundB(int b);
double d_min(double a, double b);
set<int> getRandomVertices(Graph* g, int k);
set<pair<int, int>> edgesUnion(set<pair<int, int>> a, set<pair<int, int>>& b);
set<pair<int, int>> fromGraphToSetOfEdges(Graph* g);
set<pair<int, int>> edgesMinusEdges(set<pair<int, int>>& a, set<pair<int, int>>& b);
Graph* induced_graph(Graph* g, set<int>& vertices);
Graph* subtractVertices(Graph* g, set<int>& vertices);
Graph* subtractEdges(Graph* g, set<pair<int, int>>& edges);
Graph* addIntegerToNegativeEdges(Graph* g, int e);
Graph* applyPriceFunction(Graph* g, PriceFunction p);
Graph* mergeGraphs(Graph* g1, Graph* g2);
Graph* transpose(Graph* g);
Graph* addDummySource(Graph* g);
Graph* onlyEdgesInsideSCCs(Graph* g, vector<set<int>>& SCCs);
void log(bool _cerr, int depth, string message);
void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s);
void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component);
void topoDFS(int index, vector<bool>& visited, stack<int>& s, vector<vector<int>>& sccAdj);
void printGraph(Graph* g, int depth, bool cerr);
void print_shortest_path_tree(SSSP_Result result);
bool isSubset(set<int>& a, set<int>& b);
bool checkConstantOutDegree(Graph* g);
vector<set<int>> computeSCCs(Graph* g, int depth);
vector<set<int>> topologicalSort(vector<set<int>>& SCCs, Graph* graph, vector<int>& fromVertixToSCC);
map<pair<int, int>, int> createSupportMap(Graph* g);
SSSP_Result dijkstra(Graph* g, int s);


#endif
