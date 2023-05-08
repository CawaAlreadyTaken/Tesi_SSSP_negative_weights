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

using namespace std;

extern bool terminateLDD;

int roundB(int b);
double d_min(double a, double b);
SSSP_Result dijkstra(Graph* g, int s, int INPUT_N);
set<int> getRandomVertices(Graph* g, int k, int INPUT_N);
Graph* induced_graph(Graph* g, set<int> vertices, int INPUT_N);
set<pair<int, int>> intersect(set<pair<int, int>> a, set<pair<int, int>> b);
Graph* subtractVertices(Graph* g, set<int> vertices, int INPUT_N);
Graph* subtractEdges(Graph* g, set<pair<int, int>> edges, int INPUT_N);
set<pair<int, int>> fromMatrixToSet(vector<vector<bool>> isEdge);
set<pair<int, int>> edgesMinusEdges(set<pair<int, int>> a, set<pair<int, int>> b);
bool isSubset(set<int> a, set<int> b);
Graph* addIntegerToNegativeEdges(Graph* g, int e, int INPUT_N);
Graph* applyPriceFunction(Graph* g, PriceFunction p, int INPUT_N);
Graph* mergeGraphs(Graph* g1, Graph* g2, int INPUT_N);
vector<set<int>> computeSCCs(Graph* g, int INPUT_N);
Graph* transpose(Graph* g, int INPUT_N);
void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s);
void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component);
vector<set<int>> topologicalSort(vector<set<int>> SCCs, Graph* graph, vector<int> fromVertixToSCC);
void topoDFS(int index, vector<bool>& visited, stack<int>& s, vector<vector<int>>& sccAdj);
void printGraph(Graph* g);
bool checkConstantOutDegree(Graph* g);
Graph* addDummySource(Graph* g, int INPUT_N);

#endif
