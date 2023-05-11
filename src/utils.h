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
extern int INPUT_N;

void log(bool _cerr, int depth, string message);

int roundB(int b);
double d_min(double a, double b);
SSSP_Result dijkstra(Graph* g, int s);
set<int> getRandomVertices(Graph* g, int k);
Graph* induced_graph(Graph* g, set<int> vertices);
set<pair<int, int>> edgesUnion(set<pair<int, int>> a, set<pair<int, int>> b);
Graph* subtractVertices(Graph* g, set<int> vertices);
Graph* subtractEdges(Graph* g, set<pair<int, int>> edges);
set<pair<int, int>> fromMatrixToSet(vector<vector<bool>> isEdge);
set<pair<int, int>> edgesMinusEdges(set<pair<int, int>> a, set<pair<int, int>> b);
bool isSubset(set<int> a, set<int> b);
Graph* addIntegerToNegativeEdges(Graph* g, int e);
Graph* applyPriceFunction(Graph* g, PriceFunction p);
Graph* mergeGraphs(Graph* g1, Graph* g2);
vector<set<int>> computeSCCs(Graph* g);
Graph* transpose(Graph* g);
void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s);
void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component);
vector<set<int>> topologicalSort(vector<set<int>> SCCs, Graph* graph, vector<int> fromVertixToSCC);
void topoDFS(int index, vector<bool>& visited, stack<int>& s, vector<vector<int>>& sccAdj);
void printGraph(Graph* g, int depth, bool cerr);
bool checkConstantOutDegree(Graph* g);
Graph* addDummySource(Graph* g);

void print_shortest_path_tree(SSSP_Result result);

#endif
