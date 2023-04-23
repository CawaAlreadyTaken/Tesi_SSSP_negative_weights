#ifndef UTILS_H
#define UTILS_H

#include <bits/stdc++.h> // TODO: use more specific libraries
#include "graph.h"
#include "sssp.h"

using namespace std;

int roundB(int b);
double d_min(double a, double b);
SSSP_Result dijkstra(Graph* g, int s, int INPUT_N);
set<int> getRandomVertices(Graph* g, int k, int INPUT_N);
Graph* induced_graph(Graph* g, set<int> vertices, int INPUT_N);
set<pair<int, int>> intersect(set<pair<int, int>> a, set<pair<int, int>> b);
Graph* subtractVertices(Graph* g, set<int> vertices, int INPUT_N);
Graph* subtractEdges(Graph* g, set<pair<int, int>> edges, int INPUT_N);
set<pair<int, int>> fromMatrixToSet(vector<vector<bool>> isEdge);
bool isSubset(set<int> a, set<int> b);
Graph* addIntegerToEdges(Graph* g, int e, int INPUT_N);
Graph* applyPriceFunction(Graph* g, PriceFunction p, int INPUT_N);
Graph* mergeGraphs(Graph* g1, Graph* g2, int INPUT_N);
vector<set<int>> computeSCCs(Graph* g, int INPUT_N);
Graph* transpose(Graph* g, int INPUT_N);
void DFS(Graph* g, int v, vector<bool>& visited, stack<int>& s);
void DFSaddComp(Graph* g, int v, vector<bool>& visited, set<int>& component);

#endif
