#ifndef UTILS_H
#define UTILS_H

#include <bits/stdc++.h> // TODO: use more specific libraries
#include "graph.h"

using namespace std;

int roundB(int b);
double d_min(double a, double b);
SSSP_Result dijkstra(Graph* g, int s);
set<int> getRandomVertices(Graph* g, int k);

#endif

