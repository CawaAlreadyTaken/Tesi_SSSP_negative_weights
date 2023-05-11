#!/usr/bin/env python3
# Author: Daniele Cabassi

# Create a graph with N vertices numbered from 1 to N. Start is S.
# Each vertex has out-degree D and random integral weight [-1, X].
# P is the probability for a weight of being -1, 1-P for [0, X].

import random

S = 1
N = 100
D = 6
X = 10
P = 0.005
assert D < N
assert P >= 0 and P <= 1
assert X >= 0
negative_edges_number = 0
weight_sum = 0

adj = [[] for _ in range(N+1)]
for i in range(1, N+1):
    tmp = [j for j in range(1, N+1) if j != i]
    random.shuffle(tmp)
    for _ in range(D):
        r = random.random()
        weight = -1 if r <= P else random.randint(0, X)
        if weight == -1:
            negative_edges_number += 1
        weight_sum += weight
        adj[i].append((tmp[0], weight))
        tmp = tmp[1:]

padding = "____________________________________________________\n"

fileInfo = "Fileinfo generated automatically by create_graphs.py\n\n"
fileInfo += padding
fileInfo += f"Start vertex: {S}\n"
fileInfo += padding
fileInfo += f"Total vertices number: {N}\n"
fileInfo += padding
fileInfo += f"Out-degree: {D}, total edges number: {N*D}\n"
fileInfo += padding
fileInfo += f"Negative weights probability: {P}, negative weights number: {negative_edges_number}, non-negative weights number: {N*D-negative_edges_number}\n"
fileInfo += padding
fileInfo += f"Maximum weight possible: {X}, average weight: {weight_sum/(N*D)}\n"

inputFile = f"{N} {N*D} {S}\n"
for i in range(1, N+1):
    for j, w in adj[i]:
        inputFile += f"{i} {j} {w}\n"

with open("input.info", "w") as f:
    f.write(fileInfo)

with open("input.txt", "w") as f:
    f.write(inputFile)
