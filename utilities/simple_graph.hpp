



#pragma once
#ifndef SIMPLE_GRAPH_HPP
#define SIMPLE_GRAPH_HPP

#include <vector>
#include <iostream>
#include <queue>
#include "vec3.hpp"
#include "MPI_utilities.hpp"

// We'll assume vertices are labeled 0, 1, 2, ... up to (n-1).
// 'adjList[v]' is a list of the neighbors of vertex 'v'.
using Graph = std::vector<std::vector<int>>;


void addEdge(Graph& g, int u, int v);
void makeGraph(
    Graph& g,
    const vec3* pos,
    const double* R,
    int n,
    const vec3& boxdims=vec3{0,0,0},
    bool periodic=false);
// bool isConnected(vec3* pos, double* R, int n);




#endif