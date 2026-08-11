
#include <cmath>
#include "simple_graph.hpp"

#define TOUCH_TOLERANCE 1.01


void addEdge(Graph& g, int u, int v) 
{
    g[u].push_back(v);
    g[v].push_back(u); // if undirected
}

void makeGraph(
    Graph& g,
    const vec3* pos,
    const double* R,
    int n,
    const vec3& boxdims,
    bool periodic)
{
    // Make exactly n empty vertices.
    g.clear();
    g.resize(n);

    // Loop over each unique pair once.
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            vec3 dr = pos[j] - pos[i];

            if (periodic)
            {
                // Minimum-image convention.
                dr.x -= boxdims.x * std::round(dr.x / boxdims.x);
                dr.y -= boxdims.y * std::round(dr.y / boxdims.y);
                dr.z -= boxdims.z * std::round(dr.z / boxdims.z);
            }

            const double dist = dr.norm();
            const double sum_radii = R[i] + R[j];

            if (dist < TOUCH_TOLERANCE * sum_radii)
            {
                addEdge(g, i, j);
            }
        }
    }
}





// int main()
// {
// 	// int num_particles = 3;
// 	// vec3* pos = nullptr;
// 	// double* R = nullptr;
// 	// pos = new vec3[num_particles];
// 	// R = new double[num_particles];

// 	// R[0] = 1;
// 	// R[1] = 1;
// 	// R[2] = 1;

// 	// pos[0] = {0,0,0};
// 	// pos[1] = {0,2,0};
// 	// pos[2] = {0,4,0};

// 	std::string path = '/home/kolanzl/Desktop/Visualize/V19/';

// 	std::cerr<<isConnected(pos,R,3)<<std::endl;
// }