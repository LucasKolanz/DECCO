
#include "simple_graph.hpp"
// #include "vec3.hpp"

#define TOUCH_TOLERANCE 1.01


void addEdge(Graph& g, int u, int v) 
{
    g[u].push_back(v);
    g[v].push_back(u); // if undirected
}

void makeGraph(Graph& g, vec3* pos, double* R, int n)
{
	//loop over balls
	for (int i = 0; i < n; ++i)
	{
		g.push_back({});
	}


	//loop over pairs adding edges if balls touch 
	// for (int i = 1; i < n; ++i)
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			//check if balls are touching (give some leeway if they arent exactly touching)
			const double dist = (pos[j]-pos[i]).norm();
			const double sum_radii = R[i] + R[j];
			if (TOUCH_TOLERANCE*sum_radii-dist > 0)//touching condition
			{
				addEdge(g,i,j);
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