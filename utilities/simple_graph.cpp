
#include "simple_graph.hpp"
// #include "vec3.hpp"


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
			if (1.1*sum_radii-dist > 0)//touching condition
			{
				addEdge(g,i,j);
			}
		}
	}

}

bool isConnected(vec3* pos, double* R, int n)
{
	Graph g;
	std::cerr<<"n: "<<n<<std::endl;
	makeGraph(g,pos,R,n);

	if (g.empty()) //If it is empty we have a problem somewhere
	{
		MPIsafe_print(std::cerr,"Error makeing graph in simple_graph. Graph is empty. Now exiting. . .\n");
        MPIsafe_exit(-1);
	}

    std::vector<bool> visited(g.size(), false);
    std::queue<int> q;
    
    // Start BFS from vertex 0
    visited[0] = true;
    q.push(0);
    
    int countVisited = 1;

    while (!q.empty()) 
    {
        int current = q.front();
        q.pop();
        for (int neighbor : g[current]) 
        {
            if (!visited[neighbor]) 
            {
                visited[neighbor] = true;
                q.push(neighbor);
                ++countVisited;
            }
        }
    }

    std::cerr<<"countVisited: "<<countVisited<<std::endl;
    std::cerr<<"(int)g.size(): "<<(int)g.size()<<std::endl;
    // If we've visited all vertices, the graph is connected
    return (countVisited == n);
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