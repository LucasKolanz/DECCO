
#include "simple_graph.hpp"



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
	for (int i = 1; i < n; ++i)
	{
		for (int j = 0; j < i; ++j)
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

    while (!q.empty()) {
        int current = q.front();
        q.pop();
        for (int neighbor : g[current]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
                ++countVisited;
            }
        }
    }

    // If we've visited all vertices, the graph is connected
    return (countVisited == (int)g.size());
}