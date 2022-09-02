#include <iostream>
#include <vector>
#include "vec3.hpp"

struct node
{
	int ball_index;
	// vec3 pos;
};

class graph
{
public:
	int numVerts;
	// std::vector<node> nodes;
	std::vector<std::vector<node>> adj;
	std::vector<node> nodes;

	graph() = default;

	void addEdge(node n1, node n2)
	{
		int n1ind = getNodeNum(n1);
		int n2ind = getNodeNum(n2);
		adj[n1ind].push_back(n2);
		adj[n2ind].push_back(n1);
	}

	void addNode(int ball_ind)
	{
		
	}

	int getNodeIndex(node n);
	{
		for (auto it = begin (nodes); it != end (nodes); ++it) 
		{
			if (it -> ball_index == n.ball_index)
			{
				return node
			}
		}
		return -1
	}
	~graph() {};
};

// class graph_group
// {
// public:
// 	std::vector<graph> graphs;
// 	int num_graphs = 0;
// 	graph_group() = default;
// 	~graph_group() {};
	
// 	void addGraph(graph g)
// 	{
// 		graphs.push_back(g)
// 		num_graphs++;
// 		// for (auto it = begin (graphs); it != end (graphs); ++it) 
// 		// {
//   //   		it->doSomething ();
// 		// }
// 		return;
// 	}

	

	
// 	// mergeGroups(graph_group gg);
// private:	
// };