#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include "vec3.hpp"


class graph
{
public:
	int numVerts;
	std::vector<int>* adj = nullptr;

	graph() = default;

	// graph(int num_balls)
	void init(int num_balls)
	{
		adj = new std::vector<int>[num_balls];
		numVerts = num_balls;
		std::cout<<"numVerts in init: "<<numVerts<<std::endl;
		std::cout<<"num_balls in init: "<<num_balls<<std::endl;
		return;
	}

	void addEdge(int n1, int n2)
	{
		if (std::find(adj[n1].begin(), adj[n1].end(), n2) == adj[n1].end()) 
		{
			adj[n1].push_back(n2);
		}
		if (std::find(adj[n2].begin(), adj[n2].end(), n1) == adj[n2].end()) 
		{
			adj[n2].push_back(n1);
		}
		return;
	}

	void deleteEdge(int n1, int n2)
	{
		auto find_me = std::find(adj[n1].begin(), adj[n1].end(), n2);
		if (find_me != adj[n1].end()) 
		{
			adj[n1].erase(find_me);
		}
		find_me = std::find(adj[n2].begin(), adj[n2].end(), n1);
		if (find_me != adj[n2].end()) 
		{
			adj[n2].erase(find_me);
		}
		return;
	}

	void printGraph() 
	{
		std::cout<<"numVerts in gg: "<<numVerts<<std::endl;
    	for (int d = 0; d < numVerts; ++d) 
    	{
    		std::cout << "\n Vertex "<< d << ":";
    		for (auto x : adj[d])
      		{
      			std::cout << "-> " << x;
  			}
    		std::cout<<std::endl;
  		}
  		return;
  	}

  	void printGraph(std::string file) 
	{
  		std::cout<<"numVerts in printGraph: "<<numVerts<<std::endl;
		std::string to_file = "";
    	for (int d = 0; d < numVerts; ++d) 
    	{
      		// std::cout<<"d: "<<d<<std::endl;
      		// std::cout<<"what is this: "<< "Vertex " + std::to_string(d) + ":"<<std::endl;
    		to_file += "\n Vertex " + std::to_string(d) + ":";
    		for (auto x : adj[d])
      		{
      			to_file += "-> " + std::to_string(x);
  			}
    		to_file += '\n';
  		}
  		// std::cout<<"HERERERERERER0"<<std::endl;
		std::ofstream write;
  		// std::cout<<"HERERERERERER1"<<std::endl;
		write.open(file,std::ofstream::out);
  		// std::cout<<"HERERERERERER2"<<std::endl;
		write << to_file;
  		// std::cout<<"HERERERERERER3"<<std::endl;
		write.close();
  		// std::cout<<"HERERERERERER4"<<std::endl;
  		return;
  	}


  	void printVector(std::vector<int> vec)
  	{
  		for (auto it = begin(vec); it != end(vec); ++it) 
		{
			std::cout<<*it<<std::endl;
		}	
		return;
  	}

  	std::vector<int> nNearestNeighbors(int start, int numNeighbor)
  	{
  		std::vector<int> NN;
  		std::queue<int> q;
  		int n = 0;
  		q.push(start);
  		int safety_net = 100000;
  		int i = 0;
  		while (!q.empty())
  		{
  			i++;
  			int front = q.front();
  			q.pop();
  			for (auto it = begin(adj[front]); it != end(adj[front]); ++it) 
			{	
				if (*it != front and *it != start)
				{
					if (std::find(NN.begin(), NN.end(), *it) == NN.end()) 
					{
						q.push(*it);
						// if (std::find(NN.begin(), NN.end(), start) == NN.end()) 
						// {
						NN.push_back(*it);
						n++;
						if (n >= numNeighbor)
						{
							return NN;
						}
						// }
					}
				}
			}
  		}
  		return NN;
  	}


	// int getNodeIndex(node n);
	// {
	// 	for (auto it = begin (nodes); it != end (nodes); ++it) 
	// 	{
	// 		if (it -> ball_index == n.ball_index)
	// 		{
	// 			return node
	// 		}
	// 	}
	// 	return -1
	// }
	~graph() { delete [] adj; }
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