#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <numeric>
#include "vec3.hpp"

class graph
{
public:
	bool initialized = false;
	int numVerts;
	std::vector<int> *adj = nullptr;
	int *groups = nullptr;

	// bool *visited = nullptr;
	// std::

	void resetGroups()
	{
		for (int i = 0; i < numVerts; i++)
		{
			groups[i] = -1;
		}
	}

	graph() = default;

	graph(int num_balls)
	{
		numVerts = num_balls;
		adj = new std::vector<int>[num_balls];

		//To see what group a ball belongs to
		//index into groups with ball index and 
		//the value is what group that ball belongs to.
		//Groups also acts as a visited list.
		//If group[index] = -1 then it hasnt been visited
		groups = new int[num_balls];
		findGroups();
		initialized = true;
		return;
	}

	~graph() 
	{
		delete[] adj;
		delete[] groups;
 	}

	//Copy constructor
	graph(const graph& copyme)
	{
		numVerts = copyme.numVerts;
		adj = copyme.adj;
		groups = copyme.groups;
		initialized = copyme.initialized;
	}

	//Assignment operator
	graph& operator=(const graph& src)
	{
		if (this != &src)
		{
			groups = src.groups;
			adj = src.adj;
			numVerts = src.numVerts;
			initialized = src.initialized;
		}
		return *this;
	}

	// int sumVisited()
	// {
	// 	int accum = 0;
	// 	for (int i = 0; i < numVerts; i++)
	// 	{
	// 		accum += visited[i];
	// 	}
	// 	return accum;
	// }

	int notVisited()
	{
		for (int i = 0; i < numVerts; i++)
		{
			if (groups[i] == -1)
			{
				return i;
			}
		}
		return -1;
	}

	void findGroups(int start=0)
	{
  		resetGroups();
		int curr = start;
		int group_id = 0;
  		do {
  			getWalk(curr, group_id);
  			group_id++;
  			curr = notVisited();
  			// printGroups();
  		}while(curr != -1);
  		return;
	}

	void getWalk(int start,int group_id)
	{
  		std::queue<int> q;
		q.push(start); 
  		while (!q.empty())
  		{
  			int front = q.front();
  			q.pop();
  			if (adj[front].size() == 0)
  			{
  				groups[front] = group_id;
  			}
  			else
  			{
	  			for (auto it = begin(adj[front]); it != end(adj[front]); ++it) 
				{	

					if (groups[*it] == -1) 
					{
						q.push(*it);
						groups[*it] = group_id;
					}
				}
			}
		}
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
	
	void printGroups()
	{
		std::cout<<"================Groups================"<<std::endl;
		for (int i = 0; i < numVerts; i++)
		{
			std::cout<<"(ball; group): "<<i<<"; "<<groups[i]<<std::endl;
		}
		std::cout<<"======================================"<<std::endl;
		return;
	}

	void printGraph() 
	{
		std::cout<<"================Graph================"<<std::endl;
		if (initialized)
		{
	    	for (int d = 0; d < numVerts; ++d) 
	    	{
	    		std::cout << "\n Vertex "<< d << ":";
	    		for (auto x : adj[d])
	      		{
	      			std::cout << "-> " << x;
	  			}
	    		std::cout<<std::endl;
	  		}
	  	}
	  	else
	  	{
	  		std::cout<<"ERROR: graph not initalized"<<std::endl;
	  	}
		std::cout<<"====================================="<<std::endl;
  		return;
  	}

  	void printGraph(std::string file) 
	{
		if (initialized)
		{
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
			std::ofstream write;
			write.open(file,std::ofstream::out);
			write << to_file;
			write.close();
		}
	  	else
	  	{
	  		std::cout<<"ERROR: graph not initalized"<<std::endl;
	  	}
  		return;
  	}


  	void printVector(std::vector<int> vec)
  	{
		std::cout<<"================Vector================"<<std::endl;

  		for (auto it = begin(vec); it != end(vec); ++it) 
		{
			std::cout<<*it<<", ";
		}	
		std::cout<<std::endl;
		std::cout<<"======================================"<<std::endl;
		return;
  	}

  	// void mergeGraphs(const graph src)

  	std::vector<int> nNearestNeighbors(int start, int numNeighbor)
  	{

  		//If there are no balls touching the start ball
  		//TODO: return group information for non-contact forces
  		//At the moment: return all balls besides start ball 
  		std::vector<int> NN;
  		if (adj[start].size() == 0)
  		{
  			NN.resize(numVerts);
  			std::iota (std::begin(NN),std::end(NN),0);
  			NN.erase(NN.begin() + start);
  			// std::cout<<start<<": "<<numVerts<<std::endl;
  			// printVector(NN);
  			// exit(0);
  		}
  		else
  		{
  			for (auto it = begin(adj[start]); it != end(adj[start]); ++it)
  			{
  				if (it != begin(adj[start]))
  				{
  					NN.push_back(*it);
  				}
  			}
  		// 	int n = 0;
	  	// 	std::queue<int> q;
	  	// 	q.push(start);
	  	// 	while (!q.empty())
	  	// 	{
	  	// 		int front = q.front();
	  	// 		q.pop();
	  	// 		for (auto it = begin(adj[front]); it != end(adj[front]); ++it) 
				// {	
				// 	if (*it != front and *it != start)
				// 	{
				// 		if (std::find(NN.begin(), NN.end(), *it) == NN.end()) 
				// 		{
				// 			q.push(*it);

				// 			NN.push_back(*it);
				// 			n++;
				// 			if (n >= numNeighbor)
				// 			{
				// 				return NN;
				// 			}
				// 		}
				// 	}
				// }
	  		// }
  		}
  		return NN;
  	}

};







// class graph_group
// {
// public:

// 	int groups = -1;
// 	graph *g = nullptr;

// 	graph_group() = default;

// 	void init(int *adj_matrix, int num_particles)
// 	{

// 		// g = make_graphs()
// 	}

// 	~graph_group()
// 	{
// 		delete[] g;
// 	}
	
// };