#include <iostream>
#include "graph_group.hpp"


int main(int argc, char const *argv[])
{
	graph g(10);

	g.addEdge(0,1);
	g.addEdge(1,2);
	g.addEdge(1,3);
	g.addEdge(1,4);
	g.addEdge(1,9);
	g.addEdge(3,4);
	g.addEdge(3,5);
	g.addEdge(3,6);
	g.addEdge(3,7);
	g.addEdge(5,8);
	g.addEdge(5,9);


	// std::vector<int> test = g.nNearestNeighbors(1,5);

	// g.printVector(test);
	g.printGraph();


	g.deleteEdge(1,9);
	std::cout<<"=========================================="<<std::endl;

	g.printGraph();


	return 0;
}

