#include <iostream>
#include "graph_group.hpp"

class wrapper
{
public:
	graph *g = nullptr;
	
	wrapper() = default;

	wrapper(int num)
	{
		g = new graph(num);
	}

	void addEdge()
	{
		g -> addEdge(0,1);
		g -> addEdge(1,2);
		g -> addEdge(1,3);
		g -> addEdge(1,4);
		g -> addEdge(1,9);
		g -> addEdge(3,4);
		g -> addEdge(3,5);
		g -> addEdge(3,6);
		g -> addEdge(3,7);
		g -> addEdge(5,8);
		g -> addEdge(5,9);

		g -> addEdge(10,11);
		// g -> addEdge(11,12);
		g -> addEdge(11,13);
		g -> addEdge(11,14);
		g -> addEdge(11,19);
		g -> addEdge(13,14);
		g -> addEdge(13,15);
		g -> addEdge(13,16);
		g -> addEdge(13,17);
		g -> addEdge(15,18);
		g -> addEdge(15,19);

		g -> findGroups();
	}

	~wrapper()
	{
		delete g;
	}
	
};

void test_NN()
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

	g.printVector(g.nNearestNeighbors(3,20));

}

void test_delete()
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

}

void test_groups()
{
	graph *g = nullptr;

	g = new graph(20);

	g -> addEdge(0,1);
	g -> addEdge(1,2);
	g -> addEdge(1,3);
	g -> addEdge(1,4);
	g -> addEdge(1,9);
	g -> addEdge(3,4);
	g -> addEdge(3,5);
	g -> addEdge(3,6);
	g -> addEdge(3,7);
	g -> addEdge(5,8);
	g -> addEdge(5,9);

	g -> addEdge(10,11);
	// g -> addEdge(11,12);
	g -> addEdge(11,13);
	g -> addEdge(11,14);
	g -> addEdge(11,19);
	g -> addEdge(13,14);
	g -> addEdge(13,15);
	g -> addEdge(13,16);
	g -> addEdge(13,17);
	g -> addEdge(15,18);
	g -> addEdge(15,19);

	g -> findGroups();
	g -> printGraph();
	g -> printGroups();
}

void test_wrapper_methods()
{
	std::cout<<"Start wrapper test"<<std::endl;
	wrapper w(20);
	w.addEdge();

	
	w.g -> printGraph();
	w.g -> printGroups();
	std::cout<<"END TEST"<<std::endl;
}

int main(int argc, char const *argv[])
{
	test_wrapper_methods();
	// test_NN();

	return 0;
}

