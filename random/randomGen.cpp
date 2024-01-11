#include <iostream>
#include <fstream>
#include <random>
#include "../Utils.hpp"

int main()
{

	// std::random_device rd;
	// std::mt19937 random_generator(rd());
	int n = 1000000;
	srand(static_cast<int>(time(nullptr)));
	std::ofstream G("randomGaussian.txt",std::ios::app);
	std::ofstream R("rand.txt",std::ios::app);

	for (int i = 0; i < n; i++)
	{
		G<<random_gaussian()<<std::endl;
		R<<rand()<<std::endl;
	}

	// std::cout<<random_gaussian()<<std::endl;

}