#include <iostream>
#include <fstream>
#include "../Utils.hpp"


int main()
{
	double a_max = 2;
	double sigma = 0.2;
	int nums = 100000;
	std::ofstream file;
	file.open("lognorm_output.csv");
	for (int i = 0; i < nums; i++)
	{
		if (i < nums-1)
		{
			file << lognorm_dist(a_max,sigma) << ",";
		}
		else
		{
			file << lognorm_dist(a_max,sigma);
		}
	}
	file.close();
	return 0;
}