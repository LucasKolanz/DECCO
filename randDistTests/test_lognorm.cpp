#include <iostream>
#include <fstream>
#include <math.h>
#include "../utilities/Utils.hpp"


int main()
{
	double lnSigma = 0.2; //in cm
	double a_max = std::pow(10,-5)*std::exp(-5*std::pow(lnSigma,2)/2); //in cm
	double accumulate = 0.0;
	double density = 2.25; //(in g/cm^3)
	int nums = 1000;
	std::ofstream file;
	// file.open("lognorm_output.csv");
	for (int i = 0; i < nums; i++)
	{
		accumulate += std::pow(lognorm_dist(a_max,lnSigma),3);
		if (i%100==0)
		{
		std::cerr<<i<<std::endl;
		}
	}

	std::cerr<<"constant mass: "<< (4/3) * M_PI*std::pow(std::pow(10,-5),3) * density<<std::endl;
	std::cerr<<"average lognorm mass: "<< (4/3) * M_PI*accumulate/nums * density<<std::endl;

	return 0;
}