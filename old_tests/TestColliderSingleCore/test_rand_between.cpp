#include "../Utils.hpp"
#include <iostream>


int main()
{
	// Ball_group testBG(20, true, v_custom); // dust formation constructor

	// safetyChecks()

    // int inputs = 100;
	int tries = 10;

	double rand_z;
	double rand_y;
	double radius = 1.0;
    // double input_radius[inputs];
    double output_radius[tries];
    // double avg_output_radius[inputs];


	for (int i = 0; i < tries; i++)
	{

    	rand_z = rand_between(-radius, radius);
    	rand_y = rand_between(-radius, radius);
        output_radius[i] = std::sqrt(std::pow(rand_z, 2) + std::pow(rand_y, 2));

    }

    // double row_sum;
    // for (int i = 0; i < inputs; i++)
    // {
    //     row_sum = 0;
    //     for (int j = 0; j < tries; j++)
    //     {
    //         row_sum += output_radius[i][j];
    //     }
    //     avg_output_radius[i] = row_sum/inputs;
    // }



    std::cerr<<"input radii: "<<radius<<std::endl;
    // for (int i = 0; i < inputs; i++)
    // {
    //     std::cerr<<input_radius[i]<<", ";
    // } 
    // std::cerr<<std::endl;
    std::cerr<<"output radii:"<<std::endl;
    for (int i = 0; i < tries; i++)
    {
        std::cerr<<output_radius[i]<<", ";
    } 

    std::cout<<RAND_MAX<<std::endl;

}


