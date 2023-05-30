#include <random>
#include <iostream>
#include <unistd.h>
#include "test_add.hpp"
#include "vec3.hpp"


int main()
{
	

	gen.seed(100);

	// std::normal_distribution<double> d(0, 1);
    // std::cout<<d(rand_gen)<<std::endl;
    // std::cout<<d(rand_gen)<<std::endl;
    // std::cout<<d(rand_gen)<<std::endl;
	std::cout<<rand_unit_vec3()<<std::endl;
	std::cout<<rand_unit_vec3()<<std::endl;
	std::cout<<rand_unit_vec3()<<std::endl;
	
	// rand_gen.seed(100);
	// std::normal_distribution<double> dd(0, 1);

	// std::cout<<dd(rand_gen)<<std::endl;
    // std::cout<<dd(rand_gen)<<std::endl;
    // std::cout<<dd(rand_gen)<<std::endl;
	// std::cout<<rand_unit_vec3(rand_gen)<<std::endl;
	// std::cout<<rand_unit_vec3(rand_gen)<<std::endl;
	// std::cout<<rand_unit_vec3(rand_gen)<<std::endl;
}