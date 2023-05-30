#include <random>
#include <cmath>
#include <iostream>
#include "vec3.hpp"

std::random_device rd;
std::mt19937 gen(rd());
double
random_gaussian(const double mean = 0, const double standard_deviation = 1)
{
    std::normal_distribution<double> d(mean, standard_deviation);
    return d(gen);
}


vec3
rand_unit_vec3()
{
    return vec3(random_gaussian(), random_gaussian(), random_gaussian()).normalized();
}