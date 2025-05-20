#pragma once

#ifndef UTILS_HPP
#define UTILS_HPP

#include "../external/json/single_include/nlohmann/json.hpp"
#include "linalg.hpp"
#include "vec3.hpp"

using json = nlohmann::json;

using namespace linalg;
using linalg::cross;
using linalg::normalize;
using linalg::aliases::double3;
using linalg::aliases::double3x3;

//if the attribute_key exists in input, set variable based on input[attribute_key]
template <typename T>
void set_attribute(const json& input, const std::string &attribute_key, T &variable)
{
    if (input.contains(attribute_key))
    {
        variable = input[attribute_key];
    }
    else
    {
        std::cerr<<"WARNING: attribute '"<<attribute_key<<"' does not exist."<<std::endl;
    }
}

void seed_generators(size_t seed);

bool isAllDigits(const std::string& s);
//Returns all the folders in a particular directory
std::vector<std::string> get_folders_in_directory(const std::string directory);
//@brief: Given a folder base in the form:
//          /*Global path to Spacelab_data folder*/SpaceLab_data/jobs/jobsetname{a}/N_{n}/eta_{e}/T_{t}/
//      This function finds all possible values of {a}, {n}, {e}, and/or {t}
//      and returns a random folder with this scheme. If there isn't a specified 
//      value for one of these options, then it is ignored as a possible random value.
//The order of folders after the job folder shouldn't matter to this function
//This function does not check if the job in the random folder is complete or not
std::string get_rand_projectile_folder(std::string folder);
nlohmann::json getJsonFromFolder(std::string location);
int extractNumberFromString(const std::string& s);
std::string dToSci(double value);
std::string vToSci(vec3 value);
// Convert from vec3 to double3
double3 to_double3(const vec3& vec);
// Convert from double3 to vec3
vec3 to_vec3(const double3& vec);
double random_gaussian(const double mean = 0, const double standard_deviation = 1);
// I got this from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
bool line_sphere_intersect(const vec3& position, const vec3& velocity, const vec3& center, const double radius);
// Generate a vector orthogonal to the given
double3x3 local_coordinates(const double3& x);
vec3 perpendicular_shift(const double3x3 local_basis, const double y, const double z);
void print_m33(double3x3& m33);
std::string vec_string(const double3 vec);
// Rounding
std::string rounder(double value, int digits);
// Scientific Notation
std::string scientific(const double value);
std::string scientific(const vec3 value);
// Output a nice title bar in terminal:
void titleBar(const std::string title);
// Ask a yes or no question:
bool input(const std::string& question);
double rand_between(const double min, const double max);
int rand_int_between(const int min, const int max);
// Returns a random unit vector.
vec3 rand_unit_vec3();
// // Returns a vector within the desired radius, and optionally outside an inner radius (shell).
vec3 rand_vec3(double outer_radius, double inner_radius = 0);
// @brief - returns maxwell boltzmann probability density function value
//          @param x where @param a = sqrt(K*T/m) where K is boltzmanns constant
//          T is tempurature and m is mass.
double mbdpdf(double a, double x);
// @brief - returns a velocity from the maxwell boltzmann distribution given 
//          @param a, which is the same as @param a from mbdpdf()
double max_bolt_dist(double a);
double lndpdf(double a,double sigma,double a_max);
double lognorm_dist(double a_max,double sigma);

#endif