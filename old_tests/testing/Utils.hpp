#pragma once
#include "vec3.hpp"
#include "linalg.hpp"
#include <random>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>

using namespace linalg;
using linalg::cross;
using linalg::normalize;
using linalg::aliases::double3;
using linalg::aliases::double3x3;

std::random_device rd;
std::mt19937 random_generator(rd());

// std::string rand_file = "/global/homes/l/lpkolanz/SpaceLab/random/rand.txt";
// std::string gaus_file = "/global/homes/l/lpkolanz/SpaceLab/random/randomGaussian.txt";

// std::ifstream rando(rand_file,std::ios::in);
// std::ifstream gauss(gaus_file,std::ios::in);

double get_gaus()
{
    double ret;
    std::string line;
    if (!std::getline(gauss,line))
    {
        gauss.close();
        std::ifstream gauss(gaus_file,std::ios::in);
        std::getline(gauss,line);
    }
    std::cerr<<"GAUS RAND: "<<line<<std::endl;
    std::stringstream ss;  
    ss << line;  
    ss >> ret; 
    return ret;
}

int get_rand()
{
    int ret;
    std::string line;
    if (!std::getline(rando,line))
    {
        rando.close();
        std::ifstream rando(rand_file,std::ios::in);
        std::getline(rando,line);
    }
    std::cerr<<"RAND RAND: "<<line<<std::endl;
    std::stringstream ss;  
    ss << line;  
    ss >> ret; 
    return ret;
}

// Convert from vec3 to double3
double3
to_double3(const vec3& vec)
{
    return {vec.x, vec.y, vec.z};
}

// Convert from double3 to vec3
vec3
to_vec3(const double3& vec)
{
    return {vec.x, vec.y, vec.z};
}

double
random_gaussian(const double mean = 0, const double standard_deviation = 1)
{
    std::normal_distribution<double> d(mean, standard_deviation);
    return d(random_generator);
}

// I got this from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
bool
line_sphere_intersect(const vec3& position, const vec3& velocity, const vec3& center, const double radius)
{
    vec3 direction = velocity.normalized();

    double grad = direction.dot((position - center)) * direction.dot((position - center)) -
                  ((position - center).normsquared() - radius * radius);

    // If I ever want the distance to the two intersections:
    // double d1 = -u_vec.dot(origin - center) + sqrtf(grad);
    // double d2 = -u_vec.dot(origin - center) - sqrtf(grad);
    // std::cout << '\t' << d1 << '\t' << d2 << '\n';

    if (grad < 0) {
        // std::cout << "Position"; position.print();
        // std::cout << "Direction"; direction.print();
        // std::cout << "Target"; center.print();
        // std::cout << "Radius " << radius << '\n';
        // std::cout << "Grad " << grad << '\n';
        return false;
    } else {
        // std::cout << "Intersection exists.";
        return true;
    }
}


// Generate a vector orthogonal to the given
double3x3
local_coordinates(const double3& x)
{
    const auto& x_hat = normalize(x);
    const auto& y_hat = cross(x_hat, double3(0, 0, 1));
    const auto& z_hat = cross(y_hat, x_hat);

    return double3x3(x_hat, normalize(y_hat), normalize(z_hat));
}


vec3
perpendicular_shift(const double3x3 local_basis, const double y, const double z)
{
    return to_vec3(local_basis.y * y + local_basis.z * z);

    // const auto basis33_transpose = transpose(local_basis);
    // return to_vec3(mul(basis33_transpose, double3(0, y, z)));
}

void
print_m33(double3x3& m33)
{
    for (const auto& row : m33) {
        for (const auto& el : row) { std::cout << '\t' << el; }
        std::cout << "\n";
    }
}

std::string
vec_string(const double3 vec)
{
    return '(' + std::to_string(vec.x) + ',' + std::to_string(vec.y) + ',' + std::to_string(vec.z) + ')';
}

// Rounding
inline std::string
rounder(double value, int digits)
{
    return std::to_string(value).substr(0, digits);
}

// Scientific Notation
inline std::string
scientific(double value)
{
    std::stringstream ss;
    ss << std::setprecision(0) << std::scientific << value;
    return ss.str();
}

// Output a nice title bar in terminal:
inline void
titleBar(std::string title)
{
    std::cerr << '\n';
    for (size_t i = 0; i < ((62 - title.size()) / 2); i++) { std::cerr << '='; }
    std::cerr << ' ' << title << ' ';
    for (size_t i = 0; i < ((62 - title.size()) / 2); i++) { std::cerr << '='; }
    std::cerr << "\n\n";
}

// // Print anything:
// template <typename theType>
// void print(theType value)
// {
//  std::cerr << value;
// }

// Ask a yes or no question:
inline bool
input(const std::string& question)
{
    char answer;
    std::cerr << question;
    std::cin >> answer;
    if (answer == 'y') {
        return true;
    } else {
        return false;
    }
}

// inline double
// lin_rand()

inline double
rand_between(const double min, const double max)
{
    double f = (double)rand() / RAND_MAX;
    // double f = (double)get_rand() / RAND_MAX;
    return min + f * (max - min);
}

// Returns a random unit vector.
vec3
rand_unit_vec3()
{
    return vec3(random_gaussian(), random_gaussian(), random_gaussian()).normalized();
    // return vec3(get_gaus(), get_gaus(), get_gaus()).normalized();
}

// // Returns a vector within the desired radius, and optionally outside an inner radius (shell).
vec3
rand_vec3(double outer_radius, double inner_radius = 0)
{
    double r3 = outer_radius * outer_radius * outer_radius;
    double r = cbrt(rand_between(inner_radius, r3));
    return rand_unit_vec3() * r;
}

// @brief - returns maxwell boltzmann probability density function value
//          @param x where @param a = sqrt(K*T/m) where K is boltzmanns constant
//          T is tempurature and m is mass.
double mbdpdf(double a, double x)
{
    return std::sqrt(2/M_PI)*(std::pow(x,2)/std::pow(a,3))*std::exp(-(std::pow(x,2))/(2*std::pow(a,2)));
}

// @brief - returns a velocity from the maxwell boltzmann distribution given 
//          @param a, which is the same as @param a from mbdpdf()
double max_bolt_dist(double a)
{
    double v0,Fv0,sigma,test;
    double maxVal;
    
    sigma = a*std::sqrt((3*M_PI-8)/M_PI);
    maxVal = mbdpdf(a,a*M_SQRT2);

    do
    {
        Fv0 = rand_between(0,20*M_SQRT2*a*sigma);
        v0 = Fv0/(a*M_SQRT2);
        test = rand_between(0,maxVal);
    }while(test > mbdpdf(a,v0));
    
    return v0;
}

double lndpdf(double a,double sigma,double a_max)
{
    //M_2_SQRTPI is 2/sqrt(pi)
    return M_2_SQRTPI/(a*sigma*2*M_SQRT2)*
            std::exp(-std::pow(log(a/a_max)-std::pow(sigma,2),2)/(2*std::pow(sigma,2)));
}

double lognorm_dist(double a_max,double sigma)
{
    double Fa0,a0,test,maxVal;

    maxVal = lndpdf(a_max,sigma,a_max);

    do
    {
        Fa0 = rand_between(0,20*a_max*sigma);
        a0 = Fa0/(a_max);
        test = rand_between(0,maxVal);
    }while(test > lndpdf(a0,sigma,a_max));
    
    return a0;
}