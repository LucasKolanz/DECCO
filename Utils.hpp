#pragma once
#include "vec3.hpp"
#include "linalg.hpp"
#include <random>
#include <cmath>
#include <iostream>

using namespace linalg;
using linalg::cross;
using linalg::normalize;
using linalg::aliases::double3;
using linalg::aliases::double3x3;

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

std::random_device rd;
std::mt19937 gen(rd());
double
random_gaussian(const double mean = 0, const double standard_deviation = 1)
{
    std::normal_distribution<double> d(mean, standard_deviation);
    return d(gen);
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
// 	std::cerr << value;
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
    return min + f * (max - min);
}

// Returns a random unit vector.
vec3
rand_unit_vec3()
{
    return vec3(random_gaussian(), random_gaussian(), random_gaussian()).normalized();
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