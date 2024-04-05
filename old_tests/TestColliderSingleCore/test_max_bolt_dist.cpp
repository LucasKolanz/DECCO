#define _USE_MATH_DEFINES
#include "../Utils.hpp"
#include <iostream>
#include <fstream>
// #include <random>
#include <vector>
#include <cmath>
#include <sstream>


double mbdpdf(double a, double x)
{
    return std::sqrt(2/M_PI)*(std::pow(x,2)/std::pow(a,3))*std::exp(-(std::pow(x,2))/(2*std::pow(a,2)));
}

double max_bolt_dist(double a)
{
    double v0,Fv0,sigma,test;
    double maxVal;
    
    sigma = a*std::sqrt((3*M_PI-8)/M_PI);
    maxVal = mbdpdf(a,a*M_SQRT2);
    // maxVal = 1.0;

    // std::cout<<"cpp max: "<<maxVal<<std::endl;

    do
    {
        Fv0 = rand_between(0,20*M_SQRT2*a*sigma);
        v0 = Fv0/(a*M_SQRT2);
        test = rand_between(0,maxVal);
    }while(test > mbdpdf(a,v0));
    
    return v0;
}

int main(int argc, char **argv)
{
    std::istringstream ss1(argv[1]);
    std::istringstream ss2(argv[2]);
    int tries;
    double a;
    ss1 >> tries;
    ss2 >> a;

    

    std::vector<double> output(tries);

    for (int i = 0; i < tries; i++)
    {
        output[i] = max_bolt_dist(a);
    }

    std::ofstream myfile;
    myfile.open("/home/lucas/Desktop/Research/SpaceLabTesting/SpaceLab/output/dump/test_max_bolt_dist_out.csv");
    for (int i = 0; i < tries-1; i++)
    {
        myfile << output[i] << ",";
    }
    myfile << output[tries-1];
    myfile.close();

    return 0;
}


