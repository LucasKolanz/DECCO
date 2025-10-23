#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

// Random uniform in [a,b]
double rand_between(double a, double b)
{
    return a + (b - a) * (rand() / (double)RAND_MAX);
}

double lndpdf(double a,double sigma,double a_max)
{
    //M_2_SQRTPI is 2/sqrt(pi)
    return M_2_SQRTPI/(a*sigma*2*M_SQRT2)*
            std::exp(-std::pow(log(a/a_max)-std::pow(sigma,2),2)/(2*std::pow(sigma,2)));
}

//a_max is the value with the max probability (the mode)
//sigma is the standard deviation
//MPIsafe
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

int main()
{
    std::srand((unsigned int)std::time(nullptr));  // seed RNG

    const int Nsamples = 100000;
    // const double sigma = 0.2;        // example spread
    const double sigma = 1.0;        // example spread
    // const double scaleBalls = 1e-5;   // base radius (arbitrary units)
    // const double a_max = scaleBalls * std::exp(-5.0 * sigma * sigma / 2.0);
    const double a_max = 5.0;

    std::vector<double> radii;
    radii.reserve(Nsamples);

    for (int i = 0; i < Nsamples; ++i)
    {
        std::cerr<<i<<std::endl;
        radii.push_back(lognorm_dist(a_max, sigma));
    }

    // Write to file
    std::ofstream outfile("lognorm_samples.txt");
    if (!outfile)
    {
        std::cerr << "Error: Could not open output file.\n";
        return 1;
    }

    for (double r : radii)
        outfile << r << "\n";

    outfile.close();
    std::cout << "Wrote " << Nsamples << " samples to lognorm_samples.txt\n";
    return 0;
}