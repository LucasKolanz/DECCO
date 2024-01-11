#pragma once
#include <cmath>
#include <string>
// #include <numbers>

////////////////////////////////////
// Initialization (Units are CGS) //
////////////////////////////////////

// T(K)	v(cm / s)
// 3    0.36301555459799423
// 10	0.6627726932618989
// 30	1.1479559785988718
// 100	2.0958712816717324
// 300	3.6301555459799424
// 1000	6.627726932618988

constexpr double Kb = 1.380649e-16; //in erg/K
constexpr double pi = 3.14159265358979311599796346854;


bool dynamicTime;


// using std::numbers::pi;
double G;  // Gravitational constant
double density;
double u_s;                // Coeff of sliding friction
double u_r;               // Coeff of rolling friction
double sigma;              // Poisson ratio for rolling friction.
double Y;               // Young's modulus in erg/cm3
double cor;                // Coeff of restitution
double simTimeSeconds;  // Seconds
double timeResolution;    // Seconds - This is duration between exported steps.
double fourThirdsPiRho;  // for fraction of smallest sphere radius.
double scaleBalls;                         // base radius of ball.
double maxOverlap;                           // of scaleBalls
double KEfactor;                              // Determines collision velocity based on KE/PE
double v_custom;  // Velocity cm/s
double temp;          //tempurature of simulation in Kelvin
double kConsts;
double impactParameter;  // Impact angle radians
double Ha;         // Hamaker constant for vdw force
double h_min;  // 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
double cone;  // Cone of particles ignored moving away from center of mass. Larger angle ignores more.

// Simulation Structure
int properties;  // Number of columns in simData file per ball
int genBalls;
int attempts;  // How many times to try moving every ball touching another in generator.


double spaceRange;  // Rough minimum space required
double spaceRangeIncrement;
double z0Rot;  // Cluster one z axis rotation
double y0Rot;  // Cluster one y axis rotation
double z1Rot;  // Cluster two z axis rotation
double y1Rot;  // Cluster two y axis rotation
double simTimeElapsed;

// File from which to proceed with further simulations
// std::string project_path;
// std::string output_folder;
// std::string projectileName;
// std::string targetName;
// std::string output_prefix;
