#pragma once
#include <cmath>
#include <string>
#include <numbers>

////////////////////////////////////
// Initialization (Units are CGS) //
////////////////////////////////////

constexpr bool dynamicTime = false;

// T(K)	v(cm / s)
// 3    0.36301555459799423
// 10	0.6627726932618989
// 30	1.1479559785988718
// 100	2.0958712816717324
// 300	3.6301555459799424
// 1000	6.627726932618988

using std::numbers::pi;
constexpr double G = 6.67e-8;  // Gravitational constant
constexpr double density;
constexpr double u_s;                // Coeff of sliding friction
constexpr double u_r;               // Coeff of rolling friction
constexpr double sigma;              // Poisson ratio for rolling friction.
constexpr double Y;               // Young's modulus in erg/cm3
constexpr double cor;                // Coeff of restitution
constexpr double simTimeSeconds;  // Seconds
constexpr double timeResolution;    // Seconds - This is duration between exported steps.
constexpr double fourThirdsPiRho;  // for fraction of smallest sphere radius.
constexpr double scaleBalls;                         // base radius of ball.
constexpr double maxOverlap;                           // of scaleBalls
constexpr double KEfactor;                              // Determines collision velocity based on KE/PE
constexpr double v_custom;            // Velocity cm/s
constexpr double kConsts;
constexpr double impactParameter;  // Impact angle radians
constexpr double Ha;         // Hamaker constant for vdw force
constexpr double h_min;  // 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
constexpr double cone;  // Cone of particles ignored moving away from center of mass. Larger angle ignores more.

// Simulation Structure
constexpr int properties = 11;  // Number of columns in simData file per ball
constexpr int genBalls;
constexpr int attempts;  // How many times to try moving every ball touching another in generator.

int skip;  // Steps thrown away before recording a step to the buffer. 500*.04 is every 20 seconds in sim.
int steps;

double dt;
double kin;  // Spring constant
double kout;
double spaceRange;  // Rough minimum space required
double spaceRangeIncrement;
double z0Rot;  // Cluster one z axis rotation
double y0Rot;  // Cluster one y axis rotation
double z1Rot;  // Cluster two z axis rotation
double y1Rot;  // Cluster two y axis rotation
double simTimeElapsed;

// File from which to proceed with further simulations
std::string path;
std::string output_folder;
std::string projectileName;
std::string targetName;
std::string output_prefix;
