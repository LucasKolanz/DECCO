#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
////////////////////////////////////
// Initialization (Units are CGS) //
////////////////////////////////////

constexpr bool dynamicTime = true;

constexpr double
G = 6.67e-8,   // Gravitational constant
density = 2.7, //2.7, // Typically based on some rock density
u_s = 0.9,      // Coeff of friction
u_r = 0.9,      // Coeff of friction
sigma = .29, // Poisson ratio for rolling friction.
Y = 2.0e12, // Young's modulus in erg/cm3
Ha = 0, // Hamaker constant for vdw force
h_min = 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
cor = 0.8,     // Coeff of restitution
simTimeSeconds = 18000., // Seconds
timeResolution = 20., // Seconds - This is duration between exported steps
maxOverlap = .1,
fourThirdsPiRho = 4. / 3. * M_PI * density, // for fraction of smallest sphere radius.
kConsts = fourThirdsPiRho / (maxOverlap * maxOverlap),
scaleBalls = 100, // base radius of balls
KEfactor = 0,       // Determines collision velocity based on KE/PE
v_custom = 1,
impactParameter = 0, // Impact angle radians
cone = M_PI_2; // Cone of paritlces ignored moving away from center of mass. Larger angle ignores more.

// Simulation Structure
inline const unsigned int
properties = 11, // Number of columns in simData file per ball
genBalls = 6,
attempts = 200; // How many times to try moving every ball touching another in generator.

inline size_t
skip = 0,     // Steps thrown away before recording a step to the buffer. 500*.04 is every 20 seconds in sim.
steps = 0;

inline double
dt = -1,
kin = -1,      // Spring constant
kout = -1,
spaceRange = std::pow((1. / .74 * scaleBalls * scaleBalls * scaleBalls * genBalls), 1. / 3.), // Rough minimum space required
spaceRangeIncrement = scaleBalls * 3,
z0Rot = 0,           // Cluster one z axis rotation
y0Rot = 0,           // Cluster one y axis rotation
z1Rot = 0,           // Cluster two z axis rotation
y1Rot = 0,           // Cluster two y axis rotation
simTimeElapsed = 0;

// File from which to proceed with further simulations
inline std::string
path = "/home/lucas/Desktop/Research/SpaceLabTesting/SpaceLab",
projectileName = "cohTest_",
targetName = "cohTest_",
output_prefix = "Unnamed_";
