#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include "../vector3d.hpp"
#include "sphere_physics.h"

// String buffers to hold data in memory until worth writing to file:
std::stringstream ballBuffer;
std::stringstream energyBuffer;

// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
const time_t start = time(nullptr);        // For end of program analysis
time_t startProgress; // For progress reporting (gets reset)
time_t lastWrite;     // For write control (gets reset)
bool writeStep;       // This prevents writing to file every step (which is slow).

/// @brief The ballGroup run by the main sim looper.


// Prototypes
void simOneStep(const unsigned int& Step);
[[noreturn]] void simLooper();
void safetyChecks();

Cosmos O(path, projectileName, targetName, vCustom); // Collision
//ballGroup O(path, targetName, 0); // Continue
//ballGroup O(genBalls, true, vCustom); // Generate

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int main(const int argc, char const* argv[])
{
	energyBuffer.precision(12); // Need more precision on momentum.

	safetyChecks();
	O.simInitWrite(outputPrefix);
	simLooper();

} // end main
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void simOneStep(const unsigned int& Step)
{
	// Check if this is a write step:
	if (Step % skip == 0)
	{
		writeStep = true;

		simTimeElapsed += dt * skip;

		// Progress reporting:
		float eta = ((time(nullptr) - startProgress) / static_cast<float>(skip) * static_cast<float>(steps - Step)) / 3600.f; // Hours.
		float real = (time(nullptr) - start) / 3600.f;
		float simmed = simTimeElapsed / 3600.f;
		float progress = (static_cast<float>(Step) / static_cast<float>(steps) * 100.f);
		fprintf(stderr, "Step: %u\tProgress: %2.0f%%\tETA: %5.2lf hr\tReal: %5.2f hr\tSim: %5.2f hr\tReal/Sim: %5.2\n", Step, progress, eta, real, simmed, real / simmed);
		fprintf(stdout, "Step: %u\tProgress: %2.0f%%\tETA: %5.2lf hr\tReal: %5.2f hr\tSim: %5.2f hr\tReal/Sim: %5.2\t\r", Step, progress, eta, real, simmed, real / simmed);
		fflush(stdout);
		startProgress = time(nullptr);
	}
	else
	{
		writeStep = false;
	}

	
	/// FIRST PASS - Update Kinematic Parameters:
	for (unsigned int Ball = 0; Ball < O.n; Ball++)
	{
		// Update velocity half step:
		O.g[Ball].velh = O.g[Ball].vel + .5 * O.g[Ball].acc * dt;

		// Update angular velocity half step:
		O.g[Ball].wh = O.g[Ball].w + .5 * O.g[Ball].aacc * dt;

		// Update position:
		O.g[Ball].pos += O.g[Ball].velh * dt;

		// Reinitialize acceleration to be recalculated:
		O.g[Ball].acc = { 0, 0, 0 };

		// Reinitialize angular acceleration to be recalculated:
		O.g[Ball].aacc = { 0, 0, 0 };
	}

	/// SECOND PASS - Check for collisions, apply forces and torques:
	for (unsigned int A = 1; A < O.n; A++) //cuda
	{
		/// DONT DO ANYTHING HERE. A STARTS AT 1.
		for (unsigned int B = 0; B < A; B++)
		{
			const double sumRaRb = O.g[A].R + O.g[B].R;
			vec3 rVec = O.g[B].pos - O.g[A].pos; // Start with rVec from a to b.
			const double dist = (rVec).norm();
			vec3 totalForce;

			// Check for collision between Ball and otherBall:
			double overlap = sumRaRb - dist;

			// Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
			unsigned int e = static_cast<unsigned>(A * (A - 1) * .5) + B; // a^2-a is always even, so this works.
			double oldDist = O.distances[e];

			// Check for collision between Ball and otherBall.
			if (overlap > 0)
			{
				double k;
				// Apply coefficient of restitution to balls leaving collision.
				if (dist >= oldDist)
				{
					k = kout;
				}
				else
				{
					k = kin;
				}

				// Cohesion:
				// h is the "separation" of the particles at particle radius - maxOverlap.
				// This allows particles to be touching while under vdwForce.
				const double h = maxOverlap * 1.01 - overlap;
				const double Ra = O.g[A].R;
				const double Rb = O.g[B].R;
				const double h2 = h * h;
				const double twoRah = 2 * Ra * h;
				const double twoRbh = 2 * Rb * h;
				const vec3 vdwForce =
					Ha / 6 *
					64 * Ra * Ra * Ra * Rb * Rb * Rb *
					(h + Ra + Rb) /
					(
						(h2 + twoRah + twoRbh) *
						(h2 + twoRah + twoRbh) *
						(h2 + twoRah + twoRbh + 4 * Ra * Rb) *
						(h2 + twoRah + twoRbh + 4 * Ra * Rb)
						) *
					rVec.normalized();

				// Elastic a:
				vec3 elasticForce = -k * overlap * .5 * (rVec / dist);

				// Friction a:
				vec3 dVel = O.g[B].vel - O.g[A].vel;
				vec3 frictionForce = { 0, 0, 0 };
				const vec3 relativeVelOfA = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - O.g[A].w.cross(O.g[A].R / sumRaRb * rVec) - O.g[B].w.cross(O.g[B].R / sumRaRb * rVec);
				double relativeVelMag = relativeVelOfA.norm();
				if (relativeVelMag > 1e-10) // When relative velocity is very low, dividing its vector components by its magnitude below is unstable.
				{
					frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfA / relativeVelMag);
				}

				// Torque a:
				const vec3 aTorque = (O.g[A].R / sumRaRb) * rVec.cross(frictionForce);

				// Gravity on a:
				const vec3 gravForceOnA = (G * O.g[A].m * O.g[B].m / (dist * dist)) * (rVec / dist);

				// Total forces on a:
				totalForce = gravForceOnA + elasticForce + frictionForce + vdwForce;

				// Elastic and Friction b:
				// Flip direction b -> a:
				rVec = -rVec;
				dVel = -dVel;
				elasticForce = -elasticForce;

				const vec3 relativeVelOfB = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - O.g[B].w.cross(O.g[B].R / sumRaRb * rVec) - O.g[A].w.cross(O.g[A].R / sumRaRb * rVec);
				relativeVelMag = relativeVelOfB.norm(); // todo - This should be the same as mag for A. Same speed different direction.
				if (relativeVelMag > 1e-10)
				{
					frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfB / relativeVelMag);
				}
				const vec3 bTorque = (O.g[B].R / sumRaRb) * rVec.cross(frictionForce);

				O.g[A].aacc += aTorque / O.g[A].moi;
				O.g[B].aacc += bTorque / O.g[B].moi;


				if (writeStep)
				{
					// Calculate potential energy. Important to recognize that the factor of 1/2 is not in front of K because this is for the spring potential in each ball and they are the same potential.
					O.U += -G * O.g[A].m * O.g[B].m / dist + 0.5 * k * overlap * overlap;
				}
			}
			else
			{
				// No collision: Include gravity only:
				const vec3 gravForceOnA = (G * O.g[A].m * O.g[B].m / (dist * dist)) * (rVec / dist);
				totalForce = gravForceOnA;
				if (writeStep)
				{
					O.U += -G * O.g[A].m * O.g[B].m / dist;
				}

				// For expanding overlappers:
				//O.g[A].vel = { 0,0,0 };
				//O.g[B].vel = { 0,0,0 };
			}

			// Newton's equal and opposite forces applied to acceleration of each ball:
			O.g[A].acc += totalForce / O.g[A].m;
			O.g[B].acc -= totalForce / O.g[B].m;

			// So last distance can be known for COR:
			O.distances[e] = dist;
		}
		// DONT DO ANYTHING HERE. A STARTS AT 1.
	}

	if (writeStep)
	{
		ballBuffer << '\n'; // Prepares a new line for incoming data.
	}

	// THIRD PASS - Calculate velocity for next step:
	for (unsigned int Ball = 0; Ball < O.n; Ball++)
	{
		// Velocity for next step:
		O.g[Ball].vel = O.g[Ball].velh + .5 * O.g[Ball].acc * dt;
		O.g[Ball].w = O.g[Ball].wh + .5 * O.g[Ball].aacc * dt;

		if (writeStep)
		{
			// Send positions and rotations to buffer:
			if (Ball == 0)
			{
				ballBuffer
					<< O.g[Ball].pos[0] << ','
					<< O.g[Ball].pos[1] << ','
					<< O.g[Ball].pos[2] << ','
					<< O.g[Ball].w[0] << ','
					<< O.g[Ball].w[1] << ','
					<< O.g[Ball].w[2] << ','
					<< O.g[Ball].w.norm() << ','
					<< O.g[Ball].vel.x << ','
					<< O.g[Ball].vel.y << ','
					<< O.g[Ball].vel.z << ','
					<< 0;
			}
			else
			{
				ballBuffer
					<< ',' << O.g[Ball].pos[0] << ','
					<< O.g[Ball].pos[1] << ','
					<< O.g[Ball].pos[2] << ','
					<< O.g[Ball].w[0] << ','
					<< O.g[Ball].w[1] << ','
					<< O.g[Ball].w[2] << ','
					<< O.g[Ball].w.norm() << ','
					<< O.g[Ball].vel.x << ',' <<
					O.g[Ball].vel.y << ','
					<< O.g[Ball].vel.z << ','
					<< 0;
			}

			O.T += .5 * O.g[Ball].m * O.g[Ball].vel.normsquared() + .5 * O.g[Ball].moi * O.g[Ball].w.normsquared(); // Now includes rotational kinetic energy.
			O.mom += O.g[Ball].m * O.g[Ball].vel;
			O.ang_mom += O.g[Ball].m * O.g[Ball].pos.cross(O.g[Ball].vel) + O.g[Ball].moi * O.g[Ball].w;
		}
	} // THIRD PASS END

	if (writeStep)
	{

		// Write energy to stream:
		energyBuffer << '\n'
			<< simTimeElapsed << ','
			<< O.U << ','
			<< O.T << ','
			<< O.U + O.T << ','
			<< O.mom.norm() << ','
			<< O.ang_mom.norm(); // the two zeros are bound and unbound mass

		// Reinitialize energies for next step:
		O.T = 0;
		O.U = 0;
		O.mom = { 0, 0, 0 };
		O.ang_mom = { 0, 0, 0 };
		// unboundMass = 0;
		// boundMass = massTotal;

		////////////////////////////////////////////////////////////////////
		/// Data Export ////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////
		if (time(nullptr) - lastWrite > 1800 || Step / skip % 10 == 0)
		{
			// Report vMax:
			std::cerr << "vMax = " << O.getVelMax() << " Steps recorded: " << Step / skip << '\n';
			std::cerr << "Data Write\n\n";

			// Write simData to file and clear buffer.
			std::ofstream ballWrite;
			ballWrite.open(outputPrefix + "simData.csv", std::ofstream::app);
			ballWrite << ballBuffer.rdbuf(); // Barf buffer to file.
			ballBuffer.str("");              // Empty the stream for next filling.
			ballWrite.close();

			// Write Energy data to file and clear buffer.
			std::ofstream energyWrite;
			energyWrite.open(outputPrefix + "energy.csv", std::ofstream::app);
			energyWrite << energyBuffer.rdbuf();
			energyBuffer.str(""); // Empty the stream for next filling.
			energyWrite.close();

			lastWrite = time(nullptr);
		} // Data export end
		if (dynamicTime)
		{
			O.calibrateDT(Step, false);
		}
	} // writestep end
} // Steps end


[[noreturn]] void simLooper()
{
	std::cerr << "Beginning simulation...\n";

	startProgress = time(nullptr);

	for (int Step = 1; Step < steps; Step++) // Steps start at 1 because the 0 step is initial conditions.
	{
		simOneStep(Step);
	}
	const time_t end = time(nullptr);

	std::cerr << "Simulation complete!\n"
		<< O.n << " Particles and " << steps << " Steps.\n"
		<< "Simulated time: " << steps * dt << " seconds\n"
		<< "Computation time: " << end - start << " seconds\n";
	std::cerr << "\n===============================================================\n";
	// I know the number of balls in each file and the order they were brought in, so I can effect individual clusters.
	//
	// Implement calculation of total mom vector and make it 0 mag

	exit(EXIT_SUCCESS);
} // end simLooper


void safetyChecks()
{
	titleBar("SAFETY CHECKS");

	if (O.soc <= 0)
	{
		fprintf(stderr, "\nvSOC NOT SET\n");
		exit(EXIT_FAILURE);
	}

	if (O.v_collapse <= 0)
	{
		fprintf(stderr, "\nvCollapse NOT SET\n");
		exit(EXIT_FAILURE);
	}

	if (skip == 0)
	{
		fprintf(stderr, "\nSKIP NOT SET\n");
		exit(EXIT_FAILURE);
	}

	if (kin < 0)
	{
		fprintf(stderr, "\nSPRING CONSTANT NOT SET\n");
		exit(EXIT_FAILURE);
	}

	if (dt <= 0)
	{
		fprintf(stderr, "\nDT NOT SET\n");
		exit(EXIT_FAILURE);
	}

	if (steps == 0)
	{
		fprintf(stderr, "\nSTEPS NOT SET\n");
		exit(EXIT_FAILURE);
	}

	if (O.initial_radius <= 0)
	{
		fprintf(stderr, "\nCluster initialRadius not set\n");
		exit(EXIT_FAILURE);
	}

	for (unsigned int Ball = 0; Ball < O.n; Ball++)
	{
		if (O.g[Ball].pos.norm() < vec3(1e-10, 1e-10, 1e-10).norm())
		{
			fprintf(stderr, "\nA ball position is [0,0,0]. Possibly didn't initialize balls properly.\n");
			exit(EXIT_FAILURE);
		}

		if (O.g[Ball].R <= 0)
		{
			fprintf(stderr, "\nA balls radius <= 0.\n");
			exit(EXIT_FAILURE);
		}

		if (O.g[Ball].m <= 0)
		{
			fprintf(stderr, "\nA balls mass <= 0.\n");
			exit(EXIT_FAILURE);
		}
	}
	titleBar("SAFETY PASSED");
}




//void setGuidDT(const double& vel)
//{
//	// Guidos k and dt:
//	dt = .01 * O.getRmin() / fabs(vel);
//}
//
//void setGuidK(const double& vel)
//{
//	kin = O.getMassMax() * vel * vel / (.1 * O.g[0].R * .1 * O.g[0].R);
//	kout = cor * kin;
//}


