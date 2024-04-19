#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <vector>
#include <algorithm>
#include <execution>
#include <mutex>

#include "../vector3d.hpp"
#include "../initializations.hpp"

// String buffers to hold data in memory until worth writing to file:
std::stringstream ballBuffer;
std::stringstream energyBuffer;

// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
const time_t start = time(nullptr);        // For end of program analysis
time_t startProgress; // For progress reporting (gets reset)
time_t lastWrite;     // For write control (gets reset)
bool writeStep;       // This prevents writing to file every step (which is slow).
std::mutex g_mutex;

//ballGroup O(path, projectileName, targetName, vCustom); // Collision
//ballGroup O(path, targetName, 0); // Continue
//ballGroup O(genBalls, true, vCustom); // Generate


struct P;
struct P_pair;
struct P_group;
void update_kinematics(P& P);
void compute_acceleration(P_pair& p_pair);
void compute_velocity(P& P);
void write_to_buffer(P_group& p_group);

int main()
{
	int n = 10000; // Number of particles
	P_group O(n); // Particle system
	std::vector<P_pair> pairs = make_p_pairs(O);

	std::for_each(std::execution::par_unseq, O.p_group.begin(), O.p_group.end(), update_kinematics);
	std::for_each(std::execution::par, pairs.begin(), pairs.end(), compute_acceleration);
	//std::ranges::for_each( arr, [i=0](auto &e) mutable { long_function(e,i++); } );
	if (writeStep)
	{
		ballBuffer << '\n'; // Prepares a new line for incoming data.
		write_to_buffer(O);
	}
}

struct P
{
	vector3d
		pos,
		vel,
		velh,
		w,
		wh,
		acc,
		aacc;

	double
		R = 0,
		m = 0,
		moi = 0;
};

struct P_pair
{
	P* A;
	P* B;
	double dist;
	P_pair() = default;
	P_pair(P* a, P* b, double d) : A(a), B(b), dist(d) {} // Init all
	P_pair(P* a, P* b) : A(a), B(b), dist(-1.0) {} // init pairs and set D to illogical distance
	//p_pair(const p_pair& in_pair) : A(in_pair.A), B(in_pair.B), dist(in_pair.dist) {} // Copy constructor
};

struct P_group
{
	int n; // Number of particles in the group
	std::vector<P> p_group; // Group of particles
	double U; // Potential energy
	double T; // Kinetic Energy
	vector3d mom;
	vector3d ang_mom;
	// Useful values:
	double rMin = -1;
	double rMax = -1;
	double mTotal = -1;
	double initialRadius = -1;
	double vCollapse = 0;
	double vMax = -1;
	double vMaxPrev = HUGE_VAL;
	double soc = -1;

	P_group() = default;
	P_group(int n) : n(n) {}
	P_group(std::vector<P> p_group) : p_group(p_group) {}
};

std::vector<P_pair> make_p_pairs(P_group& O)
{
	int n = O.p_group.size();
	int n_pairs = n * (n - 1) / 2;
	std::vector<P_pair> p_pairs(n_pairs); // All particle pairs
	for (size_t i = 0; i < n_pairs; i++)
	{
		// Pair Combinations [A,B] [B,C] [C,D]... [A,C] [B,D] [C,E]... ...
		int A = i % n;
		int stride = 1 + i / n; // Stride increases by 1 after each full set of pairs
		int B = (A + stride) % n;

		// Create particle* pair
		p_pairs[i] = { &O.p_group[A], &O.p_group[B] };
	}
	return p_pairs;
}

void write_to_buffer(P_group& p_group)
{
	for (size_t i = 0; i < p_group.n; i++)
	{
		P& cp = p_group.p_group[i]; // Current particle

		// Send positions and rotations to buffer:
		if (i == 0)
		{
			ballBuffer
				<< cp.pos[0] << ','
				<< cp.pos[1] << ','
				<< cp.pos[2] << ','
				<< cp.w[0] << ','
				<< cp.w[1] << ','
				<< cp.w[2] << ','
				<< cp.w.norm() << ','
				<< cp.vel.x << ','
				<< cp.vel.y << ','
				<< cp.vel.z << ','
				<< 0;
		}
		else
		{
			ballBuffer
				<< ',' << cp.pos[0] << ','
				<< cp.pos[1] << ','
				<< cp.pos[2] << ','
				<< cp.w[0] << ','
				<< cp.w[1] << ','
				<< cp.w[2] << ','
				<< cp.w.norm() << ','
				<< cp.vel.x << ','
				<< cp.vel.y << ','
				<< cp.vel.z << ','
				<< 0;
		}

		p_group.T += .5 * cp.m * cp.vel.normsquared() + .5 * cp.moi * cp.w.normsquared(); // Now includes rotational kinetic energy.
		p_group.mom += cp.m * cp.vel;
		p_group.ang_mom += cp.m * cp.pos.cross(cp.vel) + cp.moi * cp.w;
	}
}

void update_kinematics(P& P)
{
	// Update velocity half step:
	P.velh = P.vel + .5 * P.acc * dt;

	// Update angular velocity half step:
	P.wh = P.w + .5 * P.aacc * dt;

	// Update position:
	P.pos += P.velh * dt;

	// Reinitialize acceleration to be recalculated:
	P.acc = { 0, 0, 0 };

	// Reinitialize angular acceleration to be recalculated:
	P.aacc = { 0, 0, 0 };
}

void compute_acceleration(P_pair& p_pair)
{
	const double Ra = p_pair.A->R;
	const double Rb = p_pair.B->R;
	const double m_a = p_pair.A->m;
	const double m_b = p_pair.B->m;
	const double sumRaRb = Ra + Rb;
	vector3d rVec = p_pair.B->pos - p_pair.A->pos; // Start with rVec from a to b.
	const double dist = (rVec).norm();
	vector3d totalForce;

	// Check for collision between Ball and otherBall:
	double overlap = sumRaRb - dist;

	double oldDist = p_pair.dist;

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
		const double h2 = h * h;
		const double twoRah = 2 * Ra * h;
		const double twoRbh = 2 * Rb * h;
		const vector3d vdwForce =
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
		vector3d elasticForce = -k * overlap * .5 * (rVec / dist);

		// Friction a:
		vector3d dVel = p_pair.B->vel - p_pair.A->vel;
		vector3d frictionForce = { 0, 0, 0 };
		const vector3d relativeVelOfA = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - p_pair.A->w.cross(p_pair.A->R / sumRaRb * rVec) - p_pair.B->w.cross(p_pair.B->R / sumRaRb * rVec);
		double relativeVelMag = relativeVelOfA.norm();
		if (relativeVelMag > 1e-10) // When relative velocity is very low, dividing its vector components by its magnitude below is unstable.
		{
			frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfA / relativeVelMag);
		}

		// Torque a:
		const vector3d aTorque = (p_pair.A->R / sumRaRb) * rVec.cross(frictionForce);

		// Gravity on a:
		const vector3d gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);

		// Total forces on a:
		totalForce = gravForceOnA + elasticForce + frictionForce + vdwForce;

		// Elastic and Friction b:
		// Flip direction b -> a:
		rVec = -rVec;
		dVel = -dVel;
		elasticForce = -elasticForce;

		const vector3d relativeVelOfB = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - p_pair.B->w.cross(p_pair.B->R / sumRaRb * rVec) - p_pair.A->w.cross(p_pair.A->R / sumRaRb * rVec);
		relativeVelMag = relativeVelOfB.norm(); // todo - This should be the same as mag for A. Same speed different direction.
		if (relativeVelMag > 1e-10)
		{
			frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfB / relativeVelMag);
		}
		const vector3d bTorque = (p_pair.B->R / sumRaRb) * rVec.cross(frictionForce);

		{
			const std::lock_guard<std::mutex> lock(g_mutex);
			p_pair.A->aacc += aTorque / p_pair.A->moi;
		}
		{
			const std::lock_guard<std::mutex> lock(g_mutex);
			p_pair.B->aacc += bTorque / p_pair.B->moi;
		}


		if (writeStep)
		{
			// Calculate potential energy. Important to recognize that the factor of 1/2 is not in front of K because this is for the spring potential in each ball and they are the same potential.
			//O.PE += -G * pair.A->m * pair.B->m / dist + 0.5 * k * overlap * overlap;
		}
	}
	else
	{
		// No collision: Include gravity only:
		const vector3d gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);
		totalForce = gravForceOnA;
		if (writeStep)
		{
			//O.PE += -G * pair.A->m * pair.B->m / dist;
		}

		// For expanding overlappers:
		//pair.A->vel = { 0,0,0 };
		//pair.B->vel = { 0,0,0 };
	}

	// Newton's equal and opposite forces applied to acceleration of each ball:
	{
		const std::lock_guard<std::mutex> lock(g_mutex);
		p_pair.A->acc += totalForce / p_pair.A->m;
	}
	{
		const std::lock_guard<std::mutex> lock(g_mutex);
		p_pair.B->acc -= totalForce / p_pair.B->m;
	}
}

void compute_velocity(P& P)
{
	// Velocity for next step:
	P.vel = P.velh + .5 * P.acc * dt;
	P.w = P.wh + .5 * P.aacc * dt;
}
