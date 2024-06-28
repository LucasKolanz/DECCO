#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits.h>
#include <cstring>
#include <vector>
#include <mutex>
#include <thread>
#include <execution>
#include "../vector3d.hpp"

class Cosmos
{
public:
	int n; // Particle count.

	// Useful values:
	double r_min = -1;
	double r_max = -1;
	double m_total = -1;
	double initial_radius = -1;
	double v_collapse = 0;
	double v_max = -1;
	double v_max_prev = HUGE_VAL;

	vec3 mom = { 0, 0, 0 };
	vec3 ang_mom = { 0, 0, 0 }; // Can be vec3 because they only matter for writing out to file. Can process on host.

	double U = 0, T = 0;

	std::vector<Sphere> g;
	std::vector<Sphere_pair> pairs;

	Cosmos() = default;

	/// @brief For creating a new ballGroup of size n
	/// @param nBalls Number of balls to allocate.
	explicit Cosmos(const int n)
		:
		n(n),
		g(std::vector<Sphere>(n)),
		pairs(std::vector<Sphere_pair>(n* (n - 1) / 2))
	{
		make_pairs();
	}

	/// @brief For generating a new ballGroup of size nBalls
	/// @param nBalls Number of balls to allocate.
	/// @param generate Just here to get you to the right constructor. This is definitely wrong.
	/// @param customVel To condition for specific vMax.
	Cosmos(const int n, const bool generate, const double& customVel)
		:
		n(n),
		g(std::vector<Sphere>(n)),
		pairs(std::vector<Sphere_pair>(n* (n - 1) / 2))
	{

		form_cluster(n);
		calc_v_collapse();
		calibrateDT(0, customVel);
		simInitCondAndCenter();
	}


	/// @brief For continuing a sim.
	/// @param fullpath is the filename and path excluding the suffix _simData.csv, _constants.csv, etc.
	/// @param customVel To condition for specific vMax.
	explicit Cosmos(const std::string& fullpath, const double& customVel)
	{
		simContinue(path, fullpath);
		calc_v_collapse();
		calibrateDT(0, customVel);
		simInitCondAndCenter();
	}

	/// @brief For two cluster sim.
	/// @param projectileName 
	/// @param targetName 
	/// @param customVel To condition for specific vMax.
	explicit Cosmos(const std::string& path, const std::string& projectileName, const std::string& targetName, const double& customVel)
	{
		simInitTwoCluster(path, projectileName, targetName);
		calc_v_collapse();
		calibrateDT(0, customVel);
		simInitCondAndCenter();
	}




	void calibrateDT(const unsigned int& Step, const double& customSpeed = -1.)
	{
		const double dtOld = dt;

		if (customSpeed > 0.)
		{
			updateDTK(customSpeed);
			std::cerr << "CUSTOM SPEED: " << customSpeed;
		}
		else
		{
			//std::cerr << vCollapse << " <- vCollapse | Lazz Calc -> " << M_PI * M_PI * G * pow(density, 4. / 3.) * pow(mTotal, 2. / 3.) * rMax;

			v_max = getVelMax();

			std::cerr << '\n';

			// Take whichever velocity is greatest:
			std::cerr << v_collapse << " = vCollapse | vMax = " << v_max;
			if (v_max < v_collapse)
			{
				v_max = v_collapse;
			}

			if (v_max < v_max_prev)
			{
				updateDTK(v_max);
				v_max_prev = v_max;
				std::cerr << "\nk: " << kin << "\tdt: " << dt;
			}
		}

		if (Step == 0 or dtOld < 0)
		{
			steps = static_cast<unsigned>(simTimeSeconds / dt);
			std::cerr << "\tInitial Steps: " << steps;
		}
		else
		{
			steps = static_cast<unsigned>(dtOld / dt * (steps - Step) + Step);
			std::cerr << "\tSteps: " << steps;
		}

		if (timeResolution / dt > 1.)
		{
			skip = static_cast<size_t>(floor(timeResolution / dt));
			std::cerr << "\tSkip: " << skip << '\n';
		}
		else
		{
			std::cerr << "Desired time resolution is lower than dt. Setting to 1 second per skip.\n";
			skip = static_cast<size_t>(floor(1. / dt));
		}
	}

	// todo - make bigger balls favor the middle, or, smaller balls favor the outside.
	/// @brief Push balls apart until no overlaps
	void pushApart()
	{
		std::cerr << "Separating spheres - Current max overlap:\n";
		/// Using acc array as storage for accumulated position change.
		int* counter = new int[n];
		for (size_t Ball = 0; Ball < n; Ball++)
		{
			g[Ball].acc = { 0, 0, 0 };
			counter[Ball] = 0;
		}

		double overlapMax = -1;
		const double pseudoDT = r_min * .1;
		int step = 0;

		while (true)
		{
			//if (step % 10 == 0)
			//{
			//	simDataWrite("pushApart_");
			//}

			for (unsigned int A = 0; A < n; A++)
			{
				for (unsigned int B = A + 1; B < n; B++)
				{
					// Check for Ball overlap.
					vec3 rVecab = g[B].pos - g[A].pos;
					vec3 rVecba = -1 * rVecab;
					const double dist = (rVecab).norm();
					const double sumRaRb = g[A].R + g[B].R;
					const double overlap = sumRaRb - dist;

					if (overlapMax < overlap)
					{
						overlapMax = overlap;
					}

					if (overlap > 0)
					{
						g[A].acc += overlap * (rVecba / dist);
						g[B].acc += overlap * (rVecab / dist);
						counter[A] += 1;
						counter[B] += 1;
					}
				}
			}

			for (size_t Ball = 0; Ball < n; Ball++)
			{
				if (counter[Ball] > 0)
				{
					g[Ball].pos += g[Ball].acc.normalized() * pseudoDT;
					g[Ball].acc = { 0, 0, 0 };
					counter[Ball] = 0;
				}
			}

			if (overlapMax > 0)
			{
				std::cerr << overlapMax << "                        \r";
			}
			else
			{
				std::cerr << "\nSuccess!\n";
				break;
			}
			overlapMax = -1;
			step++;
		}
		delete[] counter;
	}

	void calc_v_collapse()
	{
		// Sim fall velocity onto cluster:
		// vCollapse shrinks if a ball escapes but velMax should take over at that point, unless it is ignoring far balls.
		double position = 0;
		while (position < initial_radius)
		{
			v_collapse += G * m_total / (initial_radius * initial_radius) * 0.1;
			position += v_collapse * 0.1;
		}
		v_collapse = fabs(v_collapse);
	}

	/// get max velocity
	[[nodiscard]] double getVelMax()
	{
		v_max = 0;

		// todo - make this a manual set true or false to use soc so we know if it is being used or not.
		if (soc > 0)
		{
			int counter = 0;
			for (unsigned int Ball = 0; Ball < n; Ball++)
			{
				// Only consider balls moving toward com and within 4x initial radius around it.
				// todo - this cone may be too aggressive (larger cone ignores more):
				constexpr double cone = M_PI_2 + (.5 * M_PI_2);
				const vec3 fromCOM = g[Ball].pos - getCOM();
				if (acos(g[Ball].vel.normalized().dot(fromCOM.normalized())) > cone && fromCOM.norm() < soc)
				{
					if (g[Ball].vel.norm() > v_max)
					{
						v_max = g[Ball].vel.norm();
					}
				}
				else
				{
					counter++;
				}
			}
			std::cerr << '(' << counter << " spheres ignored" << ") ";
		}
		else
		{
			for (unsigned int Ball = 0; Ball < n; Ball++)
			{
				if (g[Ball].vel.norm() > v_max)
				{
					v_max = g[Ball].vel.norm();
				}
			}

			// Is vMax for some reason unreasonably small? Don't proceed. Probably a finished sim.
			// This shouldn't apply to extremely destructive collisions because it is possible that no particles are considered, so it will keep pausing.
			if (v_max < 1e-10)
			{
				std::cerr << "\nMax velocity in system is less than 1e-10.\n";
				system("pause");
			}
		}

		return v_max;
	}

	// Kick ballGroup (give the whole thing a velocity)
	void kick(const double& vx, const double& vy, const double& vz)
	{
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].vel += {vx, vy, vz};
		}
	}


	void checkMomentum(const std::string& of) const
	{
		vec3 pTotal = { 0, 0, 0 };
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			pTotal += g[Ball].m * g[Ball].vel;
		}
		fprintf(stderr, "%s Momentum Check: %.2e, %.2e, %.2e\n", of.c_str(), pTotal.x, pTotal.y, pTotal.z);
	}

	// offset cluster
	void offset(const double& rad1, const double& rad2, const double& impactParam)
	{
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].pos.x += (rad1 + rad2) * cos(impactParam);
			g[Ball].pos.y += (rad1 + rad2) * sin(impactParam);
		}
	}

	/// Approximate the radius of the ballGroup.
	[[nodiscard]] double getRadius() const
	{
		double radius = 0;

		if (n > 1)
		{
			for (unsigned int A = 0; A < n; A++)
			{
				for (unsigned int B = A + 1; B < n; B++)
				{
					// Identify two farthest balls from each other. That is diameter of cluster.
					const double diameter = (g[A].pos - g[B].pos).norm();
					if (diameter * .5 > radius)
					{
						radius = diameter * .5;
					}
				}
			}
		}
		else
		{
			radius = g[0].R;
		}

		return radius;
	}

	// Update Potential Energy:
	void updatePE()
	{
		U = 0;

		if (n > 1) // Code below only necessary for effects between balls.
		{
			for (unsigned int A = 1; A < n; A++)
			{
				for (unsigned int B = 0; B < A; B++)
				{
					const double sumRaRb = g[A].R + g[B].R;
					const double dist = (g[A].pos - g[B].pos).norm();
					const double overlap = sumRaRb - dist;

					// Check for collision between Ball and otherBall.
					if (overlap > 0)
					{
						U += -G * g[A].m * g[B].m / dist + kin * ((sumRaRb - dist) * .5) * ((sumRaRb - dist) * .5);
					}
					else
					{
						U += -G * g[A].m * g[B].m / dist;
					}
				}
			}
		}
		else // For the case of just one ball:
		{
			U = 0;
		}
	}

	void simInitWrite(std::string& filename)
	{
		// Check if file name already exists.
		std::ifstream checkForFile;
		checkForFile.open(filename + "simData.csv", std::ifstream::in);
		// Add a counter to the file name until it isn't overwriting anything:
		int counter = 0;
		while (checkForFile.is_open())
		{
			counter++;
			checkForFile.close();
			checkForFile.open(std::to_string(counter) + '_' + filename + "simData.csv", std::ifstream::in);
		}

		if (counter > 0)
		{
			filename.insert(0, std::to_string(counter) + '_');
		}

		// Complete file names:
		std::string simDataFilename = filename + "simData.csv";
		std::string energyFilename = filename + "energy.csv";
		std::string constantsFilename = filename + "constants.csv";

		std::cerr << "New file tag: " << filename;

		// Open all file streams:
		std::ofstream energyWrite, ballWrite, constWrite;
		energyWrite.open(energyFilename, std::ofstream::app);
		ballWrite.open(simDataFilename, std::ofstream::app);
		constWrite.open(constantsFilename, std::ofstream::app);

		// Make column headers:
		energyWrite << "Time,PE,KE,E,p,L";
		ballWrite << "x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0";

		for (unsigned int Ball = 1; Ball < n; Ball++) // Start at 2nd ball because first one was just written^.
		{
			std::string thisBall = std::to_string(Ball);
			ballWrite
				<< ",x" + thisBall
				<< ",y" + thisBall
				<< ",z" + thisBall
				<< ",wx" + thisBall
				<< ",wy" + thisBall
				<< ",wz" + thisBall
				<< ",wmag" + thisBall
				<< ",vx" + thisBall
				<< ",vy" + thisBall
				<< ",vz" + thisBall
				<< ",bound" + thisBall;
		}

		// Write constant data:
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{

			constWrite
				<< g[Ball].R << ','
				<< g[Ball].m << ','
				<< g[Ball].moi
				<< '\n';
		}

		// Write energy data to buffer:
		energyBuffer
			<< '\n'
			<< simTimeElapsed << ','
			<< U << ','
			<< T << ','
			<< U + T << ','
			<< mom.norm() << ','
			<< ang_mom.norm();
		energyWrite << energyBuffer.rdbuf();
		energyBuffer.str("");

		// Reinitialize energies for next step:
		T = 0;
		U = 0;
		mom = { 0, 0, 0 };
		ang_mom = { 0, 0, 0 };

		// Send position and rotation to buffer:
		ballBuffer << '\n'; // Necessary new line after header.
		ballBuffer
			<< g[0].pos.x << ','
			<< g[0].pos.y << ','
			<< g[0].pos.z << ','
			<< g[0].w.x << ','
			<< g[0].w.y << ','
			<< g[0].w.z << ','
			<< g[0].w.norm() << ','
			<< g[0].vel.x << ','
			<< g[0].vel.y << ','
			<< g[0].vel.z << ','
			<< 0; //bound[0];
		for (unsigned int Ball = 1; Ball < n; Ball++)
		{
			ballBuffer
				<< ',' << g[Ball].pos.x << ',' // Needs comma start so the last bound doesn't have a dangling comma.
				<< g[Ball].pos.y << ','
				<< g[Ball].pos.z << ','
				<< g[Ball].w.x << ','
				<< g[Ball].w.y << ','
				<< g[Ball].w.z << ','
				<< g[Ball].w.norm() << ','
				<< g[Ball].vel.x << ','
				<< g[Ball].vel.y << ','
				<< g[Ball].vel.z << ','
				<< 0; //bound[Ball];
		}
		// Write position and rotation data to file:
		ballWrite << ballBuffer.rdbuf();
		ballBuffer.str(""); // Resets the stream buffer to blank.

		// Close Streams for user viewing:
		energyWrite.close();
		ballWrite.close();
		constWrite.close();

		std::cerr << "\nSimulating " << steps * dt / 60 / 60 << " hours.\n";
		std::cerr << "Total mass: " << m_total << '\n';
		std::cerr << "\n===============================================================\n";
	}


	[[nodiscard]] vec3 getCOM() const
	{
		if (m_total > 0)
		{
			vec3 comNumerator;
			for (unsigned int Ball = 0; Ball < n; Ball++)
			{
				comNumerator += g[Ball].m * g[Ball].pos;
			}
			vec3 com = comNumerator / m_total;
			return com;
		}
		else
		{
			std::cerr << "Mass of cluster is zero.\n";
			exit(EXIT_FAILURE);
		}
	}

	void zeroVel()
	{
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].vel = { 0, 0, 0 };
		}
	}

	void zeroAngVel()
	{
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].w = { 0, 0, 0 };
		}
	}

	void toOrigin()
	{
		const vec3 com = getCOM();

		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].pos -= com;
		}
	}

	// Set velocity of all balls such that the cluster spins:
	void comSpinner(const double& spinX, const double& spinY, const double& spinZ)
	{
		const vec3 comRot = { spinX, spinY, spinZ }; // Rotation axis and magnitude
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].vel += comRot.cross(g[Ball].pos - getCOM());
			g[Ball].w += comRot;
		}
	}

	void rotAll(const char axis, const double angle)
	{
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			g[Ball].pos = g[Ball].pos.rot(axis, angle);
			g[Ball].vel = g[Ball].vel.rot(axis, angle);
			g[Ball].w = g[Ball].w.rot(axis, angle);
		}
	}

private:
	// String buffers to hold data in memory until worth writing to file:
	std::stringstream ballBuffer;
	std::stringstream energyBuffer;

	void make_pairs()
	{
		for (size_t i = 0; i < pairs.size(); i++)
		{
			// Pair Combinations [A,B] [B,C] [C,D]... [A,C] [B,D] [C,E]... ...
			int A = i % n;
			int stride = 1 + i / n; // Stride increases by 1 after each full set of pairs
			int B = (A + stride) % n;

			// Create particle* pair
			pairs[i] = { &g[A], &g[B] };
		}
	}

	void update_kinematics(Sphere& Sphere)
	{
		// Update velocity half step:
		Sphere.velh = Sphere.vel + .5 * Sphere.acc * dt;

		// Update angular velocity half step:
		Sphere.wh = Sphere.w + .5 * Sphere.aacc * dt;

		// Update position:
		Sphere.pos += Sphere.velh * dt;

		// Reinitialize acceleration to be recalculated:
		Sphere.acc = { 0, 0, 0 };

		// Reinitialize angular acceleration to be recalculated:
		Sphere.aacc = { 0, 0, 0 };
	}


	void write_to_buffer()
	{
		for (size_t i = 0; i < n; i++)
		{
			// Send positions and rotations to buffer:
			if (i == 0)
			{
				ballBuffer
					<< g[i].pos[0] << ','
					<< g[i].pos[1] << ','
					<< g[i].pos[2] << ','
					<< g[i].w[0] << ','
					<< g[i].w[1] << ','
					<< g[i].w[2] << ','
					<< g[i].w.norm() << ','
					<< g[i].vel.x << ','
					<< g[i].vel.y << ','
					<< g[i].vel.z << ','
					<< 0;
			}
			else
			{
				ballBuffer
					<< ',' << g[i].pos[0] << ','
					<< g[i].pos[1] << ','
					<< g[i].pos[2] << ','
					<< g[i].w[0] << ','
					<< g[i].w[1] << ','
					<< g[i].w[2] << ','
					<< g[i].w.norm() << ','
					<< g[i].vel.x << ','
					<< g[i].vel.y << ','
					<< g[i].vel.z << ','
					<< 0;
			}

			T += .5 * g[i].m * g[i].vel.normsquared() + .5 * g[i].moi * g[i].w.normsquared(); // Now includes rotational kinetic energy.
			mom += g[i].m * g[i].vel;
			ang_mom += g[i].m * g[i].pos.cross(g[i].vel) + g[i].moi * g[i].w;
		}
	}

	void compute_acceleration(Sphere_pair& p_pair)
	{
		const double Ra = p_pair.A->R;
		const double Rb = p_pair.B->R;
		const double m_a = p_pair.A->m;
		const double m_b = p_pair.B->m;
		const double sumRaRb = Ra + Rb;
		vec3 rVec = p_pair.B->pos - p_pair.A->pos; // Start with rVec from a to b.
		const double dist = (rVec).norm();
		vec3 totalForce;

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
			vec3 dVel = p_pair.B->vel - p_pair.A->vel;
			vec3 frictionForce = { 0, 0, 0 };
			const vec3 relativeVelOfA = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - p_pair.A->w.cross(p_pair.A->R / sumRaRb * rVec) - p_pair.B->w.cross(p_pair.B->R / sumRaRb * rVec);
			double relativeVelMag = relativeVelOfA.norm();
			if (relativeVelMag > 1e-10) // When relative velocity is very low, dividing its vector components by its magnitude below is unstable.
			{
				frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfA / relativeVelMag);
			}

			// Torque a:
			const vec3 aTorque = (p_pair.A->R / sumRaRb) * rVec.cross(frictionForce);

			// Gravity on a:
			const vec3 gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);

			// Total forces on a:
			totalForce = gravForceOnA + elasticForce + frictionForce + vdwForce;

			// Elastic and Friction b:
			// Flip direction b -> a:
			rVec = -rVec;
			dVel = -dVel;
			elasticForce = -elasticForce;

			const vec3 relativeVelOfB = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - p_pair.B->w.cross(p_pair.B->R / sumRaRb * rVec) - p_pair.A->w.cross(p_pair.A->R / sumRaRb * rVec);
			relativeVelMag = relativeVelOfB.norm(); // todo - This should be the same as mag for A. Same speed different direction.
			if (relativeVelMag > 1e-10)
			{
				frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfB / relativeVelMag);
			}
			const vec3 bTorque = (p_pair.B->R / sumRaRb) * rVec.cross(frictionForce);

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
			const vec3 gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);
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



	[[nodiscard]] double getRmin()
	{
		r_min = g[0].R;
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			if (g[Ball].R < r_min)
			{
				r_min = g[Ball].R;
			}
		}
		return r_min;
	}

	[[nodiscard]] double getRmax()
	{
		r_max = g[0].R;
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			if (g[Ball].R > r_max)
			{
				r_max = g[Ball].R;
			}
		}
		return r_max;
	}


	[[nodiscard]] double getMassMax() const
	{
		double mMax = g[0].m;
		for (unsigned int Ball = 0; Ball < n; Ball++)
		{
			if (g[Ball].m > mMax)
			{
				mMax = g[Ball].m;
			}
		}
		return mMax;
	}





	void parseSimData(std::string line)
	{
		std::string lineElement;

		// Get number of balls in file
		unsigned int count = std::count(line.begin(), line.end(), ',') / properties + 1;
		g.resize(count);

		std::stringstream chosenLine(line); // This is the last line of the read file, containing all data for all balls at last time step

		// Get position and angular velocity data:
		for (unsigned int A = 0; A < n; A++)
		{

			for (unsigned int i = 0; i < 3; i++) // Position
			{
				std::getline(chosenLine, lineElement, ',');
				g[A].pos[i] = std::stod(lineElement);
				//std::cerr << tclus.g[A].pos[i]<<',';
			}
			for (unsigned int i = 0; i < 3; i++) // Angular Velocity
			{
				std::getline(chosenLine, lineElement, ',');
				g[A].w[i] = std::stod(lineElement);
			}
			std::getline(chosenLine, lineElement, ','); // Angular velocity magnitude skipped
			for (unsigned int i = 0; i < 3; i++)                 // velocity
			{
				std::getline(chosenLine, lineElement, ',');
				g[A].vel[i] = std::stod(lineElement);
			}
			for (unsigned int i = 0; i < properties - 10; i++) // We used 10 elements. This skips the rest.
			{
				std::getline(chosenLine, lineElement, ',');
			}
		}
	}


	/// Get previous sim constants by filename.
	void loadConsts(const std::string& path, const std::string& filename)
	{
		// Get radius, mass, moi:
		std::string constantsFilename = path + filename + "constants.csv";
		if (auto ConstStream = std::ifstream(constantsFilename, std::ifstream::in))
		{
			std::string line, lineElement;
			for (unsigned int A = 0; A < n; A++)
			{
				std::getline(ConstStream, line); // Ball line.
				std::stringstream chosenLine(line);
				std::getline(chosenLine, lineElement, ','); // Radius.
				g[A].R = std::stod(lineElement);
				std::getline(chosenLine, lineElement, ','); // Mass.
				g[A].m = std::stod(lineElement);
				std::getline(chosenLine, lineElement, ','); // Moment of inertia.
				g[A].moi = std::stod(lineElement);
			}
		}
		else
		{
			std::cerr << "Could not open constants file: " << constantsFilename << "... Existing program." << '\n';
			exit(EXIT_FAILURE);
		}
	}



	/// Get last line of previous simData by filename.
	[[nodiscard]] static std::string getLastLine(const std::string& path, const std::string& filename)
	{
		std::string simDataFilepath = path + filename + "simData.csv";

		if (auto simDataStream = std::ifstream(simDataFilepath, std::ifstream::in))
		{

			std::cerr << "\nParsing last line of data.\n";

			simDataStream.seekg(-1, std::ios_base::end); // go to one spot before the EOF

			bool keepLooping = true;
			while (keepLooping)
			{
				char ch = ' ';
				simDataStream.get(ch); // Get current byte's data

				if (static_cast<int>(simDataStream.tellg()) <= 1)
				{                           // If the data was at or before the 0th byte
					simDataStream.seekg(0); // The first line is the last line
					keepLooping = false;    // So stop there
				}
				else if (ch == '\n')
				{                        // If the data was a newline
					keepLooping = false; // Stop at the current position.
				}
				else
				{                                                // If the data was neither a newline nor at the 0 byte
					simDataStream.seekg(-2, std::ios_base::cur); // Move to the front of that data, then to the front of the data before it
				}
			}
			std::string line;
			std::getline(simDataStream, line); // Read the current line

			return line;
		}
		else
		{
			std::cerr << "Could not open simData file: " << simDataFilepath << "... Existing program." << '\n';
			exit(EXIT_FAILURE);
		}

	}




	void simDataWrite(std::string& outFilename)
	{
		// todo - for some reason I need checkForFile instead of just using ballWrite. Need to work out why.
		// Check if file name already exists. If not, initialize
		std::ifstream checkForFile;
		checkForFile.open(outFilename + "simData.csv", std::ifstream::in);
		if (checkForFile.is_open() == false)
		{
			simInitWrite(outFilename);
		}
		else
		{
			ballBuffer << '\n'; // Prepares a new line for incoming data.

			for (size_t Ball = 0; Ball < n; Ball++)
			{
				// Send positions and rotations to buffer:
				if (Ball == 0)
				{
					ballBuffer
						<< g[Ball].pos[0] << ','
						<< g[Ball].pos[1] << ','
						<< g[Ball].pos[2] << ','
						<< g[Ball].w[0] << ','
						<< g[Ball].w[1] << ','
						<< g[Ball].w[2] << ','
						<< g[Ball].w.norm() << ','
						<< g[Ball].vel.x << ','
						<< g[Ball].vel.y << ','
						<< g[Ball].vel.z << ','
						<< 0;
				}
				else
				{
					ballBuffer << ','
						<< g[Ball].pos[0] << ','
						<< g[Ball].pos[1] << ','
						<< g[Ball].pos[2] << ','
						<< g[Ball].w[0] << ','
						<< g[Ball].w[1] << ','
						<< g[Ball].w[2] << ','
						<< g[Ball].w.norm() << ','
						<< g[Ball].vel.x << ','
						<< g[Ball].vel.y << ','
						<< g[Ball].vel.z << ','
						<< 0;
				}
			}

			// Write simData to file and clear buffer.
			std::ofstream ballWrite;
			ballWrite.open(outFilename + "simData.csv", std::ofstream::app);
			ballWrite << ballBuffer.rdbuf(); // Barf buffer to file.
			ballBuffer.str("");              // Resets the stream for that balls to blank.
			ballWrite.close();
		}
		checkForFile.close();
	}



	[[nodiscard]] double getMass()
	{
		m_total = 0;
		{
			for (unsigned int Ball = 0; Ball < n; Ball++)
			{
				m_total += g[Ball].m;
			}
		}
		return m_total;
	}







	/// Make ballGroup from file data.
	void loadSim(const std::string& path, const std::string& filename)
	{
		parseSimData(getLastLine(path, filename));

		loadConsts(path, filename);

		r_min = getRmin();
		r_max = getRmax();
		m_total = getMass();
		initial_radius = getRadius();
		soc = 4 * r_max + initial_radius;

		std::cerr << "Balls: " << n << '\n';
		std::cerr << "Mass: " << m_total << '\n';
		std::cerr << "Approximate radius: " << initial_radius << " cm.\n";
	}



	void oneSizeSphere(const unsigned int nBalls)
	{

		for (unsigned int Ball = 0; Ball < nBalls; Ball++)
		{
			g[Ball].R = scaleBalls;
			g[Ball].m = density * 4. / 3. * 3.14159 * std::pow(g[Ball].R, 3);
			g[Ball].moi = .4 * g[Ball].m * g[Ball].R * g[Ball].R;
			g[Ball].w = { 0, 0, 0 };
			g[Ball].pos = rand_spherical_vec(spaceRange, spaceRange, spaceRange);
		}

		// Generate non-overlapping spherical particle field:
		int collisionDetected = 0;
		int oldCollisions = nBalls;

		for (unsigned int failed = 0; failed < attempts; failed++)
		{
			for (unsigned int A = 0; A < nBalls; A++)
			{
				for (unsigned int B = A + 1; B < nBalls; B++)
				{
					// Check for Ball overlap.
					const double dist = (g[A].pos - g[B].pos).norm();
					const double sumRaRb = g[A].R + g[B].R;
					const double overlap = dist - sumRaRb;
					if (overlap < 0)
					{
						collisionDetected += 1;
						// Move the other ball:
						g[B].pos = rand_spherical_vec(spaceRange, spaceRange, spaceRange);
					}
				}
			}
			if (collisionDetected < oldCollisions)
			{
				oldCollisions = collisionDetected;
				std::cerr << "Collisions: " << collisionDetected << "                        \r";
			}
			if (collisionDetected == 0)
			{
				std::cerr << "\nSuccess!\n";
				break;
			}
			if (failed == attempts - 1 || collisionDetected > static_cast<int>(1.5 * static_cast<double>(nBalls))) // Added the second part to speed up spatial constraint increase when there are clearly too many collisions for the space to be feasible.
			{
				std::cerr << "Failed " << spaceRange << ". Increasing range " << spaceRangeIncrement << "cm^3.\n";
				spaceRange += spaceRangeIncrement;
				failed = 0;
				for (unsigned int Ball = 0; Ball < nBalls; Ball++)
				{
					g[Ball].pos = rand_spherical_vec(spaceRange, spaceRange, spaceRange); // Each time we fail and increase range, redistribute all balls randomly so we don't end up with big balls near mid and small balls outside.
				}
			}
			collisionDetected = 0;
		}

		std::cerr << "Final spacerange: " << spaceRange << '\n';
		std::cerr << "Initial Radius: " << getRadius() << '\n';
		std::cerr << "Mass: " << m_total << '\n';
	}



	void updateDTK(const double& vel)
	{
		constexpr double kConsts = fourThirdsPiRho / (maxOverlap * maxOverlap);

		kin = kConsts * r_max * vel * vel;
		kout = cor * kin;
		dt = .01 * sqrt((fourThirdsPiRho / kin) * r_min * r_min * r_min);
	}

	void simInitCondAndCenter()
	{
		std::cerr << "==================" << '\n';
		std::cerr << "dt: " << dt << '\n';
		std::cerr << "k: " << kin << '\n';
		std::cerr << "Skip: " << skip << '\n';
		std::cerr << "Steps: " << steps << '\n';
		std::cerr << "==================" << '\n';

		toOrigin();

		checkMomentum("After Zeroing"); // Is total mom zero like it should be?

		// Compute physics between all balls. Distances, collision forces, energy totals, total mass:
		initConditions();

		// Name the file based on info above:
		outputPrefix +=
			"_k" + scientific(kin) +
			"_dt" + scientific(dt) +
			"_";
	}


	void simContinue(const std::string& path, const std::string& filename)
	{
		// Load file data:
		std::cerr << "Continuing Sim...\nFile: " << filename << '\n';

		loadSim(path, filename);

		std::cerr << '\n';
		checkMomentum("O");

		// Name the file based on info above:
		outputPrefix =
			std::to_string(n) +
			"_rho" + rounder(density, 4);
	}

	// Set's up a two cluster collision.
	void simInitTwoCluster(const std::string& path, const std::string& projectileName, const std::string& targetName)
	{
		// Load file data:
		std::cerr << "TWO CLUSTER SIM\nFile 1: " << projectileName << '\t' << "File 2: " << targetName << '\n';

		// DART PROBE
		//ballGroup projectile(1);
		//projectile.g[0].pos = { 8814, 0, 0 };
		//projectile.g[0].w = { 0, 0, 0 };
		//projectile.g[0].vel = { 0, 0, 0 };
		//projectile.g[0].R = 78.5;
		//projectile.g[0].m = 560000;
		//projectile.moi[0] = .4 * projectile.g[0].m * projectile.g[0].R * projectile.g[0].R;

		Cosmos projectile;
		projectile.loadSim(path, projectileName);
		Cosmos target;
		target.loadSim(path, targetName);

		// DO YOU WANT TO STOP EVERYTHING?
		projectile.zeroAngVel();
		projectile.zeroVel();
		target.zeroAngVel();
		target.zeroVel();


		// Calc info to determined cluster positioning and collisions velocity:
		projectile.updatePE();
		target.updatePE();

		projectile.offset(projectile.initial_radius, target.initial_radius + target.getRmax() * 2, impactParameter);

		const double PEsys = projectile.U + target.U + (-G * projectile.m_total * target.m_total / (projectile.getCOM() - target.getCOM()).norm());

		// Collision velocity calculation:
		const double mSmall = projectile.m_total;
		const double mBig = target.m_total;
		const double mTot = mBig + mSmall;
		//const double vSmall = -sqrt(2 * KEfactor * fabs(PEsys) * (mBig / (mSmall * mTot))); // Negative because small offsets right.
		const double vSmall = -vCustom; // DART probe override.
		const double vBig = -(mSmall / mBig) * vSmall; // Negative to oppose projectile.
		//const double vBig = 0; // Dymorphous override.

		if (isnan(vSmall) || isnan(vBig))
		{
			std::cerr << "A VELOCITY WAS NAN!!!!!!!!!!!!!!!!!!!!!!\n\n";
			exit(EXIT_FAILURE);
		}

		projectile.kick(vSmall, 0, 0);
		target.kick(vBig, 0, 0);

		fprintf(stderr, "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n", vBig, vSmall);

		std::cerr << '\n';
		projectile.checkMomentum("Projectile");
		target.checkMomentum("Target");

		g.insert(g.begin(), projectile.g.begin(), projectile.g.end());
		g.insert(g.end(), target.g.begin(), target.g.end());

		outputPrefix =
			projectileName + targetName +
			"T" + rounder(KEfactor, 4) +
			"_vBig" + scientific(vBig) +
			"_vSmall" + scientific(vSmall) +
			"_IP" + rounder(impactParameter * 180 / 3.14159, 2) +
			"_rho" + rounder(density, 4);
	}
};
