#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
#include <execution>
#include <thread>
#include "../dust_const.hpp"
#include "../vec3.hpp"

/// @brief Contains everything needed to describe the physical state of a sphere.
/// Radius, mass, and moment of inertia must be defined upon creation.
struct Sphere
{
	vec3
		pos,
		vel,
		velh,
		w,
		wh,
		force,
		torque;

	const double
		R,
		m,
		moi;

	Sphere(const double R, const double density)
		:
		R(R),
		m(density * 4. / 3. * M_PI * R * R * R),
		moi(.4 * m * R * R)
	{
	}
};

/// @brief Creates a reference to a pair of spheres.
struct Sphere_pair
{
private:
	Sphere& a_;
	Sphere& b_;
	double prev_overlap_;
	const double sum_Ra_Rb;

public:
	// Getters
	Sphere& a() const { return a_; }
	Sphere& b() const { return b_; }
	double get_prev_overlap() const { return prev_overlap_; }

	// Setters
	void a_add_force(const vec3& force) { a_.force += force; }
	void b_add_force(const vec3& force) { b_.force += force; }
	void update_overlap()
	{
		prev_overlap_ = sum_Ra_Rb - (a_.pos - b_.pos).norm();
	}

	// Constructors
	Sphere_pair(Sphere& a, Sphere& b) :
		a_(a),
		b_(b),
		sum_Ra_Rb(a.R + b.R)
	{
	}



};

class Cosmos
{
public:
	int n; // Particle count.
	std::vector<Sphere> spheres;

	// Useful values:
	double r_min = -1; // Smallest particle.
	double r_max = -1; // Largest particle.
	double m_total = -1; // Total mass of the cosmos.
	double initial_radius = -1; // Determined by the two most distance particles from each other.
	double v_collapse = 0; // Fastest particle at the moment the most distant matter collapses onto the rest.
	double v_max = -1; // Fastest particle.
	double v_max_prev = HUGE_VAL;

	vec3 mom = { 0, 0, 0 }; // Momentum of the cosmos.
	vec3 ang_mom = { 0, 0, 0 }; // Angular momentum of the cosmos.

	double U = 0; // Potential energy.
	double T = 0; // Kinetic energy.

	Cosmos() = default;
	Cosmos(std::vector<Sphere> spheres) : spheres(spheres) {}

	void operator()(int step)
	{
		// perform one timestep.

		// then in main all I do is for loop O(i)
	}
};

struct Rand_vec_in_sphere
{
	const double radius;

	Rand_vec_in_sphere(const double& radius) : radius(radius) {}

	void operator()(vec3& vec)
	{
		do
		{
			vec = { rand_double(radius), rand_double(radius), rand_double(radius) };
		} while (vec.norm() > radius);
	}
};

struct Collision_displacer
{
	void operator()(Sphere_pair& pair)
	{
		if (pair.get_prev_overlap() > 0)
		{
			pair.a().pos = { 1, 3, 5 };
		}
	}
};

vec3 vec_atob(const Sphere_pair& pair)
{
	return pair.a().pos - pair.b().pos;
}

double distance(const Sphere_pair& pair)
{
	vec3 vec = pair.b().pos - pair.a().pos;
	return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

void apply_grav_force(Sphere_pair& pair, double& dist)
{
	vec3 force = (G * pair.a().m * pair.b().m / (dist * dist)) * (vec_atob(pair) / dist);
	pair.a_add_force(force);
	pair.b_add_force(-force);
}

void apply_elastic_force(Sphere_pair& pair, double& k)
{
	//todo -this get_overlap is bad. I need all the force functions to be able to use overlap without recalcing for each. Maybe they need to be moved into the Sphere_pair
	const vec3 force = k * pair.get_overlap() * .5 * (vec_atob(pair).normalized());
	pair.a_add_force(-force);
	pair.b_add_force(force);
}

void apply_friction_sliding(Sphere& a, Sphere& b)
{

}

void apply_friction_rolling(Sphere& a, Sphere& b)
{

}

void apply_cohesion_force(Sphere& a, Sphere& b)
{

}

void update_sphere_kinematics(Sphere& sphere)
{
	sphere.velh = sphere.vel + .5 * sphere.force / sphere.m * dt;
	sphere.wh = sphere.w + .5 * sphere.torque / sphere.moi * dt;
	sphere.pos += sphere.velh * dt;
	sphere.force = { 0, 0, 0 };
	sphere.torque = { 0, 0, 0 };
}

std::vector<Sphere_pair> make_pairs(std::vector<Sphere> spheres)
{
	int n = spheres.size();
	if (n % 2 != 0)
	{
		std::cerr << "Sphere count not even.";
		std::exit(EXIT_FAILURE);
	}

	int combinations = n * (n - 1) / 2;

	std::vector<Sphere_pair> pairs;
	pairs.reserve(combinations);

	for (size_t i = 0; i < combinations; i++)
	{
		// Pair Combinations [A,B] [B,C] [C,D]... [A,C] [B,D] [C,E]... ...
		int A = i % n;
		int stride = 1 + i / n; // Stride increases by 1 after each full set of pairs
		int B = (A + stride) % n;

		// Create particle& pair
		pairs.emplace_back(spheres[A], spheres[B]);
	}
	return pairs;
}

void make_random_cluster(std::vector<Sphere>& spheres)
{
	// Seed for random cluster.
	const int seed = time(nullptr);
	srand(seed);

	const int n = spheres.size();

	const int smalls = std::round(static_cast<double>(n) * 27. / 31.375);
	const int mediums = std::round(static_cast<double>(n) * 27. / (8 * 31.375));
	const int larges = std::round(static_cast<double>(n) * 1. / 31.375);

	std::thread t1{ [&spheres, radius = 3 * scaleBalls] {
		std::for_each(
			std::execution::par_unseq,
			spheres.begin(),
			spheres.end(),
			Rand_vec_in_sphere(spaceRange));
	} };

	std::thread t2{ [&spheres, radius = 2 * scaleBalls, larges, mediums] {
		std::for_each(
			std::execution::par_unseq,
			spheres.begin() + larges,
			spheres.begin() + larges + mediums,
			Rand_vec_in_sphere(spaceRange));
	} };

	std::thread t3{ [&spheres, radius = scaleBalls, larges, mediums] {
		std::for_each(
			std::execution::par_unseq,
			spheres.begin() + larges + mediums,
			spheres.end(),
			Rand_vec_in_sphere(spaceRange));
	} };

	t1.join(); t2.join(); t3.join();
}

void make_no_collisions(std::vector<Sphere_pair> pairs)
{
	int collisionDetected = 0;
	int oldCollisions = pairs.size();

	for (int failed = 0; failed < attempts; failed++)
	{
		// todo start here next time
		std::for_each(std::execution::par_unseq, pairs.begin(), pairs.end(), Collision_displacer)
			// Check for Ball overlap.

			if (overlap < 0)
			{
				collisionDetected += 1;
				// Move the other ball:
				g[B].pos = rand_spherical_vec(spaceRange, spaceRange, spaceRange);
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
		if (failed == attempts - 1 || collisionDetected > static_cast<int>(1.5 * static_cast<double>(n))) // Added the second part to speed up spatial constraint increase when there are clearly too many collisions for the space to be feasible.
		{
			std::cerr << "Failed " << spaceRange << ". Increasing range " << spaceRangeIncrement << "cm^3.\n";
			spaceRange += spaceRangeIncrement;
			failed = 0;
			for (unsigned int Ball = 0; Ball < n; Ball++)
			{
				g[Ball].pos = rand_spherical_vec(spaceRange, spaceRange, spaceRange); // Each time we fail and increase range, redistribute all balls randomly so we don't end up with big balls near mid and small balls outside.
			}
		}
		collisionDetected = 0;
	}
}

