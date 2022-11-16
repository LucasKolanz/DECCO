#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <vector>
#include <algorithm>
#include <execution>

#include "../vec3.hpp"
#include "initializations.hpp"

// String buffers to hold data in memory until worth writing to file:
std::stringstream ballBuffer;
std::stringstream energyBuffer;

// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
const time_t start = time(nullptr);        // For end of program analysis
time_t startProgress; // For progress reporting (gets reset)
time_t lastWrite;     // For write control (gets reset)
bool writeStep;       // This prevents writing to file every step (which is slow).

//ballGroup O(path, projectileName, targetName, vCustom); // Collision
//ballGroup O(path, targetName, 0); // Continue
//ballGroup O(genBalls, true, vCustom); // Generate



// void update_kinematics(P& P);
// void compute_acceleration(P_pair& p_pair);
// void compute_velocity(P& P);
// void write_to_buffer(P_group& p_group);



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




int main(const int argc, char const* argv[])
{
	energyBuffer.precision(12);  // Need more precision on momentum.
    int num_balls;
    int rest = -1;
    int *restart = &rest;


    // Runtime arguments:
    if (argc > 2) 
    {
        std::stringstream s(argv[2]);
        // s << argv[2];
        s >> num_balls;
        // numThreads = atoi(argv[1]);
        // fprintf(stderr,"\nThread count set to %i.\n", numThreads);
        // projectileName = argv[2];
        // targetName = argv[3];
        // KEfactor = atof(argv[4]);
    }
    else
    {
        num_balls = 100;
    }
    // P_group O = make_group(argv[1],restart);

	// int n = 10000; // Number of particles
	std::string here = "";
	P_group O(true,121,here.c_str()); // Particle system
	std::cout<<O.getMass()<<std::endl;
	std::string file = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/restart_test1/N_1000/T_3";
	// parse_input_file(file.c_str(), O);
	// parse_input_file(argv[1])
	// std::vector<P_pair> pairs = make_p_pairs(O);

	// for (int i = *restart; i < num_balls; i++) {
 //    // for (int i = 0; i < 250; i++) {
 //        O.zeroAngVel();
 //        O.zeroVel();
 //        t.start_event("add_projectile");
 //        O = O.add_projectile();
 //        t.end_event("add_projectile");
 //        O.sim_init_write(ori_output_prefix, i);
 //        sim_looper(O);
 //        simTimeElapsed = 0;
 //    }

	// std::for_each(std::execution::par_unseq, O.p_group.begin(), O.p_group.end(), update_kinematics);
	// std::for_each(std::execution::par, pairs.begin(), pairs.end(), compute_acceleration);
	// //std::ranges::for_each( arr, [i=0](auto &e) mutable { long_function(e,i++); } );
	// if (writeStep)
	// {
	// 	ballBuffer << '\n'; // Prepares a new line for incoming data.
	// 	write_to_buffer(O);
	// }
}