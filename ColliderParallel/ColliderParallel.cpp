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

namespace fs = std::filesystem;


// String buffers to hold data in memory until worth writing to file:
// std::stringstream ballBuffer;
// std::stringstream energyBuffer;

// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
// const time_t start = time(nullptr);        // For end of program analysis
// time_t startProgress; // For progress reporting (gets reset)
// time_t lastWrite;     // For write control (gets reset)
// bool writeStep;       // This prevents writing to file every step (which is slow).

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



std::string check_restart(std::string folder,int* restart)
{
    std::string file;
    // int tot_count = 0;
    // int file_count = 0;
    int largest_file_index = -1;
    int file_index = -100;
    std::string largest_index_name;
    for (const auto & entry : fs::directory_iterator(folder))
    {
        file = entry.path();
        size_t pos = file.find_last_of("/");
        file = file.erase(0,pos+1);
        // tot_count++;
        if (file.substr(file.size()-4,file.size()) == ".csv")
        {
            size_t pos = file.find('_');
            file_index = stoi(file.substr(0,pos));
            if (file[pos+2] != '_')
            {
                file_index = 0;
            }
            if (file_index > largest_file_index)
            {
                largest_file_index = file_index;
                largest_index_name = file;
            }
        }
    }
    *restart = largest_file_index;
    if (*restart != -1)
    {
        size_t start,end;
        start = largest_index_name.find('_');
        end = largest_index_name.find_last_of('_');
        //Delete most recent save file as this is likely only partially 
        //complete if we are restarting

        std::string remove_file;

        if (*restart == 0)
        {
            remove_file = largest_index_name.substr(0,end+1);
        }
        else
        {
            remove_file = std::to_string(*restart) + largest_index_name.substr(start,end-start+1);
        }

        std::string file1 = folder + remove_file + "constants.csv";
        std::string file2 = folder + remove_file + "energy.csv";
        std::string file3 = folder + remove_file + "simData.csv";
        int status1 = remove(file1.c_str());
        int status2 = remove(file2.c_str());
        int status3 = remove(file3.c_str());

        if (status1 != 0)
        {
            std::cout<<"File: "<<file1<<" could not be removed, now exiting with failure."<<std::endl;
            exit(EXIT_FAILURE);
        }
        else if (status2 != 0)
        {
            std::cout<<"File: "<<file2<<" could not be removed, now exiting with failure."<<std::endl;
            exit(EXIT_FAILURE);
        }
        else if (status3 != 0)
        {
            std::cout<<"File: "<<file3<<" could not be removed, now exiting with failure."<<std::endl;
            exit(EXIT_FAILURE);
        }

        return largest_index_name.substr(start,end-start+1);
    }
    else
    {
        return "";
    }
}

//@brief sets Ball_group object based on the need for a restart or not
P_group make_group(const char *argv1,int* restart)
{
    P_group O;
    //See if run has already been started
    std::string filename = check_restart(argv1,restart);
    if (*restart > -1) //Restart is necessary unless only first write has happended so far
    {
        if (*restart > 1)
        {//TESTED
            (*restart)--;
            // filename = std::to_string(*restart) + filename;
            filename = filename.substr(1,filename.length());
            O = P_group(argv1,filename,*restart);
        }
        else if (*restart == 1) //restart from first write (different naming convension for first write)
        {//TESTED
            (*restart)--;
            filename = filename.substr(1,filename.length());
            // exit(EXIT_SUCCESS);
            O = P_group(argv1,filename,*restart);
        }
        else //if restart is 0, need to rerun whole thing
        {//TESTED
            O = P_group(true, argv1); // Generate new group
        }

    }
    else // Make new ball group
    {
        *restart = 0;
        O = P_group(true, argv1); // Generate new group
    }
    O.energyBuffer.precision(12);  // Need more precision on momentum.
    return O;
}

int main(const int argc, char const* argv[])
{
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
 
	P_group O = make_group(argv[1],restart); // Particle system

	for (int i = *restart; i < num_balls; i++) {
    // for (int i = 0; i < 250; i++) {
        O.zeroAngVel();
        O.zeroVel();
        // t.start_event("add_projectile");
        O.add_projectile();
        O.make_pairs();
        // t.end_event("add_projectile");
        O.sim_init_write(O.s_location, i);
        O.sim_looper();
        O.simTimeElapsed = 0;
    }

	// std::for_each(std::execution::par_unseq, O.p_group.begin(), O.p_group.end(), update_kinematics);
	// std::for_each(std::execution::par, pairs.begin(), pairs.end(), compute_acceleration);
	// //std::ranges::for_each( arr, [i=0](auto &e) mutable { long_function(e,i++); } );
	// if (writeStep)
	// {
	// 	ballBuffer << '\n'; // Prepares a new line for incoming data.
	// 	write_to_buffer(O);
	// }
}

