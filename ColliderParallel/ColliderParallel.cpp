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
#include "ballSim.hpp"
// #include "../timing/timing.hpp"

namespace fs = std::filesystem;

void safetyChecks(P_group &O);
std::string check_restart(std::string folder,int* restart);
P_group make_group(const char *argv1,int* restart);
// P_group make_group(const char *argv1,int* restart,timey *t);

int main(const int argc, char const* argv[])
{
    int num_balls;
    int rest = -1;
    int *restart = &rest;
    timey main_t;
    // timey* main_t = &t;
    main_t.start_event("MAIN:Whole Thing");


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
        //This will be overwritten if the parameter is given in the input file
        num_balls = 100;
    }
 
	P_group O = make_group(argv[1],restart); // Particle system
	safetyChecks(O);

    if (*restart < num_balls)
    { 
    	for (int i = *restart; i < num_balls; i++) 
        {
            main_t.start_event("MAIN:Zeroing");
            // O.zeroAngVel();
            // O.zeroVel();
            main_t.end_event("MAIN:Zeroing");
            main_t.start_event("MAIN:Add Projectile");
            O.add_projectile();
            main_t.end_event("MAIN:Add Projectile");
            main_t.start_event("MAIN:Sim init write");
            O.sim_init_write(O.s_location, i);
            main_t.end_event("MAIN:Sim init write");
            main_t.start_event("MAIN:Sim Looper");
            O.sim_looper();
            main_t.end_event("MAIN:Sim Looper");
            O.simTimeElapsed = 0;
        }
    }
    else
    {
        std::cerr<<"Simulation already complete, now exiting."<<std::endl;
    }
    main_t.end_event("MAIN:Whole Thing");
    main_t.print_events();
    main_t.save_events(std::string(argv[1])+"timing.txt");
    O.write_timing();
}

/// @brief makes sure important parameters have been set and wont error the sim
void safetyChecks(P_group &O)
{
    titleBar("SAFETY CHECKS");

    if (O.soc <= 0) {
        fprintf(stderr, "\nSOC NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.vCollapse <= 0) {
        fprintf(stderr, "\nvCollapse NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.skip == 0) {
        fprintf(stderr, "\nSKIP NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.kin < 0) {
        fprintf(stderr, "\nSPRING CONSTANT NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.dt <= 0) {
        fprintf(stderr, "\nDT NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.steps == 0) {
        fprintf(stderr, "\nSTEPS NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.initialRadius <= 0) {
        fprintf(stderr, "\nCluster initialRadius not set\n");
        exit(EXIT_FAILURE);
    }

    for (int Ball = 0; Ball < O.num_particles; Ball++) {
        if (O.p_group[Ball].pos.norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA ball position is [0,0,0]. Possibly didn't initialize balls properly.\n");
            exit(EXIT_FAILURE);
        }

        if (O.p_group[Ball].R <= 0) {
            fprintf(stderr, "\nA balls radius <= 0.\n");
            exit(EXIT_FAILURE);
        }

        if (O.p_group[Ball].m <= 0) {
            fprintf(stderr, "\nA balls mass <= 0.\n");
            exit(EXIT_FAILURE);
        }
    }
    titleBar("SAFETY PASSED");
}


/// @brief checks if a simulation needs to be restarted or needs to be a new sim
/// @param folder is the path to the folder the sim is in
/// @param restart will be set to a value to indicate the need for a restart.
///       After running this function... 
///         if restart = -2 -> Finished sim
///         if restart = -1 -> Need to start new sim
///         if restart = 0  -> Only one step has completed and we need to restart a new sim
///         if restart > 0  -> Continue the sim starting from this step
/// @returns the name of the file that needs to be restarted, or an empty string if no restart is needed
std::string check_restart(std::string folder,int* restart)
{
    std::string file;
    int largest_file_index = -1;
    int file_index = -100;
    std::string largest_index_name;
    *restart = -1;
    for (const auto & entry : fs::directory_iterator(folder))
    {
        if (file.substr(0,file.size()-4) == "timing")
        {
            *restart = -2;
            return "";
        }

        file = entry.path();
        size_t pos = file.find_last_of("/");
        file = file.erase(0,pos+1);
        // tot_count++;
        if (file.substr(file.size()-4,file.size()) == ".csv")
        {
            std::string full_file = file;
            size_t pos = file.find('_');
            file = file.substr(0,pos);
            if (file.length() < 4 || file.substr(file.size()-4,file.size()) != ".csv")
            {
                file_index = stoi(file.substr(0,pos));
            }
            else
            {
                file_index = 0;
                // file = file;
            }

            if (file_index > largest_file_index)
            {
                largest_file_index = file_index;
                largest_index_name = full_file;
            }

        }
    }

    *restart = largest_file_index;
    if (*restart != -1)
    {
        size_t start,end;
        start = largest_index_name.find('_');
        end = largest_index_name.find_last_of('_');
        // //Delete most recent save file as this is likely only partially 
        // //complete if we are restarting

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

        // std::cout<< largest_index_name <<std::endl;
        // std::cout<<start<<", "<<end<<std::endl;
        if (*restart == 0)
        {
            return largest_index_name;
        }
        else
        {   
            return largest_index_name.substr(start,end-start+1);
        }
    }
    else
    {
        return "";
    }
}

//@brief Calls correct P_group constructor based on the need for a restart or not
P_group make_group(const char *argv1,int* restart)
{
    P_group O;

    //See if run has already been started
    std::string filename = check_restart(argv1,restart);

    if (*restart > 0) //Restart is necessary unless only first write has happended so far
    {
        O = P_group(argv1,*restart);
    }
    else if (*restart > -2) // Make new ball group
    {
        *restart = 0;
        O = P_group(true, argv1); // Generate new group
    }
    else
    {
        std::cerr<<"Simulation already complete.\n";
    }
    O.energyBuffer.precision(12);  // Need more precision on momentum.
    return O;
}
