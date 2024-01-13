#include "ball_group.hpp"
// #include "../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
#include <iomanip>
// #include <filesystem>
// namespace fs = std::filesystem;


// String buffers to hold data in memory until worth writing to file:
// std::stringstream ballBuffer;
// std::stringstream energyBuffer;
std::stringstream contactBuffer;
extern const int bufferlines;


// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
bool contact = false;
bool inital_contact = true;


// Prototypes
void
safetyChecks(Ball_group &O);
int 
check_restart(std::string folder);
Ball_group 
make_group(std::string argv1);
inline int 
twoDtoOneD(const int row, const int col, const int width);
void 
BPCA(std::string path, int num_balls);
void 
collider(std::string path, std::string projectileName,std::string targetName);
/// @brief The ballGroup run by the main sim looper.
// Ball_group O(output_folder, projectileName, targetName, v_custom); // Collision
// Ball_group O(path, targetName, 0);  // Continue
// std::cerr<<"genBalls: "<<genBalls<<std::endl;
// Ball_group O(20, true, v_custom); // Generate
// Ball_group O(genBalls, true, v_custom); // Generate
timey t;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int
main(const int argc, char const* argv[])
{
    t.start_event("WholeThing");
    // energyBuffer.precision(12);  // Need more precision on momentum.

    //make dummy ball group to read input file
    std::string location;
    Ball_group dummy(1);
    if (argc == 2)
    {
        location = std::string(argv[1]);
    }
    else
    {
        location = "";
    }
    dummy.parse_input_file(location);
    // O.zeroAngVel();
    // O.pushApart();

    // Normal sim:
    // O.sim_init_write(output_prefix);
    // sim_looper();
    BPCA(dummy.output_folder.c_str(),dummy.N);
    // collider(argv[1],dummy.projectileName,dummy.targetName);

    // collider(argv[1],projTarget,projTarget);
    
    t.end_event("WholeThing");
    t.print_events();
    t.save_events(dummy.output_folder + "timing.txt");
}  // end main
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void collider(std::string path, std::string projectileName, std::string targetName)
{
    t.start_event("collider");
    Ball_group O = Ball_group(path,std::string(projectileName),std::string(targetName));
    safetyChecks(O);
    O.sim_init_write();
    sim_looper(O,O.start_step);
    t.end_event("collider");
    O.freeMemory();
    return;
}

void BPCA(std::string path, int num_balls)
{
    int rest = -1;
    Ball_group O = Ball_group(path);  
    safetyChecks(O);
    if  (O.mid_sim_restart)
    {
        O.sim_looper(O.start_step);
    }
    // exit(0);
    // Add projectile: For dust formation BPCA
    for (int i = O.start_index; i < num_balls; i++) {
    // for (int i = 0; i < 250; i++) {
        // O.zeroAngVel();
        // O.zeroVel();
        contact = false;
        inital_contact = true;

        // t.start_event("add_projectile");
        O = O.add_projectile();
        // t.end_event("add_projectile");
        O.sim_init_write(i);
        O.sim_looper(1);
        simTimeElapsed = 0;
    }
    // O.freeMemory();
    return;
}




void
safetyChecks(Ball_group &O)
{
    titleBar("SAFETY CHECKS");

    if (O.soc <= 0) {
        fprintf(stderr, "\nSOC NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.v_collapse <= 0) {
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

    if (O.initial_radius <= 0) {
        fprintf(stderr, "\nCluster initialRadius not set\n");
        exit(EXIT_FAILURE);
    }


    for (int Ball = 0; Ball < O.num_particles; Ball++) {
        if (O.pos[Ball].norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA ball position is [0,0,0]. Possibly didn't initialize balls properly.\n");
            exit(EXIT_FAILURE);
        }

        if (O.R[Ball] <= 0) {
            fprintf(stderr, "\nA balls radius <= 0.\n");
            exit(EXIT_FAILURE);
        }

        if (O.m[Ball] <= 0) {
            fprintf(stderr, "\nA balls mass <= 0.\n");
            exit(EXIT_FAILURE);
        }
    }
    titleBar("SAFETY PASSED");
}


// void setGuidDT(const double& vel)
//{
//	// Guidos k and dt:
//	dt = .01 * O.getRmin() / fabs(vel);
//}
//
// void setGuidK(const double& vel)
//{
//	kin = O.getMassMax() * vel * vel / (.1 * O.R[0] * .1 * O.R[0]);
//	kout = cor * kin;
//}

inline int twoDtoOneD(const int row, const int col, const int width)
{
    return width * row + col;
}