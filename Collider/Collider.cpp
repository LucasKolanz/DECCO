//TODO: Restarting a relax job double copies the constant file into the RELAXconstants file, doubling its size and throwing off
//      Porosity_FD.py. The double copy happens for RELAXenergy and RELAXsimData as well, but Porosity_FD only looks at the last line, 
//      which is still good data. Should probably just delete all RELAX files that already exist for a RELAX restart. 



// This file was originally written for SpaceLab/DECCO take care of initialization and calling the correct functions from ball_group.hpp to 
// carry out Soft Sphere, Discrete Element simulations 

// Authors: Job Guidos and Lucas Kolanz




#include "ball_group.hpp"
// #include "../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
#include <iomanip>
#include <omp.h>

#ifdef MPI_ENABLE
    #include <mpi.h>
#endif


//How many lines will be saved to a buffer before being saved
extern const int bufferlines;


// Prototypes
void
safetyChecks(Ball_group &O);
void 
runAggregation(std::string path, int num_balls);
void 
runRelax(std::string path);
void 
runCollider(std::string path);
timey t;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int
main(int argc, char* argv[])
{
        // MPI Initialization
    int world_rank, world_size;
    #ifdef MPI_ENABLE
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        if (world_rank == 0)
        {
            std::cerr<<"MPI enabled"<<std::endl;
        }
    #else
        world_rank = 0;
        world_size = 1;
    #endif

    if (world_rank == 0)
    {
        t.start_event("WholeThing");
        std::cerr<<"=========================================Start Simulation========================================="<<std::endl;
    }
    //Verify we have all the nodes we asked for
    fprintf(
        stderr,
        "Hello from rank %d\n",
        world_rank);
    fflush(stderr);
    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
    #endif


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

    //verify OpenMP threads
    std::cerr<<"Max of "<<omp_get_max_threads()<<" threads on rank "<<world_rank<<".\n";


    std::string radiiDist;
    if (dummy.attrs.radiiDistribution == logNorm)
    {
        radiiDist = "logNormal";
    }
    else if (dummy.attrs.radiiDistribution == constant)
    {
        radiiDist = "constant";
    }

    //verify total time and frequency of writes
    std::string message(
        "simTimeSeconds: "+dToSci(dummy.attrs.simTimeSeconds) + '\n' +
        "timeResolution: "+dToSci(dummy.attrs.timeResolution) + '\n' +
        "Using "+radiiDist+" particle radii distribution\n");
    MPIsafe_print(std::cerr,message);

    if (dummy.attrs.typeSim == collider)
    {
        // MPIsafe_print(std::cerr,"COLLIDER NOT IMPLIMENTED\n");
        // MPIsafe_exit(-1);
        #ifdef MPI_ENABLE
            MPI_Barrier(MPI_COMM_WORLD);
        #endif
        runCollider(argv[1]);
    }
    else if (dummy.attrs.typeSim == BPCA || dummy.attrs.typeSim == BCCA)
    {
        if (dummy.attrs.N >=  0)
        {
            #ifdef MPI_ENABLE
                MPI_Barrier(MPI_COMM_WORLD);
            #endif
            runAggregation(dummy.attrs.output_folder,dummy.attrs.N);
        }
        else
        {
            MPIsafe_print(std::cerr,"ERROR: if simType is BPCA, N >= 0 must be true.\n");
        }
    }
    else if (dummy.attrs.typeSim == relax)
    {
        #ifdef MPI_ENABLE
            MPI_Barrier(MPI_COMM_WORLD);
        #endif
        runRelax(dummy.attrs.output_folder);
    }
    else
    {
        MPIsafe_print(std::cerr,"ERROR: input file needs to specify a simulation type (simType).\n");
    }

    if (world_rank == 0)
    {
        std::cerr<<"=========================================Finish Simulation========================================="<<std::endl;
        t.end_event("WholeThing");
        t.print_events();
        t.save_events(dummy.attrs.output_folder + "timing.txt");
    }  // end main


    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    #endif
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//Initailizes and carries out Ball_group to collide two preexisting aggregates
void runCollider(std::string path)
{
    // t.start_event("collider");
    Ball_group O = Ball_group(path);
    std::cerr<<"Ball_group initiated"<<std::endl;
    safetyChecks(O);
    std::cerr<<"safety checked"<<std::endl;
    O.sim_init_write();
    std::cerr<<"ssim init wrote"<<std::endl;
    O.sim_looper(O.attrs.start_step);
    std::cerr<<"looper finished"<<std::endl;
    // t.end_event("collider");
    O.freeMemory();
    return;
}

//Initializes and carries out BPCA aggregate growth
void runAggregation(std::string path, int num_balls)
{
    int world_rank = getRank();

    Ball_group O = Ball_group(path);  
    safetyChecks(O);
    std::string message;
    if  (O.attrs.mid_sim_restart)
    {
        message = "Asking for "+std::to_string(O.get_num_threads())+" threads.\n";
        MPIsafe_print(std::cerr,message);
        O.sim_looper(O.attrs.start_step);
    }

    int increment = -1;

    // Add projectile: For dust formation BPCA or BCCA
    for (int i = O.attrs.start_index; i < num_balls; i+=increment) 
    {
        // t.start_event("add_projectile");
        O = O.add_projectile(O.attrs.typeSim);
        increment = O.attrs.num_particles;
        // t.end_event("add_projectile");
        if (world_rank == 0)
        {
            std::cerr<<"I: "<<i<<std::endl;
            O.sim_init_write(i);
            std::cerr<<"Asking for "<<O.get_num_threads()<<" threads."<<std::endl;
        }

        O.sim_looper(1);
        O.attrs.simTimeElapsed = 0;

        if (increment == -1)
        {
            MPIsafe_print(std::cerr,"ERROR: increment is -1");
            MPIsafe_exit(EXIT_FAILURE);
        }
    }
    // O.freeMemory();
    return;
}

//Initalizes and carries out a relaxation run
void runRelax(std::string path)
{
    int world_rank = getRank();
    Ball_group O = Ball_group(path);  
    safetyChecks(O);

    if (world_rank == 0)
    {
        O.sim_init_write(O.attrs.relax_index);
    }
    O.sim_looper(1);
}

// // Function to calculate the closest power of 2 to a given number.
// int closestPowerOf2(double number) 
// {
//     int lower = pow(2, floor(log2(number))); // Next lower power of 2
//     int higher = pow(2, ceil(log2(number))); // Next higher power of 2

//     // Return the power of 2 that is closer to the number
//     return (number - lower < higher - number) ? lower : higher;
// }


//with a known slope and intercept, givin N, the number of particles, what is the 
//optimum number of threads. The function then chooses the power of 2 that is closest
//to this optimum
int get_num_threads(Ball_group &O)
{
    int N = O.attrs.num_particles;
    // //This is from speed tests on COSINE
    // double slope = ;
    // double intercept = ;

    // double interpolatedValue = slope * n + intercept; // Linear interpolation
    // return std::min(closestPowerOf2(interpolatedValue),O.attrs.MAXOMPthreads);        // Find the closest power of 2

    //I could only test up to 16 threads so far. Not enough data for linear interp
    int threads;
    if (N < 0)
    {
        std::cerr<<"ERROR: negative number of particles."<<std::endl;
        exit(-1);
    }
    else if (N < 80)
    {
        threads = 1;
    }
    else if (N < 100)
    {
        threads = 2;
    }
    else if (N < 175)
    {
        threads = 24;
    }

    if (threads > O.attrs.MAXOMPthreads)
    {
        threads = O.attrs.MAXOMPthreads;
    }
    return threads;
}






void
safetyChecks(Ball_group &O) //Should be ready to call sim_looper
{
    #ifdef MPI_ENABLE
        std::cerr<<std::flush;
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
    if (getRank() == 0)
    {
        titleBar("SAFETY CHECKS");
        std::cerr<<std::flush;
    }
    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
    


    if (O.attrs.output_folder == "")
    {
        fprintf(stderr, "\noutput_folder NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.data_directory == "")
    {
        fprintf(stderr, "\ndata_directory NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.soc <= 0) {
        fprintf(stderr, "\nSOC NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.v_collapse <= 0) {
        fprintf(stderr, "\nvCollapse NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.skip <= 1) {
        fprintf(stderr, "\nSKIP NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.radiiDistribution != constant && O.attrs.radiiDistribution != logNorm) {
        fprintf(stderr, "\nradiiDistribution NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    } 

    if (O.attrs.typeSim != BPCA && O.attrs.typeSim != BCCA && O.attrs.typeSim != collider && O.attrs.typeSim != relax) {
        fprintf(stderr, "\ntypeSim NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    } 

    if (O.attrs.kin < 0) {
        fprintf(stderr, "\nSPRING CONSTANT IN NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.kout < 0) {
        fprintf(stderr, "\nSPRING CONSTANT OUT NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.density < 0) {
        fprintf(stderr, "\ndensity NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.u_s < 0) {
        fprintf(stderr, "\nu_s NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.u_r < 0) {
        fprintf(stderr, "\nu_r NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.cor < 0) {
        fprintf(stderr, "\ncor NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.simTimeSeconds < 0) {
        fprintf(stderr, "\nsimTimeSeconds NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.timeResolution < 0) {
        fprintf(stderr, "\ntimeResolution NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.fourThirdsPiRho < 0) {
        fprintf(stderr, "\nfourThirdsPiRho NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.scaleBalls < 0) {
        fprintf(stderr, "\nscaleBalls NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.maxOverlap < 0) {
        fprintf(stderr, "\nmaxOverlap NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.kConsts < 0) {
        fprintf(stderr, "\nkConsts NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.impactParameter < 0) {
        fprintf(stderr, "\nimpactParameter NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.Ha < 0) {
        fprintf(stderr, "\nHa NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.genBalls < 0) {
        fprintf(stderr, "\ngenBalls NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.typeSim == BPCA && O.attrs.N < 0) {
        fprintf(stderr, "\nN NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.typeSim == relax && O.attrs.relax_index < 0) {
        fprintf(stderr, "\nrestart_index NOT SET and in relax mode\n");
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.attempts < 0) {
        fprintf(stderr, "\nattempts NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.spaceRange < 0) {
        fprintf(stderr, "\nspaceRange NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.spaceRangeIncrement < 0) {
        fprintf(stderr, "\nspaceRangeIncrement NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.h_min < 0) {
        fprintf(stderr, "\nh_min NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.properties < 0) {
        fprintf(stderr, "\nproperties NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.dt <= 0) {
        fprintf(stderr, "\nDT NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.steps == 0) {
        fprintf(stderr, "\nSTEPS NOT SET for rank %1d\n",getRank());
        MPIsafe_exit(EXIT_FAILURE);
    }

    if (O.attrs.initial_radius <= 0) {
        fprintf(stderr, "\nCluster initialRadius not set\n");
        MPIsafe_exit(EXIT_FAILURE);
    }


    for (int Ball = 0; Ball < O.attrs.num_particles; Ball++) {
        if (O.pos[Ball].norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA ball position is [0,0,0]. Possibly didn't initialize positions properly for rank %1d\n",getRank());
            MPIsafe_exit(EXIT_FAILURE);
        }


        if (O.acc[Ball].norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA balls acc is [0,0,0]. Possibly didn't initialize acceleration properly for rank %1d\n",getRank());
            MPIsafe_exit(EXIT_FAILURE);
        }

        if (O.R[Ball] <= 0) {
            fprintf(stderr, "\nA balls radius <= 0 for rank %1d\n",getRank());
            MPIsafe_exit(EXIT_FAILURE);
        }

        if (O.m[Ball] <= 0) {
            fprintf(stderr, "\nA balls mass <= 0 for rank %1d\n",getRank());
            MPIsafe_exit(EXIT_FAILURE);
        }

        if (O.moi[Ball] <= 0) {
            fprintf(stderr, "\nA balls moi <= 0 for rank %1d\n",getRank());
            MPIsafe_exit(EXIT_FAILURE);
        }
    }

    fprintf(stderr, "SAFETY PASSED rank %1d\n",getRank());
    fflush(stderr);

    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
    #endif
}


// void setGuidDT(const double& vel)
//{
//  // Guidos k and dt:
//  dt = .01 * O.getRmin() / fabs(vel);
//}
//
// void setGuidK(const double& vel)
//{
//  kin = O.getMassMax() * vel * vel / (.1 * O.R[0] * .1 * O.R[0]);
//  kout = cor * kin;
//}

// inline int twoDtoOneD(const int row, const int col, const int width)
// {
//     return width * row + col;
// }