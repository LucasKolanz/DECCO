#include "ball_group.hpp"
// #include "../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
#include <iomanip>

#ifdef MPI_ENABLE
    #include <mpi.h>
#endif
// #include <filesystem>
// namespace fs = std::filesystem;


// String buffers to hold data in memory until worth writing to file:
// std::stringstream ballBuffer;
// std::stringstream energyBuffer;
// std::stringstream contactBuffer;
extern const int bufferlines;


// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
// const time_t start = time(nullptr);  // For end of program analysis
// time_t startProgress;                // For progress reporting (gets reset)
// time_t lastWrite;                    // For write control (gets reset)
// bool writeStep;                      // This prevents writing to file every step (which is slow).
// bool contact = false;
// bool inital_contact = true;


// Prototypes
void
sim_looper(Ball_group &O,unsigned long long start_step);
void
safetyChecks(Ball_group &O);
void 
BPCA(std::string path, int num_balls);
void 
relax(std::string path);
void 
collider(std::string path, std::string projectileName,std::string targetName);
int get_num_threads(Ball_group &O);
timey t;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int
main(int argc, char* argv[])
{
    t.start_event("WholeThing");
        // MPI Initialization
    int world_rank, world_size;
    #ifdef MPI_ENABLE
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    #else
        world_rank = 0;
        world_size = 1;
    #endif

    std::cerr<<"=========================================Start Simulation========================================="<<std::endl;
    //Verify we have all the nodes we asked for
    fprintf(
        stderr,
        "Hello from rank %d\n",
        world_rank);
    fflush(stderr);


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

    //verify OpenMP threads
    std::cerr<<"Max of "<<omp_get_max_threads()<<" threads on this machine."<<std::endl;

    //verify total time and frequency of writes
    std::cerr<<"simTimeSeconds: "<<dummy.attrs.simTimeSeconds<<std::endl;
    std::cerr<<"timeResolution: "<<dummy.attrs.timeResolution<<std::endl;

    std::string radiiDist;
    if (dummy.attrs.radiiDistribution == dummy.attrs.logNorm)
    {
        radiiDist = "logNormal";
    }
    else if (dummy.attrs.radiiDistribution == dummy.attrs.constant)
    {
        radiiDist = "constant";
    }

    std::cerr<<"Using "<<radiiDist<<" particle radii distribution"<<std::endl;

    if (dummy.attrs.typeSim == dummy.attrs.collider)
    {
        std::cerr<<"COLLIDER NOT IMPLIMENTED"<<std::endl;
        exit(-1);
        // #ifdef MPI_ENABLE
        //     MPI_Barrier(MPI_COMM_WORLD);
        // #endif
        // collider(argv[1],dummy.projectileName,dummy.targetName);
    }
    else if (dummy.attrs.typeSim == dummy.attrs.BPCA)
    {
        if (dummy.attrs.N >=  0)
        {
            #ifdef MPI_ENABLE
                MPI_Barrier(MPI_COMM_WORLD);
            #endif
            BPCA(dummy.attrs.output_folder,dummy.attrs.N);
        }
        else
        {
            std::cerr<<"ERROR: if simType is BPCA, N >= 0 must be true."<<std::endl;
        }
    }
    else if (dummy.attrs.typeSim == dummy.attrs.relax)
    {
        #ifdef MPI_ENABLE
            MPI_Barrier(MPI_COMM_WORLD);
        #endif
        relax(dummy.attrs.output_folder);
    }
    else
    {
        std::cerr<<"ERROR: input file needs to specify a simulation type (simType)."<<std::endl;
    }

    
    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    #endif

    std::cerr<<"=========================================Finish Simulation========================================="<<std::endl;

    t.end_event("WholeThing");
    t.print_events();
    t.save_events(dummy.attrs.output_folder + "timing.txt");
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
    sim_looper(O,O.attrs.start_step);
    t.end_event("collider");
    O.freeMemory();
    return;
}

void BPCA(std::string path, int num_balls)
{
    int world_rank = getRank();

    // int rest = -1;
    Ball_group O = Ball_group(path);  
    safetyChecks(O);
    if  (O.attrs.mid_sim_restart)
    {
        std::cerr<<"Asking for "<<get_num_threads(O)<<" threads."<<std::endl;
        sim_looper(O,O.attrs.start_step);
    }

    // Add projectile: For dust formation BPCA
    for (int i = O.attrs.start_index; i < num_balls; i++) {
        std::cerr<<"I: "<<i<<std::endl;
        // t.start_event("add_projectile");
        O = O.add_projectile();
        // t.end_event("add_projectile");
        if (world_rank == 0)
        {
            O.sim_init_write(i);
        }
        std::cerr<<"Asking for "<<get_num_threads(O)<<" threads."<<std::endl;
        sim_looper(O,1);
        O.attrs.simTimeElapsed = 0;
    }
    // O.freeMemory();
    return;
}

//Load file index, zero out velocity and angular velocity and run to let aggregate relax
void relax(std::string path)
{
    int world_rank = getRank();
    Ball_group O = Ball_group(path);  
    safetyChecks(O);

    if (world_rank == 0)
    {
        O.sim_init_write(O.attrs.relax_index);
    }
    sim_looper(O,1);
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
        threads = 32;
    }

    if (threads > O.attrs.MAXOMPthreads)
    {
        threads = O.attrs.MAXOMPthreads;
    }
    return threads;
}


void
sim_looper(Ball_group &O,unsigned long long start_step=1)
{
    int world_rank = getRank();

    O.attrs.num_writes = 0;
    unsigned long long Step;
    bool writeStep = false;

    if (world_rank == 0)
    {   
        std::cerr << "Beginning simulation...\n";
        O.attrs.startProgress = time(nullptr);

        std::cerr<<"start step: "<<start_step<<std::endl;
        std::cerr<<"Stepping through "<<O.attrs.steps<<" steps"<<std::endl;
        std::cerr<<"Simulating "<<O.attrs.simTimeSeconds<<" seconds per sim."<<std::endl;
        std::cerr<<"Writing out every "<<O.attrs.timeResolution<<" seconds."<<std::endl;
        std::cerr<<"For a total of "<<O.attrs.simTimeSeconds/O.attrs.timeResolution<<" timesteps saved per sim."<<std::endl;
    }


    //Set the number of threads to be appropriate
    O.attrs.OMPthreads = get_num_threads(O);

    for (Step = start_step; Step < O.attrs.steps; Step++)  // Steps start at 1 for non-restart because the 0 step is initial conditions.
    {
        // simTimeElapsed += dt; //New code #1
        // Check if this is a write step:
        if (Step % O.attrs.skip == 0) {
            t.start_event("writeProgressReport");
            writeStep = true;
            // std::cerr<<"Write step "<<Step<<std::endl;

            /////////////////////// Original code #1
            O.attrs.simTimeElapsed += O.attrs.dt * O.attrs.skip;
            ///////////////////////

            if (world_rank == 0)
            {
                // Progress reporting:
                float eta = ((time(nullptr) - O.attrs.startProgress) / static_cast<float>(O.attrs.skip) *
                             static_cast<float>(O.attrs.steps - Step)) /
                            3600.f;  // Hours.
                float real = (time(nullptr) - O.attrs.start) / 3600.f;
                float simmed = static_cast<float>(O.attrs.simTimeElapsed / 3600.f);
                float progress = (static_cast<float>(Step) / static_cast<float>(O.attrs.steps) * 100.f);
                fprintf(
                    stderr,
                    "%llu\t%2.0f%%\tETA: %5.2lf\tReal: %5.2f\tSim: %5.2f hrs\tR/S: %5.2f\n",
                    Step,
                    progress,
                    eta,
                    real,
                    simmed,
                    real / simmed);
                // fprintf(stdout, "%u\t%2.0f%%\tETA: %5.2lf\tReal: %5.2f\tSim: %5.2f hrs\tR/S: %5.2f\n", Step,
                // progress, eta, real, simmed, real / simmed);
                fflush(stdout);
                O.attrs.startProgress = time(nullptr);
            }
            t.end_event("writeProgressReport");
        } else {
            writeStep = O.attrs.debug;
        }

        // Physics integration step:
        O.sim_one_step(writeStep);

        if (writeStep) {
            // t.start_event("writeStep");
            // Write energy to stream:
            ////////////////////////////////////
            //TURN THIS ON FOR REAL RUNS!!!
            // O.energyBuffer = std::vector<double> (data->getWidth("energy"));
            // std::cerr<<"start,num_writes: "<<start<<','<<O.num_writes<<std::endl;
            if (world_rank == 0)
            {    
                int start = O.data->getWidth("energy")*(O.attrs.num_writes-1);
                O.energyBuffer[start] = O.attrs.simTimeElapsed;
                O.energyBuffer[start+1] = O.PE;
                O.energyBuffer[start+2] = O.KE;
                O.energyBuffer[start+3] = O.PE+O.KE;
                O.energyBuffer[start+4] = O.mom.norm();
                O.energyBuffer[start+5] = O.ang_mom.norm();

                if (Step / O.attrs.skip % 10 == 0) 
                {

                    std::cerr << "vMax = " << O.getVelMax() << " Steps recorded: " << Step / O.attrs.skip << '\n';
                    std::cerr << "Data Write to "<<O.data->getFileName()<<"\n";
                    
                    O.data->Write(O.ballBuffer,"simData",bufferlines);

                    O.ballBuffer.clear();
                    O.ballBuffer = std::vector<double>(O.data->getWidth("simData")*bufferlines);
                    O.data->Write(O.energyBuffer,"energy");
                    O.energyBuffer.clear();
                    O.energyBuffer = std::vector<double>(O.data->getWidth("energy")*bufferlines);

                    O.attrs.num_writes = 0;

                }  // Data export end
                
                O.attrs.lastWrite = time(nullptr);
            }
            
            // Reinitialize energies for next step:
            O.KE = 0;
            O.PE = 0;
            O.mom = {0, 0, 0};
            O.ang_mom = {0, 0, 0};

            if (O.attrs.dynamicTime) { O.calibrate_dt(Step, false); }
            // t.end_event("writeStep");
        }  // writestep end
    }


    if (world_rank == 0)
    {
        const time_t end = time(nullptr);

        std::cerr << "Simulation complete! \n"
                  << O.attrs.num_particles << " Particles and " << Step << '/' << O.attrs.steps << " Steps.\n"
                  << "Simulated time: " << O.attrs.steps * O.attrs.dt << " seconds\n"
                  << "Computation time: " << end - O.attrs.start << " seconds\n";
        std::cerr << "\n===============================================================\n";
    }


}  // end simLooper


void
safetyChecks(Ball_group &O) //Should be ready to call sim_looper
{
    titleBar("SAFETY CHECKS");


    if (O.attrs.output_folder == "")
    {
        fprintf(stderr, "\noutput_folder NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.data_directory == "")
    {
        fprintf(stderr, "\ndata_directory NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.soc <= 0) {
        fprintf(stderr, "\nSOC NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.v_collapse <= 0) {
        fprintf(stderr, "\nvCollapse NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.skip <= 1) {
        fprintf(stderr, "\nSKIP NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.radiiDistribution != O.attrs.constant && O.attrs.radiiDistribution != O.attrs.logNorm) {
        fprintf(stderr, "\nradiiDistribution NOT SET\n");
        exit(EXIT_FAILURE);
    } 

    if (O.attrs.typeSim != O.attrs.BPCA && O.attrs.typeSim != O.attrs.collider && O.attrs.typeSim != O.attrs.relax) {
        fprintf(stderr, "\ntypeSim NOT SET\n");
        exit(EXIT_FAILURE);
    } 

    if (O.attrs.kin < 0) {
        fprintf(stderr, "\nSPRING CONSTANT IN NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.kout < 0) {
        fprintf(stderr, "\nSPRING CONSTANT OUT NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.density < 0) {
        fprintf(stderr, "\ndensity NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.u_s < 0) {
        fprintf(stderr, "\nu_s NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.u_r < 0) {
        fprintf(stderr, "\nu_r NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.cor < 0) {
        fprintf(stderr, "\ncor NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.simTimeSeconds < 0) {
        fprintf(stderr, "\nsimTimeSeconds NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.timeResolution < 0) {
        fprintf(stderr, "\ntimeResolution NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.fourThirdsPiRho < 0) {
        fprintf(stderr, "\nfourThirdsPiRho NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.scaleBalls < 0) {
        fprintf(stderr, "\nscaleBalls NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.maxOverlap < 0) {
        fprintf(stderr, "\nmaxOverlap NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.kConsts < 0) {
        fprintf(stderr, "\nkConsts NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.impactParameter < 0) {
        fprintf(stderr, "\nimpactParameter NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.Ha < 0) {
        fprintf(stderr, "\nHa NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.genBalls < 0) {
        fprintf(stderr, "\ngenBalls NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.typeSim == O.attrs.BPCA && O.attrs.N < 0) {
        fprintf(stderr, "\nN NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.typeSim == O.attrs.relax && O.attrs.relax_index < 0) {
        fprintf(stderr, "\nrestart_index NOT SET and in relax mode\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.attempts < 0) {
        fprintf(stderr, "\nattempts NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.spaceRange < 0) {
        fprintf(stderr, "\nspaceRange NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.spaceRangeIncrement < 0) {
        fprintf(stderr, "\nspaceRangeIncrement NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.h_min < 0) {
        fprintf(stderr, "\nh_min NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.properties < 0) {
        fprintf(stderr, "\nproperties NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.dt <= 0) {
        fprintf(stderr, "\nDT NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.steps == 0) {
        fprintf(stderr, "\nSTEPS NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.initial_radius <= 0) {
        fprintf(stderr, "\nCluster initialRadius not set\n");
        exit(EXIT_FAILURE);
    }


    for (int Ball = 0; Ball < O.attrs.num_particles; Ball++) {
        if (O.pos[Ball].norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA ball position is [0,0,0]. Possibly didn't initialize positions properly.\n");
            exit(EXIT_FAILURE);
        }


        if (O.acc[Ball].norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA balls acc is [0,0,0]. Possibly didn't initialize acceleration properly.\n");
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

        if (O.moi[Ball] <= 0) {
            fprintf(stderr, "\nA balls moi <= 0.\n");
            exit(EXIT_FAILURE);
        }
    }

    titleBar("SAFETY PASSED");
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

inline int twoDtoOneD(const int row, const int col, const int width)
{
    return width * row + col;
}