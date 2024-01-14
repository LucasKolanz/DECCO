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
Ball_group 
make_group(std::string argv1);
void
sim_looper(Ball_group &O,int start_step);

void 
BPCA(std::string path, int num_balls);
void 
collider(std::string path, std::string projectileName,std::string targetName);

timey t;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int
main(int argc, char* argv[])
{
    t.start_event("WholeThing");
    // energyBuffer.precision(12);  // Need more precision on momentum.

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

    //Verify we have all the nodes we asked for
    fprintf(
        stderr,
        "Hello from rank %d\n",
        world_rank);
    fflush(stderr);

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


    //Do whichever type of sim they asked for
    if (dummy.attrs.typeSim == dummy.attrs.collider)
    {
        std::cerr<<"ERROR: collider not implimented yet."<<std::endl;
        exit(EXIT_FAILURE);
        // #ifdef MPI_ENABLE
        //     MPI_Barrier(MPI_COMM_WORLD);
        // #endif
        // collider(argv[1],dummy.projectileName,dummy.targetName);
    }
    else if (dummy.attrs.typeSim == dummy.attrs.BPCA)
    {
        if (dummy.attrs.N >= 0)
        {
            #ifdef MPI_ENABLE
                MPI_Barrier(MPI_COMM_WORLD);
            #endif
            
            BPCA(dummy.attrs.output_folder.c_str(),dummy.attrs.N);
        }
        else
        {
            std::cerr<<"ERROR: if simType is BPCA, N >= 0 must be true."<<std::endl;
        }
    }
    else
    {
        std::cerr<<"ERROR: input file needs to specify a simulation type (simType)."<<std::endl;
    }





    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    #endif
    
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
    int rest = -1;
    Ball_group O = Ball_group(path);  
    safetyChecks(O);
    if  (O.attrs.mid_sim_restart)
    {
        sim_looper(O,O.attrs.start_step);
    }

    // Add projectile: For dust formation BPCA
    for (int i = O.attrs.start_index; i < num_balls; i++) 
    {
        // contact = false;
        // inital_contact = true;

        // t.start_event("add_projectile");
        O = O.add_projectile();
        // t.end_event("add_projectile");
        if (world_rank == 0)
        {
            O.sim_init_write(i);
        }
        sim_looper(O,1);
        O.attrs.simTimeElapsed = 0;
    }
    // O.freeMemory();
    return;
}




void
safetyChecks(Ball_group &O)
{
    titleBar("SAFETY CHECKS");

    if (O.attrs.soc <= 0) {
        fprintf(stderr, "\nSOC NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.v_collapse <= 0) {
        fprintf(stderr, "\nvCollapse NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.skip == 0) {
        fprintf(stderr, "\nSKIP NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (O.attrs.kin < 0) {
        fprintf(stderr, "\nSPRING CONSTANT NOT SET\n");
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


void sim_looper(Ball_group &O,int start_step=1)
{
    bool writeStep = false;
    O.attrs.num_writes = 0;
    std::cerr << "Beginning simulation...\n";

    std::cerr<<"start step: "<<start_step<<std::endl;

    O.attrs.startProgress = time(nullptr);

    std::cerr<<"Stepping through "<<O.attrs.steps<<" steps"<<std::endl;

    int Step;

    for (Step = start_step; Step < O.attrs.steps; Step++)  // Steps start at 1 for non-restart because the 0 step is initial conditions.
    {
        // simTimeElapsed += dt; //New code #1
        // Check if this is a write step:
        if (Step % O.attrs.skip == 0) {
            // t.start_event("writeProgressReport");
            writeStep = true;
            // std::cerr<<"Write step "<<Step<<std::endl;

            /////////////////////// Original code #1
            O.attrs.simTimeElapsed += O.attrs.dt * O.attrs.skip;
            ///////////////////////

            // Progress reporting:
            float eta = ((time(nullptr) - O.attrs.startProgress) / static_cast<float>(O.attrs.skip) *
                         static_cast<float>(O.attrs.steps - Step)) /
                        3600.f;  // Hours.
            float real = (time(nullptr) - O.attrs.start) / 3600.f;
            float simmed = static_cast<float>(O.attrs.simTimeElapsed / 3600.f);
            float progress = (static_cast<float>(Step) / static_cast<float>(O.attrs.steps) * 100.f);
            fprintf(
                stderr,
                "%u\t%2.0f%%\tETA: %5.2lf\tReal: %5.2f\tSim: %5.2f hrs\tR/S: %5.2f\n",
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
            // t.end_event("writeProgressReport");
        } else {
            writeStep = O.attrs.debug;
        }

        // Physics integration step:
        ///////////
        // if (write_all)
        // {
        //     zeroSaveVals();
        // }
        ///////////
        O.sim_one_step_single_core(writeStep);

        if (writeStep) {
            // t.start_event("writeStep");
            // Write energy to stream:
            ////////////////////////////////////
            //TURN THIS ON FOR REAL RUNS!!!
            int start = O.data->getWidth("energy")*(O.attrs.num_writes-1);
            O.energyBuffer[start] = O.attrs.simTimeElapsed;
            O.energyBuffer[start+1] = O.PE;
            O.energyBuffer[start+2] = O.KE;
            O.energyBuffer[start+3] = O.PE+O.KE;
            O.energyBuffer[start+4] = O.mom.norm();
            O.energyBuffer[start+5] = O.ang_mom.norm();



            // Reinitialize energies for next step:
            O.KE = 0;
            O.PE = 0;
            O.mom = {0, 0, 0};
            O.ang_mom = {0, 0, 0};


            // Data Export. Exports every 10 writeSteps (10 new lines of data) and also if the last write was
            // a long time ago.
            // if (time(nullptr) - lastWrite > 1800 || Step / skip % 10 == 0) {
            if (Step / O.attrs.skip % 10 == 0) {
                // Report vMax:

                std::cerr << "vMax = " << O.getVelMax() << " Steps recorded: " << Step / O.attrs.skip << '\n';
                std::cerr << "Data Write to "<<O.attrs.output_folder<<"\n";
                // std::cerr<<"output_prefix: "<<output_prefix<<std::endl;
                
                O.data->Write(O.ballBuffer,"simData",bufferlines);
                O.ballBuffer.clear();
                O.ballBuffer = std::vector<double>(O.data->getWidth("simData")*bufferlines);
                O.data->Write(O.energyBuffer,"energy");
                O.energyBuffer.clear();
                O.energyBuffer = std::vector<double>(O.data->getWidth("energy")*bufferlines);

                O.attrs.num_writes = 0;
                O.attrs.lastWrite = time(nullptr);

                // if (num_particles > 5)
                // {
                //     std::cerr<<"EXITING, step: "<<Step<<std::endl;
                //     exit(0);
                // }

            }  // Data export end


            if (O.attrs.dynamicTime) { O.calibrate_dt(Step, false); }
            // t.end_event("writeStep");
        }  // writestep end
    }

    const time_t end = time(nullptr);

    std::cerr << "Simulation complete! \n"
              << O.attrs.num_particles << " Particles and " << Step << '/' << O.attrs.steps << " Steps.\n"
              << "Simulated time: " << O.attrs.steps * O.attrs.dt << " seconds\n"
              << "Computation time: " << end - O.attrs.start << " seconds\n";
    std::cerr << "\n===============================================================\n";


}  // end simLooper