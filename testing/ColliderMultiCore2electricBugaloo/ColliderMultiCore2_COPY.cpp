#include "../ball_group.hpp"
// #include "../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
#include <iomanip>
#include <filesystem>
namespace fs = std::filesystem;


// String buffers to hold data in memory until worth writing to file:
std::stringstream ballBuffer;
std::stringstream energyBuffer;
std::stringstream contactBuffer;


// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
const time_t start = time(nullptr);  // For end of program analysis
time_t startProgress;                // For progress reporting (gets reset)
time_t lastWrite;                    // For write control (gets reset)
bool writeStep;                      // This prevents writing to file every step (which is slow).
bool contact = false;
bool inital_contact = true;


// Prototypes
void
sim_one_step(const bool write_step, Ball_group &O);
void
sim_looper(Ball_group &O);
void
safetyChecks(Ball_group &O);
std::string 
check_restart(std::string folder,int* restart);
Ball_group 
make_group(const char *argv1,int* restart);
inline int 
twoDtoOneD(const int row, const int col, const int width);
void 
BPCA(const char *path, int num_balls);
void 
collider(const char *path, std::string projectileName,std::string targetName);
/// @brief The ballGroup run by the main sim looper.
// Ball_group O(output_folder, projectileName, targetName, v_custom); // Collision
// Ball_group O(path, targetName, 0);  // Continue
// std::cerr<<"genBalls: "<<genBalls<<std::endl;
// Ball_group O(20, true, v_custom); // Generate
// Ball_group O(genBalls, true, v_custom); // Generate
// timey t;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int
main(const int argc, char const* argv[])
{
    // t.start_event("WholeThing");
    energyBuffer.precision(12);  // Need more precision on momentum.
    int num_balls;
    
    
    
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

    // Ball_group dummy(1);
    // dummy.parse_input_file(argv[1]);
    // O.zeroAngVel();
    // O.pushApart();

    // Normal sim:
    // O.sim_init_write(output_prefix);
    // sim_looper();
    BPCA(argv[1],num_balls);
    // std::string proj = "15_2_R2e-05_v0e+00_cor0.63_mu0.1_rho2.25_k3e-15_Ha5e-12_dt5e-10_";
    // std::string targ = "15_2_R2e-05_v0e+00_cor0.63_mu0.1_rho2.25_k3e-15_Ha5e-12_dt5e-10_"; 
    // collider(argv[1],proj,targ);

    // collider(argv[1],projTarget,projTarget);
    
    // t.end_event("WholeThing");
    // t.print_events();
    // t.save_events(output_folder + "timing.txt");
}  // end main
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void collider(const char *path, std::string projectileName, std::string targetName)
{
    // t.start_event("collider");
    Ball_group O = Ball_group(std::string(path),std::string(projectileName),std::string(targetName));
    safetyChecks(O);
    O.sim_init_write(output_prefix);
    sim_looper(O);
    // t.end_event("collider");
    // O.freeMemory();
    return;
}

void BPCA(const char *path, int num_balls)
{
    int rest = -1;
    int *restart = &rest;
    Ball_group O = make_group(path,restart);    
    safetyChecks(O);
    // Add projectile: For dust formation BPCA
    std::string ori_output_prefix = output_prefix;
    for (int i = *restart; i < num_balls; i++) {
    // for (int i = 0; i < 250; i++) {
        // O.zeroAngVel();
        // O.zeroVel();
        contact = false;
        inital_contact = true;
        // t.start_event("add_projectile");
        O = O.add_projectile();
        // t.end_event("add_projectile");
        O.sim_init_write(ori_output_prefix, i);
        sim_looper(O);
        simTimeElapsed = 0;
    }
    // O.freeMemory();
    return;
}


//@brief sets Ball_group object based on the need for a restart or not
Ball_group make_group(const char *argv1,int* restart)
{
    Ball_group O;
    
    //See if run has already been started
    std::string filename = check_restart(argv1,restart);
    if (*restart > -1) //Restart is necessary unless only first write has happended so far
    {
        if (*restart > 1)
        {//TESTED
            (*restart)--;
            // filename = std::to_string(*restart) + filename;
            filename = filename.substr(1,filename.length());
            O = Ball_group(argv1,filename,v_custom,*restart);
        }
        else if (*restart == 1) //restart from first write (different naming convension for first write)
        {//TESTED
            (*restart)--;
            filename = filename.substr(1,filename.length());
            // exit(EXIT_SUCCESS);
            O = Ball_group(argv1,filename,v_custom,*restart);
        }
        else //if restart is 0, need to rerun whole thing
        {//TESTED
            O = Ball_group(true, v_custom, argv1); // Generate new group
        }

    }
    else if (*restart == -1) // Make new ball group
    {
        *restart = 0;
        O = Ball_group(true, v_custom, argv1); // Generate new group
    }
    else
    {
        std::cerr<<"Simulation already complete.\n";
    }
    return O;
}

// @brief checks if this is new job or restart
std::string check_restart(std::string folder,int* restart)
{
    std::string file;
    // int tot_count = 0;
    // int file_count = 0;
    int largest_file_index = -1;
    int file_index;
    std::string largest_index_name;
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
            // file_count++;
            if (file[3] == '_')
            {
                file_index = stoi(file.substr(0,file.find("_")));
            }
            else if (file[1] == '_' and file[3] != '_')
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




void
sim_looper(Ball_group &O)
{
    std::cerr << "Beginning simulation...\n";

    startProgress = time(nullptr);

    std::cerr<<"Stepping through "<<steps<<" steps"<<std::endl;

    int maxThreads;

    // Get the maximum number of threads
    #pragma omp parallel
    {
        #pragma omp single
        maxThreads = omp_get_max_threads();
    }

    std::cerr<<"Max threads: "<<maxThreads<<std::endl;

    for (int Step = 1; Step < steps; Step++)  // Steps start at 1 because the 0 step is initial conditions.
    {
        // simTimeElapsed += dt; //New code #1
        // Check if this is a write step:
        if (Step % skip == 0) {
            // t.start_event("writeProgressReport");
            writeStep = true;

            /////////////////////// Original code #1
            simTimeElapsed += dt * skip;
            ///////////////////////

            // Progress reporting:
            float eta = ((time(nullptr) - startProgress) / static_cast<float>(skip) *
                         static_cast<float>(steps - Step)) /
                        3600.f;  // Hours.
            float real = (time(nullptr) - start) / 3600.f;
            float simmed = static_cast<float>(simTimeElapsed / 3600.f);
            float progress = (static_cast<float>(Step) / static_cast<float>(steps) * 100.f);
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
            startProgress = time(nullptr);
            // t.end_event("writeProgressReport");
        } else {
            writeStep = O.debug;
        }

        // Physics integration step:
        ///////////
        if (O.write_all)
        {
            O.zeroSaveVals();
        }
        ///////////
        O.sim_one_step(writeStep);

        ///////////////////////////////////////
        //TAKE THIS OUT FOR REAL RUNS
        // energyBuffer << '\n'
        //              << simTimeElapsed << ',' << O.PE << ',' << O.KE << ',' << O.PE + O.KE << ','
        //              << O.mom.norm() << ','
        //              << O.ang_mom.norm();  // the two zeros are bound and unbound mass

        // // Reinitialize energies for next step:
        // O.KE = 0;
        // O.PE = 0;
        // O.mom = {0, 0, 0};
        // O.ang_mom = {0, 0, 0};
        // // unboundMass = 0;
        // // boundMass = massTotal;
        // if (Step % 10 == 0)
        // {
        //     std::ofstream energyWrite;
        //     energyWrite.open(output_folder + output_prefix + "energy.csv", std::ofstream::app);
        //     energyWrite << energyBuffer.rdbuf();
        //     energyBuffer.str("");  // Empty the stream for next filling.
        //     energyWrite.close();
        // }
        ///////////////////////////////////////
        //////////////////////////////////////
        if (contact and inital_contact)
        {
            inital_contact = false;
            std::ofstream contact_write;
            contact_write.open(output_folder + "contact.csv", std::ofstream::app);
            contact_write << contactBuffer.rdbuf();  // Barf buffer to file.
            contactBuffer.str("");               // Empty the stream for next filling.
            contact_write.close();
        }
        //////////////////////////////////////
        if (O.write_all)
        {
            //////////////////////////////////////
            std::ofstream vdwWrite;
            vdwWrite.open(output_folder + output_prefix + "vdwData.csv", std::ofstream::app);
            vdwWrite<<O.vdwForce[0][0]<<','<<O.vdwForce[0][1]<<','<<O.vdwForce[0][0];
            for (int i = 1; i < O.num_particles; ++i)
            {
                vdwWrite<<','<<O.vdwForce[i][0]<<','<<O.vdwForce[i][1]<<','<<O.vdwForce[i][2];
            }
            vdwWrite << '\n';  
            vdwWrite.close();

            std::ofstream distWrite;
            distWrite.open(output_folder + output_prefix + "distData.csv", std::ofstream::app);
            int dist_length = (O.num_particles * O.num_particles / 2) - (O.num_particles / 2);
            distWrite<<O.distances[0];
            for (int i = 1; i < dist_length; ++i)
            {
                distWrite<<','<<O.distances[i];
            }
            distWrite << '\n';  
            distWrite.close();

            // std::ofstream elasticWrite;
            // elasticWrite.open(output_folder + output_prefix + "elasticData.csv", std::ofstream::app);
            // elasticWrite<<O.elasticForce[0][0]<<','<<O.elasticForce[0][1]<<','<<O.elasticForce[0][2];
            // for (int i = 1; i < O.num_particles; ++i)
            // {
            //     elasticWrite<<','<<O.elasticForce[i][0]<<','<<O.elasticForce[i][1]<<','<<O.elasticForce[i][2];
            // }
            // elasticWrite << '\n';  
            // elasticWrite.close();

            // std::ofstream slideWrite;
            // slideWrite.open(output_folder + output_prefix + "slideData.csv", std::ofstream::app);
            // slideWrite<<O.slideForce[0][0]<<','<<O.slideForce[0][1]<<','<<O.slideForce[0][2];
            // for (int i = 1; i < O.num_particles*O.num_particles; ++i)
            // {
            //     slideWrite<<','<<O.slideForce[i][0]<<','<<O.slideForce[i][1]<<','<<O.slideForce[i][2];
            // }
            // slideWrite << '\n';  
            // slideWrite.close();

            // std::ofstream rollWrite;
            // rollWrite.open(output_folder + output_prefix + "rollData.csv", std::ofstream::app);
            // rollWrite<<O.rollForce[0][0]<<','<<O.rollForce[0][1]<<','<<O.rollForce[0][2];
            // for (int i = 1; i < O.num_particles*O.num_particles; ++i)
            // {
            //     rollWrite<<','<<O.rollForce[i][0]<<','<<O.rollForce[i][1]<<','<<O.rollForce[i][2];
            // }
            // rollWrite << '\n';  
            // rollWrite.close();

            // std::ofstream torqueWrite;
            // torqueWrite.open(output_folder + output_prefix + "torqueData.csv", std::ofstream::app);
            // torqueWrite<<O.torqueForce[0][0]<<','<<O.torqueForce[0][1]<<','<<O.torqueForce[0][2];
            // for (int i = 1; i < O.num_particles*O.num_particles; ++i)
            // {
            //     torqueWrite<<','<<O.torqueForce[i][0]<<','<<O.torqueForce[i][1]<<','<<O.torqueForce[i][2];
            // }
            // torqueWrite << '\n';  
            // torqueWrite.close();

            // std::ofstream wWrite;
            // wWrite.open(output_folder + output_prefix + "wData.csv", std::ofstream::app);
            // wWrite<<O.w[0][0]<<','<<O.w[0][1]<<','<<O.w[0][2];
            // for (int i = 1; i < O.num_particles; ++i)
            // {
            //     wWrite<<','<<O.w[i][0]<<','<<O.w[i][1]<<','<<O.w[i][2];
            // }
            // wWrite << '\n';  
            // wWrite.close();

            // std::ofstream posWrite;
            // posWrite.open(output_folder + output_prefix + "posData.csv", std::ofstream::app);
            // posWrite<<O.pos[0][0]<<','<<O.pos[0][1]<<','<<O.pos[0][2];
            // for (int i = 1; i < O.num_particles; ++i)
            // {
            //     posWrite<<','<<O.pos[i][0]<<','<<O.pos[i][1]<<','<<O.pos[i][2];
            // }
            // posWrite << '\n';  
            // posWrite.close();

            // std::ofstream velWrite;
            // velWrite.open(output_folder + output_prefix + "velData.csv", std::ofstream::app);
            // velWrite<<O.vel[0][0]<<','<<O.vel[0][1]<<','<<O.vel[0][2];
            // for (int i = 1; i < O.num_particles; ++i)
            // {
            //     velWrite<<','<<O.vel[i][0]<<','<<O.vel[i][1]<<','<<O.vel[i][2];
            // }
            // velWrite << '\n';  
            // velWrite.close();
        }
        //////////////////////////////////////
        

        if (writeStep) {
            // t.start_event("writeStep");
            // Write energy to stream:
            ////////////////////////////////////
            //TURN THIS ON FOR REAL RUNS!!!
            energyBuffer << '\n'
                         << simTimeElapsed << ',' << O.PE << ',' << O.KE << ',' << O.PE + O.KE << ','
                         << O.mom.norm() << ','
                         << O.ang_mom.norm();  // the two zeros are bound and unbound mass

            // Reinitialize energies for next step:
            O.KE = 0;
            O.PE = 0;
            O.mom = {0, 0, 0};
            O.ang_mom = {0, 0, 0};
            //unboundMass = 0;
            //boundMass = massTotal;
            ////////////////////////////////////

            // Data Export. Exports every 10 writeSteps (10 new lines of data) and also if the last write was
            // a long time ago.
            // if (time(nullptr) - lastWrite > 1800 || Step / skip % 10 == 0) {
            if (Step / skip % 10 == 0) {
                // Report vMax:

                std::cerr << "vMax = " << O.getVelMax() << " Steps recorded: " << Step / skip << '\n';
                std::cerr << "Data Write to "<<output_folder<<"\n";
                // std::cerr<<"output_prefix: "<<output_prefix<<std::endl;


                // Write simData to file and clear buffer.
                std::ofstream ballWrite;
                ballWrite.open(output_folder + output_prefix + "simData.csv", std::ofstream::app);
                ballWrite << ballBuffer.rdbuf();  // Barf buffer to file.
                ballBuffer.str("");               // Empty the stream for next filling.
                ballWrite.close();

                // Write Energy data to file and clear buffer.
                ////////////////////////////////////////////
                //TURN ON FOR REAL SIM
                std::ofstream energyWrite;
                energyWrite.open(output_folder + output_prefix + "energy.csv", std::ofstream::app);
                energyWrite << energyBuffer.rdbuf();
                energyBuffer.str("");  // Empty the stream for next filling.
                energyWrite.close();
                // ////////////////////////////////////////////

                lastWrite = time(nullptr);
            }  // Data export end


            if (dynamicTime) { O.calibrate_dt(Step, false); }
            // t.end_event("writeStep");
        }  // writestep end
    }

    if (true)
    {
        for (int i = 0; i < O.num_particles; i++)
        {
            std::cerr<<"===================================="<<std::endl;
            std::cerr<<O.pos[i]<<std::endl;
            std::cerr<<O.vel[i]<<std::endl;
            std::cerr<<"===================================="<<std::endl;
        }
    }

    const time_t end = time(nullptr);

    std::cerr << "Simulation complete!\n"
              << O.num_particles << " Particles and " << steps << " Steps.\n"
              << "Simulated time: " << steps * dt << " seconds\n"
              << "Computation time: " << end - start << " seconds\n";
    std::cerr << "\n===============================================================\n";
}  // end simLooper


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

    if (skip == 0) {
        fprintf(stderr, "\nSKIP NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (kin < 0) {
        fprintf(stderr, "\nSPRING CONSTANT NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (dt <= 0) {
        fprintf(stderr, "\nDT NOT SET\n");
        exit(EXIT_FAILURE);
    }

    if (steps == 0) {
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