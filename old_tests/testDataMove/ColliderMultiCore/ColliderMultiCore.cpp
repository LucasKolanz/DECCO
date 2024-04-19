#include "ball_group_multi_core.hpp"
// #include "../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include <sstream>
#include <string>
#include <iomanip>
// #include <filesystem>
#include <string.h>
// namespace fs = std::filesystem;


// String buffers to hold data in memory until worth writing to file:
// std::stringstream O.ballBuffer;
// std::stringstream energyBuffer;
// std::stringstream contactBuffer;


// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
// const time_t start = time(nullptr);  // For end of program analysis
// #pragma acc enter data copyin(start)
// bool writeStep;                      // This prevents writing to file every step (which is slow).
// bool contact = false;
// bool inital_contact = true;


// Prototypes
// void
// sim_one_step(const bool write_step, Ball_group &O);
// void 
// bufferBarf(Ball_group &O);
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
int 
determine_index(std::string s, char del);
bool 
isNumber(const std::string& str);
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
main(int argc, char* argv[])
{

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
    
    // t.start_event("WholeThing");
    
    
    
    // Runtime arguments:
    // if (argc > 2) 
    // {
    //     std::stringstream s(argv[2]);
    //     // s << argv[2];
    //     s >> num_balls;
    //     // numThreads = atoi(argv[1]);
    //     // fprintf(stderr,"\nThread count set to %i.\n", numThreads);
    //     // projectileName = argv[2];
    //     // targetName = argv[3];
    //     // KEfactor = atof(argv[4]);
    // }
    // else
    // {
    //     num_balls = 100;
    // }

    // O.zeroAngVel();
    // O.pushApart();

    // Normal sim:
    // O.sim_init_write(output_prefix);
    Ball_group dummy(1);
    dummy.parse_input_file(argv[1]);
    if (dummy.typeSim == dummy.collider)
    {
        #ifdef MPI_ENABLE
            MPI_Barrier(MPI_COMM_WORLD);
        #endif
        collider(argv[1],dummy.projectileName,dummy.targetName);
    }
    else if (dummy.typeSim == dummy.BPCA)
    {
        if (dummy.total_balls_to_add >= 0)
        {
            #ifdef MPI_ENABLE
                MPI_Barrier(MPI_COMM_WORLD);
            #endif
            BPCA(argv[1],dummy.total_balls_to_add);
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
    // sim_looper();

    // collider(argv[1],projTarget,projTarget);
    // #pragma acc exit data delete(start)
    // t.end_event("WholeThing");
    // t.print_events();
    // t.save_events(output_folder + "timing.txt");
    #ifdef MPI_ENABLE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    #endif
}  // end main
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void collider(const char *path, std::string projectileName, std::string targetName)
{

    int world_rank = getRank();

    // t.start_event("collider");
    Ball_group O = Ball_group(std::string(path),std::string(projectileName),std::string(targetName));
    if (world_rank == 0)
    {
        safetyChecks(O);
        O.sim_init_write(output_prefix);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    O.sim_looper();
    // t.end_event("collider");
    // O.freeMemory();
    return;
}

void BPCA(const char *path,int num_balls)
{
    int world_rank = getRank();

    int rest = -1;
    int *restart = &rest;
    Ball_group O = make_group(path,restart);    
    std::cout<<path<<std::endl;
    safetyChecks(O);
    // Add projectile: For dust formation BPCA
    std::string ori_output_prefix = output_prefix;
    for (int i = *restart; i < num_balls; i++) {
    // for (int i = 0; i < 250; i++) {
        // O.zeroAngVel();
        // O.zeroVel();
        // contact = false;
        // inital_contact = true;
        // t.start_event("add_projectile");
        O = O.add_projectile();
        // t.end_event("add_projectile");
        if (world_rank == 0)
        {
            O.sim_init_write(ori_output_prefix, i);
        }
        O.sim_looper();
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
            /////////////
            (*restart)--;
            /////////////
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

bool isNumber(const std::string& str)
{
    for (char const &c : str) {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}

int determine_index(std::string s, char del)
{
    std::stringstream ss(s);
    std::string word;
    std::string firstWord;
    int count = 0;
    while (!ss.eof()) {
        getline(ss, word, del);
        // std::cerr<<word<<std::endl;
        if (!isNumber(word))
        {
            if (count == 1)
            {
                return 0;
            }
            else
            {
                return stoi(firstWord);
            }
        }
        else if (count == 0)
        {
            firstWord = word;
        }
        count++;
    }
    return -1;
}

// @brief checks if this is new job or restart
std::string check_restart(std::string folder,int* restart)
{
    return "0";
    // int world_rank = getRank();
    // std::string file;
    // // int tot_count = 0;
    // // int file_count = 0;
    // int largest_file_index = -1;
    // int file_index;
    // std::string largest_index_name;
    // for (const auto & entry : fs::directory_iterator(folder))
    // {
    //     file = entry.path();
    //     if (file.substr(0,file.size()-4) == "timing")
    //     {
    //         *restart = -2;
    //         return "";
    //     }

    //     size_t pos = file.find_last_of("/");
    //     file = file.erase(0,pos+1);
    //     // std::cerr<<"file1: "<<file<<std::endl;
    //     // tot_count++;
    //     if (file.substr(file.size()-4,file.size()) == ".csv")
    //     {
    //         // file_count++;
    //         // std::cerr<<"file1: "<<file<<std::endl;
            
    //         // if (file[3] == '_')
    //         // {
    //         //     std::cerr<<stoi(file.substr(0,file.find("_")))<<std::endl;
    //         //     std::cerr<<file.substr(0,file.find("_"))<<std::endl;
    //         //     file_index = stoi(file.substr(0,file.find("_")));
    //         // }
    //         // else if (file[1] == '_' and file[3] != '_')
    //         // {
    //         //     file_index = 0;
    //         // }
    //         file_index = determine_index(file,'_');
    //         if (file_index > largest_file_index)
    //         {
    //             largest_file_index = file_index;
    //             largest_index_name = file;
    //         }
    //     }
    // }
    // *restart = largest_file_index;
    // // std::cerr<<largest_file_index<<std::endl;
    // if (*restart != -1)
    // {
    //     size_t sta,end;
    //     sta = largest_index_name.find('_');
    //     end = largest_index_name.find_last_of('_');
    //     //Delete most recent save file as this is likely only partially 
    //     //complete if we are restarting

    //     std::string remove_file;

    //     if (*restart == 0)
    //     {
    //         remove_file = largest_index_name.substr(0,end+1);
    //     }
    //     else
    //     {
    //         remove_file = std::to_string(*restart) + largest_index_name.substr(sta,end-sta+1);
    //     }

    //     std::string file1 = folder + remove_file + "constants.csv";
    //     std::string file2 = folder + remove_file + "energy.csv";
    //     std::string file3 = folder + remove_file + "simData.csv";
        
    //     #ifdef MPI_ENABLE
    //         MPI_Barrier(MPI_COMM_WORLD);
    //     #endif
    //     if (world_rank == 0)
    //     {
    //         int status1 = remove(file1.c_str());
    //         int status2 = remove(file2.c_str());
    //         int status3 = remove(file3.c_str());

    //         if (status1 != 0)
    //         {
    //             std::cout<<"File: "<<file1<<" could not be removed, now exiting with failure."<<std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //         else if (status2 != 0)
    //         {
    //             std::cout<<"File: "<<file2<<" could not be removed, now exiting with failure."<<std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //         else if (status3 != 0)
    //         {
    //             std::cout<<"File: "<<file3<<" could not be removed, now exiting with failure."<<std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //     }

    //     return largest_index_name.substr(sta,end-sta+1);
    // }
    // else
    // {
    //     return "";
    // }
}




// void bufferBarf(Ball_group &O)
// {
//     // Write simData to file and clear buffer.
//     std::ofstream ballWrite;
//     ballWrite.open(output_folder + output_prefix + "simData.csv", std::ofstream::app);
//     ballWrite << O.ballBuffer.rdbuf();  // Barf buffer to file.
//     O.ballBuffer.str("");               // Empty the stream for next filling.
//     ballWrite.close();

//     // Write Energy data to file and clear buffer.
//     ////////////////////////////////////////////
//     //TURN ON FOR REAL SIM
//     std::ofstream energyWrite;
//     energyWrite.open(output_folder + output_prefix + "energy.csv", std::ofstream::app);
//     energyWrite << O.energyBuffer.rdbuf();
//     O.energyBuffer.str("");  // Empty the stream for next filling.
//     energyWrite.close();
//     // ////////////////////////////////////////////
// }


void
safetyChecks(Ball_group &O)
{
    int world_rank = getRank();

    if (world_rank == 0)
    {
        titleBar("SAFETY CHECKS");
    }

    if (O.soc <= 0) {
        fprintf(stderr, "\nSOC NOT SET, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }

    if (O.v_collapse <= 0) {
        fprintf(stderr, "\nvCollapse NOT SET, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }

    if (skip == 0) {
        fprintf(stderr, "\nSKIP NOT SET, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }

    if (kin < 0) {
        fprintf(stderr, "\nSPRING CONSTANT NOT SET, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }

    if (dt <= 0) {
        fprintf(stderr, "\nDT NOT SET, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }

    if (steps == 0) {
        fprintf(stderr, "\nSTEPS NOT SET, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }

    if (O.initial_radius <= 0) {
        fprintf(stderr, "\nCluster initialRadius not set, rank %d\n",world_rank);
        exit(EXIT_FAILURE);
    }


    for (int Ball = 0; Ball < O.num_particles; Ball++) {
        if (O.pos[Ball].norm() < vec3(1e-10, 1e-10, 1e-10).norm()) {
            fprintf(stderr, "\nA ball position is [0,0,0]. Possibly didn't initialize balls properly, rank %d\n",world_rank);
            exit(EXIT_FAILURE);
        }

        if (O.R[Ball] <= 0) {
            fprintf(stderr, "\nA balls radius <= 0, rank %d\n",world_rank);
            exit(EXIT_FAILURE);
        }

        if (O.m[Ball] <= 0) {
            fprintf(stderr, "\nA balls mass <= 0, rank %d\n",world_rank);
            exit(EXIT_FAILURE);
        }
    }
    if (world_rank == 0)
    {
        titleBar("SAFETY PASSED");
    }
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