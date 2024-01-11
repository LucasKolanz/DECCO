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
// std::stringstream O.ballBuffer;
// std::stringstream energyBuffer;
// std::stringstream contactBuffer;


// These are used within simOneStep to keep track of time.
// They need to survive outside its scope, and I don't want to have to pass them all.
const time_t start = time(nullptr);  // For end of program analysis
// bool writeStep;                      // This prevents writing to file every step (which is slow).
// bool contact = false;
// bool inital_contact = true;


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
    // collider(argv[1],dummy.projectileName,dummy.targetName);

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
        // contact = false;
        // inital_contact = true;
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
        size_t sta,end;
        sta = largest_index_name.find('_');
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
            remove_file = std::to_string(*restart) + largest_index_name.substr(sta,end-sta+1);
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

        return largest_index_name.substr(sta,end-sta+1);
    }
    else
    {
        return "";
    }
}

// void
// sim_one_step(const bool write_step, Ball_group &O)
// {
//     /// FIRST PASS - Update Kinematic Parameters:
//     // t.start_event("UpdateKinPar");
//     for (int Ball = 0; Ball < O.num_particles; Ball++) {
//         // Update velocity half step:
//         O.velh[Ball] = O.vel[Ball] + .5 * O.acc[Ball] * dt;

//         // Update angular velocity half step:
//         O.wh[Ball] = O.w[Ball] + .5 * O.aacc[Ball] * dt;

//         // Update position:
//         O.pos[Ball] += O.velh[Ball] * dt;

//         // Reinitialize acceleration to be recalculated:
//         O.acc[Ball] = {0, 0, 0};

//         // Reinitialize angular acceleration to be recalculated:
//         O.aacc[Ball] = {0, 0, 0};
//     }
//     // t.end_event("UpdateKinPar");

//     // std::ofstream accWrite, aaccWrite;
//     // accWrite.open(output_folder+"accWrite_"+std::to_string(O.num_particles)+".txt",std::ios::app);
//     // aaccWrite.open(output_folder+"aaccWrite_"+std::to_string(O.num_particles)+".txt",std::ios::app);

//     /// SECOND PASS - Check for collisions, apply forces and torques:
//     // t.start_event("CalcForces/loopApplicablepairs");
//     for (int A = 1; A < O.num_particles; A++)  
//     {
//         /// DONT DO ANYTHING HERE. A STARTS AT 1.
//         for (int B = 0; B < A; B++) {
//             const double sumRaRb = O.R[A] + O.R[B];
//             const vec3 rVecab = O.pos[B] - O.pos[A];  // Vector from a to b.
//             const vec3 rVecba = -rVecab;
//             const double dist = (rVecab).norm();

//             //////////////////////
//             // const double grav_scale = 3.0e21;
//             //////////////////////

//             // Check for collision between Ball and otherBall:
//             double overlap = sumRaRb - dist;

//             vec3 totalForceOnA{0, 0, 0};

//             // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
//             int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
//             double oldDist = O.distances[e];
//             /////////////////////////////
//             // double inoutT;
//             /////////////////////////////
//             // Check for collision between Ball and otherBall.
//             if (overlap > 0) {

//                 // if (!contact && A == O.num_particles-1)
//                 // {
//                 //     // std::cout<<"CONTACT MADE"<<std::endl;
//                 //     contact = true;
//                 //     contactBuffer<<A<<','<<simTimeElapsed<<'\n';
//                 // }

//                 double k;
//                 if (dist >= oldDist) {
//                     k = kout;
//                 } else {
//                     k = kin;
//                 }

//                 // Cohesion (in contact) h must always be h_min:
//                 // constexpr double h = h_min;
//                 const double h = h_min;
//                 const double Ra = O.R[A];
//                 const double Rb = O.R[B];
//                 const double h2 = h * h;
//                 // constexpr double h2 = h * h;
//                 const double twoRah = 2 * Ra * h;
//                 const double twoRbh = 2 * Rb * h;
//                 const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
//                                          ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
//                                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
//                                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
//                                          rVecab.normalized();
                
//                 // std::cout<<vdwForceOnA[0]<<std::endl;  
//                 // const vec3 vdwForceOnA = {0,0,0}; // FOR TESTING
//                 /////////////////////////////
//                 if (O.write_all)
//                 {
//                     O.vdwForce[A] += vdwForceOnA;
//                     O.vdwForce[B] -= vdwForceOnA;
//                 }
//                 /////////////////////////////
//                 // Elastic force:
//                 // vec3 elasticForceOnA{0, 0, 0};
//                 // if (std::fabs(overlap) > 1e-6)
//                 // {
//                 //     elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
//                 // }
//                 const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
//                 ///////////////////////////////
//                 // O.elasticForce[A] += elasticForceOnA;
//                 // O.elasticForce[B] -= elasticForceOnA;
//                 ///////////////////////////////
//                 ///////////////////////////////
//                 ///////material parameters for silicate composite from Reissl 2023
//                 // const double Estar = 1e5*169; //in Pa
//                 // const double nu2 = 0.27*0.27; // nu squared (unitless)
//                 // const double prevoverlap = sumRaRb - oldDist;
//                 // const double rij = sqrt(std::pow(Ra,2)-std::pow((Ra-overlap/2),2));
//                 // const double Tvis = 15e-12; //Viscoelastic timescale (15ps)
//                 // // const double Tvis = 5e-12; //Viscoelastic timescale (5ps)
//                 // const vec3 viscoelaticforceOnA = -(2*Estar/nu2) * 
//                 //                                  ((overlap - prevoverlap)/dt) * 
//                 //                                  rij * Tvis * (rVecab / dist);
//                 const vec3 viscoelaticforceOnA = {0,0,0};
//                 ///////////////////////////////

//                 // Gravity force:
//                 // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] * grav_scale / (dist * dist)) * (rVecab / dist); //SCALE MASS
//                 const vec3 gravForceOnA = {0,0,0};
//                 // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

//                 // Sliding and Rolling Friction:
//                 vec3 slideForceOnA{0, 0, 0};
//                 vec3 rollForceA{0, 0, 0};
//                 vec3 torqueA{0, 0, 0};
//                 vec3 torqueB{0, 0, 0};

//                 // Shared terms:
//                 const double elastic_force_A_mag = elasticForceOnA.norm();
//                 const vec3 r_a = rVecab * O.R[A] / sumRaRb;  // Center to contact point
//                 const vec3 r_b = rVecba * O.R[B] / sumRaRb;
//                 const vec3 w_diff = O.w[A] - O.w[B];

//                 // Sliding friction terms:
//                 const vec3 d_vel = O.vel[B] - O.vel[A];
//                 const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
//                                            O.w[A].cross(r_a) - O.w[B].cross(r_a);

//                 // Compute sliding friction force:
//                 const double rel_vel_mag = frame_A_vel_B.norm();
//                 // if (rel_vel_mag > 1e-20)  // Divide by zero protection.
//                 // if (rel_vel_mag > 1e-8)  // Divide by zero protection.
//                 ////////////////////////////////////////// CALC THIS AT INITIALIZATION for all combos os Ra,Rb
//                 // const double u_scale = O.calc_VDW_force_mag(Ra,Rb,O.h_min_physical)/
//                 //                         vdwForceOnA.norm();         //Friction coefficient scale factor
//                 //////////////////////////////////////////
//                 if (rel_vel_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
//                 {
//                     // slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
//                     // In the frame of A, B applies force in the direction of B's velocity.
//                     ///////////////////////////////////
//                     // if (O.mu_scale)
//                     // {
//                     //     if (O.u_scale[e]*u_s > O.max_mu)
//                     //     {
//                     //         slideForceOnA = O.max_mu * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
//                     //     }
//                     //     else
//                     //     {
//                     //         slideForceOnA = O.u_scale[e] * u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
//                     //     }
//                     // }
//                     // else
//                     // {
//                         slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
//                     // }
//                     ///////////////////////////////////
//                 }
//                 //////////////////////////////////////
//                 // O.slideForce[A] += slideForceOnA;
//                 // O.slideForce[B] -= slideForceOnA;
//                 //////////////////////////////////////


//                 // Compute rolling friction force:
//                 const double w_diff_mag = w_diff.norm();
//                 // if (w_diff_mag > 1e-20)  // Divide by zero protection.
//                 // if (w_diff_mag > 1e-8)  // Divide by zero protection.
//                 if (w_diff_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
//                 {
//                     // rollForceA = 
//                     //     -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
//                     //     (w_diff).cross(r_a).norm();
//                     /////////////////////////////////////
//                     // if (O.mu_scale)
//                     // {
//                     //     if (O.u_scale[e]*u_r > O.max_mu)
//                     //     {
//                     //         rollForceA = 
//                     //             -O.max_mu * elastic_force_A_mag * (w_diff).cross(r_a) / 
//                     //             (w_diff).cross(r_a).norm();
//                     //     }
//                     //     else
//                     //     {
//                     //         rollForceA = 
//                     //             -O.u_scale[e] * u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
//                     //             (w_diff).cross(r_a).norm();
//                     //     }
//                     // }
//                     // else
//                     // {
//                         rollForceA = 
//                             -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
//                             (w_diff).cross(r_a).norm();
//                     // }
//                     /////////////////////////////////////
//                 }
//                 //////////////////////////////////////
//                 // O.rollForce[A] += rollForceA;
//                 // O.rollForce[B] -= rollForceA;
//                 //////////////////////////////////////

//                 /////////////////////////////////
//                 // if (A == O.num_particles-1)
//                 // {
//                 //     // O.slidDir[B] += (frame_A_vel_B / rel_vel_mag);
//                 //     // O.rollDir[B] += (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
//                 //     O.slidB3[B] += slideForceOnA;
//                 //     O.rollB3[B] += rollForceA;
//                 //     // O.distB3[B] += dist;
//                 //     // O.inout[B] += inoutT;
//                 // }
//                 // O.slidFric[A] += slideForceOnA;
//                 // O.slidFric[B] -= slideForceOnA;
//                 /////////////////////////////////
//                 /////////////////////////////////
//                 // O.rollFric[A] += rollForceA;
//                 // O.rollFric[B] -= rollForceA;
//                 /////////////////////////////////

//                 // Total forces on a:
//                 // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;
//                 ////////////////////////////////
//                 totalForceOnA = viscoelaticforceOnA + gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;
//                 ////////////////////////////////

//                 // Total torque a and b:
//                 torqueA = r_a.cross(slideForceOnA + rollForceA);
//                 torqueB = r_b.cross(-slideForceOnA + rollForceA); // original code

//                 // aaccWrite<<"["<<A<<';'<<B<<";("<<torqueA<<");("<<torqueB<<")],";
//                 //////////////////////////////////////
//                 // torqueB = r_b.cross(slideForceOnA + rollForceA); // test code
//                 //////////////////////////////////////
//                 //////////////////////////////////////
//                 // O.torqueForce[twoDtoOneD(A,B,O.num_particles)] = torqueA;
//                 // O.torqueForce[twoDtoOneD(B,A,O.num_particles)] = torqueB;
//                 //////////////////////////////////////

//                 O.aacc[A] += torqueA / O.moi[A];
//                 O.aacc[B] += torqueB / O.moi[B];

//                 if (write_step) {
//                     // No factor of 1/2. Includes both spheres:
//                     // O.PE += -G * O.m[A] * O.m[B] * grav_scale / dist + 0.5 * k * overlap * overlap;
//                     // O.PE += -G * O.m[A] * O.m[B] / dist + 0.5 * k * overlap * overlap;

//                     // Van Der Waals + elastic:
//                     const double diffRaRb = O.R[A] - O.R[B];
//                     const double z = sumRaRb + h;
//                     const double two_RaRb = 2 * O.R[A] * O.R[B];
//                     const double denom_sum = z * z - (sumRaRb * sumRaRb);
//                     const double denom_diff = z * z - (diffRaRb * diffRaRb);
//                     const double U_vdw =
//                         -Ha / 6 *
//                         (two_RaRb / denom_sum + two_RaRb / denom_diff + 
//                         log(denom_sum / denom_diff));
//                     O.PE += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
//                 }
//             } else  // Non-contact forces:
//             {

//                 // No collision: Include gravity and vdw:
//                 // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] * grav_scale / (dist * dist)) * (rVecab / dist);
//                 const vec3 gravForceOnA = {0.0,0.0,0.0};
//                 // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
//                 double h = std::fabs(overlap);
//                 if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
//                 {
//                     h = h_min;
//                 }
//                 const double Ra = O.R[A];
//                 const double Rb = O.R[B];
//                 const double h2 = h * h;
//                 const double twoRah = 2 * Ra * h;
//                 const double twoRbh = 2 * Rb * h;
//                 const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
//                                          ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
//                                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
//                                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
//                                          rVecab.normalized();
//                 // const vec3 vdwForceOnA = {0.0,0.0,0.0};
//                 /////////////////////////////
//                 if (O.write_all)
//                 {
//                     O.vdwForce[A] += vdwForceOnA;
//                     O.vdwForce[B] -= vdwForceOnA;
//                 }
//                 /////////////////////////////
//                 /////////////////////////////
//                 totalForceOnA = vdwForceOnA + gravForceOnA;
//                 // totalForceOnA = vdwForceOnA;
//                 // totalForceOnA = gravForceOnA;
//                 /////////////////////////////
//                 if (write_step) {
//                     // O.PE += -G * O.m[A] * O.m[B] * grav_scale / dist; // Gravitational

//                     const double diffRaRb = O.R[A] - O.R[B];
//                     const double z = sumRaRb + h;
//                     const double two_RaRb = 2 * O.R[A] * O.R[B];
//                     const double denom_sum = z * z - (sumRaRb * sumRaRb);
//                     const double denom_diff = z * z - (diffRaRb * diffRaRb);
//                     const double U_vdw =
//                         -Ha / 6 *
//                         (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
//                     O.PE += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
//                 }

//                 // todo this is part of push_apart. Not great like this.
//                 // For pushing apart overlappers:
//                 // O.vel[A] = { 0,0,0 };
//                 // O.vel[B] = { 0,0,0 };
//             }

//             // Newton's equal and opposite forces applied to acceleration of each ball:
//             O.acc[A] += totalForceOnA / O.m[A];
//             O.acc[B] -= totalForceOnA / O.m[B];

//             // accWrite<<"["<<A<<';'<<B<<";("<<totalForceOnA<<")],";

//             // So last distance can be known for COR:
//             O.distances[e] = dist;
//             //////////////////////////
//             // if (A == O.num_particles-1)
//             // {
//             //     O.distB3[B] += dist;
//             // }
//             //////////////////////////
//         }
//         // DONT DO ANYTHING HERE. A STARTS AT 1.
//     }

//     // aaccWrite<<'\n';
//     // accWrite<<'\n';

//     //////////////////////////////
//     // Write out sliding/rolling fric
//     // std::ofstream fricWrite;
//     // fricWrite.open(output_folder + output_prefix + "fricData.csv", std::ofstream::app);
//     // fricWrite << '\n';
//     // for (int i = 0; i < O.num_particles; ++i)
//     // {
//     //     fricWrite << O.rollFric[i] << "," << O.slidFric[i];
//     //     if (i != O.num_particles - 1)
//     //     {
//     //         fricWrite << ',';
//     //     } 
//     // }
//     // fricWrite.close();

//     // std::ofstream fricB3Write;
//     // fricB3Write.open(output_folder + output_prefix + "fricB3Data.csv", std::ofstream::app);
//     // fricB3Write << '\n';
//     // for (int i = 0; i < O.num_particles; ++i)
//     // {
//     //     fricB3Write << O.rollB3[i] << "," << O.slidB3[i];
//     //     if (i != O.num_particles - 1)
//     //     {
//     //         fricB3Write << ',';
//     //     } 
//     // }
//     // fricB3Write.close();

//     // std::ofstream dirB3Write;
//     // dirB3Write.open(output_folder + output_prefix + "dirB3Data.csv", std::ofstream::app);
//     // dirB3Write << '\n';
//     // for (int i = 0; i < O.num_particles; ++i)
//     // {
//     //     dirB3Write << O.rollDir[i] << "," << O.slidDir[i];
//     //     if (i != O.num_particles - 1)
//     //     {
//     //         dirB3Write << ',';
//     //     } 
//     // }
//     // dirB3Write.close();

//     // std::ofstream distB3Write;
//     // distB3Write.open(output_folder + output_prefix + "distB3Data.csv", std::ofstream::app);
//     // distB3Write << '\n';
//     // for (int i = 0; i < O.num_particles-1; ++i)
//     // {
//     //     distB3Write << O.distB3[i];
//     //     if (i != O.num_particles - 2)
//     //     {
//     //         distB3Write << ',';
//     //     } 
//     // }
//     // distB3Write.close();

//     // std::ofstream inoutB3Write;
//     // inoutB3Write.open(output_folder + output_prefix + "inoutB3Data.csv", std::ofstream::app);
//     // inoutB3Write << '\n';
//     // for (int i = 0; i < O.num_particles-1; ++i)
//     // {
//     //     inoutB3Write << O.inout[i];
//     //     if (i != O.num_particles - 2)
//     //     {
//     //         inoutB3Write << ',';
//     //     } 
//     // }
//     // inoutB3Write.close();
    
//     //////////////////////////////

//     // t.end_event("CalcForces/loopApplicablepairs");

//     if (write_step) {
//         O.ballBuffer << '\n';  // Prepares a new line for incoming data.
//         // std::cerr<<"Writing "<<O.num_particles<<" balls"<<std::endl;
//     }

//     // THIRD PASS - Calculate velocity for next step:
//     // t.start_event("CalcVelocityforNextStep");
//     for (int Ball = 0; Ball < O.num_particles; Ball++) {
//         // Velocity for next step:
//         O.vel[Ball] = O.velh[Ball] + .5 * O.acc[Ball] * dt;
//         O.w[Ball] = O.wh[Ball] + .5 * O.aacc[Ball] * dt;

//         /////////////////////////////////
//         // if (true) {
//         /////////////////////////////////
//         if (write_step) {
//             // Send positions and rotations to buffer:
//             // std::cerr<<"Write ball "<<Ball<<std::endl;
//             if (Ball == 0) {
//                 O.ballBuffer << O.pos[Ball][0] << ',' << O.pos[Ball][1] << ',' << O.pos[Ball][2] << ','
//                            << O.w[Ball][0] << ',' << O.w[Ball][1] << ',' << O.w[Ball][2] << ','
//                            << O.w[Ball].norm() << ',' << O.vel[Ball].x << ',' << O.vel[Ball].y << ','
//                            << O.vel[Ball].z << ',' << 0;
//             } else {
//                 O.ballBuffer << ',' << O.pos[Ball][0] << ',' << O.pos[Ball][1] << ',' << O.pos[Ball][2] << ','
//                            << O.w[Ball][0] << ',' << O.w[Ball][1] << ',' << O.w[Ball][2] << ','
//                            << O.w[Ball].norm() << ',' << O.vel[Ball].x << ',' << O.vel[Ball].y << ','
//                            << O.vel[Ball].z << ',' << 0;
//             }

//             O.KE += .5 * O.m[Ball] * O.vel[Ball].normsquared() +
//                     .5 * O.moi[Ball] * O.w[Ball].normsquared();  // Now includes rotational kinetic energy.
//             O.mom += O.m[Ball] * O.vel[Ball];
//             O.ang_mom += O.m[Ball] * O.pos[Ball].cross(O.vel[Ball]) + O.moi[Ball] * O.w[Ball];
//         }
//     }  // THIRD PASS END
//     // t.end_event("CalcVelocityforNextStep");
// }  // one Step end


void
sim_looper(Ball_group &O)
{
    std::cerr << "Beginning simulation...\n";

    // startProgress = ;

    std::cerr<<"Stepping through "<<steps<<" steps"<<std::endl;

    bool writeStep;
    time_t startProgress = time(nullptr);                // For progress reporting (gets reset)
    time_t lastWrite;                    // For write control (gets reset)

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
        // if (O.write_all)
        // {
        //     O.zeroSaveVals();
        // }
        ///////////
        // sim_one_step(writeStep,O);
        O.sim_one_step(writeStep);

        if (writeStep) {
            // t.start_event("writeStep");
            // Write energy to stream:
            ////////////////////////////////////
            //TURN THIS ON FOR REAL RUNS!!!
            O.energyBuffer << '\n'
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
                ballWrite << O.ballBuffer.rdbuf();  // Barf buffer to file.
                O.ballBuffer.str("");               // Empty the stream for next filling.
                ballWrite.close();

                // Write Energy data to file and clear buffer.
                ////////////////////////////////////////////
                //TURN ON FOR REAL SIM
                std::ofstream energyWrite;
                energyWrite.open(output_folder + output_prefix + "energy.csv", std::ofstream::app);
                energyWrite << O.energyBuffer.rdbuf();
                O.energyBuffer.str("");  // Empty the stream for next filling.
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