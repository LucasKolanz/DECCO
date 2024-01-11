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
sim_looper(Ball_group &O,int start_step);
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
        sim_looper(O,O.start_step);
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

        sim_looper(O,1);

        simTimeElapsed = 0;
    }
    // O.freeMemory();
    return;
}


void
sim_one_step(const bool write_step, Ball_group &O)
{
    /// FIRST PASS - Update Kinematic Parameters:
    t.start_event("UpdateKinPar");
    for (int Ball = 0; Ball < O.num_particles; Ball++) {
        // Update velocity half step:
        O.velh[Ball] = O.vel[Ball] + .5 * O.acc[Ball] * O.dt;

        // Update angular velocity half step:
        O.wh[Ball] = O.w[Ball] + .5 * O.aacc[Ball] * O.dt;

        // Update position:
        O.pos[Ball] += O.velh[Ball] * O.dt;

        // Reinitialize acceleration to be recalculated:
        O.acc[Ball] = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        O.aacc[Ball] = {0, 0, 0};
    }
    t.end_event("UpdateKinPar");

    // std::ofstream accWrite, aaccWrite;
    // accWrite.open(output_folder+"accWrite_"+std::to_string(O.num_particles)+".txt",std::ios::app);
    // aaccWrite.open(output_folder+"aaccWrite_"+std::to_string(O.num_particles)+".txt",std::ios::app);

    /// SECOND PASS - Check for collisions, apply forces and torques:
    t.start_event("CalcForces/loopApplicablepairs");
    for (int A = 1; A < O.num_particles; A++)  
    {
        /// DONT DO ANYTHING HERE. A STARTS AT 1.
        for (int B = 0; B < A; B++) {
            const double sumRaRb = O.R[A] + O.R[B];
            const vec3 rVecab = O.pos[B] - O.pos[A];  // Vector from a to b.
            const vec3 rVecba = -rVecab;
            const double dist = (rVecab).norm();

            //////////////////////
            // const double grav_scale = 3.0e21;
            //////////////////////

            // Check for collision between Ball and otherBall:
            double overlap = sumRaRb - dist;

            vec3 totalForceOnA{0, 0, 0};

            // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
            int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
            double oldDist = O.distances[e];
            /////////////////////////////
            // double inoutT;
            /////////////////////////////
            // Check for collision between Ball and otherBall.
            if (overlap > 0) {

                // if (!contact && A == O.num_particles-1)
                // {
                //     // std::cout<<"CONTACT MADE"<<std::endl;
                //     contact = true;
                //     contactBuffer<<A<<','<<simTimeElapsed<<'\n';
                // }

                double k;
                if (dist >= oldDist) {
                    k = O.kout;
                } else {
                    k = O.kin;
                }

                // Cohesion (in contact) h must always be h_min:
                // constexpr double h = h_min;
                const double h = h_min;
                const double Ra = O.R[A];
                const double Rb = O.R[B];
                const double h2 = h * h;
                // constexpr double h2 = h * h;
                const double twoRah = 2 * Ra * h;
                const double twoRbh = 2 * Rb * h;
                const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                         ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                         rVecab.normalized();
                
                // std::cout<<vdwForceOnA[0]<<std::endl;  
                // const vec3 vdwForceOnA = {0,0,0}; // FOR TESTING
                /////////////////////////////
                if (O.write_all)
                {
                    O.vdwForce[A] += vdwForceOnA;
                    O.vdwForce[B] -= vdwForceOnA;
                }
                /////////////////////////////
                // Elastic force:
                // vec3 elasticForceOnA{0, 0, 0};
                // if (std::fabs(overlap) > 1e-6)
                // {
                //     elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
                // }
                const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
                ///////////////////////////////
                // O.elasticForce[A] += elasticForceOnA;
                // O.elasticForce[B] -= elasticForceOnA;
                ///////////////////////////////
                ///////////////////////////////
                ///////material parameters for silicate composite from Reissl 2023
                // const double Estar = 1e5*169; //in Pa
                // const double nu2 = 0.27*0.27; // nu squared (unitless)
                // const double prevoverlap = sumRaRb - oldDist;
                // const double rij = sqrt(std::pow(Ra,2)-std::pow((Ra-overlap/2),2));
                // const double Tvis = 15e-12; //Viscoelastic timescale (15ps)
                // // const double Tvis = 5e-12; //Viscoelastic timescale (5ps)
                // const vec3 viscoelaticforceOnA = -(2*Estar/nu2) * 
                //                                  ((overlap - prevoverlap)/dt) * 
                //                                  rij * Tvis * (rVecab / dist);
                const vec3 viscoelaticforceOnA = {0,0,0};
                ///////////////////////////////

                // Gravity force:
                // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] * grav_scale / (dist * dist)) * (rVecab / dist); //SCALE MASS
                const vec3 gravForceOnA = {0,0,0};
                // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

                // Sliding and Rolling Friction:
                vec3 slideForceOnA{0, 0, 0};
                vec3 rollForceA{0, 0, 0};
                vec3 torqueA{0, 0, 0};
                vec3 torqueB{0, 0, 0};

                // Shared terms:
                const double elastic_force_A_mag = elasticForceOnA.norm();
                const vec3 r_a = rVecab * O.R[A] / sumRaRb;  // Center to contact point
                const vec3 r_b = rVecba * O.R[B] / sumRaRb;
                const vec3 w_diff = O.w[A] - O.w[B];

                // Sliding friction terms:
                const vec3 d_vel = O.vel[B] - O.vel[A];
                const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
                                           O.w[A].cross(r_a) - O.w[B].cross(r_a);

                // Compute sliding friction force:
                const double rel_vel_mag = frame_A_vel_B.norm();
                // if (rel_vel_mag > 1e-20)  // Divide by zero protection.
                // if (rel_vel_mag > 1e-8)  // Divide by zero protection.
                ////////////////////////////////////////// CALC THIS AT INITIALIZATION for all combos os Ra,Rb
                // const double u_scale = O.calc_VDW_force_mag(Ra,Rb,O.h_min_physical)/
                //                         vdwForceOnA.norm();         //Friction coefficient scale factor
                //////////////////////////////////////////
                if (rel_vel_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
                {
                    // slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                    // In the frame of A, B applies force in the direction of B's velocity.
                    ///////////////////////////////////
                    // if (O.mu_scale)
                    // {
                    //     if (O.u_scale[e]*u_s > O.max_mu)
                    //     {
                    //         slideForceOnA = O.max_mu * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                    //     }
                    //     else
                    //     {
                    //         slideForceOnA = O.u_scale[e] * u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                    //     }
                    // }
                    // else
                    // {
                        slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                    // }
                    ///////////////////////////////////
                }
                //////////////////////////////////////
                // O.slideForce[A] += slideForceOnA;
                // O.slideForce[B] -= slideForceOnA;
                //////////////////////////////////////


                // Compute rolling friction force:
                const double w_diff_mag = w_diff.norm();
                // if (w_diff_mag > 1e-20)  // Divide by zero protection.
                // if (w_diff_mag > 1e-8)  // Divide by zero protection.
                if (w_diff_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
                {
                    // rollForceA = 
                    //     -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
                    //     (w_diff).cross(r_a).norm();
                    /////////////////////////////////////
                    // if (O.mu_scale)
                    // {
                    //     if (O.u_scale[e]*u_r > O.max_mu)
                    //     {
                    //         rollForceA = 
                    //             -O.max_mu * elastic_force_A_mag * (w_diff).cross(r_a) / 
                    //             (w_diff).cross(r_a).norm();
                    //     }
                    //     else
                    //     {
                    //         rollForceA = 
                    //             -O.u_scale[e] * u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
                    //             (w_diff).cross(r_a).norm();
                    //     }
                    // }
                    // else
                    // {
                        rollForceA = 
                            -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
                            (w_diff).cross(r_a).norm();
                    // }
                    /////////////////////////////////////
                }


                // Total forces on a:
                // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;
                ////////////////////////////////
                totalForceOnA = viscoelaticforceOnA + gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;
                ////////////////////////////////

                // Total torque a and b:
                torqueA = r_a.cross(slideForceOnA + rollForceA);
                torqueB = r_b.cross(-slideForceOnA + rollForceA); // original code



                O.aacc[A] += torqueA / O.moi[A];
                O.aacc[B] += torqueB / O.moi[B];

                if (write_step) {
                    // No factor of 1/2. Includes both spheres:
                    // O.PE += -G * O.m[A] * O.m[B] * grav_scale / dist + 0.5 * k * overlap * overlap;
                    // O.PE += -G * O.m[A] * O.m[B] / dist + 0.5 * k * overlap * overlap;

                    // Van Der Waals + elastic:
                    const double diffRaRb = O.R[A] - O.R[B];
                    const double z = sumRaRb + h;
                    const double two_RaRb = 2 * O.R[A] * O.R[B];
                    const double denom_sum = z * z - (sumRaRb * sumRaRb);
                    const double denom_diff = z * z - (diffRaRb * diffRaRb);
                    const double U_vdw =
                        -Ha / 6 *
                        (two_RaRb / denom_sum + two_RaRb / denom_diff + 
                        log(denom_sum / denom_diff));
                    O.PE += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
                }
            } else  // Non-contact forces:
            {

                // No collision: Include gravity and vdw:
                // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] * grav_scale / (dist * dist)) * (rVecab / dist);
                const vec3 gravForceOnA = {0.0,0.0,0.0};
                // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
                double h = std::fabs(overlap);
                if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
                {
                    h = h_min;
                }
                const double Ra = O.R[A];
                const double Rb = O.R[B];
                const double h2 = h * h;
                const double twoRah = 2 * Ra * h;
                const double twoRbh = 2 * Rb * h;
                const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                         ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                         rVecab.normalized();
                // const vec3 vdwForceOnA = {0.0,0.0,0.0};
                /////////////////////////////
                if (O.write_all)
                {
                    O.vdwForce[A] += vdwForceOnA;
                    O.vdwForce[B] -= vdwForceOnA;
                }
                /////////////////////////////
                /////////////////////////////
                totalForceOnA = vdwForceOnA + gravForceOnA;
                // totalForceOnA = vdwForceOnA;
                // totalForceOnA = gravForceOnA;
                /////////////////////////////
                if (write_step) {
                    // O.PE += -G * O.m[A] * O.m[B] * grav_scale / dist; // Gravitational

                    const double diffRaRb = O.R[A] - O.R[B];
                    const double z = sumRaRb + h;
                    const double two_RaRb = 2 * O.R[A] * O.R[B];
                    const double denom_sum = z * z - (sumRaRb * sumRaRb);
                    const double denom_diff = z * z - (diffRaRb * diffRaRb);
                    const double U_vdw =
                        -Ha / 6 *
                        (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                    O.PE += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
                }

                // todo this is part of push_apart. Not great like this.
                // For pushing apart overlappers:
                // O.vel[A] = { 0,0,0 };
                // O.vel[B] = { 0,0,0 };
            }

            // Newton's equal and opposite forces applied to acceleration of each ball:
            O.acc[A] += totalForceOnA / O.m[A];
            O.acc[B] -= totalForceOnA / O.m[B];


            // So last distance can be known for COR:
            O.distances[e] = dist;

        }
        // DONT DO ANYTHING HERE. A STARTS AT 1.
    }

    

    t.end_event("CalcForces/loopApplicablepairs");

    // if (write_step) {
    //     ballBuffer << '\n';  // Prepares a new line for incoming data.
    //     // std::cerr<<"Writing "<<O.num_particles<<" balls"<<std::endl;
    // }

    // THIRD PASS - Calculate velocity for next step:
    t.start_event("CalcVelocityforNextStep");
    for (int Ball = 0; Ball < O.num_particles; Ball++) 
    {
        // Velocity for next step:
        O.vel[Ball] = O.velh[Ball] + .5 * O.acc[Ball] * O.dt;
        O.w[Ball] = O.wh[Ball] + .5 * O.aacc[Ball] * O.dt;

        /////////////////////////////////
        // if (true) {
        /////////////////////////////////
        if (write_step) 
        {
            // Send positions and rotations to buffer:
            // std::cerr<<"Write ball "<<Ball<<std::endl;
            
            int start = O.data->getWidth("simData")*O.num_writes+Ball*O.data->getSingleWidth("simData");
            // std::cout<<"O.data->getSingleWidth(simData)"<<O.data->getSingleWidth("simData")<<std::endl;
            // std::cout<<"O.data->getWidth(simData)"<<O.data->getWidth("simData")<<std::endl;
            // std::cout<<"O.num_writes: "<<O.num_writes<<"\tball: "<<Ball<<"\tstart: "<<start<<std::endl;
            O.ballBuffer[start] = O.pos[Ball][0];
            O.ballBuffer[start+1] = O.pos[Ball][1];
            O.ballBuffer[start+2] = O.pos[Ball][2];
            O.ballBuffer[start+3] = O.w[Ball][0];
            O.ballBuffer[start+4] = O.w[Ball][1];
            O.ballBuffer[start+5] = O.w[Ball][2];
            O.ballBuffer[start+6] = O.w[Ball].norm();
            O.ballBuffer[start+7] = O.vel[Ball][0];
            O.ballBuffer[start+8] = O.vel[Ball][1];
            O.ballBuffer[start+9] = O.vel[Ball][2];
            O.ballBuffer[start+10] = 0;

            // O.ballBuffer.push_back(O.pos[Ball][0]);
            // O.ballBuffer.push_back(O.pos[Ball][1]);
            // O.ballBuffer.push_back(O.pos[Ball][2]);
            // O.ballBuffer.push_back(O.w[Ball][0]);
            // O.ballBuffer.push_back(O.w[Ball][1]);
            // O.ballBuffer.push_back(O.w[Ball][2]);
            // O.ballBuffer.push_back(O.w[Ball].norm());
            // O.ballBuffer.push_back(O.vel[Ball][0]);
            // O.ballBuffer.push_back(O.vel[Ball][1]);
            // O.ballBuffer.push_back(O.vel[Ball][2]);
            // O.ballBuffer.push_back(0);

            // if (Ball == 0) {
            //     ballBuffer << O.pos[Ball][0] << ',' << O.pos[Ball][1] << ',' << O.pos[Ball][2] << ','
            //                << O.w[Ball][0] << ',' << O.w[Ball][1] << ',' << O.w[Ball][2] << ','
            //                << O.w[Ball].norm() << ',' << O.vel[Ball].x << ',' << O.vel[Ball].y << ','
            //                << O.vel[Ball].z << ',' << 0;
            // } else {
            //     ballBuffer << ',' << O.pos[Ball][0] << ',' << O.pos[Ball][1] << ',' << O.pos[Ball][2] << ','
            //                << O.w[Ball][0] << ',' << O.w[Ball][1] << ',' << O.w[Ball][2] << ','
            //                << O.w[Ball].norm() << ',' << O.vel[Ball].x << ',' << O.vel[Ball].y << ','
            //                << O.vel[Ball].z << ',' << 0;
            // }

            
            O.KE += .5 * O.m[Ball] * O.vel[Ball].normsquared() +
                    .5 * O.moi[Ball] * O.w[Ball].normsquared();  // Now includes rotational kinetic energy.
            O.mom += O.m[Ball] * O.vel[Ball];
            O.ang_mom += O.m[Ball] * O.pos[Ball].cross(O.vel[Ball]) + O.moi[Ball] * O.w[Ball];
        }
    }  // THIRD PASS END
    if (writeStep)
    {
        O.num_writes ++;
    }
    t.end_event("CalcVelocityforNextStep");
}  // one Step end


void
sim_looper(Ball_group &O,int start_step=1)
{
    O.num_writes = 0;
    std::cerr << "Beginning simulation...\n";

    std::cerr<<"start step: "<<start_step<<std::endl;

    startProgress = time(nullptr);

    std::cerr<<"Stepping through "<<O.steps<<" steps"<<std::endl;

    for (int Step = start_step; Step < O.steps; Step++)  // Steps start at 1 for non-restart because the 0 step is initial conditions.
    {
        // simTimeElapsed += dt; //New code #1
        // Check if this is a write step:
        if (Step % O.skip == 0) {
            t.start_event("writeProgressReport");
            writeStep = true;

            /////////////////////// Original code #1
            simTimeElapsed += O.dt * O.skip;
            ///////////////////////

            // Progress reporting:
            float eta = ((time(nullptr) - startProgress) / static_cast<float>(O.skip) *
                         static_cast<float>(O.steps - Step)) /
                        3600.f;  // Hours.
            float real = (time(nullptr) - start) / 3600.f;
            float simmed = static_cast<float>(simTimeElapsed / 3600.f);
            float progress = (static_cast<float>(Step) / static_cast<float>(O.steps) * 100.f);
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
            t.end_event("writeProgressReport");
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
        sim_one_step(writeStep,O);

        if (writeStep) {
            // t.start_event("writeStep");
            // Write energy to stream:
            ////////////////////////////////////
            //TURN THIS ON FOR REAL RUNS!!!
            // O.energyBuffer = std::vector<double> (data->getWidth("energy"));
            int start = O.data->getWidth("energy")*(O.num_writes-1);
            // std::cerr<<"start,num_writes: "<<start<<','<<O.num_writes<<std::endl;
            O.energyBuffer[start] = simTimeElapsed;
            O.energyBuffer[start+1] = O.PE;
            O.energyBuffer[start+2] = O.KE;
            O.energyBuffer[start+3] = O.PE+O.KE;
            O.energyBuffer[start+4] = O.mom.norm();
            O.energyBuffer[start+5] = O.ang_mom.norm();

            
            // O.energyBuffer.push_back(simTimeElapsed);
            // O.energyBuffer.push_back(O.PE);
            // O.energyBuffer.push_back(O.KE);
            // O.energyBuffer.push_back(O.PE+O.KE);
            // O.energyBuffer.push_back(O.mom.norm());
            // O.energyBuffer.push_back(O.ang_mom.norm());


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
            if (Step / O.skip % 10 == 0) {
                // Report vMax:

                std::cerr << "vMax = " << O.getVelMax() << " Steps recorded: " << Step / O.skip << '\n';
                std::cerr << "Data Write to "<<O.output_folder<<"\n";
                // std::cerr<<"output_prefix: "<<output_prefix<<std::endl;
                
                O.data->Write(O.ballBuffer,"simData",bufferlines);
                // if (O.num_particles > 5)
                // {
                //     for (int i = 0; i < O.ballBuffer.size(); i++)
                //     {
                //         std::cerr<<O.ballBuffer[i]<<", ";
                //     }
                //     std::cerr<<std::endl;
                // }
                O.ballBuffer.clear();
                O.ballBuffer = std::vector<double>(O.data->getWidth("simData")*bufferlines);
                O.data->Write(O.energyBuffer,"energy");
                O.energyBuffer.clear();
                O.energyBuffer = std::vector<double>(O.data->getWidth("energy")*bufferlines);

                O.num_writes = 0;
                lastWrite = time(nullptr);

                // if (O.num_particles > 5)
                // {
                //     std::cerr<<"EXITING, step: "<<Step<<std::endl;
                //     exit(0);
                // }

            }  // Data export end


            if (dynamicTime) { O.calibrate_dt(Step, false); }
            // t.end_event("writeStep");
        }  // writestep end
    }

    // if (true)
    // {
    //     for (int i = 0; i < O.num_particles; i++)
    //     {
    //         std::cerr<<"===================================="<<std::endl;
    //         std::cerr<<O.pos[i]<<std::endl;
    //         std::cerr<<O.vel[i]<<std::endl;
    //         std::cerr<<"===================================="<<std::endl;
    //     }
    // }

    const time_t end = time(nullptr);

    std::cerr << "Simulation complete!\n"
              << O.num_particles << " Particles and " << O.steps << " Steps.\n"
              << "Simulated time: " << O.steps * O.dt << " seconds\n"
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