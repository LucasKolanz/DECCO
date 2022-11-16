#include <iostream>
#include <cmath>
#include <omp.h>
#include <fstream>
#include "../ball_group.hpp"
#include "../vec3.hpp"
#include "../timing/timing.hpp"
#include "../grid_group/grid_group.hpp"
namespace fs = std::filesystem;

std::string 
check_restart(std::string folder,int* restart);

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
    else // Make new ball group
    {
        *restart = 0;
        O = Ball_group(true, v_custom, argv1); // Generate new group
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




void parseSimData(std::string line,vec3 *pos,vec3 *vel,vec3 *w,int num_particles)
{
    std::string lineElement;
    // Get number of balls in file
    int count = 269;
    int properties = 11;
    
    // int count = std::count(line.begin(), line.end(), ',') / properties + 1;
    // allocate_group(count);
    std::stringstream chosenLine(line);  // This is the last line of the read file, containing all data
                                         // for all balls at last time step
    // Get position and angular velocity data:
    for (int A = 0; A < num_particles; A++) {
        for (int i = 0; i < 3; i++)  // Position
        {
            std::getline(chosenLine, lineElement, ',');
            pos[A][i] = std::stod(lineElement);
        }
        for (int i = 0; i < 3; i++)  // Angular Velocity
        {
            std::getline(chosenLine, lineElement, ',');
            w[A][i] = std::stod(lineElement);
        }
        std::getline(chosenLine, lineElement, ',');  // Angular velocity magnitude skipped
        for (int i = 0; i < 3; i++)                  // velocity
        {
            std::getline(chosenLine, lineElement, ',');
            vel[A][i] = std::stod(lineElement);
        }
        for (int i = 0; i < properties - 10; i++)  // We used 10 elements. This skips the rest.
        {
            std::getline(chosenLine, lineElement, ',');
        }
    }
}

/// Get last line of previous simData by filename.
[[nodiscard]] static std::string getLastLine(const std::string& filename)
{
    std::string simDataFilepath = filename;
    if (auto simDataStream = std::ifstream(simDataFilepath, std::ifstream::in)) {
        std::cerr << "\nParsing last line of data.\n";

        simDataStream.seekg(-1, std::ios_base::end);  // go to one spot before the EOF

        bool keepLooping = true;
        while (keepLooping) {
            char ch = ' ';
            simDataStream.get(ch);  // Get current byte's data

            if (static_cast<int>(simDataStream.tellg()) <=
                1) {                     // If the data was at or before the 0th byte
                simDataStream.seekg(0);  // The first line is the last line
                keepLooping = false;     // So stop there
            } else if (ch == '\n') {     // If the data was a newline
                keepLooping = false;     // Stop at the current position.
            } else {                     // If the data was neither a newline nor at the 0 byte
                simDataStream.seekg(-2, std::ios_base::cur);  // Move to the front of that data, then to
                                                              // the front of the data before it
            }
        }
        std::string line;
        std::getline(simDataStream, line);  // Read the current line
        return line;
    } else {
        std::cerr << "Could not open simData file: " << simDataFilepath << "... Existing program."
                  << '\n';
        exit(EXIT_FAILURE);
    }
    return "-1";
}

// inline void ball_interactions(int A, int B, vec3 *acc,vec3 *aacc,Ball_group& O, bool write_step)
// {
//     const double sumRaRb = O.R[A] + O.R[B];
//     const vec3 rVecab = O.pos[B] - O.pos[A];  // Vector from a to b.
//     const vec3 rVecba = -rVecab;
//     const double dist = (rVecab).norm();

//     // Check for collision between Ball and otherBall:
//     double overlap = sumRaRb - dist;

//     vec3 totalForceOnA{0, 0, 0};

//     // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
//     int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
//     double oldDist = O.distances[e];

//     // Check for collision between Ball and otherBall.
//     if (overlap > 0) {
//         double k;
//         // Apply coefficient of restitution to balls leaving collision.
//         if (dist >= oldDist) {
//             k = kout;
//         } else {
//             k = kin;
//         }

//         // Cohesion (in contact) h must always be h_min:
//         // constexpr double h = h_min;
//         const double h = h_min;
//         const double Ra = O.R[A];
//         const double Rb = O.R[B];
//         const double h2 = h * h;
//         // constexpr double h2 = h * h;
//         const double twoRah = 2 * Ra * h;
//         const double twoRbh = 2 * Rb * h;
//         const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
//                                  ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
//                                                    (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
//                                                    (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
//                                  rVecab.normalized();

//         // Elastic force:
//         const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

//         // Gravity force:
//         // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

//         // Sliding and Rolling Friction:
//         vec3 slideForceOnA{0, 0, 0};
//         vec3 rollForceA{0, 0, 0};
//         vec3 torqueA{0, 0, 0};
//         vec3 torqueB{0, 0, 0};

//         // Shared terms:
//         const double elastic_force_A_mag = elasticForceOnA.norm();
//         const vec3 r_a = rVecab * O.R[A] / sumRaRb;  // Center to contact point
//         const vec3 r_b = rVecba * O.R[B] / sumRaRb;
//         const vec3 w_diff = O.w[A] - O.w[B];

//         // Sliding friction terms:
//         const vec3 d_vel = O.vel[B] - O.vel[A];
//         const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
//                                    O.w[A].cross(r_a) - O.w[B].cross(r_a);

//         // Compute sliding friction force:
//         const double rel_vel_mag = frame_A_vel_B.norm();
//         if (rel_vel_mag > 1e-13)  // Divide by zero protection.
//         {
//             // In the frame of A, B applies force in the direction of B's velocity.
//             slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
//         }

//         // Compute rolling friction force:
//         const double w_diff_mag = w_diff.norm();
//         if (w_diff_mag > 1e-13)  // Divide by zero protection.
//         {
//             rollForceA =
//                 -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
//         }

//         // Total forces on a:
//         //Took out gravity force
//         totalForceOnA = elasticForceOnA + slideForceOnA + vdwForceOnA;
//         // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;

//         // Total torque a and b:
//         torqueA = r_a.cross(slideForceOnA + rollForceA);
//         torqueB = r_b.cross(-slideForceOnA + rollForceA);

//         aacc[A] += torqueA / O.moi[A];
//         aacc[B] += torqueB / O.moi[B];

//         if (write_step) {
//             // No factor of 1/2. Includes both spheres:
//             // O.PE += -G * O.m[A] * O.m[B] / dist + 0.5 * k * overlap * overlap;

//             // Van Der Waals + elastic:
//             const double diffRaRb = O.R[A] - O.R[B];
//             const double z = sumRaRb + h;
//             const double two_RaRb = 2 * O.R[A] * O.R[B];
//             const double denom_sum = z * z - (sumRaRb * sumRaRb);
//             const double denom_diff = z * z - (diffRaRb * diffRaRb);
//             const double U_vdw =
//                 -Ha / 6 *
//                 (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
//             O.PE += U_vdw + 0.5 * k * overlap * overlap;
//         }
//     } else  // Non-contact forces:
//     {
//         // No collision: Include gravity and vdw:
//         // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

//         // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
//         double h = std::fabs(overlap);
//         if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
//         {
//             h = h_min;
//         }
//         const double Ra = O.R[A];
//         const double Rb = O.R[B];
//         const double h2 = h * h;
//         const double twoRah = 2 * Ra * h;
//         const double twoRbh = 2 * Rb * h;
//         const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
//                                  ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
//                                                    (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
//                                                    (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
//                                  rVecab.normalized();

//         totalForceOnA = vdwForceOnA;  // +gravForceOnA;
//         if (write_step) {
//             // O.PE += -G * O.m[A] * O.m[B] / dist; // Gravitational

//             const double diffRaRb = O.R[A] - O.R[B];
//             const double z = sumRaRb + h;
//             const double two_RaRb = 2 * O.R[A] * O.R[B];
//             const double denom_sum = z * z - (sumRaRb * sumRaRb);
//             const double denom_diff = z * z - (diffRaRb * diffRaRb);
//             const double U_vdw =
//                 -Ha / 6 *
//                 (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
//             O.PE += U_vdw;  // Van Der Waals
//         }

//         // todo this is part of push_apart. Not great like this.
//         // For pushing apart overlappers:
//         // O.vel[A] = { 0,0,0 };
//         // O.vel[B] = { 0,0,0 };
//     }

//     // Newton's equal and opposite forces applied to acceleration of each ball:
//     acc[A] += totalForceOnA / O.m[A];
//     acc[B] -= totalForceOnA / O.m[B];

//     // So last distance can be known for COR:
//     O.distances[e] = dist;
// }

inline void ball_interactions(int A, int B,Ball_group& bg, bool write_step)
{
    const double sumRaRb = bg.R[A] + bg.R[B];
    const vec3 rVecab = bg.pos[B] - bg.pos[A];  // Vector from a to b.
    const vec3 rVecba = -rVecab;
    const double dist = (rVecab).norm();

    // Check for collision between Ball and otherBall:
    double overlap = sumRaRb - dist;

    vec3 totalForceOnA{0, 0, 0};

    // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
    int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
    double oldDist = bg.distances[e];

    // Check for collision between Ball and otherBall.
    if (overlap > 0) {
        double k;
        // Apply coefficient of restitution to balls leaving collision.
        if (dist >= oldDist) {
            k = kout;
        } else {
            k = kin;
        }

        // Cohesion (in contact) h must always be h_min:
        // constexpr double h = h_min;
        const double h = h_min;
        const double Ra = bg.R[A];
        const double Rb = bg.R[B];
        const double h2 = h * h;
        // constexpr double h2 = h * h;
        const double twoRah = 2 * Ra * h;
        const double twoRbh = 2 * Rb * h;
        const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                 ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                 rVecab.normalized();

        // Elastic force:
        const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

        // Gravity force:
        // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

        // Sliding and Rolling Friction:
        vec3 slideForceOnA{0, 0, 0};
        vec3 rollForceA{0, 0, 0};
        vec3 torqueA{0, 0, 0};
        vec3 torqueB{0, 0, 0};

        // Shared terms:
        const double elastic_force_A_mag = elasticForceOnA.norm();
        const vec3 r_a = rVecab * bg.R[A] / sumRaRb;  // Center to contact point
        const vec3 r_b = rVecba * bg.R[B] / sumRaRb;
        const vec3 w_diff = bg.w[A] - bg.w[B];

        // Sliding friction terms:
        const vec3 d_vel = bg.vel[B] - bg.vel[A];
        const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
                                   bg.w[A].cross(r_a) - bg.w[B].cross(r_a);

        // Compute sliding friction force:
        const double rel_vel_mag = frame_A_vel_B.norm();
        if (rel_vel_mag > 1e-13)  // Divide by zero protection.
        {
            // In the frame of A, B applies force in the direction of B's velocity.
            slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
        }

        // Compute rolling friction force:
        const double w_diff_mag = w_diff.norm();
        if (w_diff_mag > 1e-13)  // Divide by zero protection.
        {
            rollForceA =
                -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
        }

        // Total forces on a:
        //Took out gravity force
        totalForceOnA = elasticForceOnA + slideForceOnA + vdwForceOnA;
        // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;

        // Total torque a and b:
        torqueA = r_a.cross(slideForceOnA + rollForceA);
        torqueB = r_b.cross(-slideForceOnA + rollForceA);

        bg.aacc[A] += torqueA / bg.moi[A];
        bg.aacc[B] += torqueB / bg.moi[B];

        if (write_step) {
            // No factor of 1/2. Includes both spheres:
            // bg.PE += -G * bg.m[A] * bg.m[B] / dist + 0.5 * k * overlap * overlap;

            // Van Der Waals + elastic:
            const double diffRaRb = bg.R[A] - bg.R[B];
            const double z = sumRaRb + h;
            const double two_RaRb = 2 * bg.R[A] * bg.R[B];
            const double denom_sum = z * z - (sumRaRb * sumRaRb);
            const double denom_diff = z * z - (diffRaRb * diffRaRb);
            const double U_vdw =
                -Ha / 6 *
                (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
            bg.PE += U_vdw + 0.5 * k * overlap * overlap;
        }
    } else  // Non-contact forces:
    {
        // No collision: Include gravity and vdw:
        // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

        // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
        double h = std::fabs(overlap);
        if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
        {
            h = h_min;
        }
        const double Ra = bg.R[A];
        const double Rb = bg.R[B];
        const double h2 = h * h;
        const double twoRah = 2 * Ra * h;
        const double twoRbh = 2 * Rb * h;
        const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                 ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                 rVecab.normalized();

        totalForceOnA = vdwForceOnA;  // +gravForceOnA;
        if (write_step) {
            // O.PE += -G * O.m[A] * O.m[B] / dist; // Gravitational

            const double diffRaRb = bg.R[A] - bg.R[B];
            const double z = sumRaRb + h;
            const double two_RaRb = 2 * bg.R[A] * bg.R[B];
            const double denom_sum = z * z - (sumRaRb * sumRaRb);
            const double denom_diff = z * z - (diffRaRb * diffRaRb);
            const double U_vdw =
                -Ha / 6 *
                (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
            bg.PE += U_vdw;  // Van Der Waals
        }

        // todo this is part of push_apart. Not great like this.
        // For pushing apart overlappers:
        // bg.vel[A] = { 0,0,0 };
        // bg.vel[B] = { 0,0,0 };
    }

    // Newton's equal and opposite forces applied to acceleration of each ball:
    bg.acc[A] += totalForceOnA / bg.m[A];
    bg.acc[B] -= totalForceOnA / bg.m[B];

    // So last distance can be known for COR:
    bg.distances[e] = dist;
}


void serial_accel(const bool write_step, Ball_group& P)
{
	double dt = 1e-5;
    // #pragma omp parallel for default(none) private(O) shared(dt) schedule(static,256)
    for (int Ball = 0; Ball < P.num_particles; Ball++) {
        // Update velocity half step:
        P.velh[Ball] = P.vel[Ball] + .5 * P.acc[Ball] * dt;

        // Update angular velocity half step:
        P.wh[Ball] = P.w[Ball] + .5 * P.aacc[Ball] * dt;

        // Update position:
        P.pos[Ball] += P.velh[Ball] * dt;

        // Reinitialize acceleration to be recalculated:
        P.acc[Ball] = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        P.aacc[Ball] = {0, 0, 0};
    }

    grid g(P.num_particles, P.R[0], P.gridSize, P.pos, P.tolerance);

    // std::cout<<"num_particles: "<<P.num_particles<<'\n';
    // std::cout<<"P.pos: "<<P.pos[10]<<'\n';
    // std::cout<<"gridSize: "<<P.gridSize<<'\n';
    // std::cout<<"tolerance: "<<P.tolerance<<'\n';
    // std::cout<<"radius: "<<P.R[0]<<'\n';
    // std::cout<<"m_total: "<<P.m_total<<'\n';

    // #pragma omp for 
    for (int A = 0; A < P.num_particles; ++A)  // cuda
    {

        // t.start_event("getNearestNeighbors");
        std::vector<int> nearest_neighbors = g.getBalls(A);
        // t.end_event("getNearestNeighbors");
        // g.printVector(nearest_neighbors);
        // g.printVector(g.getBalls(A));
        // t.start_event("loopApplicablePairs");
        // #pragma omp parallel default(none) private(A) shared(nearest_neighbors,P,write_step)
        // {
        // #pragma omp for 
        for (int B : nearest_neighbors)
        {
            if (A != B)
            {
            	// ball_interactions(A,B,P,write_step);
            	const double sumRaRb = P.R[A] + P.R[B];
			    const vec3 rVecab = P.pos[B] - P.pos[A];  // Vector from a to b.
			    const vec3 rVecba = -rVecab;
			    const double dist = (rVecab).norm();

			    // Check for collision between Ball and otherBall:
			    double overlap = sumRaRb - dist;

			    vec3 totalForceOnA{0, 0, 0};

			    // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
			    int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
			    double oldDist = P.distances[e];

			    // Check for collision between Ball and otherBall.
			    if (overlap > 0) {
			        double k;
			        // Apply coefficient of restitution to balls leaving collision.
			        if (dist >= oldDist) {
			            k = kout;
			        } else {
			            k = kin;
			        }

			        // Cohesion (in contact) h must always be h_min:
			        // constexpr double h = h_min;
			        const double h = h_min;
			        const double Ra = P.R[A];
			        const double Rb = P.R[B];
			        const double h2 = h * h;
			        // constexpr double h2 = h * h;
			        const double twoRah = 2 * Ra * h;
			        const double twoRbh = 2 * Rb * h;
			        const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
			                                 ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
			                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
			                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
			                                 rVecab.normalized();

			        // Elastic force:
			        const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

			        // Gravity force:
			        // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

			        // Sliding and Rolling Friction:
			        vec3 slideForceOnA{0, 0, 0};
			        vec3 rollForceA{0, 0, 0};
			        vec3 torqueA{0, 0, 0};
			        vec3 torqueB{0, 0, 0};

			        // Shared terms:
			        const double elastic_force_A_mag = elasticForceOnA.norm();
			        const vec3 r_a = rVecab * P.R[A] / sumRaRb;  // Center to contact point
			        const vec3 r_b = rVecba * P.R[B] / sumRaRb;
			        const vec3 w_diff = P.w[A] - P.w[B];

			        // Sliding friction terms:
			        const vec3 d_vel = P.vel[B] - P.vel[A];
			        const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
			                                   P.w[A].cross(r_a) - P.w[B].cross(r_a);

			        // Compute sliding friction force:
			        const double rel_vel_mag = frame_A_vel_B.norm();
			        if (rel_vel_mag > 1e-13)  // Divide by zero protection.
			        {
			            // In the frame of A, B applies force in the direction of B's velocity.
			            slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
			        }

			        // Compute rolling friction force:
			        const double w_diff_mag = w_diff.norm();
			        if (w_diff_mag > 1e-13)  // Divide by zero protection.
			        {
			            rollForceA =
			                -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
			        }

			        // Total forces on a:
			        //Took out gravity force
			        totalForceOnA = elasticForceOnA + slideForceOnA + vdwForceOnA;
			        // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;

			        // Total torque a and b:
			        torqueA = r_a.cross(slideForceOnA + rollForceA);
			        torqueB = r_b.cross(-slideForceOnA + rollForceA);

			        P.aacc[A] += torqueA / P.moi[A];
			        P.aacc[B] += torqueB / P.moi[B];

			        if (write_step) {
			            // No factor of 1/2. Includes both spheres:
			            // P.PE += -G * P.m[A] * P.m[B] / dist + 0.5 * k * overlap * overlap;

			            // Van Der Waals + elastic:
			            const double diffRaRb = P.R[A] - P.R[B];
			            const double z = sumRaRb + h;
			            const double two_RaRb = 2 * P.R[A] * P.R[B];
			            const double denom_sum = z * z - (sumRaRb * sumRaRb);
			            const double denom_diff = z * z - (diffRaRb * diffRaRb);
			            const double U_vdw =
			                -Ha / 6 *
			                (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
			            P.PE += U_vdw + 0.5 * k * overlap * overlap;
			        }
			    } else  // Non-contact forces:
			    {
			        // No collision: Include gravity and vdw:
			        // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

			        // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
			        double h = std::fabs(overlap);
			        if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
			        {
			            h = h_min;
			        }
			        const double Ra = P.R[A];
			        const double Rb = P.R[B];
			        const double h2 = h * h;
			        const double twoRah = 2 * Ra * h;
			        const double twoRbh = 2 * Rb * h;
			        const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
			                                 ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
			                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
			                                                   (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
			                                 rVecab.normalized();

			        totalForceOnA = vdwForceOnA;  // +gravForceOnA;
			        if (write_step) {
			            // O.PE += -G * O.m[A] * O.m[B] / dist; // Gravitational

			            const double diffRaRb = P.R[A] - P.R[B];
			            const double z = sumRaRb + h;
			            const double two_RaRb = 2 * P.R[A] * P.R[B];
			            const double denom_sum = z * z - (sumRaRb * sumRaRb);
			            const double denom_diff = z * z - (diffRaRb * diffRaRb);
			            const double U_vdw =
			                -Ha / 6 *
			                (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
			            P.PE += U_vdw;  // Van Der Waals
			        }

			        // todo this is part of push_apart. Not great like this.
			        // For pushing apart overlappers:
			        // P.vel[A] = { 0,0,0 };
			        // P.vel[B] = { 0,0,0 };
			    }

			    // Newton's equal and opposite forces applied to acceleration of each ball:
			    P.acc[A] += totalForceOnA / P.m[A];
			    P.acc[B] -= totalForceOnA / P.m[B];

			    // So last distance can be known for COR:
			    P.distances[e] = dist;
            }
        }
        // }
        // t.end_event("loopApplicablePairs");  
    }
}

bool
check_arrs(vec3* arr1,vec3* arr2,int size,bool show=false)
{
	for (int i = 0; i < size; ++i)
	{
		if (show)
		{
			std::cout<<arr1[i]<<": "<<arr2[i]<<'\n';
		}
		else
		{
			if (!(arr1[i] == arr2[i]))
			{
				return false;
			}
		}
	}
	return true;
}

void addvec3s(vec3 &a,vec3 &b)
{
	a = a + b;
}

void 
parallel_test1()
{
	#pragma omp declare reduction(addvec3s: vec3: addvec3s(omp_out,omp_in)) 
	int balls = 1000000;
	vec3* pos = new vec3[balls];
    #pragma omp parallel for default(none) shared(balls,std::cout) reduction(addvec3s:pos[0:balls])
	for (int Ball = 0; Ball < balls; Ball++)
	{
		pos[Ball] = {Ball*1.0,(Ball*1.0)*17.0,(Ball*1.0)/17.0};
	}	
}

void addptrs(vec3 &a,vec3 &b)
{
	a = a + b;
}

void  
parallel_test3()
{
	#pragma omp declare reduction(addptrs: vec3: addptrs(omp_out,omp_in)) 
	int balls = 1000;
	// vec3 pos[balls];
	// double pos = 0;
	vec3* pos = new vec3[balls];
    #pragma omp parallel for default(none) shared(balls) reduction(addptrs:pos[:balls])
	for (int Ball = 0; Ball < balls; Ball++)
	{
		pos[Ball] = {Ball*1.0,(Ball*1.0)*17.0,(Ball*1.0)/17.0};
	}	
}

void 
parallel_test2()
{
	int balls = 1000000000;
	vec3* pos = new vec3[balls];
    #pragma omp parallel for default(none) shared(balls,pos) 
	for (int Ball = 0; Ball < balls; Ball++)
	{
		pos[Ball] = {Ball*1.0,(Ball*1.0)*17.0,(Ball*1.0)/17.0};
	}	
}




void 
serial_test1()
{
	int balls = 1000000000;
	vec3* pos = new vec3[balls];
	for (int Ball = 0; Ball < balls; Ball++)
	{
		pos[Ball] = {Ball*1.0,(Ball*1.0)*17.0,(Ball*1.0)/17.0};
	}	
}

vec3*
parallel_accel(const bool write_step, Ball_group& O)
{
	// vec3 *pos = new vec3[O.num_particles];
	// vec3 *pos = O.pos;
	// for (int Ball = 0; Ball < O.num_particles; Ball++)
	// {
	// 	pos[Ball] = O.pos[Ball];
	// }

	double dt = 1e-5;
	// #pragma omp declare reduction(addptrs: vec3: addptrs(omp_out,omp_in))
    for (int Ball = 0; Ball < O.num_particles; Ball++) {
        // Update velocity half step:
        O.velh[Ball] = O.vel[Ball] + .5 * O.acc[Ball] * dt;

        // Update angular velocity half step:
        O.wh[Ball] = O.w[Ball] + .5 * O.aacc[Ball] * dt;

        // Update position:
        O.pos[Ball] += O.velh[Ball] * dt;

        // Reinitialize acceleration to be recalculated:
        O.acc[Ball] = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        O.aacc[Ball] = {0, 0, 0};
    }

    grid g(O.num_particles, O.R[0], O.gridSize, O.pos, O.tolerance);

    
    vec3 *acc = new vec3[O.num_particles]({0,0,0});
    vec3 *aacc = new vec3[O.num_particles]({0,0,0});
    // for (int i = 0; i < )
    // #pragma omp for 
    // #pragma omp parallel for default(none) shared(dt,O) reduction(addptrs:pos[0:O.num_particles])
	#pragma omp declare reduction(addptrs: vec3: addptrs(omp_out,omp_in))
    #pragma omp parallel for default(none) shared(std::cout,write_step,O,g,h_min,kin,kout,u_s,u_r,Ha) reduction(addptrs:acc[0:O.num_particles],aacc[0:O.num_particles])
    for (int A = 0; A < O.num_particles; ++A)  // cuda
    {
        // t.start_event("getNearestNeighbors");
        std::vector<int> nearest_neighbors = g.getBalls(A);
        // t.end_event("getNearestNeighbors");
        // g.printVector(nearest_neighbors);
        // g.printVector(g.getBalls(A));
        // t.start_event("loopApplicablePairs");
        // #pragma omp parallel default(none) private(A) shared(nearest_neighbors,O,write_step)
        // {
        // #pragma omp for 
        for (int B : nearest_neighbors)
        {
            if (A != B)
            {
            	const double sumRaRb = O.R[A] + O.R[B];
			    const vec3 rVecab = O.pos[B] - O.pos[A];  // Vector from a to b.
			    const vec3 rVecba = -rVecab;
			    const double dist = (rVecab).norm();

			    // Check for collision between Ball and otherBall:
			    double overlap = sumRaRb - dist;

			    vec3 totalForceOnA{0, 0, 0};

			    // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
			    int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
			    double oldDist = O.distances[e];

			    // Check for collision between Ball and otherBall.
			    if (overlap > 0) {
			        double k;
			        // Apply coefficient of restitution to balls leaving collision.
			        if (dist >= oldDist) {
			            k = kout;
			        } else {
			            k = kin;
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

			        // Elastic force:
			        const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

			        // Gravity force:
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
			        if (rel_vel_mag > 1e-13)  // Divide by zero protection.
			        {
			            // In the frame of A, B applies force in the direction of B's velocity.
			            slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
			        }

			        // Compute rolling friction force:
			        const double w_diff_mag = w_diff.norm();
			        if (w_diff_mag > 1e-13)  // Divide by zero protection.
			        {
			            rollForceA =
			                -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
			        }

			        // Total forces on a:
			        //Took out gravity force
			        totalForceOnA = elasticForceOnA + slideForceOnA + vdwForceOnA;
			        // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;

			        // Total torque a and b:
			        torqueA = r_a.cross(slideForceOnA + rollForceA);
			        torqueB = r_b.cross(-slideForceOnA + rollForceA);

			        aacc[A] += torqueA / O.moi[A];
			        aacc[B] += torqueB / O.moi[B];

			        if (write_step) {
			            // No factor of 1/2. Includes both spheres:
			            // O.PE += -G * O.m[A] * O.m[B] / dist + 0.5 * k * overlap * overlap;

			            // Van Der Waals + elastic:
			            const double diffRaRb = O.R[A] - O.R[B];
			            const double z = sumRaRb + h;
			            const double two_RaRb = 2 * O.R[A] * O.R[B];
			            const double denom_sum = z * z - (sumRaRb * sumRaRb);
			            const double denom_diff = z * z - (diffRaRb * diffRaRb);
			            const double U_vdw =
			                -Ha / 6 *
			                (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
			            O.PE += U_vdw + 0.5 * k * overlap * overlap;
			        }
			    } else  // Non-contact forces:
			    {
			        // No collision: Include gravity and vdw:
			        // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] / (dist * dist)) * (rVecab / dist);

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

			        totalForceOnA = vdwForceOnA;  // +gravForceOnA;
			        if (write_step) {
			            // O.PE += -G * O.m[A] * O.m[B] / dist; // Gravitational

			            const double diffRaRb = O.R[A] - O.R[B];
			            const double z = sumRaRb + h;
			            const double two_RaRb = 2 * O.R[A] * O.R[B];
			            const double denom_sum = z * z - (sumRaRb * sumRaRb);
			            const double denom_diff = z * z - (diffRaRb * diffRaRb);
			            const double U_vdw =
			                -Ha / 6 *
			                (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
			            O.PE += U_vdw;  // Van Der Waals
			        }

			        // todo this is part of push_apart. Not great like this.
			        // For pushing apart overlappers:
			        // O.vel[A] = { 0,0,0 };
			        // O.vel[B] = { 0,0,0 };
			    }

			    // Newton's equal and opposite forces applied to acceleration of each ball:
			    acc[A] += totalForceOnA / O.m[A];
			    acc[B] -= totalForceOnA / O.m[B];

			    // So last distance can be known for COR:
			    O.distances[e] = dist;
            }
        }
        // }
        // t.end_event("loopApplicablePairs");  
    }
    // #pragma omp parallel for 
    for (int i = 0; i < O.num_particles; ++i)
    {

    	O.acc[i] = acc[i];
    	O.aacc[i] = aacc[i];
    }
    delete[] acc;
    delete[] aacc;
    return O.acc;
}

void test_value_setting()
{
	double radius = 1e-5;
	double grid_size = 1e-4;
	double tol = 3*radius;
	int num_particles = 230;
	int times_extra_particles = 20;
	std::string folder = "/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab_stable/SpaceLab/jobs/large_aggregate_optO3_1/N_1000/T_3/";
	std::string file = "228_2_R4e-05_v4e-01_cor0.63_mu0.1_rho2.25_k4e+00_Ha5e-12_dt5e-10_";
	// vec3 *pos = nullptr;
	// pos = new vec3[200];
	// vec3 *velh = nullptr;
	// velh = new vec3[200];
	// vec3 *vel = nullptr;
	// vel = new vec3[200];
	// vec3 *acc = nullptr;
	// acc = new vec3[200];
	// vec3 *w = nullptr;
	// w = new vec3[200];
	// vec3 *wh = nullptr;
	// wh = new vec3[200];
	// vec3 *aacc = nullptr;
	// aacc = new vec3[200];
	// aacc = new vec3[200];
	// parseSimData(getLastLine(file),pos,vel,w,200);
	int partial_step = 0;
	Ball_group O = Ball_group(num_particles*times_extra_particles,folder.c_str(),file.c_str(),true);
	vec3 *par_pos = new vec3[O.num_particles];
	vec3 *seq_pos = new vec3[O.num_particles];
    parseSimData(getLastLine(folder+file+"simData.csv"),O.pos,O.vel,O.w,num_particles*times_extra_particles);
	for (int Ball = 0; Ball < O.num_particles; Ball++) {
        // Update velocity half step:
        par_pos[Ball] = O.pos[Ball];
        partial_step = Ball%(O.num_particles/times_extra_particles);
        O.pos[Ball] = O.pos[partial_step]*partial_step;
        O.w[Ball] = O.w[partial_step];
        O.vel[Ball] = O.vel[partial_step];
    }


	
	// for (int i = 0; i < O.num_particles; ++i)
	// {
 //    	par[i] = O.acc[i];
	// }


	Ball_group P = Ball_group(num_particles*times_extra_particles,folder.c_str(),file.c_str(),true);
	parseSimData(getLastLine(folder+file+"simData.csv"),P.pos,P.vel,P.w,num_particles*times_extra_particles);
	partial_step = 0;
	for (int Ball = 0; Ball < P.num_particles; Ball++) {
        seq_pos[Ball] = P.pos[Ball];
        // Update velocity half step:
        partial_step = Ball%(P.num_particles/times_extra_particles);

        // P.velh[Ball] = {0,0,0};
        // P.wh[Ball] = {0,0,0};

        // Update angular velocity half step:
        P.pos[Ball] = P.pos[partial_step]*partial_step;
        P.w[Ball] = P.w[partial_step];
        P.vel[Ball] = P.vel[partial_step];

        // Update position:

        // // Reinitialize acceleration to be recalculated:
        // P.acc[Ball] = {0, 0, 0};

        // // Reinitialize angular acceleration to be recalculated:
        // P.aacc[Ball] = {0, 0, 0};
    }
	std::cout<<"The pos arrays are the same: "<<check_arrs(par_pos,seq_pos,num_particles)<<'\n';
	std::cout<<"The vel arrays are the same: "<<check_arrs(O.vel,P.vel,num_particles)<<'\n';
	std::cout<<"The w   arrays are the same: "<<check_arrs(O.w,P.w,num_particles)<<'\n';

	// double dt = 1e-5;
	timey t;
	double time = omp_get_wtime();
	t.start_event("parallel_time");
	parallel_accel(false,O);
	t.end_event("parallel_time");
	// timey t;
	std::cout<<"parallel time: "<<omp_get_wtime()-time<<'\n';

	time = omp_get_wtime();
	t.start_event("serial_time");
    serial_accel(false,P);
	t.end_event("serial_time");
	std::cout<<"sequential time: "<<omp_get_wtime()-time<<'\n';
	t.print_events();
	std::cout<<"The arrays are the same: "<<check_arrs(P.acc,O.acc,num_particles,false)<<'\n';
	// wrapper w(200,file,2e-5);
	
}



int main(int argc, char const *argv[])
{

	test_value_setting();


	// test_value_setting();
	// timey t;
	// t.start_event("c++11");
	// parallel_test3();
	// t.end_event("c++11");

	// t.start_event("reduction");
	// parallel_accel();
	// t.end_event("reduction");

	// t.start_event("serial");
	// serial_accel();
	// t.end_event("serial");

	// t.print_events();

	return 0;
}

