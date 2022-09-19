#pragma once
#include "dust_const_init.hpp"
#include "graph_group.hpp"
// #include "dust_const.hpp"
#include "json/single_include/nlohmann/json.hpp"
// #include "json/single_include/nlohmann/json.hpp"
#include "vec3.hpp"
#include "linalg.hpp"
#include "Utils.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits.h>
#include <cstring>
#include <typeinfo>

using std::numbers::pi;
using json = nlohmann::json;

/// @brief Facilitates the concept of a group of balls with physical properties.
class Ball_group
{
public:
    int num_particles = 0;
    int num_groups = -1;
    int num_particles_added = 0;

    // Useful values:
    double r_min = -1;
    double r_max = -1;
    double m_total = -1;
    double initial_radius = -1;
    double v_collapse = 0;
    double v_max = -1;
    double v_max_prev = HUGE_VAL;
    double soc = -1;

    vec3 mom = {0, 0, 0};
    vec3 ang_mom = {
        0,
        0,
        0};  // Can be vec3 because they only matter for writing out to file. Can process on host.

    double PE = 0, KE = 0;

    double* distances = nullptr;

    vec3* pos = nullptr;
    vec3* vel = nullptr;
    vec3* velh = nullptr;  ///< Velocity half step for integration purposes.
    vec3* acc = nullptr;
    vec3* w = nullptr;
    vec3* wh = nullptr;  ///< Angular velocity half step for integration purposes.
    vec3* aacc = nullptr;
    double* R = nullptr;    ///< Radius
    double* m = nullptr;    ///< Mass
    double* moi = nullptr;  ///< Moment of inertia

    graph *g = nullptr;


    // double get_soc()
    // {
    //     return soc;
    // }

    void parse_input_file(char const* location)
    {
        std::string s_location(location);
        std::string json_file = s_location + "input.json";
        std::ifstream ifs(json_file);
        json inputs = json::parse(ifs);

        dynamicTime = inputs["dynamicTime"];
        G = inputs["G"];
        density = inputs["density"];
        u_s = inputs["u_s"];
        u_r = inputs["u_r"];
        sigma = inputs["sigma"];
        Y = inputs["Y"];
        cor = inputs["cor"];
        simTimeSeconds = inputs["simTimeSeconds"];
        timeResolution = inputs["timeResolution"];
        fourThirdsPiRho = 4. / 3. * pi * density;
        scaleBalls = inputs["scaleBalls"];
        maxOverlap = inputs["maxOverlap"];
        KEfactor = inputs["KEfactor"];
        if (inputs["v_custom"] == std::string("default"))
        {
            v_custom = 0.36301555459799423;
        }
        else
        {
            v_custom = inputs["v_custom"];
        }
        temp = inputs["temp"]; // this will modify v_custom in oneSizeSphere
        double temp_kConst = inputs["kConsts"];
        kConsts = temp_kConst * (fourThirdsPiRho / (maxOverlap * maxOverlap));
        impactParameter = inputs["impactParameter"];
        Ha = inputs["Ha"];
        double temp_h_min = inputs["h_min"];
        h_min = temp_h_min * scaleBalls;
        if (inputs["cone"] == std::string("default"))
        {
            cone = pi/2;
        }
        else
        {
            cone = inputs["cone"];
        }
        properties = inputs["properties"];
        genBalls = inputs["genBalls"];
        attempts = inputs["attempts"];
        skip = inputs["skip"];
        steps = inputs["steps"];
        dt = inputs["dt"];
        kin = inputs["kin"];
        kout = inputs["kout"];
        if (inputs["spaceRange"] == std::string("default"))
        {
            spaceRange = 4 * std::pow(
                            (1. / .74 * scaleBalls * scaleBalls * scaleBalls * genBalls),
                            1. / 3.); 
        }
        else
        {
            spaceRange = inputs["spaceRange"];
        }
        if (inputs["spaceRangeIncrement"] == std::string("default"))
        {
            spaceRangeIncrement = scaleBalls * 3;
        }
        else
        {
            spaceRangeIncrement = inputs["spaceRangeIncrement"];
        }
        z0Rot = inputs["z0Rot"];
        y0Rot = inputs["y0Rot"];
        z1Rot = inputs["z1Rot"];
        y1Rot = inputs["y1Rot"];
        simTimeElapsed = inputs["simTimeElapsed"];
        project_path = inputs["project_path"];
        if (project_path == std::string("default"))
        {
            project_path = s_location;
        }
        output_folder = inputs["output_folder"];
        if (output_folder == std::string("default"))
        {
            output_folder = s_location;
        }
        projectileName = inputs["projectileName"];
        targetName = inputs["targetName"];
        output_prefix = inputs["output_prefix"];
        if (output_prefix == std::string("default"))
        {
            output_prefix = "";
        }
    }

    Ball_group() = default;

    /// @brief For creating a new ballGroup of size nBalls
    /// @param nBalls Number of balls to allocate.
    explicit Ball_group(const int nBalls)
    {
        allocate_group(nBalls);

        for (size_t i = 0; i < nBalls; i++) {
            R[i] = 1;
            m[i] = 1;
            moi[i] = calc_moi(R[i], m[i]);
        }
    }

    /// @brief For generating a new ballGroup of size nBalls
    /// @param nBalls Number of balls to allocate.
    /// @param generate Just here to get you to the right constructor. This is definitely wrong.
    /// @param customVel To condition for specific vMax.
    Ball_group(const bool generate, const double& customVel, const char* path)
    {
        parse_input_file(path);
        generate_ball_field(genBalls);

        // Hack - Override and creation just 2 balls position and velocity.
        pos[0] = {0, 1.101e-5, 0};
        pos[1] = {0, -1.101e-5, 0};
        vel[0] = {0, 0, 0};
        vel[1] = {0, 0, 0};

        m_total = getMass();
        calc_v_collapse();
        calibrate_dt(0, customVel);
        simInit_cond_and_center(true);
    }

    /// @brief For continuing a sim.
    /// @param fullpath is the filename and path excluding the suffix _simData.csv, _constants.csv, etc.
    /// @param customVel To condition for specific vMax.
    explicit Ball_group(const std::string& path, const std::string& filename, const double& customVel, int start_file_index=0)
    {
        parse_input_file(path.c_str());
        sim_continue(path, filename,start_file_index);
        calc_v_collapse();
        calibrate_dt(0, customVel);
        simInit_cond_and_center(false);
    }

    /// @brief For two cluster sim.
    /// @param projectileName
    /// @param targetName
    /// @param customVel To condition for specific vMax.
    explicit Ball_group(
        const std::string& path,
        const std::string& projectileName,
        const std::string& targetName,
        const double& customVel)
    {
        std::cerr<<path<<std::endl;
        sim_init_two_cluster(path, projectileName, targetName);
        calc_v_collapse();
        calibrate_dt(0, customVel);
        simInit_cond_and_center(true);
    }

    Ball_group& operator=(const Ball_group& rhs)
    {
        num_particles = rhs.num_particles;
        num_particles_added = rhs.num_particles_added;

        // Useful values:
        r_min = rhs.r_min;
        r_max = rhs.r_max;
        m_total = rhs.m_total;
        initial_radius = rhs.initial_radius;
        v_collapse = rhs.v_collapse;
        v_max = rhs.v_max;
        v_max_prev = rhs.v_max_prev;
        soc = rhs.soc;

        mom = rhs.mom;
        ang_mom = rhs.ang_mom;  // Can be vec3 because they only matter for writing out to file. Can process
                                // on host.

        PE = rhs.PE;

        distances = rhs.distances;

        pos = rhs.pos;
        vel = rhs.vel;
        velh = rhs.velh;  ///< Velocity half step for integration purposes.
        acc = rhs.acc;
        w = rhs.w;
        wh = rhs.wh;  ///< Angular velocity half step for integration purposes.
        aacc = rhs.aacc;
        R = rhs.R;      ///< Radius
        m = rhs.m;      ///< Mass
        moi = rhs.moi;  ///< Moment of inertia

        g = rhs.g; ///graph class

        return *this;
    }

    Ball_group(const Ball_group& rhs)
    {
        num_particles = rhs.num_particles;
        num_particles_added = rhs.num_particles_added;

        // Useful values:
        r_min = rhs.r_min;
        r_max = rhs.r_max;
        m_total = rhs.m_total;
        initial_radius = rhs.initial_radius;
        v_collapse = rhs.v_collapse;
        v_max = rhs.v_max;
        v_max_prev = rhs.v_max_prev;
        soc = rhs.soc;

        mom = rhs.mom;
        ang_mom = rhs.ang_mom;  // Can be vec3 because they only matter for writing out to file. Can process
                                // on host.

        PE = rhs.PE;

        distances = rhs.distances;

        pos = rhs.pos;
        vel = rhs.vel;
        velh = rhs.velh;  ///< Velocity half step for integration purposes.
        acc = rhs.acc;
        w = rhs.w;
        wh = rhs.wh;  ///< Angular velocity half step for integration purposes.
        aacc = rhs.aacc;
        R = rhs.R;      ///< Radius
        m = rhs.m;      ///< Mass
        moi = rhs.moi;  ///< Moment of inertia

        g = rhs.g; ///graph class
    }

    void calc_helpfuls()
    {
        r_min = getRmin();
        r_max = getRmax();
        m_total = getMass();
        initial_radius = get_radius(getCOM());
        soc = 4 * r_max + initial_radius;
    }

    void calibrate_dt(int const Step, const double& customSpeed = -1.)
    {
        const double dtOld = dt;

        if (customSpeed > 0.) {
            updateDTK(customSpeed);
            std::cerr << "CUSTOM SPEED: " << customSpeed;
        } else {
            // std::cerr << vCollapse << " <- vCollapse | Lazz Calc -> " << M_PI * M_PI * G * pow(density, 4.
            // / 3.) * pow(mTotal, 2. / 3.) * rMax;

            v_max = getVelMax();

            std::cerr << '\n';

            // Take whichever velocity is greatest:
            std::cerr << v_collapse << " = vCollapse | vMax = " << v_max;
            if (v_max < v_collapse) { v_max = v_collapse; }

            if (v_max < v_max_prev) {
                updateDTK(v_max);
                v_max_prev = v_max;
                std::cerr << "\nk: " << kin << "\tdt: " << dt;
            }
        }

        if (Step == 0 or dtOld < 0) {
            steps = static_cast<int>(simTimeSeconds / dt);
            std::cerr << "\tInitial Steps: " << steps << '\n';
        } else {
            steps = static_cast<int>(dtOld / dt) * (steps - Step) + Step;
            std::cerr << "\tSteps: " << steps;
        }

        if (timeResolution / dt > 1.) {
            skip = static_cast<int>(floor(timeResolution / dt));
            std::cerr << "\tSkip: " << skip << '\n';
        } else {
            std::cerr << "Desired time resolution is lower than dt. Setting to 1 second per skip.\n";
            skip = static_cast<int>(floor(1. / dt));
        }
    }

    // todo - make bigger balls favor the middle, or, smaller balls favor the outside.
    /// @brief Push balls apart until no overlaps
    void pushApart() const
    {
        std::cerr << "Separating spheres - Current max overlap:\n";
        /// Using acc array as storage for accumulated position change.
        int* counter = new int[num_particles];
        for (int Ball = 0; Ball < num_particles; Ball++) {
            acc[Ball] = {0, 0, 0};
            counter[Ball] = 0;
        }

        double overlapMax = -1;
        const double pseudoDT = r_min * .1;
        int step = 0;

        while (true) {
            // if (step % 10 == 0)
            //{
            //	simDataWrite("pushApart_");
            //}

            for (int A = 0; A < num_particles; A++) {
                for (int B = A + 1; B < num_particles; B++) {
                    // Check for Ball overlap.
                    vec3 rVecab = pos[B] - pos[A];
                    vec3 rVecba = -1 * rVecab;
                    const double dist = (rVecab).norm();
                    const double sumRaRb = R[A] + R[B];
                    const double overlap = sumRaRb - dist;

                    if (overlapMax < overlap) { overlapMax = overlap; }

                    if (overlap > 0) {
                        acc[A] += overlap * (rVecba / dist);
                        acc[B] += overlap * (rVecab / dist);
                        counter[A] += 1;
                        counter[B] += 1;
                    }
                }
            }

            for (int Ball = 0; Ball < num_particles; Ball++) {
                if (counter[Ball] > 0) {
                    pos[Ball] += acc[Ball].normalized() * pseudoDT;
                    acc[Ball] = {0, 0, 0};
                    counter[Ball] = 0;
                }
            }

            if (overlapMax > 0) {
                std::cerr << overlapMax << "                        \r";
            } else {
                std::cerr << "\nSuccess!\n";
                break;
            }
            overlapMax = -1;
            step++;
        }
        delete[] counter;
    }

    void calc_v_collapse()
    {
        // Sim fall velocity onto cluster:
        // vCollapse shrinks if a ball escapes but velMax should take over at that point, unless it is
        // ignoring far balls.
        double position = 0;
        while (position < initial_radius) {
            // todo - include vdw!!!
            v_collapse += G * m_total / (initial_radius * initial_radius) * 0.1;
            position += v_collapse * 0.1;
        }
        v_collapse = fabs(v_collapse);
    }

    /// get max velocity
    [[nodiscard]] double getVelMax()
    {
        v_max = 0;

        // todo - make this a manual set true or false to use soc so we know if it is being used or not.
        if (soc > 0) {
            int counter = 0;
            for (int Ball = 0; Ball < num_particles; Ball++) {
                // Only consider balls moving toward com and within 4x initial radius around it.
                const vec3 fromCOM = pos[Ball] - getCOM();
                if (acos(vel[Ball].normalized().dot(fromCOM.normalized())) > cone && fromCOM.norm() < soc) {
                    if (vel[Ball].norm() > v_max) { v_max = vel[Ball].norm(); }
                } else {
                    counter++;
                }
            }
            std::cerr << '(' << counter << " spheres ignored"
                      << ") ";
        } else {
            for (int Ball = 0; Ball < num_particles; Ball++) {
                if (vel[Ball].norm() > v_max) { v_max = vel[Ball].norm(); }
            }

            // Is vMax for some reason unreasonably small? Don't proceed. Probably a finished sim.
            // This shouldn't apply to extremely destructive collisions because it is possible that no
            // particles are considered, so it will keep pausing.
            if (v_max < 1e-10) {
                std::cerr << "\nMax velocity in system is less than 1e-10.\n";
                system("pause");
            }
        }

        return v_max;
    }

    // Kick ballGroup (give the whole thing a velocity)
    void kick(const vec3& vec) const
    {
        for (int Ball = 0; Ball < num_particles; Ball++) { vel[Ball] += vec; }
    }


    vec3 calc_momentum(const std::string& of = "") const
    {
        vec3 pTotal = {0, 0, 0};
        for (int Ball = 0; Ball < num_particles; Ball++) { pTotal += m[Ball] * vel[Ball]; }
        // fprintf(stderr, "%s Momentum Check: %.2e, %.2e, %.2e\n", of.c_str(), pTotal.x, pTotal.y, pTotal.z);
        return pTotal;
    }

    // offset cluster
    void offset(const double& rad1, const double& rad2, const double& impactParam) const
    {
        for (int Ball = 0; Ball < num_particles; Ball++) {
            pos[Ball].x += (rad1 + rad2) * cos(impactParam);
            pos[Ball].y += (rad1 + rad2) * sin(impactParam);
        }
    }

    /// Approximate the radius of the ballGroup.
    [[nodiscard]] double get_radius(const vec3& center) const
    {
        double radius = 0;

        if (num_particles > 1) {
            for (size_t i = 0; i < num_particles; i++) {
                const auto this_radius = (pos[i] - center).norm();
                if (this_radius > radius) radius = this_radius;
            }
        } else {
            radius = R[0];
        }

        return radius;
    }

    // Update Gravitational Potential Energy:
    void updateGPE()
    {
        PE = 0;

        if (num_particles > 1)  // Code below only necessary for effects between balls.
        {
            for (int A = 1; A < num_particles; A++) {
                for (int B = 0; B < A; B++) {
                    const double sumRaRb = R[A] + R[B];
                    const double dist = (pos[A] - pos[B]).norm();
                    const double overlap = sumRaRb - dist;

                    // Check for collision between Ball and otherBall.
                    if (overlap > 0) {
                        PE +=
                            -G * m[A] * m[B] / dist + kin * ((sumRaRb - dist) * .5) * ((sumRaRb - dist) * .5);
                    } else {
                        PE += -G * m[A] * m[B] / dist;
                    }
                }
            }
        } else  // For the case of just one ball:
        {
            PE = 0;
        }
    }

    void sim_init_write(std::string filename, int counter = 0)
    {
        // Create string for file name identifying spin combination negative is 2, positive is 1 on each axis.
        // std::string spinCombo = "";
        // for ( int i = 0; i < 3; i++)
        //{
        //	if (spins[i] < 0) { spinCombo += "2"; }
        //	else if (spins[i] > 0) { spinCombo += "1"; }
        //	else { spinCombo += "0"; }
        //}

        // todo - filename is now a copy and this works. Need to consider how old way worked for
        // compatibility. What happens without setting output_prefix = filename? Check if file name already
        // exists.
        std::ifstream checkForFile;
        checkForFile.open(output_folder + filename + "simData.csv", std::ifstream::in);
        // Add a counter to the file name until it isn't overwriting anything:
        while (checkForFile.is_open()) {
            counter++;
            checkForFile.close();
            checkForFile.open(output_folder + std::to_string(counter) + '_' + filename + "simData.csv", std::ifstream::in);
        }

        if (counter > 0) { filename.insert(0, std::to_string(counter) + '_'); }

        output_prefix = filename;

        // Complete file names:
        std::string simDataFilename = output_folder + filename + "simData.csv";
        std::string energyFilename = output_folder + filename + "energy.csv";
        std::string constantsFilename = output_folder + filename + "constants.csv";

        std::cerr << "New file tag: " << filename;

        // Open all file streams:
        std::ofstream energyWrite, ballWrite, constWrite;
        energyWrite.open(energyFilename, std::ofstream::app);
        ballWrite.open(simDataFilename, std::ofstream::app);
        constWrite.open(constantsFilename, std::ofstream::app);

        // Make column headers:
        energyWrite << "Time,PE,KE,E,p,L";
        ballWrite << "x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0";

        for (int Ball = 1; Ball < num_particles;
             Ball++)  // Start at 2nd ball because first one was just written^.
        {
            std::string thisBall = std::to_string(Ball);
            ballWrite << ",x" + thisBall << ",y" + thisBall << ",z" + thisBall << ",wx" + thisBall
                      << ",wy" + thisBall << ",wz" + thisBall << ",wmag" + thisBall << ",vx" + thisBall
                      << ",vy" + thisBall << ",vz" + thisBall << ",bound" + thisBall;
        }

        // Write constant data:
        for (int Ball = 0; Ball < num_particles; Ball++) {
            constWrite << R[Ball] << ',' << m[Ball] << ',' << moi[Ball] << '\n';
        }

        // Write energy data to buffer:
        energyBuffer << '\n'
                     << simTimeElapsed << ',' << PE << ',' << KE << ',' << PE + KE << ',' << mom.norm() << ','
                     << ang_mom.norm();
        energyWrite << energyBuffer.rdbuf();
        energyBuffer.str("");

        // Reinitialize energies for next step:
        KE = 0;
        PE = 0;
        mom = {0, 0, 0};
        ang_mom = {0, 0, 0};

        // Send position and rotation to buffer:
        ballBuffer << '\n';  // Necessary new line after header.
        ballBuffer << pos[0].x << ',' << pos[0].y << ',' << pos[0].z << ',' << w[0].x << ',' << w[0].y << ','
                   << w[0].z << ',' << w[0].norm() << ',' << vel[0].x << ',' << vel[0].y << ',' << vel[0].z
                   << ',' << 0;  // bound[0];
        for (int Ball = 1; Ball < num_particles; Ball++) {
            ballBuffer << ',' << pos[Ball].x
                       << ','  // Needs comma start so the last bound doesn't have a dangling comma.
                       << pos[Ball].y << ',' << pos[Ball].z << ',' << w[Ball].x << ',' << w[Ball].y << ','
                       << w[Ball].z << ',' << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y
                       << ',' << vel[Ball].z << ',' << 0;  // bound[Ball];
        }
        // Write position and rotation data to file:
        ballWrite << ballBuffer.rdbuf();
        ballBuffer.str("");  // Resets the stream buffer to blank.

        // Close Streams for user viewing:
        energyWrite.close();
        ballWrite.close();
        constWrite.close();

        std::cerr << "\nSimulating " << steps * dt / 60 / 60 << " hours.\n";
        std::cerr << "Total mass: " << m_total << '\n';
        std::cerr << "\n===============================================================\n";
    }


    [[nodiscard]] vec3 getCOM() const
    {
        if (m_total > 0) {
            vec3 comNumerator;
            for (int Ball = 0; Ball < num_particles; Ball++) { comNumerator += m[Ball] * pos[Ball]; }
            vec3 com = comNumerator / m_total;
            return com;
        } else {
            std::cerr << "Mass of cluster is zero.\n";
            exit(EXIT_FAILURE);
        }
    }

    void zeroVel() const
    {
        for (int Ball = 0; Ball < num_particles; Ball++) { vel[Ball] = {0, 0, 0}; }
    }

    void zeroAngVel() const
    {
        for (int Ball = 0; Ball < num_particles; Ball++) { w[Ball] = {0, 0, 0}; }
    }

    void to_origin() const
    {
        const vec3 com = getCOM();

        for (int Ball = 0; Ball < num_particles; Ball++) { pos[Ball] -= com; }
    }

    // Set velocity of all balls such that the cluster spins:
    void comSpinner(const double& spinX, const double& spinY, const double& spinZ) const
    {
        const vec3 comRot = {spinX, spinY, spinZ};  // Rotation axis and magnitude
        for (int Ball = 0; Ball < num_particles; Ball++) {
            vel[Ball] += comRot.cross(pos[Ball] - getCOM());
            w[Ball] += comRot;
        }
    }

    void rotAll(const char axis, const double angle) const
    {
        for (int Ball = 0; Ball < num_particles; Ball++) {
            pos[Ball] = pos[Ball].rot(axis, angle);
            vel[Ball] = vel[Ball].rot(axis, angle);
            w[Ball] = w[Ball].rot(axis, angle);
        }
    }

    double calc_mass(const double& radius, const double& density)
    {
        return density * 4. / 3. * 3.14159 * std::pow(radius, 3);
    }

    double calc_moi(const double& radius, const double& mass) { return .4 * mass * radius * radius; }

    Ball_group spawn_particles(const int count)
    {
        // Load file data:
        std::cerr << "Add Particle\n";

        // Random particle to origin
        Ball_group projectile(count);
        // Particle random position at twice radius of target:
        // We want the farthest from origin since we are offsetting form origin. Not com.
        const auto cluster_radius = 3;

        const vec3 projectile_direction = rand_vec3(1).normalized();
        projectile.pos[0] = projectile_direction * (cluster_radius + scaleBalls * 4);
        projectile.w[0] = {0, 0, 0};
        // Velocity toward origin:
        projectile.vel[0] = -v_custom * projectile_direction;
        projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
        projectile.m[0] = density * 4. / 3. * pi * std::pow(projectile.R[0], 3);
        projectile.moi[0] = calc_moi(projectile.R[0], projectile.m[0]);


        const double3x3 local_coords = local_coordinates(to_double3(projectile_direction));
        // to_vec3(local_coords.y).print();
        // to_vec3(local_coords.z).print();
        // projectile.pos[0].print();
        for (int i = 1; i < projectile.num_particles - 3; i++) {
            const auto rand_y = rand_between(-cluster_radius, cluster_radius);
            const auto rand_z = rand_between(-cluster_radius, cluster_radius);
            // projectile.pos[i] = projectile.pos[0] + perpendicular_shift(local_coords, rand_y, rand_z);
            projectile.pos[i] = projectile.pos[0] + perpendicular_shift(local_coords, rand_y, rand_z);
            // std::cout << rand_y << '\t' << to_vec3(local_coords.y * rand_y) <<'\t'<< projectile.pos[i] <<
            // '\n';
        }
        projectile.pos[projectile.num_particles - 3] = projectile_direction * 2;
        projectile.pos[projectile.num_particles - 2] = projectile_direction * 4;
        projectile.pos[projectile.num_particles - 1] = projectile_direction * 6;

        Ball_group new_group{projectile.num_particles + num_particles};

        new_group.merge_ball_group(*this);
        new_group.merge_ball_group(projectile);
        // new_group.init_graph();

        // new_group.calibrate_dt(0, 1);
        // new_group.init_conditions();

        // new_group.to_origin();

        return new_group;
    }

    //@brief returns new position of particle after it is given random offset
    //@param local_coords is plane perpendicular to direction of projectile
    //@param projectile_pos is projectile's position before offset is applied
    //@param projectile_vel is projectile's velocity
    //@param projectile_rad is projectile's radius
    vec3 dust_agglomeration_offset(
        const double3x3 local_coords,
        vec3 projectile_pos,
        vec3 projectile_vel,
        const double projectile_rad)
    {
        const auto cluster_radius = get_radius(vec3(0, 0, 0));
        bool intersect = false;
        int count = 0;
        vec3 new_position = vec3(0,0,0);
        do {
            const auto rand_y = rand_between(-cluster_radius, cluster_radius);
            const auto rand_z = rand_between(-cluster_radius, cluster_radius);
            auto test_pos = projectile_pos + perpendicular_shift(local_coords, rand_y, rand_z);

            count++;
            for (size_t i = 0; i < num_particles; i++) {
                // Check that velocity intersects one of the spheres:
                if (line_sphere_intersect(test_pos, projectile_vel, pos[i], R[i] + projectile_rad)) {
                    new_position = test_pos;
                    intersect = true;
                    break;
                }
            }
        } while (!intersect);
        return new_position;
    }

    // @brief returns new ball group consisting of one particle
    //        where particle is given initial conditions
    //        including an random offset linearly dependant on radius 
    Ball_group dust_agglomeration_particle_init()
    {
        // Random particle to origin
        Ball_group projectile(1);

        // Particle random position at twice radius of target:
        // We want the farthest from origin since we are offsetting form origin. Not com.
        const auto cluster_radius = get_radius(vec3(0, 0, 0));

        const vec3 projectile_direction = rand_unit_vec3();
        projectile.pos[0] = projectile_direction * (cluster_radius + scaleBalls * 4);
        projectile.R[0] = scaleBalls;  // rand_between(1,3)*1e-5;
        projectile.w[0] = {0, 0, 0};
        projectile.m[0] = density * 4. / 3. * pi * std::pow(projectile.R[0], 3);
        // Velocity toward origin:
        if (temp > 0)
        {
            double a = std::sqrt(Kb*temp/projectile.m[0]);
            v_custom = max_bolt_dist(a); 
        }
        projectile.vel[0] = -v_custom * projectile_direction;
        // projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
        projectile.moi[0] = calc_moi(projectile.R[0], projectile.m[0]);

        const double3x3 local_coords = local_coordinates(to_double3(projectile_direction));

        projectile.pos[0] = dust_agglomeration_offset(local_coords,projectile.pos[0],projectile.vel[0],projectile.R[0]);

        
        return projectile;
    } 

    

    // Uses previous O as target and adds one particle to hit it:
    Ball_group add_projectile()
    {
        // Load file data:
        std::cerr << "Add Particle\n";

        
        Ball_group projectile = dust_agglomeration_particle_init();

        
        // Collision velocity calculation:
        const vec3 p_target{calc_momentum("p_target")};
        const vec3 p_projectile{projectile.calc_momentum("p_particle")};
        const vec3 p_total{p_target + p_projectile};
        const double m_target{getMass()};
        const double m_projectile{projectile.getMass()};
        const double m_total{m_target + m_projectile};
        const vec3 v_com = p_total / m_total;

        // Negate total system momentum:
        projectile.kick(-v_com);
        kick(-v_com);

        fprintf(
            stderr,
            "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n",
            vel[0].norm(),
            projectile.vel[0].norm());

        std::cerr << '\n';
        projectile.calc_momentum("Projectile");
        calc_momentum("Target");

        Ball_group new_group{projectile.num_particles + num_particles};

        new_group.merge_ball_group(*this);
        new_group.merge_ball_group(projectile);

        // Hack - Calibrate to vel = 1 so we don't have to reform the pair. Probly fine?
        new_group.calibrate_dt(0, v_custom);
        new_group.calibrate_dt(0, 1);
        new_group.init_conditions();
        // new_group.g.printGraph();

        new_group.to_origin();

        // new_group.g.printGraph();
        return new_group;
    }

    /// @brief Add another ballGroup into this one.
    /// @param src The ballGroup to be added.
    void merge_ball_group(const Ball_group& src)
    {
        // Copy incoming data to the end of the currently loaded data.
        std::memcpy(
            &distances[num_particles_added], src.distances, sizeof(src.distances[0]) * src.num_particles);
        std::memcpy(&pos[num_particles_added], src.pos, sizeof(src.pos[0]) * src.num_particles);
        std::memcpy(&vel[num_particles_added], src.vel, sizeof(src.vel[0]) * src.num_particles);
        std::memcpy(&velh[num_particles_added], src.velh, sizeof(src.velh[0]) * src.num_particles);
        std::memcpy(&acc[num_particles_added], src.acc, sizeof(src.acc[0]) * src.num_particles);
        std::memcpy(&w[num_particles_added], src.w, sizeof(src.w[0]) * src.num_particles);
        std::memcpy(&wh[num_particles_added], src.wh, sizeof(src.wh[0]) * src.num_particles);
        std::memcpy(&aacc[num_particles_added], src.aacc, sizeof(src.aacc[0]) * src.num_particles);
        std::memcpy(&R[num_particles_added], src.R, sizeof(src.R[0]) * src.num_particles);
        std::memcpy(&m[num_particles_added], src.m, sizeof(src.m[0]) * src.num_particles);
        std::memcpy(&moi[num_particles_added], src.moi, sizeof(src.moi[0]) * src.num_particles);

        //TODO 
        // g = new graph(num_particles_added);

        // Keep track of now loaded ball set to start next set after it:
        num_particles_added += src.num_particles;
        calc_helpfuls();
    }

    void printGraph(std::string file)
    {
        std::cout<<"File: "<<file<<std::endl;
        g -> printGraph(file);
    }

    void printGraph()
    {
        // std::cout<<"File: "<<file<<std::endl;
        g -> printGraph();
    }

private:
    // String buffers to hold data in memory until worth writing to file:
    std::stringstream ballBuffer;
    std::stringstream energyBuffer;

    // void init_graph()
    // {
    //     for (int A = 1; A < num_particles; A++) 
    //     {
    //         for (int B = 0; B < A; B++) 
    //         {
    //             const double sumRaRb = R[A] + R[B];
    //             const double dist = (pos[A] - pos[B]).norm();
    //             const double overlap = sumRaRb - dist;

    //             // Check for collision between Ball and otherBall.
    //             // Collision exists
    //             if (overlap >= 0) 
    //             {
    //                 g.addEdge(A,B);  
    //             } 
    //         }
    //     }
    // }

    /// Allocate balls
    void allocate_group(const int nBalls)
    {
        num_particles = nBalls;

        try {
            distances = new double[(num_particles * num_particles / 2) - (num_particles / 2)];

            g = new graph(num_particles);
            
            pos = new vec3[num_particles];
            vel = new vec3[num_particles];
            velh = new vec3[num_particles];
            acc = new vec3[num_particles];
            w = new vec3[num_particles];
            wh = new vec3[num_particles];
            aacc = new vec3[num_particles];
            R = new double[num_particles];
            m = new double[num_particles];
            moi = new double[num_particles];


            // init(num_particles);
        } catch (const std::exception& e) {
            std::cerr << "Failed trying to allocate group. " << e.what() << '\n';
        }
    }


    /// @brief Deallocate arrays to recover memory.
    void freeMemory() const
    {
        delete[] distances;
        delete[] pos;
        delete[] vel;
        delete[] velh;
        delete[] acc;
        delete[] w;
        delete[] wh;
        delete[] aacc;
        delete[] R;
        delete[] m;
        delete[] moi;

        delete g;

    }


    // Initialize accelerations and energy calculations:
    void init_conditions()
    {
        std::cout<<"INIT COND CALLED"<<std::endl;
        // SECOND PASS - Check for collisions, apply forces and torques:
        for (int A = 1; A < num_particles; A++)  // cuda
        {
            // DONT DO ANYTHING HERE. A STARTS AT 1.
            for (int B = 0; B < A; B++) {
                const double sumRaRb = R[A] + R[B];
                const vec3 rVecab = pos[B] - pos[A];  // Vector from a to b.
                const vec3 rVecba = -rVecab;
                const double dist = (rVecab).norm();

                // Check for collision between Ball and otherBall:
                double overlap = sumRaRb - dist;

                vec3 totalForceOnA{0, 0, 0};

                // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
                int e = static_cast<int>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
                double oldDist = distances[e];

                // Check for collision between Ball and otherBall.
                if (overlap > 0) {
                    g -> addEdge(A,B);
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
                    const double Ra = R[A];
                    const double Rb = R[B];
                    const double h2 = h * h;
                    // constexpr double h2 = h * h;
                    const double twoRah = 2 * Ra * h;
                    const double twoRbh = 2 * Rb * h;
                    const vec3 vdwForceOnA =
                        Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                        ((h + Ra + Rb) /
                         ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                          (h2 + twoRah + twoRbh + 4 * Ra * Rb) * (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                        rVecab.normalized();

                    // Elastic force:
                    const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

                    // Gravity force:
                    const vec3 gravForceOnA = (G * m[A] * m[B] / (dist * dist)) * (rVecab / dist);

                    // Sliding and Rolling Friction:
                    vec3 slideForceOnA{0, 0, 0};
                    vec3 rollForceA{0, 0, 0};
                    vec3 torqueA{0, 0, 0};
                    vec3 torqueB{0, 0, 0};

                    // Shared terms:
                    const double elastic_force_A_mag = elasticForceOnA.norm();
                    const vec3 r_a = rVecab * R[A] / sumRaRb;  // Center to contact point
                    const vec3 r_b = rVecba * R[B] / sumRaRb;
                    const vec3 w_diff = w[A] - w[B];

                    // Sliding friction terms:
                    const vec3 d_vel = vel[B] - vel[A];
                    const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
                                               w[A].cross(r_a) - w[B].cross(r_a);

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
                    totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;

                    // Total torque a and b:
                    torqueA = r_a.cross(slideForceOnA + rollForceA);
                    torqueB = r_b.cross(-slideForceOnA + rollForceA);

                    aacc[A] += torqueA / moi[A];
                    aacc[B] += torqueB / moi[B];


                    // No factor of 1/2. Includes both spheres:
                    // PE += -G * m[A] * m[B] / dist + 0.5 * k * overlap * overlap;

                    // Van Der Waals + elastic:
                    const double diffRaRb = R[A] - R[B];
                    const double z = sumRaRb + h;
                    const double two_RaRb = 2 * R[A] * R[B];
                    const double denom_sum = z * z - (sumRaRb * sumRaRb);
                    const double denom_diff = z * z - (diffRaRb * diffRaRb);
                    const double U_vdw =
                        -Ha / 6 *
                        (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                    PE += U_vdw + 0.5 * k * overlap * overlap;

                } else  // Non-contact forces:
                {
                    // No collision: Include gravity and vdw:
                    // const vec3 gravForceOnA = (G * m[A] * m[B] / (dist * dist)) * (rVecab / dist);

                    // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic
                    // cancellation:
                    double h = std::fabs(overlap);
                    if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
                    {
                        h = h_min;
                    }
                    const double Ra = R[A];
                    const double Rb = R[B];
                    const double h2 = h * h;
                    const double twoRah = 2 * Ra * h;
                    const double twoRbh = 2 * Rb * h;
                    const vec3 vdwForceOnA =
                        Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                        ((h + Ra + Rb) /
                         ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                          (h2 + twoRah + twoRbh + 4 * Ra * Rb) * (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                        rVecab.normalized();

                    totalForceOnA = vdwForceOnA;  // +gravForceOnA;

                    // PE += -G * m[A] * m[B] / dist; // Gravitational

                    const double diffRaRb = R[A] - R[B];
                    const double z = sumRaRb + h;
                    const double two_RaRb = 2 * R[A] * R[B];
                    const double denom_sum = z * z - (sumRaRb * sumRaRb);
                    const double denom_diff = z * z - (diffRaRb * diffRaRb);
                    const double U_vdw =
                        -Ha / 6 *
                        (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                    PE += U_vdw;  // Van Der Waals


                    // todo this is part of push_apart. Not great like this.
                    // For pushing apart overlappers:
                    // vel[A] = { 0,0,0 };
                    // vel[B] = { 0,0,0 };
                }

                // Newton's equal and opposite forces applied to acceleration of each ball:
                acc[A] += totalForceOnA / m[A];
                acc[B] -= totalForceOnA / m[B];

                // So last distance can be known for COR:
                distances[e] = dist;
            }
            // DONT DO ANYTHING HERE. A STARTS AT 1.
        }

        // Calc energy:
        for (int Ball = 0; Ball < num_particles; Ball++) {
            KE += .5 * m[Ball] * vel[Ball].dot(vel[Ball]) + .5 * moi[Ball] * w[Ball].dot(w[Ball]);
            mom += m[Ball] * vel[Ball];
            ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
        }
    }


    [[nodiscard]] double getRmin()
    {
        r_min = R[0];
        for (int Ball = 1; Ball < num_particles; Ball++) {
            if (R[Ball] < r_min) { r_min = R[Ball]; }
        }
        return r_min;
    }

    [[nodiscard]] double getRmax()
    {
        r_max = R[0];
        for (int Ball = 0; Ball < num_particles; Ball++) {
            if (R[Ball] > r_max) { r_max = R[Ball]; }
        }
        return r_max;
    }


    [[nodiscard]] double getMassMax() const
    {
        double mMax = m[0];
        for (int Ball = 0; Ball < num_particles; Ball++) {
            if (m[Ball] > mMax) { mMax = m[Ball]; }
        }
        return mMax;
    }


    void parseSimData(std::string line)
    {
        std::string lineElement;
        // Get number of balls in file
        int count = 54 / properties;
        if (num_particles > 0)
        {
            count = num_particles;
        }
        // int count = std::count(line.begin(), line.end(), ',') / properties + 1;
        allocate_group(count);
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


    /// Get previous sim constants by filename.
    void loadConsts(const std::string& path, const std::string& filename)
    {
        // Get radius, mass, moi:
        std::string constantsFilename = path + filename + "constants.csv";
        if (auto ConstStream = std::ifstream(constantsFilename, std::ifstream::in)) {
            std::string line, lineElement;
            for (int A = 0; A < num_particles; A++) {
                std::getline(ConstStream, line);  // Ball line.
                std::stringstream chosenLine(line);
                std::getline(chosenLine, lineElement, ',');  // Radius.
                R[A] = std::stod(lineElement);
                std::getline(chosenLine, lineElement, ',');  // Mass.
                m[A] = std::stod(lineElement);
                std::getline(chosenLine, lineElement, ',');  // Moment of inertia.
                moi[A] = std::stod(lineElement);
            }
        } else {
            std::cerr << "Could not open constants file: " << constantsFilename << "... Existing program."
                      << '\n';
            exit(EXIT_FAILURE);
        }
    }


    /// Get last line of previous simData by filename.
    [[nodiscard]] static std::string getLastLine(const std::string& path, const std::string& filename)
    {
        std::string simDataFilepath = path + filename + "simData.csv";
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


    void simDataWrite(std::string& outFilename)
    {
        // todo - for some reason I need checkForFile instead of just using ballWrite. Need to work out why.
        // Check if file name already exists. If not, initialize
        std::ifstream checkForFile;
        checkForFile.open(output_folder + outFilename + "simData.csv", std::ifstream::in);
        if (checkForFile.is_open() == false) {
            sim_init_write(outFilename);
        } else {
            ballBuffer << '\n';  // Prepares a new line for incoming data.

            for (int Ball = 0; Ball < num_particles; Ball++) {
                // Send positions and rotations to buffer:
                if (Ball == 0) {
                    ballBuffer << pos[Ball][0] << ',' << pos[Ball][1] << ',' << pos[Ball][2] << ','
                               << w[Ball][0] << ',' << w[Ball][1] << ',' << w[Ball][2] << ','
                               << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y << ','
                               << vel[Ball].z << ',' << 0;
                } else {
                    ballBuffer << ',' << pos[Ball][0] << ',' << pos[Ball][1] << ',' << pos[Ball][2] << ','
                               << w[Ball][0] << ',' << w[Ball][1] << ',' << w[Ball][2] << ','
                               << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y << ','
                               << vel[Ball].z << ',' << 0;
                }
            }

            // Write simData to file and clear buffer.
            std::ofstream ballWrite;
            ballWrite.open(output_folder + outFilename + "simData.csv", std::ofstream::app);
            ballWrite << ballBuffer.rdbuf();  // Barf buffer to file.
            ballBuffer.str("");               // Resets the stream for that balls to blank.
            ballWrite.close();
        }
        checkForFile.close();
    }


    [[nodiscard]] double getMass()
    {
        m_total = 0;
        {
            for (int Ball = 0; Ball < num_particles; Ball++) { m_total += m[Ball]; }
        }
        return m_total;
    }

    void threeSizeSphere(const int nBalls)
    {
        // Make nBalls of 3 sizes in CGS with ratios such that the mass is distributed evenly among the 3
        // sizes (less large nBalls than small nBalls).
        const int smalls = static_cast<int>(std::round(
            static_cast<double>(nBalls) * 27. /
            31.375));  // Just here for reference. Whatever nBalls are left will be smalls.
        const int mediums = static_cast<int>(std::round(static_cast<double>(nBalls) * 27. / (8 * 31.375)));
        const int larges = static_cast<int>(std::round(static_cast<double>(nBalls) * 1. / 31.375));


        for (int Ball = 0; Ball < larges; Ball++) {
            // Below comment maintains asteroid radius while increasing particle count.
            // std::pow(1. / (double)nBalls, 1. / 3.) * 3. * scaleBalls;

            R[Ball] = 3. * scaleBalls;
            m[Ball] = density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
            moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
            w[Ball] = {0, 0, 0};
            pos[Ball] = rand_vec3(spaceRange);
        }

        for (int Ball = larges; Ball < (larges + mediums); Ball++) {
            R[Ball] = 2. * scaleBalls;  // std::pow(1. / (double)nBalls, 1. / 3.) * 2. * scaleBalls;
            m[Ball] = density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
            moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
            w[Ball] = {0, 0, 0};
            pos[Ball] = rand_vec3(spaceRange);
        }
        for (int Ball = (larges + mediums); Ball < nBalls; Ball++) {
            R[Ball] = 1. * scaleBalls;  // std::pow(1. / (double)nBalls, 1. / 3.) * 1. * scaleBalls;
            m[Ball] = density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
            moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
            w[Ball] = {0, 0, 0};
            pos[Ball] = rand_vec3(spaceRange);
        }

        m_total = 0;
        for (int i = 0; i < nBalls; i++)
        {
            m_total += m[i];
            std::cerr<<"Ball "<<i<<"\tmass is "<<m[i]<<"\t"<<"radius is "<<R[i]<<std::endl;
        }

        std::cerr << "Smalls: " << smalls << " Mediums: " << mediums << " Larges: " << larges << '\n';

        // Generate non-overlapping spherical particle field:
        int collisionDetected = 0;
        int oldCollisions = nBalls;

        for (int failed = 0; failed < attempts; failed++) {
            for (int A = 0; A < nBalls; A++) {
                for (int B = A + 1; B < nBalls; B++) {
                    // Check for Ball overlap.
                    const double dist = (pos[A] - pos[B]).norm();
                    const double sumRaRb = R[A] + R[B];
                    const double overlap = dist - sumRaRb;
                    if (overlap < 0) {
                        collisionDetected += 1;
                        // Move the other ball:
                        pos[B] = rand_vec3(spaceRange);
                    }
                }
            }
            if (collisionDetected < oldCollisions) {
                oldCollisions = collisionDetected;
                std::cerr << "Collisions: " << collisionDetected << "                        \r";
            }
            if (collisionDetected == 0) {
                std::cerr << "\nSuccess!\n";
                break;
            }
            if (failed == attempts - 1 ||
                collisionDetected >
                    static_cast<int>(
                        1.5 *
                        static_cast<double>(
                            nBalls)))  // Added the second part to speed up spatial constraint increase when
                                       // there are clearly too many collisions for the space to be feasible.
            {
                std::cerr << "Failed " << spaceRange << ". Increasing range " << spaceRangeIncrement
                          << "cm^3.\n";
                spaceRange += spaceRangeIncrement;
                failed = 0;
                for (int Ball = 0; Ball < nBalls; Ball++) {
                    pos[Ball] = rand_vec3(
                        spaceRange);  // Each time we fail and increase range, redistribute all balls randomly
                                      // so we don't end up with big balls near mid and small balls outside.
                }
            }
            collisionDetected = 0;
        }

        std::cerr << "Final spacerange: " << spaceRange << '\n';
        std::cerr << "m_total: " << m_total << '\n';
        std::cerr << "Initial Radius: " << get_radius(getCOM()) << '\n';
        std::cerr << "Mass: " << getMass() << '\n';
    }

    void generate_ball_field(const int nBalls)
    {
        std::cerr << "CLUSTER FORMATION\n";

        allocate_group(nBalls);

        // Create new random number set.
        //		const int seedSave = static_cast<int>(time(nullptr));
        srand(0);  // srand(seedSave);

        oneSizeSphere(nBalls);
        
        calc_helpfuls();
        // threeSizeSphere(nBalls);

        output_prefix = std::to_string(nBalls) + "_R" + scientific(get_radius(getCOM())) + "_v" +
                        scientific(v_custom) + "_cor" + rounder(sqrtf(cor), 4) + "_mu" + rounder(u_s, 3) +
                        "_rho" + rounder(density, 4);
    }

    /// Make ballGroup from file data.
    void loadSim(const std::string& path, const std::string& filename)
    {
        parseSimData(getLastLine(path, filename));
        loadConsts(path, filename);

        calc_helpfuls();

        std::cerr << "Balls: " << num_particles << '\n';
        std::cerr << "Mass: " << m_total << '\n';
        std::cerr << "Approximate radius: " << initial_radius << " cm.\n";
    }


    void oneSizeSphere(const int nBalls)
    {

        for (int Ball = 0; Ball < nBalls; Ball++) {
            R[Ball] = scaleBalls;
            m[Ball] = density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
            moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
            w[Ball] = {0, 0, 0};
            pos[Ball] = rand_vec3(spaceRange);
        }

        m_total = getMass();

        // Generate non-overlapping spherical particle field:
        int collisionDetected = 0;
        int oldCollisions = nBalls;

        for (int failed = 0; failed < attempts; failed++) {
            for (int A = 0; A < nBalls; A++) {
                for (int B = A + 1; B < nBalls; B++) {
                    // Check for Ball overlap.
                    const double dist = (pos[A] - pos[B]).norm();
                    const double sumRaRb = R[A] + R[B];
                    const double overlap = dist - sumRaRb;
                    if (overlap < 0) {
                        collisionDetected += 1;
                        // Move the other ball:
                        pos[B] = rand_vec3(spaceRange);
                    }
                }
            }
            if (collisionDetected < oldCollisions) {
                oldCollisions = collisionDetected;
                std::cerr << "Collisions: " << collisionDetected << "                        \r";
            }
            if (collisionDetected == 0) {
                std::cerr << "\nSuccess!\n";
                break;
            }
            if (failed == attempts - 1 ||
                collisionDetected >
                    static_cast<int>(
                        1.5 *
                        static_cast<double>(
                            nBalls)))  // Added the second part to speed up spatial constraint increase when
                                       // there are clearly too many collisions for the space to be feasible.
            {
                std::cerr << "Failed " << spaceRange << ". Increasing range " << spaceRangeIncrement
                          << "cm^3.\n";
                spaceRange += spaceRangeIncrement;
                failed = 0;
                for (int Ball = 0; Ball < nBalls; Ball++) {
                    pos[Ball] = rand_vec3(
                        spaceRange);  // Each time we fail and increase range, redistribute all balls randomly
                                      // so we don't end up with big balls near mid and small balls outside.
                }
            }
            collisionDetected = 0;
        }

        std::cerr << "Final spacerange: " << spaceRange << '\n';
        std::cerr << "Initial Radius: " << get_radius(getCOM()) << '\n';
        std::cerr << "Mass: " << m_total << '\n';
    }


    void updateDTK(const double& velocity)
    {
        calc_helpfuls();
        kin = kConsts * r_max * velocity * velocity;
        kout = cor * kin;
        const double h2 = h_min * h_min;
        const double four_R_min = 4 * r_min * h_min;
        const double vdw_force_max = Ha / 6 * 64 * r_min * r_min * r_min * r_min * r_min * r_min *
                                     ((h_min + r_min + r_min) / ((h2 + four_R_min) * (h2 + four_R_min) *
                                                                 (h2 + four_R_min + 4 * r_min * r_min) *
                                                                 (h2 + four_R_min + 4 * r_min * r_min)));
        // todo is it rmin*rmin or rmin*rmax
        const double elastic_force_max = kin * maxOverlap * r_min;
        const double regime = (vdw_force_max > elastic_force_max) ? vdw_force_max : elastic_force_max;
        const double regime_adjust = regime / (maxOverlap * r_min);
        dt = .01 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min);
    }

    void simInit_cond_and_center(bool add_prefix)
    {
        std::cerr << "==================" << '\n';
        std::cerr << "dt: " << dt << '\n';
        std::cerr << "k: " << kin << '\n';
        std::cerr << "Skip: " << skip << '\n';
        std::cerr << "Steps: " << steps << '\n';
        std::cerr << "==================" << '\n';

        to_origin();

        calc_momentum("After Zeroing");  // Is total mom zero like it should be?

        // Compute physics between all balls. Distances, collision forces, energy totals, total mass:
        init_conditions();

        // Name the file based on info above:
        if (add_prefix)
        {   
            output_prefix += "_k" + scientific(kin) + "_Ha" + scientific(Ha) + "_dt" + scientific(dt) + "_";
        }
    }


    void sim_continue(const std::string& path, const std::string& filename, int start_file_index=0)
    {
        // Load file data:
        num_particles = 3 + start_file_index;
        if (start_file_index == 0)
        {
            std::cerr << "Continuing Sim...\nFile: " << filename << '\n';
            loadSim(path, filename);
        }
        else
        {
            std::cerr << "Continuing Sim...\nFile: " << start_file_index << filename << '\n';
            loadSim(path, std::to_string(start_file_index) + "_" + filename);
        }



        std::cerr << '\n';
        calc_momentum("O");

        // Name the file based on info above:
        // output_prefix = std::to_string(num_particles) + "_rho" + rounder(density, 4);
        output_prefix = filename;
    }


    // Set's up a two cluster collision.
    void sim_init_two_cluster(
        const std::string& path,
        const std::string& projectileName,
        const std::string& targetName)
    {
        // Load file data:
        std::cerr << "TWO CLUSTER SIM\nFile 1: " << projectileName << '\t' << "File 2: " << targetName
                  << '\n';

        // DART PROBE
        // ballGroup projectile(1);
        // projectile.pos[0] = { 8814, 0, 0 };
        // projectile.w[0] = { 0, 0, 0 };
        // projectile.vel[0] = { 0, 0, 0 };
        // projectile.R[0] = 78.5;
        // projectile.m[0] = 560000;
        // projectile.moi[0] = .4 * projectile.m[0] * projectile.R[0] * projectile.R[0];


        Ball_group projectile;
        projectile.loadSim(path, projectileName);
        Ball_group target;
        target.loadSim(path, targetName);

        // DO YOU WANT TO STOP EVERYTHING?
        // projectile.zeroAngVel();
        // projectile.zeroVel();
        // target.zeroAngVel();
        // target.zeroVel();


        // Calc info to determined cluster positioning and collisions velocity:
        projectile.updateGPE();
        target.updateGPE();

        projectile.offset(
            projectile.initial_radius, target.initial_radius + target.getRmax() * 2, impactParameter);

        //		const double PEsys = projectile.PE + target.PE + (-G * projectile.mTotal * target.mTotal /
        //(projectile.getCOM() - target.getCOM()).norm());

        // Collision velocity calculation:
        const double mSmall = projectile.m_total;
        const double mBig = target.m_total;
        //		const double mTot = mBig + mSmall;
        // const double vSmall = -sqrt(2 * KEfactor * fabs(PEsys) * (mBig / (mSmall * mTot))); // Negative
        // because small offsets right.
        const double vSmall = -v_custom;                // DART probe override.
        const double vBig = -(mSmall / mBig) * vSmall;  // Negative to oppose projectile.
        // const double vBig = 0; // Dymorphous override.

        if (std::isnan(vSmall) || std::isnan(vBig)) {
            std::cerr << "A VELOCITY WAS NAN!!!!!!!!!!!!!!!!!!!!!!\n\n";
            exit(EXIT_FAILURE);
        }

        projectile.kick(vec3(vSmall, 0, 0));
        target.kick(vec3(vBig, 0, 0));

        fprintf(stderr, "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n", vBig, vSmall);

        std::cerr << '\n';
        projectile.calc_momentum("Projectile");
        target.calc_momentum("Target");

        allocate_group(projectile.num_particles + target.num_particles);

        merge_ball_group(target);
        merge_ball_group(projectile);  // projectile second so smallest ball at end and largest ball at front
                                       // for dt/k calcs.

        output_prefix = projectileName + targetName + "T" + rounder(KEfactor, 4) + "_vBig" +
                        scientific(vBig) + "_vSmall" + scientific(vSmall) + "_IP" +
                        rounder(impactParameter * 180 / 3.14159, 2) + "_rho" + rounder(density, 4);
    }
};