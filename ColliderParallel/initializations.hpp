#include "../json/single_include/nlohmann/json.hpp"
#include "../Utils.hpp"
#include <mutex>

std::mutex g_mutex;

using json = nlohmann::json;
using std::numbers::pi;

struct P;
struct P_pair;
struct P_group;

struct P
{
    vec3
        pos,
        vel,
        velh,
        w,
        wh,
        acc,
        aacc;

    double
        R = 0,
        m = 0,
        moi = 0;
};

struct P_pair
{
    P* A;
    P* B;
    double dist;
    P_pair() = default;
    P_pair(P* a, P* b, double d) : A(a), B(b), dist(d) {} // Init all
    P_pair(P* a, P* b) : A(a), B(b), dist(-1.0) {} // init pairs and set D to illogical distance
    //p_pair(const p_pair& in_pair) : A(in_pair.A), B(in_pair.B), dist(in_pair.dist) {} // Copy constructor
};

struct P_group
{
    int n; // Number of particles in the group
    std::vector<P> p_group; // Group of particles
    double U; // Potential energy
    double T; // Kinetic Energy
    vec3 mom;
    vec3 ang_mom;
    // Useful values:
    double rMin = -1;
    double rMax = -1;
    double mTotal = -1;
    double initialRadius = -1;
    double vCollapse = 0;
    double vMax = -1;
    double vMaxPrev = HUGE_VAL;
    double soc = -1;
    double gridSize = -1.0;
    double PE = 0.0, KE = 0.0;
    double density = -1.0;
    int seed = -1;
    std::string s_location;
    std::string projectileName;
    std::string targetName;
    std::string output_prefix;

    double u_s = -1.0, u_r = -1.0;
    double cor  =-1.0;
    double simTimeSeconds = -1.0;
    double timeResolution = -1.0;
    double fourThirdsPiRho = -1.0;
    double scaleBalls = -1.0; //ball radius
    double maxOverlap = -1.0;
    double v_custom = -1.0;
    double temp = -1.0;
    double kConsts = -1.0;
    double impactParameter = -1.0;
    double Ha = -1.0;
    double h_min = -1.0;
    double cone = -1.0;
    double dt = -1.0;
    double kin = -1.0;
    double kout = -1.0;
    double spaceRange = -1.0;
    double spaceRangeIncrement = -1.0;
    double simTimeElapsed = -1.0;
    double G = -1.0;

    int properties = -1;
    int genBalls = -1;
    int attempts = -1;
    int steps = -1;
    int skip = -1;
    int num_particles = -1;

    void parse_input_file(char const* location);
    void oneSizeSphere(const int nBalls);
    [[nodiscard]] vec3 getCOM() const;
    [[nodiscard]] double get_radius(const vec3& center) const;
    void calibrate_dt(int const Step, const double& customSpeed = -1.);
    void updateDTK(const double& velocity);
    void calc_v_collapse();
    [[nodiscard]] double getVelMax();
    [[nodiscard]] double getRmax();
    [[nodiscard]] double getRmin();
    void calc_helpfuls();


    P_group() = default;
    P_group(int n) : n(n) {}
    P_group(std::vector<P> p_group) : p_group(p_group) {}

    P_group(const bool generate, const double& customVel, const char* path)
    {
        parse_input_file(path);//should be first in constructor
        num_particles = genBalls;
        generate_ball_field(genBalls);

        // calc_v_collapse();
        calibrate_dt(0, customVel);
        // simInit_cond_and_center(true);
    }

    // [[nodiscard]] double getMass()
    double getMass()
    {
        auto lambda = [&](double sum, const P &b){return sum + b.m; };
        
        return std::accumulate(p_group.begin(), p_group.end(), 0.0, lambda);        
    }

    void generate_ball_field(const int nBalls)
    {
        std::cerr << "CLUSTER FORMATION\n";

        // allocate_group(nBalls);

        // Create new random number set.
        //      const int seedSave = static_cast<int>(time(nullptr));
        srand(seed);  // srand(seedSave);

        oneSizeSphere(nBalls);
        
        
        calc_helpfuls();
        // threeSizeSphere(nBalls);

        // output_prefix = std::to_string(nBalls) + "_R" + scientific(get_radius(getCOM())) + "_v" +
        //                 scientific(v_custom) + "_cor" + rounder(sqrtf(cor), 4) + "_mu" + rounder(u_s, 3) +
        //                 "_rho" + rounder(density, 4);
    }

    void update_kinematics(P& P)
    {
        // Update velocity half step:
        P.velh = P.vel + .5 * P.acc * dt;

        // Update angular velocity half step:
        P.wh = P.w + .5 * P.aacc * dt;

        // Update position:
        P.pos += P.velh * dt;

        // Reinitialize acceleration to be recalculated:
        P.acc = { 0, 0, 0 };

        // Reinitialize angular acceleration to be recalculated:
        P.aacc = { 0, 0, 0 };
    }

    void compute_velocity(P& P)
    {
        // Velocity for next step:
        P.vel = P.velh + .5 * P.acc * dt;
        P.w = P.wh + .5 * P.aacc * dt;
    }

    void compute_acceleration(P_pair& p_pair, bool writeStep)
    {
        const double Ra = p_pair.A->R;
        const double Rb = p_pair.B->R;
        const vec3 rVecab = p_pair.B -> pos - p_pair.A -> pos;
        const vec3 rVecba = -rVecab;
        const double m_a = p_pair.A->m;
        const double m_b = p_pair.B->m;
        const double sumRaRb = Ra + Rb;
        vec3 rVec = p_pair.B->pos - p_pair.A->pos; // Start with rVec from a to b.
        const double dist = (rVec).norm();
        vec3 totalForce;

        // Cohesion:
        // h is the "separation" of the particles at particle radius - maxOverlap.
        // This allows particles to be touching while under vdwForce.
        // const double h = maxOverlap * 1.01 - overlap;
        const double h = h_min;
        const double h2 = h * h;
        const double twoRah = 2 * Ra * h;
        const double twoRbh = 2 * Rb * h;
        const vec3 vdwForce =
            Ha / 6 *
            64 * Ra * Ra * Ra * Rb * Rb * Rb *
            (h + Ra + Rb) /
            (
                (h2 + twoRah + twoRbh) *
                (h2 + twoRah + twoRbh) *
                (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                (h2 + twoRah + twoRbh + 4 * Ra * Rb)
                ) *
            rVec.normalized();



        // Check for collision between Ball and otherBall:
        double overlap = sumRaRb - dist;

        double oldDist = p_pair.dist;

        // Check for collision between Ball and otherBall.
        if (overlap > 0)
        {
            double k;
            // Apply coefficient of restitution to balls leaving collision.
            if (dist >= oldDist)
            {
                k = kout;
            }
            else
            {
                k = kin;
            }

            

            // Elastic a:
            vec3 elasticForce = -k * overlap * .5 * (rVec / dist);
            const double elastic_force_A_mag = elasticForce.norm();
            // Friction a:
            vec3 dVel = p_pair.B->vel - p_pair.A->vel;
            vec3 frictionForce = { 0, 0, 0 };
            const vec3 r_a = rVecab * p_pair.A -> R / sumRaRb;  // Center to contact point
            const vec3 r_b = rVecba * p_pair.B -> R / sumRaRb;
            const vec3 w_diff = p_pair.A -> w - p_pair.B -> w;
            const double w_diff_mag = w_diff.norm();
            const vec3 relativeVelOfA = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - 
                                        p_pair.A->w.cross(p_pair.A->R / sumRaRb * rVec) - 
                                        p_pair.B->w.cross(p_pair.B->R / sumRaRb * rVec);
            double relativeVelMag = relativeVelOfA.norm();
            
            vec3 slideForceOnA, rollForceA;
            if (relativeVelMag > 1e-10) // When relative velocity is very low, dividing its vector components by its magnitude below is unstable.
            {
                slideForceOnA = u_s * (elasticForce.norm()) *
                                (relativeVelOfA / relativeVelMag);
                // frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) *
                //              (relativeVelOfA / relativeVelMag);
            }

            if (w_diff_mag > 1e-13)  // Divide by zero protection.
            {
                rollForceA =
                    -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
            }

            // Torque a:
            const vec3 aTorque = (p_pair.A->R) * rVec.cross(slideForceOnA + rollForceA);
            const vec3 bTorque = (p_pair.B->R) * rVec.cross(-slideForceOnA + rollForceA);
            // const vec3 aTorque = (p_pair.A->R / sumRaRb) * rVec.cross(slideForceOnA + rollForceA);
            // const vec3 bTorque = (p_pair.B->R / sumRaRb) * rVec.cross(-slideForceOnA + rollForceA);

            // Gravity on a:
            // const vec3 gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);

            // Total forces on a:
            totalForce = elasticForce + slideForceOnA + vdwForce;
            // totalForce = gravForceOnA + elasticForce + frictionForce + vdwForce;

            // Elastic and Friction b:
            // Flip direction b -> a:
            // rVec = -rVec;
            // dVel = -dVel;
            // elasticForce = -elasticForce;

            // const vec3 relativeVelOfB = dVel - dVel.dot(rVec) * (rVec / (dist * dist)) - p_pair.B->w.cross(p_pair.B->R / sumRaRb * rVec) - p_pair.A->w.cross(p_pair.A->R / sumRaRb * rVec);
            // relativeVelMag = relativeVelOfB.norm(); // todo - This should be the same as mag for A. Same speed different direction.
            // // if (relativeVelMag > 1e-10)
            // // {
            // //     frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) * (relativeVelOfB / relativeVelMag);
            // // }
            // const vec3 bTorque = (p_pair.B->R / sumRaRb) * rVec.cross(frictionForce);

            {
                const std::lock_guard<std::mutex> lock(g_mutex);
                p_pair.A->aacc += aTorque / p_pair.A->moi;
            }
            {
                const std::lock_guard<std::mutex> lock(g_mutex);
                p_pair.B->aacc += bTorque / p_pair.B->moi;
            }


            if (writeStep)
            {
                // Calculate potential energy. Important to recognize that the factor of 1/2 is not in front of K because this is for the spring potential in each ball and they are the same potential.
                //O.PE += -G * pair.A->m * pair.B->m / dist + 0.5 * k * overlap * overlap;
            }
        }
        else
        {
            // No collision: Include gravity only:
            // const vec3 gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);
            totalForce = vdwForce;
            // totalForce = gravForceOnA;
            if (writeStep)
            {
                //O.PE += -G * pair.A->m * pair.B->m / dist;
            }

        //  // For expanding overlappers:
        //  //pair.A->vel = { 0,0,0 };
        //  //pair.B->vel = { 0,0,0 };
        }

        // Newton's equal and opposite forces applied to acceleration of each ball:
        {
            const std::lock_guard<std::mutex> lock(g_mutex);
            p_pair.A->acc += totalForce / p_pair.A->m;
        }
        {
            const std::lock_guard<std::mutex> lock(g_mutex);
            p_pair.B->acc -= totalForce / p_pair.B->m;
        }
    }

    
};

// //@brief sets Ball_group object based on the need for a restart or not
// P_group make_group(const char *argv1,int* restart)
// {
//     P_group O;
    
//     //See if run has already been started
//     std::string filename = check_restart(argv1,restart);
//     if (*restart > -1) //Restart is necessary unless only first write has happended so far
//     {
//         if (*restart > 1)
//         {//TESTED
//             (*restart)--;
//             // filename = std::to_string(*restart) + filename;
//             filename = filename.substr(1,filename.length());
//             O = Ball_group(argv1,filename,v_custom,*restart);
//         }
//         else if (*restart == 1) //restart from first write (different naming convension for first write)
//         {//TESTED
//             (*restart)--;
//             filename = filename.substr(1,filename.length());
//             // exit(EXIT_SUCCESS);
//             O = Ball_group(argv1,filename,v_custom,*restart);
//         }
//         else //if restart is 0, need to rerun whole thing
//         {//TESTED
//             O = Ball_group(true, v_custom, argv1); // Generate new group
//         }

//     }
//     else // Make new ball group
//     {
//         *restart = 0;
//         O = Ball_group(true, v_custom, argv1); // Generate new group
//     }
//     return O;
// }

void P_group::oneSizeSphere(const int nBalls)
{

    for (int Ball = 0; Ball < nBalls; ++Ball) 
    {
        P new_p;
        new_p.R = scaleBalls;
        new_p.m = density * 4. / 3. * 3.14159 * std::pow(new_p.R, 3);
        new_p.moi = .4 * new_p.m * new_p.R * new_p.R;
        new_p.w = {0,0,0};
        new_p.pos = rand_vec3(spaceRange);
        p_group.push_back(new_p);
    }

    mTotal = getMass();

    // Generate non-overlapping spherical particle field:
    int collisionDetected = 0;
    int oldCollisions = nBalls;

    for (int failed = 0; failed < attempts; failed++) {
        for (int A = 0; A < nBalls; A++) {
            for (int B = A + 1; B < nBalls; B++) {
                // Check for Ball overlap.
                const double dist = (p_group[A].pos - p_group[B].pos).norm();
                const double sumRaRb = p_group[A].R + p_group[B].R;
                const double overlap = dist - sumRaRb;
                if (overlap < 0) {
                    collisionDetected += 1;
                    // Move the other ball:
                    p_group[B].pos = rand_vec3(spaceRange);
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
                p_group[Ball].pos = rand_vec3(
                    spaceRange);  // Each time we fail and increase range, redistribute all balls randomly
                                  // so we don't end up with big balls near mid and small balls outside.
            }
        }
        collisionDetected = 0;
    }

    std::cerr << "Final spacerange: " << spaceRange << '\n';
    std::cerr << "Initial Radius: " << get_radius(getCOM()) << '\n';
    std::cerr << "Mass: " << mTotal << '\n';
}

void P_group::parse_input_file(char const* location)
{
    std::string loc(location);
    s_location = loc;
    std::string json_file = loc + "input.json";
    // std::cout<<json_file<<std::endl;
    std::ifstream ifs(json_file);
    json inputs = json::parse(ifs);

    if (inputs["seed"] == std::string("default"))
    {
        seed = static_cast<int>(time(nullptr));
    }
    else
    {
        seed = static_cast<int>(inputs["seed"]);
    }

    // dynamicTime = inputs["dynamicTime"];

    gridSize = inputs["gridSize"];
    // tolerance = inputs["gridTolerance"];
    
    G = inputs["G"];
    density = inputs["density"];
    u_s = inputs["u_s"];
    u_r = inputs["u_r"];
    // sigma = inputs["sigma"];
    // Y = inputs["Y"];
    cor = inputs["cor"];
    simTimeSeconds = inputs["simTimeSeconds"];
    timeResolution = inputs["timeResolution"];
    fourThirdsPiRho = 4. / 3. * pi * density;
    scaleBalls = inputs["scaleBalls"];
    maxOverlap = inputs["maxOverlap"];
    // KEfactor = inputs["KEfactor"];
    if (inputs["v_custom"] != std::string("default"))
    {
        v_custom = inputs["v_custom"];
        // v_custom = 0.36301555459799423;
    }
    // else
    // {
    // }
    temp = inputs["temp"]; // this will modify v_custom in oneSizeSphere
    double temp_kConst = inputs["kConsts"];
    kConsts = temp_kConst * (fourThirdsPiRho / (maxOverlap * maxOverlap));
    impactParameter = inputs["impactParameter"];
    Ha = inputs["Ha"];
    // std::cout<<"in Ha: "<<Ha<<'\n';
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
    // z0Rot = inputs["z0Rot"];
    // y0Rot = inputs["y0Rot"];
    // z1Rot = inputs["z1Rot"];
    // y1Rot = inputs["y1Rot"];
    simTimeElapsed = inputs["simTimeElapsed"];
    // project_path = inputs["project_path"];
    // if (project_path == std::string("default"))
    // {
    //     project_path = s_location;
    // }
    // output_folder = inputs["output_folder"];
    // out_folder = inputs["output_folder"];
    // if (output_folder == std::string("default"))
    // {
    //     output_folder = s_location;
    //     out_folder = s_location;
    // }
    projectileName = inputs["projectileName"];
    targetName = inputs["targetName"];
    output_prefix = inputs["output_prefix"];
    if (output_prefix == std::string("default"))
    {
        output_prefix = "";
    }
}

/// Approximate the radius of the ballGroup.
[[nodiscard]] double P_group::get_radius(const vec3& center) const
{
    double radius = 0;

    if (num_particles > 1) {
        for (size_t i = 0; i < num_particles; i++) {
            const auto this_radius = (p_group[i].pos - center).norm();
            if (this_radius > radius) radius = this_radius;
        }
    } else {
        radius = p_group[0].R;
    }

    return radius;
}

[[nodiscard]] vec3 P_group::getCOM() const
{
    if (mTotal > 0) {
        vec3 comNumerator;
        for (int Ball = 0; Ball < num_particles; ++Ball) { comNumerator += p_group[Ball].m * p_group[Ball].pos; }
        vec3 com = comNumerator / mTotal;
        return com;
    } else {
        std::cerr << "Mass of cluster is zero.\n";
        exit(EXIT_FAILURE);
    }
}

void P_group::calc_v_collapse()
{
    // Sim fall velocity onto cluster:
    // vCollapse shrinks if a ball escapes but velMax should take over at that point, unless it is
    // ignoring far balls.
    double position = 0;
    while (position < initialRadius) {
        // todo - include vdw!!!
        vCollapse += G * mTotal / (initialRadius * initialRadius) * 0.1;
        position += vCollapse * 0.1;
    }
    vCollapse = fabs(vCollapse);
}

void P_group::calibrate_dt(int const Step, const double& customSpeed)
{
    const double dtOld = dt;

    if (customSpeed > 0.) {
        updateDTK(customSpeed);
        std::cerr << "CUSTOM SPEED: " << customSpeed;
    } else {
        // std::cerr << vCollapse << " <- vCollapse | Lazz Calc -> " << M_PI * M_PI * G * pow(density, 4.
        // / 3.) * pow(mTotal, 2. / 3.) * rMax;

        vMax = getVelMax();

        std::cerr << '\n';

        // Take whichever velocity is greatest:
        std::cerr << vCollapse << " = vCollapse | vMax = " << vMax;
        if (vMax < vCollapse) { vMax = vCollapse; }

        if (vMax < vMaxPrev) {
            updateDTK(vMax);
            vMaxPrev = vMax;
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

void P_group::updateDTK(const double& velocity)
{
    calc_helpfuls();
    kin = kConsts * rMax * velocity * velocity;
    kout = cor * kin;
    const double h2 = h_min * h_min;
    const double four_R_min = 4 * rMin * h_min;
    const double vdw_force_max = Ha / 6 * 64 * rMin * rMin * rMin * rMin * rMin * rMin *
                                 ((h_min + rMin + rMin) / ((h2 + four_R_min) * (h2 + four_R_min) *
                                                             (h2 + four_R_min + 4 * rMin * rMin) *
                                                             (h2 + four_R_min + 4 * rMin * rMin)));
    // todo is it rmin*rmin or rmin*rmax
    const double elastic_force_max = kin * maxOverlap * rMin;
    const double regime = (vdw_force_max > elastic_force_max) ? vdw_force_max : elastic_force_max;
    const double regime_adjust = regime / (maxOverlap * rMin);
    dt = .01 * sqrt((fourThirdsPiRho / regime_adjust) * rMin * rMin * rMin);
}

/// get max velocity
[[nodiscard]] double P_group::getVelMax()
{
    vMax = 0;

    // todo - make this a manual set true or false to use soc so we know if it is being used or not.
    if (soc > 0) {
        int counter = 0;
        for (int Ball = 0; Ball < num_particles; Ball++) {
            // Only consider balls moving toward com and within 4x initial radius around it.
            const vec3 fromCOM = p_group[Ball].pos - getCOM();
            if (acos(p_group[Ball].vel.normalized().dot(fromCOM.normalized())) > cone && fromCOM.norm() < soc) {
                if (p_group[Ball].vel.norm() > vMax) { vMax = p_group[Ball].vel.norm(); }
            } else {
                counter++;
            }
        }
        std::cerr << '(' << counter << " spheres ignored"
                  << ") ";
    } else {
        for (int Ball = 0; Ball < num_particles; Ball++) {
            if (p_group[Ball].vel.norm() > vMax) { vMax = p_group[Ball].vel.norm(); }
        }

        // Is vMax for some reason unreasonably small? Don't proceed. Probably a finished sim.
        // This shouldn't apply to extremely destructive collisions because it is possible that no
        // particles are considered, so it will keep pausing.
        if (vMax < 1e-10) {
            std::cerr << "\nMax velocity in system is less than 1e-10, probably a finished sim.\n";
            std::cerr << "now exiting sim\n";
            exit(EXIT_FAILURE);
            // system("pause");
        }
    }

    return vMax;
}

[[nodiscard]] double P_group::getRmin()
{
    rMin = p_group[0].R;
    for (int Ball = 1; Ball < num_particles; Ball++) {
        if (p_group[Ball].R < rMin) { rMin = p_group[Ball].R; }
    }
    return rMin;
}

[[nodiscard]] double P_group::getRmax()
{
    rMax = p_group[0].R;
    for (int Ball = 0; Ball < num_particles; Ball++) {
        if (p_group[Ball].R > rMax) { rMax = p_group[Ball].R; }
    }
    return rMax;
}

void P_group::calc_helpfuls()
{
    rMin = getRmin();
    rMax = getRmax();
    mTotal = getMass();
    initialRadius = get_radius(getCOM());
    soc = 4 * rMax + initialRadius;
}