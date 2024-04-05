#pragma once
#include "../dust_const_init.hpp"
// #include "dust_const.hpp"
#include "../../external/json/single_include/nlohmann/json.hpp"
#include "../vec3.hpp"
#include "../linalg.hpp"
#include "../Utils.hpp"
#include "../../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits.h>
#include <cstring>
#include <typeinfo>
#include <random>
#include <omp.h>
#include <mpi.h>
#include <openacc.h>

// using std::numbers::pi;
const double pi = 3.14159265358979311599796346854;
using json = nlohmann::json;
const time_t start = time(nullptr);  // For end of program analysis

int getSize()
{
    int world_size;
    #ifdef MPI_ENABLE
        MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    #else
        world_size = 1;
    #endif
    return world_size;
}

int getRank()
{
    int world_rank;
    #ifdef MPI_ENABLE
        MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    #else
        world_rank = 0;
    #endif
    return world_rank;
}


/// @brief Facilitates the concept of a group of balls with physical properties.
class Ball_group
{
public:

    std::stringstream ballBuffer;
    std::stringstream energyBuffer;
    std::stringstream contactBuffer;

    double radiiFraction = -1;
    bool debug = false;
    bool write_all = false;

    std::string out_folder;
    int num_particles = 0;
    int num_particles_added = 0;
    int total_balls_to_add = 0;

    int OMPthreads = 8;
    double update_time = 0.0;

    std::string targetName;
    std::string projectileName;

    int seed = -1;
    int output_width = -1;
    enum distributions {constant, logNorm};
    distributions radiiDistribution;
    enum simType {BPCA, collider};
    simType typeSim;
    double lnSigma = 0.2; //sigma for log normal distribution 

    // Useful values:
    double r_min = -1;
    double r_max = -1;
    double m_total = -1;
    double initial_radius = -1;
    double v_collapse = 0;
    double v_max = -1;
    double v_max_prev = HUGE_VAL;
    double soc = -1;

    int world_rank;
    int world_size;
    int num_pairs;

    /////////////////////////////////
    const double h_min_physical = 2.1e-8; //prolly should make this a parameter/calculation
    const double max_mu = 0.5; // Make another parameter
    bool mu_scale = false;
    /////////////////////////////////

    vec3 mom = {0, 0, 0};
    vec3 ang_mom = {
        0,
        0,
        0};  // Can be vec3 because they only matter for writing out to file. Can process on host.

    double KE = 0;

    double PE = 0.0;

    double* distances = nullptr;

    vec3* pos = nullptr;
    vec3* vel = nullptr;
    vec3* velh = nullptr;  ///< Velocity half step for integration purposes.
    vec3* acc = nullptr;
    vec3* w = nullptr;
    vec3* wh = nullptr;  ///< Angular velocity half step for integration purposes.
    vec3* aacc = nullptr;
    vec3* accsq = nullptr;
    vec3* aaccsq = nullptr;
    double* R = nullptr;    ///< Radius
    double* m = nullptr;    ///< Mass
    double* moi = nullptr;  ///< Moment of inertia
    double* u_scale = nullptr; ///ADD TO COPY CONSTRUCTOR, ETC
    
    int ndev;
    int thegpu;
    ///////////////////////////
    // if (write_all)
    // {
    // vec3* slidDir = nullptr;
    // vec3* rollDir = nullptr;
    // double* inout = nullptr;
    // vec3* slidB3 = nullptr;
    // vec3* rollB3 = nullptr;
    // double* distB3 = nullptr;
    // vec3* slidFric = nullptr;
    // vec3* rollFric = nullptr;
    //////////
    //all force pointers
    vec3* vdwForce = nullptr;
    // vec3* elasticForce = nullptr;
    // vec3* slideForce = nullptr;
    // vec3* rollForce = nullptr;
    // vec3* torqueForce = nullptr;
    // }
    ///////////////////////////

    Ball_group() = default;

    explicit Ball_group(const int nBalls);
    explicit Ball_group(const std::string& path,const std::string& filename,const double& customVel, int start_file_index);
    explicit Ball_group(const std::string& path,const std::string& projectileName,const std::string& targetName,const double& customVel);
    Ball_group(const bool generate,const double& customVel,const char* path);
    Ball_group(const Ball_group& rhs);
    Ball_group& operator=(const Ball_group& rhs);
    void parse_input_file(char const* location);
    inline double calc_VDW_force_mag(const double Ra, const double Rb, const double h);
    void calc_mu_scale_factor();
    void zeroSaveVals();
    void calibrate_dt(int const Step, const double& customSpeed);
    void pushApart() const;
    void calc_v_collapse();
    [[nodiscard]] double getVelMax();
    void calc_helpfuls();
    double get_soc();    
    void kick(const vec3& vec) const;
    vec3 calc_momentum(const std::string& of) const;
    void offset(const double& rad1, const double& rad2, const double& impactParam) const;
    [[nodiscard]] double get_radius(const vec3& center) const;
    void updateGPE();
    void sim_init_write(std::string filename, int counter);
    [[nodiscard]] vec3 getCOM() const;
    void zeroVel() const;
    void zeroAngVel() const;
    void to_origin() const;
    void comSpinner(const double& spinX, const double& spinY, const double& spinZ) const;
    void rotAll(const char axis, const double angle) const;
    double calc_mass(const double& radius, const double& density);
    double calc_moi(const double& radius, const double& mass);
    Ball_group spawn_particles(const int count);
    vec3 dust_agglomeration_offset(const double3x3 local_coords,vec3 projectile_pos,vec3 projectile_vel,const double projectile_rad);
    Ball_group dust_agglomeration_particle_init();
    Ball_group add_projectile();
    void merge_ball_group(const Ball_group& src);

    void sim_one_step(const bool write_step);
    void sim_looper();
    void bufferBarf();

    // String buffers to hold data in memory until worth writing to file:
    // std::stringstream ballBuffer;
    // std::stringstream energyBuffer;

    void allocate_group(const int nBalls);
    void freeMemory() const;
    void init_conditions();
    [[nodiscard]] double getRmin();
    [[nodiscard]] double getRmax();
    [[nodiscard]] double getMassMax() const;
    void parseSimData(std::string line);
    void loadConsts(const std::string& path, const std::string& filename);
    [[nodiscard]] std::string getLastLine(const std::string& path, const std::string& filename);
    void simDataWrite(std::string& outFilename);
    [[nodiscard]] double getMass();
    void threeSizeSphere(const int nBalls);
    void generate_ball_field(const int nBalls);
    void loadSim(const std::string& path, const std::string& filename);
    void distSizeSphere(const int nBalls);
    void oneSizeSphere(const int nBalls);
    void placeBalls(const int nBalls);
    void updateDTK(const double& velocity);
    void simInit_cond_and_center(bool add_prefix);
    void sim_continue(const std::string& path, const std::string& filename, int start_file_index);
    void sim_init_two_cluster(const std::string& path,const std::string& projectileName,const std::string& targetName);
private:
};

/// @brief For creating a new ballGroup of size nBalls
/// @param nBalls Number of balls to allocate.
Ball_group::Ball_group(const int nBalls)
{
    energyBuffer.precision(12);  // Need more precision on momentum.
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
Ball_group::Ball_group(const bool generate, const double& customVel, const char* path)
{
    energyBuffer.precision(12);  // Need more precision on momentum.
    parse_input_file(path);
    generate_ball_field(genBalls);
    // Hack - Override and creation just 2 balls position and velocity.
    pos[0] = {0, 1.101e-5, 0};
    vel[0] = {0, 0, 0};
    if (genBalls > 1)
    {
        pos[1] = {0, -1.101e-5, 0};
        vel[1] = {0, 0, 0};
    }

    if (mu_scale)
    {
        calc_mu_scale_factor();
    }
    std::cerr<<initial_radius<<std::endl;

    m_total = getMass();
    calc_v_collapse();
    // std::cerr<<"INIT VCUSTOM "<<v_custom<<std::endl;
    calibrate_dt(0, v_custom);
    simInit_cond_and_center(true);
    // std::cerr<<pos[0]<<std::endl;
    // std::cerr<<vel[0]<<std::endl;

    // std::cerr<<pos[1]<<std::endl;
    // std::cerr<<vel[1]<<std::endl;
}

/// @brief For continuing a sim.
/// @param fullpath is the filename and path excluding the suffix _simData.csv, _constants.csv, etc.
/// @param customVel To condition for specific vMax.
Ball_group::Ball_group(const std::string& path, const std::string& filename, const double& customVel,int start_file_index=0)
{
    energyBuffer.precision(12);  // Need more precision on momentum.
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
Ball_group::Ball_group(
    const std::string& path,
    const std::string& projectileName,
    const std::string& targetName,
    const double& customVel=-1.)
{
    parse_input_file(path.c_str());
    // std::cerr<<path<<std::endl;
    sim_init_two_cluster(path, projectileName, targetName);
    calc_v_collapse();
    if (customVel > 0){calibrate_dt(0, customVel);}
    else {calibrate_dt(0, v_custom);}
    simInit_cond_and_center(true);
}

Ball_group& Ball_group::operator=(const Ball_group& rhs)
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
    KE = rhs.KE;

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

    radiiDistribution = rhs.radiiDistribution;
    radiiFraction = rhs.radiiFraction;


    /////////////////////////////////////////////////////////////
    // dynamicTime = rhs.dynamicTime;

    // G = rhs.G;  // Gravitational constant
    // density = rhs.density;
    // u_s = rhs.u_s;       // Coeff of sliding friction
    // u_r = rhs.u_r;               // Coeff of rolling friction
    // sigma = rhs.sigma;              // Poisson ratio for rolling friction.
    // Y = rhs.Y;               // Young's modulus in erg/cm3
    // cor = rhs.cor;                // Coeff of restitution
    // simTimeSeconds = rhs.simTimeSeconds;  // Seconds
    // timeResolution = rhs.timeResolution;    // Seconds - This is duration between exported steps.
    // fourThirdsPiRho = rhs.fourThirdsPiRho;  // for fraction of smallest sphere radius.
    // scaleBalls = rhs.scaleBalls;                         // base radius of ball.
    // maxOverlap = rhs.maxOverlap;                           // of scaleBalls
    // KEfactor = rhs.KEfactor;                              // Determines collision velocity based on KE/PE
    // v_custom = rhs.v_custom;  // Velocity cm/s
    // temp = rhs.temp;          //tempurature of simulation in Kelvin
    // kConsts = rhs.kConsts;
    // impactParameter = rhs.impactParameter;  // Impact angle radians
    // Ha = rhs.Ha;         // Hamaker constant for vdw force
    // h_min = rhs.h_min;  // 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
    // cone = rhs.cone;  // Cone of particles ignored moving away from center of mass. Larger angle ignores more.

    // // Simulation Structure
    // properties = rhs.properties;  // Number of columns in simData file per ball
    // genBalls = rhs.genBalls;
    // attempts = rhs.attempts;  // How many times to try moving every ball touching another in generator.

    // skip = rhs.skip;  // Steps thrown away before recording a step to the buffer. 500*.04 is every 20 seconds in sim.
    // steps = rhs.steps;

    // dt = rhs.dt;
    // kin = rhs.kin;  // Spring constant
    // kout = rhs.kout;
    // spaceRange = rhs.spaceRange;  // Rough minimum space required
    // spaceRangeIncrement = rhs.spaceRangeIncrement;
    // z0Rot = rhs.z0Rot;  // Cluster one z axis rotation
    // y0Rot = rhs.y0Rot;  // Cluster one y axis rotation
    // z1Rot = rhs.z1Rot;  // Cluster two z axis rotation
    // y1Rot = rhs.y1Rot;  // Cluster two y axis rotation
    // simTimeElapsed = rhs.simTimeElapsed;

    // // File from which to proceed with further simulations
    // project_path = rhs.project_path;
    // output_folder = rhs.output_folder;
    // // projectileName;
    // // targetName;
    // output_prefix = rhs.output_prefix;
    /////////////////////////////////////////////////////////////

    /////////////////////////////////////
    if (write_all)
    {
        // slidDir = rhs.slidDir;
        // rollDir = rhs.rollDir;
        // inout = rhs.inout;
        // distB3 = rhs.distB3;
        // slidB3 = rhs.slidB3;
        // rollB3 = rhs.rollB3;
        // slidFric = rhs.slidFric;
        // rollFric = rhs.rollFric;
        //////////
        vdwForce = rhs.vdwForce;
        // elasticForce = rhs.elasticForce;
        // slideForce = rhs.slideForce;
        // rollForce = rhs.rollForce;
        // torqueForce = rhs.torqueForce;
    }
    /////////////////////////////////////

    return *this;
}

Ball_group::Ball_group(const Ball_group& rhs)
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
    KE = rhs.KE;

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

    radiiDistribution = rhs.radiiDistribution;
    radiiFraction = rhs.radiiFraction;


    /////////////////////////////////////////////////////////////
    // dynamicTime = rhs.dynamicTime;

    // G = rhs.G;  // Gravitational constant
    // density = rhs.density;
    // u_s = rhs.u_s;       // Coeff of sliding friction
    // u_r = rhs.u_r;               // Coeff of rolling friction
    // sigma = rhs.sigma;              // Poisson ratio for rolling friction.
    // Y = rhs.Y;               // Young's modulus in erg/cm3
    // cor = rhs.cor;                // Coeff of restitution
    // simTimeSeconds = rhs.simTimeSeconds;  // Seconds
    // timeResolution = rhs.timeResolution;    // Seconds - This is duration between exported steps.
    // fourThirdsPiRho = rhs.fourThirdsPiRho;  // for fraction of smallest sphere radius.
    // scaleBalls = rhs.scaleBalls;                         // base radius of ball.
    // maxOverlap = rhs.maxOverlap;                           // of scaleBalls
    // KEfactor = rhs.KEfactor;                              // Determines collision velocity based on KE/PE
    // v_custom = rhs.v_custom;  // Velocity cm/s
    // temp = rhs.temp;          //tempurature of simulation in Kelvin
    // kConsts = rhs.kConsts;
    // impactParameter = rhs.impactParameter;  // Impact angle radians
    // Ha = rhs.Ha;         // Hamaker constant for vdw force
    // h_min = rhs.h_min;  // 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
    // cone = rhs.cone;  // Cone of particles ignored moving away from center of mass. Larger angle ignores more.

    // // Simulation Structure
    // properties = rhs.properties;  // Number of columns in simData file per ball
    // genBalls = rhs.genBalls;
    // attempts = rhs.attempts;  // How many times to try moving every ball touching another in generator.

    // skip = rhs.skip;  // Steps thrown away before recording a step to the buffer. 500*.04 is every 20 seconds in sim.
    // steps = rhs.steps;

    // dt = rhs.dt;
    // kin = rhs.kin;  // Spring constant
    // kout = rhs.kout;
    // spaceRange = rhs.spaceRange;  // Rough minimum space required
    // spaceRangeIncrement = rhs.spaceRangeIncrement;
    // z0Rot = rhs.z0Rot;  // Cluster one z axis rotation
    // y0Rot = rhs.y0Rot;  // Cluster one y axis rotation
    // z1Rot = rhs.z1Rot;  // Cluster two z axis rotation
    // y1Rot = rhs.y1Rot;  // Cluster two y axis rotation
    // simTimeElapsed = rhs.simTimeElapsed;

    // // File from which to proceed with further simulations
    // project_path = rhs.project_path;
    // output_folder = rhs.output_folder;
    // // projectileName;
    // // targetName;
    // output_prefix = rhs.output_prefix;
    /////////////////////////////////////////////////////////////


    /////////////////////////////////////
    // slidDir = rhs.slidDir;
    // rollDir = rhs.rollDir;
    // inout = rhs.inout;
    // distB3 = rhs.distB3;
    // slidB3 = rhs.slidB3;
    // rollB3 = rhs.rollB3;
    // slidFric = rhs.slidFric;
    // rollFric = rhs.rollFric;
    //////////
    vdwForce = rhs.vdwForce;
    // elasticForce = rhs.elasticForce;
    // slideForce = rhs.slideForce;
    // rollForce = rhs.rollForce;
    // torqueForce = rhs.torqueForce;
    /////////////////////////////////////
}



void Ball_group::parse_input_file(char const* location)
{
    std::string s_location(location);
    std::string json_file = s_location + "input.json";
    std::ifstream ifs(json_file);
    // std::cerr<<json_file<<std::endl;
    //// CANNOT USE json::parse() IF YOU RDBUF TOO
    // std::cerr<<ifs.rdbuf()<<std::endl;
    json inputs = json::parse(ifs);

    
    if (world_rank == 0)
    {
        if (inputs["seed"] == std::string("default"))
        {
            seed = static_cast<int>(time(nullptr));
        }
        else
        {
            seed = static_cast<int>(inputs["seed"]);
        }

    }
    MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

    srand(seed);
    random_generator.seed(seed);
    // std::cout<<"SEED: "<<seed<<" for rank: "<<world_rank<<std::endl;

    OMPthreads = inputs["OMPthreads"];
    total_balls_to_add = inputs["N"];
    std::string temp_typeSim = inputs["simType"];
    if (temp_typeSim == "BPCA")
    {
        typeSim = BPCA;
    }
    else
    {
        typeSim = collider;
    }
    std::string temp_distribution = inputs["radiiDistribution"];
    if (temp_distribution == "logNormal")
    {
        radiiDistribution = logNorm;
    }
    else
    {
        radiiDistribution = constant;
    }
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
    out_folder = inputs["output_folder"];
    if (output_folder == std::string("default"))
    {
        output_folder = s_location;
        out_folder = s_location;
    }
    projectileName = inputs["projectileName"];
    targetName = inputs["targetName"];
    output_prefix = inputs["output_prefix"];
    if (output_prefix == std::string("default"))
    {
        output_prefix = "";
    }

    radiiFraction = inputs["radiiFraction"];

    output_width = num_particles;
}

// @brief calculates the vdw force
inline double Ball_group::calc_VDW_force_mag(const double Ra,const double Rb,const double h)
{
    const double h2 = h * h;
    // constexpr double h2 = h * h;
    const double twoRah = 2 * Ra * h;
    const double twoRbh = 2 * Rb * h;
    return Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
             ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                               (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                               (h2 + twoRah + twoRbh + 4 * Ra * Rb)));

}

// @breif calculates the mu scaling factor for all pairs of particle sizes
void Ball_group::calc_mu_scale_factor()
{
    int e;
    for (int A = 1; A < num_particles; ++A)
    {
        for (int B = 0; B < A; ++B)
        {   
            e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
            u_scale[e] = calc_VDW_force_mag(R[A],R[B],h_min_physical)/
                                calc_VDW_force_mag(R[A],R[B],h_min);  
        }
    }
}

////////////////////////////////////
void Ball_group::zeroSaveVals()
{
    int size = num_particles*num_particles;
    for (int i = 0; i < size; ++i)
    {
        vdwForce[i] = {0,0,0};
        // elasticForce[i] = {0,0,0};
        // slideForce[i] = {0,0,0};
        // rollForce[i] = {0,0,0};
        // torqueForce[i] = {0,0,0};
    }
    // for (int i = 0; i < num_particles; ++i)
    // {
    //     // if (i < num_particles)
    //     // {
    //     // distB3[i] = 0.0;
    //     // }
    //     // slidDir[i] = {0,0,0};
    //     // rollDir[i] = {0,0,0};
    //     // inout[i] = 0.0;
    //     // slidB3[i] = {0,0,0};
    //     // rollB3[i] = {0,0,0};
    //     // // slidFric[i] = {0,0,0};
    //     // rollFric[i] = {0,0,0};
    // }
}
////////////////////////////////////

void Ball_group::calibrate_dt(int const Step, const double& customSpeed = -1.)
{

    const double dtOld = dt;

    if (customSpeed > 0.) {
        updateDTK(customSpeed);
        if (world_rank == 0)
        {
            std::cerr << "CUSTOM SPEED: " << customSpeed;
        }
    } else {
        // std::cerr << vCollapse << " <- vCollapse | Lazz Calc -> " << M_PI * M_PI * G * pow(density, 4.
        // / 3.) * pow(mTotal, 2. / 3.) * rMax;

        v_max = getVelMax();


        // Take whichever velocity is greatest:
        if (world_rank == 0)
        {
            std::cerr << '\n';
            std::cerr << v_collapse << " = vCollapse | vMax = " << v_max;
        }
        if (v_max < v_collapse) { v_max = v_collapse; }

        if (v_max < v_max_prev) {
            updateDTK(v_max);
            v_max_prev = v_max;
            if (world_rank == 0)
            {
                std::cerr << "\nk: " << kin << "\tdt: " << dt;
            }
        }
    }

    if (Step == 0 or dtOld < 0) {
        steps = static_cast<int>(simTimeSeconds / dt);
        if (world_rank == 0)
        {
            std::cerr << "\tInitial Steps: " << steps << '\n';
        }
    } else {
        steps = static_cast<int>(dtOld / dt) * (steps - Step) + Step;
        if(world_rank == 0)
        {
            std::cerr << "\tSteps: " << steps;
        }
    }

    if (timeResolution / dt > 1.) {
        skip = static_cast<int>(floor(timeResolution / dt));
        if (world_rank == 0)
        {
            std::cerr << "\tSkip: " << skip << '\n';
        }
    } else {
        if (world_rank == 0)
        {
            std::cerr << "Desired time resolution is lower than dt. Setting to 1 second per skip.\n";
        }
        skip = static_cast<int>(floor(1. / dt));
    }
}

// todo - make bigger balls favor the middle, or, smaller balls favor the outside.
/// @brief Push balls apart until no overlaps
void Ball_group::pushApart() const
{

    if (world_rank == 0)
    {
        std::cerr << "Separating spheres - Current max overlap:\n";
    }
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
        //  simDataWrite("pushApart_");
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
            if (world_rank == 0)
            {
                std::cerr << overlapMax << "                        \r";
            }
        } else {
            if (world_rank == 0)
            {
                std::cerr << "\nSuccess!\n";
            }
            break;
        }
        overlapMax = -1;
        step++;
    }
    delete[] counter;
}

void Ball_group::calc_v_collapse()
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
[[nodiscard]] double Ball_group::getVelMax()
{
    // int world_rank = getRank();

    v_max = 0;

    // todo - make this a manual set true or false to use soc so we know if it is being used or not.
    if (soc > 0) {
        int counter = 0;
        for (int Ball = 0; Ball < num_particles; Ball++) {
            if (vel[Ball].norm() > v_max) 
            { 
                v_max = vel[Ball].norm();
            }
            /////////////////SECTION COMMENTED FOR ACCURACY TESTS
            // Only consider balls moving toward com and within 4x initial radius around it.
            // const vec3 fromCOM = pos[Ball] - getCOM();
            // if (acos(vel[Ball].normalized().dot(fromCOM.normalized())) > cone && fromCOM.norm() < soc) {
            //     if (vel[Ball].norm() > v_max) { v_max = vel[Ball].norm(); }
            // } else {
            //     counter++;
            // }
        }
        if (world_rank == 0)
        {
            std::cerr << '(' << counter << " spheres ignored"
                      << ") ";
        }
    } else {
        for (int Ball = 0; Ball < num_particles; Ball++) {

            if (vel[Ball].norm() > v_max) 
            { 
                v_max = vel[Ball].norm();
            }
        }

        // Is vMax for some reason unreasonably small? Don't proceed. Probably a finished sim.
        // This shouldn't apply to extremely destructive collisions because it is possible that no
        // particles are considered, so it will keep pausing.
        if (v_max < 1e-10) {
            if (world_rank == 0)
            {
                std::cerr << "\nMax velocity in system is less than 1e-10.\n";
            }
            system("pause");
        }
    }

    return v_max;
}

double Ball_group::get_soc()
{
    return soc;
}

void Ball_group::calc_helpfuls()
{
    r_min = getRmin();
    r_max = getRmax();
    m_total = getMass();
    initial_radius = get_radius(getCOM());
    soc = 4 * r_max + initial_radius;
    // soc = -1;
}   

// Kick ballGroup (give the whole thing a velocity)
void Ball_group::kick(const vec3& vec) const
{
    for (int Ball = 0; Ball < num_particles; Ball++) { vel[Ball] += vec; }
}


vec3 Ball_group::calc_momentum(const std::string& of = "") const
{
    vec3 pTotal = {0, 0, 0};
    for (int Ball = 0; Ball < num_particles; Ball++) { pTotal += m[Ball] * vel[Ball]; }
    // fprintf(stderr, "%s Momentum Check: %.2e, %.2e, %.2e\n", of.c_str(), pTotal.x, pTotal.y, pTotal.z);
    return pTotal;
}

// offset cluster
void Ball_group::offset(const double& rad1, const double& rad2, const double& impactParam) const
{
    for (int Ball = 0; Ball < num_particles; Ball++) {
        pos[Ball].x += (rad1 + rad2) * cos(impactParam);
        pos[Ball].y += (rad1 + rad2) * sin(impactParam);
    }
}

/// Approximate the radius of the ballGroup.
[[nodiscard]] double Ball_group::get_radius(const vec3& center) const
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
void Ball_group::updateGPE()
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

//No need to wrap writes in if (world_rank == 0) because only world_rank 0 calls function
void Ball_group::sim_init_write(std::string filename, int counter = 0)
{
    // Create string for file name identifying spin combination negative is 2, positive is 1 on each axis.
    // std::string spinCombo = "";
    // for ( int i = 0; i < 3; i++)
    //{
    //  if (spins[i] < 0) { spinCombo += "2"; }
    //  else if (spins[i] > 0) { spinCombo += "1"; }
    //  else { spinCombo += "0"; }
    //}


    // todo - filename is now a copy and this works. Need to consider how old way worked for
    // compatibility. What happens without setting output_prefix = filename? Check if file name already
    // exists.
    std::ifstream checkForFile;
    std::string filenum;
    if (counter != 0)
    {
        filenum = std::to_string(counter) + '_';
    }
    else
    {
        filenum = "";
    }
    checkForFile.open(output_folder + filenum + filename + "simData.csv", std::ifstream::in);
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
    std::string timeFileName = output_folder + "time.csv";

    std::cerr << "New file tag: " << filename;

    // Open all file streams:
    std::ofstream energyWrite, ballWrite, constWrite,timeWrite;
    energyWrite.open(energyFilename, std::ofstream::app);
    ballWrite.open(simDataFilename, std::ofstream::app);
    constWrite.open(constantsFilename, std::ofstream::app);
    timeWrite.open(timeFileName, std::ofstream::app);

    // Make column headers:
    energyWrite << "Time,PE,KE,E,p,L";
    ballWrite << "x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0";
    timeWrite << "ball,update time\n";
    // std::cout << "x0,y0,z0,wx0,wy0,wz0,wmag0,vx0,vy0,vz0,bound0";
    // std::cout<<simDataFilename<<std::endl;
    // std::cout<<num_particles<<std::endl;

    for (int Ball = 1; Ball < num_particles;
         Ball++)  // Start at 2nd ball because first one was just written^.
    {
        // std::cout<<Ball<<','<<num_particles<<std::endl;
        std::string thisBall = std::to_string(Ball);
        ballWrite << ",x" + thisBall << ",y" + thisBall << ",z" + thisBall << ",wx" + thisBall
                  << ",wy" + thisBall << ",wz" + thisBall << ",wmag" + thisBall << ",vx" + thisBall
                  << ",vy" + thisBall << ",vz" + thisBall << ",bound" + thisBall;
        // std::cout << ",x" + thisBall << ",y" + thisBall << ",z" + thisBall << ",wx" + thisBall
        //           << ",wy" + thisBall << ",wz" + thisBall << ",wmag" + thisBall << ",vx" + thisBall
        //           << ",vy" + thisBall << ",vz" + thisBall << ",bound" + thisBall;

    }

    // Write constant data:
    for (int Ball = 0; Ball < num_particles; Ball++) {
        constWrite << R[Ball] << ',' << m[Ball] << ',' << moi[Ball] << '\n';
    }

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

[[nodiscard]] vec3 Ball_group::getCOM() const
{
    // int world_rank = getRank();

    if (m_total > 0) {
        vec3 comNumerator;
        for (int Ball = 0; Ball < num_particles; Ball++) { comNumerator += m[Ball] * pos[Ball]; }
        vec3 com = comNumerator / m_total;
        return com;
    } else {
        if (world_rank == 0)
        {
            std::cerr << "Mass of cluster is zero.\n";
        }
        exit(EXIT_FAILURE);
    }
}

void Ball_group::zeroVel() const
{
    for (int Ball = 0; Ball < num_particles; Ball++) { vel[Ball] = {0, 0, 0}; }
}

void Ball_group::zeroAngVel() const
{
    for (int Ball = 0; Ball < num_particles; Ball++) { w[Ball] = {0, 0, 0}; }
}

void Ball_group::to_origin() const
{
    const vec3 com = getCOM();

    for (int Ball = 0; Ball < num_particles; Ball++) { pos[Ball] -= com; }
}

// Set velocity of all balls such that the cluster spins:
void Ball_group::comSpinner(const double& spinX, const double& spinY, const double& spinZ) const
{
    const vec3 comRot = {spinX, spinY, spinZ};  // Rotation axis and magnitude
    for (int Ball = 0; Ball < num_particles; Ball++) {
        vel[Ball] += comRot.cross(pos[Ball] - getCOM());
        w[Ball] += comRot;
    }
}

void Ball_group::rotAll(const char axis, const double angle) const
{
    for (int Ball = 0; Ball < num_particles; Ball++) {
        pos[Ball] = pos[Ball].rot(axis, angle);
        vel[Ball] = vel[Ball].rot(axis, angle);
        w[Ball] = w[Ball].rot(axis, angle);
    }
}

double Ball_group::calc_mass(const double& radius, const double& density)
{
    return density * 4. / 3. * 3.14159 * std::pow(radius, 3);
}

double Ball_group::calc_moi(const double& radius, const double& mass) { return .4 * mass * radius * radius; }

//Not used anywhere at the moment
// Ball_group Ball_group::spawn_particles(const int count)
// {
//     // Load file data:
//     std::cerr << "Add Particle\n";

//     // Random particle to origin
//     Ball_group projectile(count);
//     // Particle random position at twice radius of target:
//     // We want the farthest from origin since we are offsetting form origin. Not com.
//     const auto cluster_radius = 3;

//     const vec3 projectile_direction = rand_vec3(1).normalized();
//     projectile.pos[0] = projectile_direction * (cluster_radius + scaleBalls * 4);
//     projectile.w[0] = {0, 0, 0};
//     // Velocity toward origin:
//     projectile.vel[0] = -v_custom * projectile_direction;
//     projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
//     projectile.m[0] = density * 4. / 3. * pi * std::pow(projectile.R[0], 3);
//     projectile.moi[0] = calc_moi(projectile.R[0], projectile.m[0]);

//     const vec3 local_coords[3] = local_coordinates(projectile_direction);
//     // to_vec3(local_coords.y).print();
//     // to_vec3(local_coords.z).print();
//     // projectile.pos[0].print();
//     for (int i = 1; i < projectile.num_particles - 3; i++) {
//         const auto rand_y = rand_between(-cluster_radius, cluster_radius);
//         const auto rand_z = rand_between(-cluster_radius, cluster_radius);
//         // projectile.pos[i] = projectile.pos[0] + perpendicular_shift(local_coords, rand_y, rand_z);
//         projectile.pos[i] = projectile.pos[0] + perpendicular_shift(local_coords, rand_y, rand_z);
//         // std::cout << rand_y << '\t' << to_vec3(local_coords.y * rand_y) <<'\t'<< projectile.pos[i] <<
//         // '\n';
//     }
//     projectile.pos[projectile.num_particles - 3] = projectile_direction * 2;
//     projectile.pos[projectile.num_particles - 2] = projectile_direction * 4;
//     projectile.pos[projectile.num_particles - 1] = projectile_direction * 6;

//     Ball_group new_group{projectile.num_particles + num_particles};

//     new_group.merge_ball_group(*this);
//     new_group.merge_ball_group(projectile);

//     // new_group.calibrate_dt(0, 1);
//     // new_group.init_conditions();

//     // new_group.to_origin();
//     return new_group;
// }

//@brief returns new position of particle after it is given random offset
//@param local_coords is plane perpendicular to direction of projectile
//@param projectile_pos is projectile's position before offset is applied
//@param projectile_vel is projectile's velocity
//@param projectile_rad is projectile's radius
vec3 Ball_group::dust_agglomeration_offset(
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
Ball_group Ball_group::dust_agglomeration_particle_init()
{
    // Random particle to origin
    Ball_group projectile(1);
    projectile.radiiDistribution = radiiDistribution;
    projectile.radiiFraction = radiiFraction;
    // Particle random position at twice radius of target:
    // We want the farthest from origin since we are offsetting form origin. Not com.
    const auto cluster_radius = get_radius(vec3(0, 0, 0));

    const vec3 projectile_direction = rand_unit_vec3();
    projectile.pos[0] = projectile_direction * (cluster_radius + scaleBalls * 4);
    if (radiiDistribution == constant)
    {
        // std::cout<<"radiiFraction: "<<radiiFraction<<std::endl;
        projectile.R[0] = scaleBalls;  //MAKE BOTH VERSIONS SAME
        // projectile.R[0] = scaleBalls/radiiFraction;  //limit of 1.4// rand_between(1,3)*1e-5;
        // std::cout<<"(constant) Particle added with radius of "<<projectile.R[0]<<std::endl;
    }
    else
    {
        projectile.R[0] = lognorm_dist(scaleBalls*std::exp(-5*std::pow(lnSigma,2)/2),lnSigma);
        // std::cout<<"(lognorm) Particle added with radius of "<<projectile.R[0]<<std::endl;
    }
    projectile.w[0] = {0, 0, 0};
    projectile.m[0] = density * 4. / 3. * pi * std::pow(projectile.R[0], 3);
    // Velocity toward origin:
    if (temp > 0)
    {
        double a = std::sqrt(Kb*temp/projectile.m[0]);
        v_custom = max_bolt_dist(a); 
        std::cerr<<"v_custom set to "<<v_custom<< "cm/s based on a temp of "
                <<temp<<" degrees K."<<std::endl; 
    }
    projectile.vel[0] = -v_custom * projectile_direction;

    
    // projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
    projectile.moi[0] = calc_moi(projectile.R[0], projectile.m[0]);

    //////////////////////////////////
    if (write_all)
    {
        // projectile.slidDir[0] = {0,0,0};
        // projectile.rollDir[0] = {0,0,0};
        // projectile.inout[0] = 0.0;
        // projectile.distB3[0] = 0.0;
        // projectile.slidB3[0] = {0,0,0};
        // projectile.rollB3[0] = {0,0,0};
        // projectile.slidFric[0] = {0,0,0};
        // projectile.rollFric[0] = {0,0,0};
        //////////
        projectile.vdwForce[0] = {0,0,0};
        // projectile.elasticForce[0][0] = {0,0,0};
        // projectile.slideForce[0] = {0,0,0};
        // projectile.rollForce[0] = {0,0,0};
        // projectile.torqueForce[0] = {0,0,0};
    }
    //////////////////////////////////

    const double3x3 local_coords = local_coordinates(to_double3(projectile_direction));
    
    projectile.pos[0] = dust_agglomeration_offset(local_coords,projectile.pos[0],projectile.vel[0],projectile.R[0]);
    // std::cerr<<"pos, dir: "<<projectile.pos[0]<<", "<<projectile_direction<<std::endl;
    //////////////////////////////////
    //TURN ON above LINE AND OFF REST FOR REAL SIM
    // if (num_particles == 3)
    // {
    //     projectile.vel[0] = {0,0,0};
    //     projectile.pos[0] = {0.5e-5,0.5e-5,1.2e-5};

    // }
    // else
    // {
    //     projectile.pos[0] = dust_agglomeration_offset(local_coords,projectile.pos[0],projectile.vel[0],projectile.R[0]);
    // }
    //////////////////////////////////

    
    return projectile;
}


// Uses previous O as target and adds one particle to hit it:
Ball_group Ball_group::add_projectile()
{
    // int world_rank = getRank();

    // Load file data:
    if (world_rank == 0)
    {
        std::cerr << "Add Particle\n";
    }

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

    if (world_rank == 0)
    {
        fprintf(
            stderr,
            "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n",
            vel[0].norm(),
            projectile.vel[0].norm());

        std::cerr << '\n';
    }
    projectile.calc_momentum("Projectile");
    calc_momentum("Target");
    Ball_group new_group{projectile.num_particles + num_particles};

    // std::cerr<<"EXISTING dt pre set: "<<dt<<std::endl;
    // std::cerr<<"NEW dt pre set: "<<new_group.dt<<std::endl;
    // std::cerr<<"NEW dt after set: "<<new_group.dt<<std::endl;

    new_group.merge_ball_group(*this);
    new_group.merge_ball_group(projectile);
    // std::cout<<"radiiDistribution in add_projectile(4): "<<new_group.radiiDistribution<<std::endl;

    // Hack - if v_custom is less than 1 there are problems if dt is calibrated to this
    //        if v_custom is greater than 1 you need to calibrate dt to that v_custom
    if (v_custom < 1)
    {
        new_group.calibrate_dt(0, 1);
    }
    else
    {
        new_group.calibrate_dt(0, v_custom);
    }
    new_group.init_conditions();

    new_group.to_origin();
   
    return new_group;
}

/// @brief Add another ballGroup into this one.
/// @param src The ballGroup to be added.
void Ball_group::merge_ball_group(const Ball_group& src)
{
    // Copy incoming data to the end of the currently loaded data.
    // std::memcpy(
    //     &distances[num_particles_added], src.distances, sizeof(src.distances[0]) * (src.num_particles*src.num_particles-src.num_particles)/2);
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
    //////////////////////////////////////
    if (write_all)
    {
        // std::memcpy(&distB3[num_particles_added], src.distB3, sizeof(src.distB3[0]) * src.num_particles);
        // std::memcpy(&inout[num_particles_added], src.inout, sizeof(src.inout[0]) * src.num_particles);
        // std::memcpy(&slidDir[num_particles_added], src.slidDir, sizeof(src.slidDir[0]) * src.num_particles);
        // std::memcpy(&rollDir[num_particles_added], src.rollDir, sizeof(src.rollDir[0]) * src.num_particles);
        // std::memcpy(&slidB3[num_particles_added], src.slidB3, sizeof(src.slidB3[0]) * src.num_particles);
        // std::memcpy(&rollB3[num_particles_added], src.rollB3, sizeof(src.rollB3[0]) * src.num_particles);
        // std::memcpy(&slidFric[num_particles_added], src.slidFric, sizeof(src.slidFric[0]) * src.num_particles);
        // std::memcpy(&rollFric[num_particles_added], src.rollFric, sizeof(src.rollFric[0]) * src.num_particles);
        //////////
        std::memcpy(&vdwForce[num_particles_added], src.vdwForce, sizeof(src.vdwForce[0]) * src.num_particles);
        // std::memcpy(&elasticForce[num_particles_added], src.elasticForce, sizeof(src.elasticForce[0]) * src.num_particles);
        // std::memcpy(&slideForce[num_particles_added], src.slideForce, sizeof(src.slideForce[0]) * src.num_particles);
        // std::memcpy(&rollForce[num_particles_added], src.rollForce, sizeof(src.rollForce[0]) * src.num_particles);
        // std::memcpy(&torqueForce[num_particles_added], src.torqueForce, sizeof(src.torqueForce[0]) * src.num_particles);
    }
    ////////////////////////////////////////

    // Keep track of now loaded ball set to start next set after it:
    num_particles_added += src.num_particles;
    radiiDistribution = src.radiiDistribution;
    radiiFraction = src.radiiFraction;

    //////////////////////////////////////////////////////
    // dynamicTime = src.dynamicTime;

    // G = src.G;  // Gravitational constant
    // density = src.density;
    // u_s = src.u_s;       // Coeff of sliding friction
    // u_r = src.u_r;               // Coeff of rolling friction
    // sigma = src.sigma;              // Poisson ratio for rolling friction.
    // Y = src.Y;               // Young's modulus in erg/cm3
    // cor = src.cor;                // Coeff of restitution
    // simTimeSeconds = src.simTimeSeconds;  // Seconds
    // timeResolution = src.timeResolution;    // Seconds - This is duration between exported steps.
    // fourThirdsPiRho = src.fourThirdsPiRho;  // for fraction of smallest sphere radius.
    // scaleBalls = src.scaleBalls;                         // base radius of ball.
    // maxOverlap = src.maxOverlap;                           // of scaleBalls
    // KEfactor = src.KEfactor;                              // Determines collision velocity based on KE/PE
    // v_custom = src.v_custom;  // Velocity cm/s
    // temp = src.temp;          //tempurature of simulation in Kelvin
    // kConsts = src.kConsts;
    // impactParameter = src.impactParameter;  // Impact angle radians
    // Ha = src.Ha;         // Hamaker constant for vdw force
    // h_min = src.h_min;  // 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
    // cone = src.cone;  // Cone of particles ignored moving away from center of mass. Larger angle ignores more.

    // // Simulation Structure
    // properties = src.properties;  // Number of columns in simData file per ball
    // genBalls = src.genBalls;
    // attempts = src.attempts;  // How many times to try moving every ball touching another in generator.

    // skip = src.skip;  // Steps thrown away before recording a step to the buffer. 500*.04 is every 20 seconds in sim.
    // steps = src.steps;

    // dt = src.dt;
    // kin = src.kin;  // Spring constant
    // kout = src.kout;
    // spaceRange = src.spaceRange;  // Rough minimum space required
    // spaceRangeIncrement = src.spaceRangeIncrement;
    // z0Rot = src.z0Rot;  // Cluster one z axis rotation
    // y0Rot = src.y0Rot;  // Cluster one y axis rotation
    // z1Rot = src.z1Rot;  // Cluster two z axis rotation
    // y1Rot = src.y1Rot;  // Cluster two y axis rotation
    // simTimeElapsed = src.simTimeElapsed;

    // // File from which to proceed with further simulations
    // project_path = src.project_path;
    // output_folder = src.output_folder;
    // // projectileName;
    // // targetName;
    // output_prefix = src.output_prefix;
    //////////////////////////////////////////////////////
    

    calc_helpfuls();
}

/// Allocate balls
void Ball_group::allocate_group(const int nBalls)
{
    num_particles = nBalls;
    // int world_rank = getRank();
    try {
        distances = new double[(num_particles * num_particles / 2) - (num_particles / 2)];
        
        accsq = new vec3[num_particles*num_particles];
        aaccsq = new vec3[num_particles*num_particles];

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

        // /////////////////////////
        if (write_all)
        {
            // inout = new double[num_particles-1];
            // distB3 = new double[num_particles-1];
            // slidDir = new vec3[num_particles];
            // rollDir = new vec3[num_particles];
            // slidB3 = new vec3[num_particles];
            // rollB3 = new vec3[num_particles];
            // slidFric = new vec3[num_particles];
            // rollFric = new vec3[num_particles];
            //////////
            vdwForce = new vec3[num_particles*num_particles];
            // elasticForce = new vec3[num_particles*num_particles];
            // slideForce = new vec3[num_particles*num_particles];
            // rollForce = new vec3[num_particles*num_particles];
            // torqueForce = new vec3[num_particles*num_particles];
        }
        // /////////////////////////
    } catch (const std::exception& e) {
        std::cerr <<"Rank "<<world_rank<< " failed trying to allocate group. " << e.what() << '\n';
    }
}


/// @brief Deallocate arrays to recover memory.
void Ball_group::freeMemory() const
{
    delete[] distances;
    delete[] pos;
    delete[] vel;
    delete[] velh;
    delete[] acc;
    delete[] accsq;
    delete[] w;
    delete[] wh;
    delete[] aacc;
    delete[] aaccsq;
    delete[] R;
    delete[] m;
    delete[] moi;
    /////////////////////
    if (write_all)
    {
        // delete[] slidDir;
        // delete[] rollDir;
        // delete[] inout;
        // delete[] distB3;
        // delete[] slidB3;
        // delete[] rollB3;
        // delete[] slidFric;
        // delete[] rollFric;
        //////////
        delete[] vdwForce;
        // delete[] elasticForce;
        // delete[] slideForce;
        // delete[] rollForce;
        // delete[] torqueForce;
    }
    /////////////////////
}


// Initialize accelerations and energy calculations:
void Ball_group::init_conditions()
{
    // int world_rank = getRank();

    update_time = 0.0;
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
            // double oldDist = distances[e];

            // Check for collision between Ball and otherBall.
            if (overlap > 0) {
                double k;
                k = kin;
                // Apply coefficient of restitution to balls leaving collision.
                // if (dist >= oldDist) {
                //     k = kout;
                // } else {
                //     k = kin;
                // }

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

    #ifdef MPI_ENABLE
        double local_PE = PE;
        PE = 0.0;
        MPI_Reduce(&local_PE,&PE,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    #endif

    // Calc energy:
    if (world_rank == 0)
    {
        for (int Ball = 0; Ball < num_particles; Ball++) {
            KE += .5 * m[Ball] * vel[Ball].dot(vel[Ball]) + .5 * moi[Ball] * w[Ball].dot(w[Ball]);
            mom += m[Ball] * vel[Ball];
            ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
        }
    }
}

[[nodiscard]] double Ball_group::getRmin()
{
    r_min = R[0];
    for (int Ball = 1; Ball < num_particles; Ball++) {
        if (R[Ball] < r_min) { r_min = R[Ball]; }
    }
    return r_min;
}

[[nodiscard]] double Ball_group::getRmax()
{
    r_max = R[0];
    for (int Ball = 0; Ball < num_particles; Ball++) {
        if (R[Ball] > r_max) { r_max = R[Ball]; }
    }
    return r_max;
}


[[nodiscard]] double Ball_group::getMassMax() const
{
    double mMax = m[0];
    for (int Ball = 0; Ball < num_particles; Ball++) {
        if (m[Ball] > mMax) { mMax = m[Ball]; }
    }
    return mMax;
}


void Ball_group::parseSimData(std::string line)
{
    std::string lineElement;
    // Get number of balls in file
    int count = 54 / properties;
    num_particles = static_cast<int>((static_cast<int>(std::count(line.begin(), line.end(), ',')) + 1)/11);
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
void Ball_group::loadConsts(const std::string& path, const std::string& filename)
{

    // int world_rank = getRank();
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
        if (world_rank == 0)
        {
            std::cerr << "Could not open constants file: " << constantsFilename << "... Existing program."
                      << '\n';
        }
        exit(EXIT_FAILURE);
    }
}


//This used to be [[nodiscard]] static std::string ... but wont compile outside the actual class definition
/// Get last line of previous simData by filename.
[[nodiscard]] std::string Ball_group::getLastLine(const std::string& path, const std::string& filename)
{
    // int world_rank = getRank();
    std::string simDataFilepath = path + filename + "simData.csv";
    if (auto simDataStream = std::ifstream(simDataFilepath, std::ifstream::in)) {
        if (world_rank == 0)
        {
            std::cerr << "\nParsing last line of data.\n";
        }

        simDataStream.seekg(-1, std::ios_base::end);  // go to 
         // spot before the EOF

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
        if (world_rank == 0)
        {
            std::cerr << "Could not open simData file: " << simDataFilepath << "... Existing program."
                      << '\n';
        }
        exit(EXIT_FAILURE);
    }
}

void Ball_group::simDataWrite(std::string& outFilename)
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


[[nodiscard]] double Ball_group::getMass()
{
    m_total = 0;
    {
        for (int Ball = 0; Ball < num_particles; Ball++) { m_total += m[Ball]; }
    }
    return m_total;
}

void Ball_group::threeSizeSphere(const int nBalls)
{
    // int world_rank = getRank();

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
        if (world_rank == 0)
        {
            std::cerr<<"Ball "<<i<<"\tmass is "<<m[i]<<"\t"<<"radius is "<<R[i]<<std::endl;
        }
    }

    if (world_rank == 0)
    {
        std::cerr << "Smalls: " << smalls << " Mediums: " << mediums << " Larges: " << larges << '\n';
    }

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
            if (world_rank == 0)
            {
                std::cerr << "Collisions: " << collisionDetected << "                        \r";
            }
        }
        if (collisionDetected == 0) {
            if (world_rank == 0)
            {
                std::cerr << "\nSuccess!\n";
            }
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
            if (world_rank == 0)
            {
                std::cerr << "Failed " << spaceRange << ". Increasing range " << spaceRangeIncrement
                          << "cm^3.\n";
            }
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

    if (world_rank == 0)
    {
        std::cerr << "Final spacerange: " << spaceRange << '\n';
        std::cerr << "m_total: " << m_total << '\n';
        std::cerr << "Initial Radius: " << get_radius(getCOM()) << '\n';
        std::cerr << "Mass: " << getMass() << '\n';
    }
}

void Ball_group::generate_ball_field(const int nBalls)
{
    // int world_rank = getRank();
    if (world_rank == 0)
    {
        std::cerr << "CLUSTER FORMATION\n";
    }

    allocate_group(nBalls);

    // Create new random number set.
        //This should be d
         // in parse_input_file
    // const int seedSave = static_cast<int>(time(nullptr));
    // srand(seedSave);
    if (radiiDistribution == constant)
    {
        oneSizeSphere(nBalls);
    }
    else
    {
        distSizeSphere(nBalls);
    }
    
    calc_helpfuls();
    // threeSizeSphere(nBalls);

    output_prefix = std::to_string(nBalls) + "_R" + scientific(get_radius(getCOM())) + "_v" +
                    scientific(v_custom) + "_cor" + rounder(sqrtf(cor), 4) + "_mu" + rounder(u_s, 3) +
                    "_rho" + rounder(density, 4);
}

/// Make ballGroup from file data.
void Ball_group::loadSim(const std::string& path, const std::string& filename)
{
    // int world_rank = getRank();
    parseSimData(getLastLine(path, filename));
    loadConsts(path, filename);

    calc_helpfuls();

    if (world_rank == 0)
    {        
        std::cerr << "Balls: " << num_particles << '\n';
        std::cerr << "Mass: " << m_total << '\n';
        std::cerr << "Approximate radius: " << initial_radius << " cm.\n";
    }
}

void Ball_group::distSizeSphere(const int nBalls)
{
    for (int Ball = 0; Ball < nBalls; Ball++) {
        R[Ball] = lognorm_dist(scaleBalls*std::exp(-5*std::pow(lnSigma,2)/2),lnSigma);
        m[Ball] = density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(spaceRange);
    }

    m_total = getMass();

    placeBalls(nBalls);
}

void Ball_group::oneSizeSphere(const int nBalls)
{
    for (int Ball = 0; Ball < nBalls; Ball++) {
        R[Ball] = scaleBalls;
        m[Ball] = density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(spaceRange);
        ////////////////////////////
        // if (Ball < nBalls-1)
        // {
        //     inout[Ball] = 0.0;
        //     distB3[Ball] = 0.0;
        // }
        // slidDir[Ball] = {0,0,0};
        // rollDir[Ball] = {0,0,0};
        // slidB3[Ball] = {0,0,0};
        // rollB3[Ball] = {0,0,0};
        // slidFric[Ball] = {0,0,0};
        // rollFric[Ball] = {0,0,0};
        ////////////////////////////
    }

    m_total = getMass();

    placeBalls(nBalls);
}

void Ball_group::placeBalls(const int nBalls)
{
    // int world_rank = getRank();

    // Generate non-overlapping spherical particle field:
    int collisionDetected = 0;
    int oldCollisions = nBalls;

    if (nBalls == 1)
    {
        pos[0] = {0,1e-5,0};
    }

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
            if (world_rank == 0)
            {
                std::cerr << "Collisions: " << collisionDetected << "                        \r";
            }
        }
        if (collisionDetected == 0) {
            if (world_rank == 0)
            {
                std::cerr << "\nSuccess!\n";
            }
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
            if (world_rank == 0)
            {
                std::cerr << "Failed " << spaceRange << ". Increasing range " << spaceRangeIncrement
                          << "cm^3.\n";
            }
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
    if (world_rank == 0)
    {
        std::cerr << "Final spacerange: " << spaceRange << '\n';
        std::cerr << "Initial Radius: " << get_radius(getCOM()) << '\n';
        std::cerr << "Mass: " << m_total << '\n';
    }
}

void Ball_group::updateDTK(const double& velocity)
{
    // int world_rank = getRank();

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

    // dt = .02 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min);
    dt = .01 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min); //NORMAL ONE
    // dt = .005 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min);
    // dt = .0025 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min); //NORMAL ONE
    if (world_rank == 0)
    {
        std::cerr << "==================" << '\n';
        std::cerr << "dt set to: " << dt << '\n';
        std::cerr << "kin set to: " << kin << '\n';
        std::cerr << "kout set to: " << kout << '\n';
        std::cerr << "h_min set to: " << h_min << '\n';
        std::cerr << "Ha set to: " << Ha << '\n';
        std::cerr << "u_s set to: " << u_s << '\n';
        std::cerr << "u_r set to: " << u_r << '\n';
        if (vdw_force_max > elastic_force_max)
        {
            std::cerr << "In the vdw regime."<<std::endl;
        }
        else
        {
            std::cerr << "In the elastic regime."<<std::endl;
        }
        std::cerr << "==================" << '\n';
    }
}


void Ball_group::simInit_cond_and_center(bool add_prefix)
{
    // int world_rank = getRank();

    PE = 0.0;

    if (world_rank == 0)
    {
        std::cerr << "==================" << '\n';
        std::cerr << "dt: " << dt << '\n';
        std::cerr << "k: " << kin << '\n';
        std::cerr << "Skip: " << skip << '\n';
        std::cerr << "Steps: " << steps << '\n';
        std::cerr << "==================" << '\n';
    }

    if (num_particles > 1)
    {
        to_origin();
    }

    calc_momentum("After Zeroing");  // Is total mom zero like it should be?

    // Compute physics between all balls. Distances, collision forces, energy totals, total mass:
    init_conditions();

    // Name the file based on info above:
    if (add_prefix)
    {   
        output_prefix += "_k" + scientific(kin) + "_Ha" + scientific(Ha) + "_dt" + scientific(dt) + "_";
    }
}


void Ball_group::sim_continue(const std::string& path, const std::string& filename, int start_file_index=0)
{
    // int world_rank = getRank();
    // Load file data:
    num_particles = 3 + start_file_index;
    if (start_file_index == 0)
    {
        if (world_rank == 0)
        {
            std::cerr << "Continuing Sim...\nFile: " << filename << '\n';
        }
        loadSim(path, filename);
    }
    else
    {
        if (world_rank == 0)
        {
            std::cerr << "Continuing Sim...\nFile: " << start_file_index << '_' << filename << '\n';
        }
        loadSim(path, std::to_string(start_file_index) + "_" + filename);
    }



    std::cerr << '\n';
    calc_momentum("O");

    // Name the file based on info above:
    // output_prefix = std::to_string(num_particles) + "_rho" + rounder(density, 4);
    output_prefix = filename;
}


// Set's up a two cluster collision.
void Ball_group::sim_init_two_cluster(
    const std::string& path,
    const std::string& projectileName,
    const std::string& targetName)
{
    // int world_rank = getRank();

    // Load file data:
    if (world_rank == 0)
    {
        std::cerr << "TWO CLUSTER SIM\nFile 1: " << projectileName << '\t' << "File 2: " << targetName
                  << '\n';
    }

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

    num_particles = projectile.num_particles + target.num_particles;
    
    if (world_rank == 0)
    {
        std::cerr<<"Total number of particles in sim: "<<num_particles<<std::endl;
    }

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

    //      const double PEsys = projectile.PE + target.PE + (-G * projectile.mTotal * target.mTotal /
    //(projectile.getCOM() - target.getCOM()).norm());

    // Collision velocity calculation:
    const double mSmall = projectile.m_total;
    const double mBig = target.m_total;
    //      const double mTot = mBig + mSmall;
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

    if (world_rank == 0)
    {
        fprintf(stderr, "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n", vBig, vSmall);
    }

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

void Ball_group::bufferBarf()
{
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
}

void Ball_group::sim_looper()
{
    world_rank = getRank();
    world_size = getSize();
    ndev=acc_get_num_devices(acc_device_nvidia);
    thegpu=world_rank%ndev;
    acc_set_device_num(thegpu,acc_device_nvidia);
    std::cout<<ndev<<" devices"<<std::endl;
    std::cout<<thegpu<<" thegpu"<<std::endl;
    std::cout<<world_rank<<" world_rank"<<std::endl;

    if (world_rank == 0)
    {    
        std::cerr << "Beginning simulation...\n";

        // startProgress = ;

        std::cerr<<"Stepping through "<<steps<<" steps"<<std::endl;
    }

    bool write_step;
    time_t startProgress = time(nullptr);                // For progress reporting (gets reset)
    time_t lastWrite;                    // For write control (gets reset)
    num_pairs = (num_particles*num_particles-num_particles)*0.5;

    #pragma acc enter data copyin(this) 
    #pragma acc enter data copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) 
    #pragma acc enter data copyin(accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],acc[0:num_particles],aacc[0:num_particles],velh[0:num_particles],wh[0:num_particles]) 
    #pragma acc enter data copyin(dt,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,world_rank,world_size)
    
    for (int Step = 1; Step < steps; Step++)  // Steps start at 1 because the 0 step is initial conditions.
    {
        // simTimeElapsed += dt; //New code #1
        // Check if this is a write step:
        if (Step % skip == 0) {
            // t.start_event("writeProgressReport");
            write_step = true;

            /////////////////////// Original code #1
            simTimeElapsed += dt * skip;
            ///////////////////////

            // Progress reporting:
            if (world_rank == 0)
            {
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
            }
            startProgress = time(nullptr);
            // t.end_event("writeProgressReport");
        } else {
            write_step = debug;
            // write_step = true;
        }

        // Physics integration step:
        ///////////
        // if (write_all)
        // {
        //     zeroSaveVals();
        // }
        ///////////
        // sim_one_step(write_step,O);
        sim_one_step(write_step);
        // std::cerr<<"STEP: "<<Step<<std::endl;

        if (write_step) 
        {
            if (world_rank == 0)
            {
                // t.start_event("write_step");
                // Write energy to stream:
                ////////////////////////////////////
                //TURN THIS ON FOR REAL RUNS!!!
                energyBuffer << '\n'
                             << simTimeElapsed << ',' << PE << ',' << KE << ',' << PE + KE << ','
                             << mom.norm() << ','
                             << ang_mom.norm();  // the two zeros are bound and unbound mass


                // Data Export. Exports every 10 write_steps (10 new lines of data) and also if the last write was
                // a long time ago.
                // if (time(nullptr) - lastWrite > 1800 || Step / skip % 10 == 0) {
                if (Step / skip % 10 == 0 || Step == steps-1) {
                    // Report vMax:

                    std::cerr << "vMax = " << getVelMax() << " Steps recorded: " << Step / skip << '\n';
                    std::cerr << "Data Write to "<<output_folder<<"\n";
                    // std::cerr<<"output_prefix: "<<output_prefix<<std::endl;


                    bufferBarf();

                    lastWrite = time(nullptr);
                }  // Data export end


                if (dynamicTime) { calibrate_dt(Step, false); }
            }
            // Reinitialize energies for next step:
            KE = 0;
            PE = 0;
            mom = {0, 0, 0};
            ang_mom = {0, 0, 0};
            //unboundMass = 0;
            //boundMass = massTotal;
            ////////////////////////////////////
                // t.end_event("write_step");
        }  // write_step end

    }
    // #pragma acc exit data delete(this)
    // #pragma acc exit data delete(acc[0:num_particles],aacc[0:num_particles],PE[0:1])
    // #pragma acc exit data delete(velh[0:num_particles],wh[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs])
    // #pragma acc exit data delete(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    #pragma acc exit data delete(accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],acc[0:num_particles],aacc[0:num_particles])
    #pragma acc exit data delete(m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs])
    #pragma acc exit data delete(dt,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,world_rank,world_size)
    #pragma acc exit data delete(this)

    if (world_rank == 0)
    {
        std::ofstream timeWrite;
        timeWrite.open("time.csv", std::ofstream::app);
        timeWrite << num_particles << ',' <<update_time << std::endl;
    

        // if (true)
        // {
        //     for (int i = 0; i < num_particles; i++)
        //     {
        //         std::cerr<<"===================================="<<std::endl;
        //         std::cerr<<pos[i]<<std::endl;
        //         std::cerr<<vel[i]<<std::endl;
        //         std::cerr<<"===================================="<<std::endl;
        //     }
        // }

        const time_t end = time(nullptr);

        std::cerr << "Simulation complete!\n"
                  << num_particles << " Particles and " << steps << " Steps.\n"
                  << "Simulated time: " << steps * dt << " seconds\n"
                  << "Computation time: " << end - start << " seconds\n";
        std::cerr << "\n===============================================================\n";
    }
    if (!ballBuffer.str().empty())
    {
        ballBuffer<<'\n';
        energyBuffer<<'\n';
        bufferBarf();
    }

}  // end simLooper

void Ball_group::sim_one_step(const bool writeStep)
{
    /// FIRST PASS - Update Kinematic Parameters:
    #pragma acc parallel loop gang worker present(this,velh[0:num_particles],vel[0:num_particles],acc[0:num_particles],dt,wh[0:num_particles],w[0:num_particles],aacc[0:num_particles],pos[0:num_particles],num_particles)
    for (int Ball = 0; Ball < num_particles; Ball++) {
        // Update velocity half step:
        velh[Ball] = vel[Ball] + .5 * acc[Ball] * dt;

        // Update angular velocity half step:
        wh[Ball] = w[Ball] + .5 * aacc[Ball] * dt;

        // Update position:
        pos[Ball] += velh[Ball] * dt;

        // Reinitialize acceleration to be recalculated:
        acc[Ball] = {0.0,0.0,0.0};

        // Reinitialize angular acceleration to be recalculated:
        aacc[Ball] = {0.0,0.0,0.0};
    }

    #pragma acc parallel loop gang worker num_gangs(108) present(this,num_particles,accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles])
    for (int i = 0; i < num_particles*num_particles; ++i)
    {
        accsq[i] = {0.0,0.0,0.0};
        aaccsq[i] = {0.0,0.0,0.0};
    }


    double pe = 0.0;
    #pragma acc enter data copyin(pe)
    #pragma acc enter data copyin(writeStep)

    double t0 = omp_get_wtime();
    // std::cerr<<"IN simonestep"<<std::endl;

    #pragma acc parallel loop gang worker num_gangs(108) num_workers(256) reduction(+:pe) present(pe,this,accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],m[0:num_particles],moi[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs],num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,writeStep,world_rank,world_size)
    for (int pc = world_rank+1; pc <= num_pairs; pc += world_size)
    {
        // if (writeStep)
        // {
        //     pe += 1;
        // }
        
        double pd = (double)pc;
        pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
        pd -= 0.00001;
        int A = (int)pd;
        int B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);

        const double sumRaRb =  R[A] + R[B];
        const vec3 rVecab = pos[B] - pos[A];  // Vector from a to b.
        const vec3 rVecba = -rVecab;
        const double dist = (rVecab).norm();

        // Check for collision between Ball and otherBall:
        double overlap = sumRaRb - dist;

        vec3 totalForceOnA{0, 0, 0};

        int e = pc-1;
        double oldDist = distances[e];

        // if (true) {
        if (overlap > 0) {

            double k;
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

            // ==========================================
            // Test new vdw force equation with less division
            const double d1 = h2 + twoRah + twoRbh;
            const double d2 = d1 + 4 * Ra * Rb;
            const double numer = 64*Ha*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
            const double denomrecip = 1/(6*d1*d1*d2*d2);
            const vec3 vdwForceOnA = (numer*denomrecip)*rVecab.normalized();
            // ==========================================
            // const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
            //                          ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
            //                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
            //                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
            //                          rVecab.normalized();
            

            const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

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
            // const vec3 gravForceOnA = (G * m[A] * m[B] * grav_scale / (dist * dist)) * (rVecab / dist); //SCALE MASS
            const vec3 gravForceOnA = {0,0,0};
            // const vec3 gravForceOnA = (G * m[A] * m[B] / (dist * dist)) * (rVecab / dist);

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

            if (rel_vel_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
            {
                slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
            }



            // Compute rolling friction force:
            const double w_diff_mag = w_diff.norm();
            // if (w_diff_mag > 1e-20)  // Divide by zero protection.
            // if (w_diff_mag > 1e-8)  // Divide by zero protection.
            if (w_diff_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
            {
                rollForceA = 
                        -u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
                        (w_diff).cross(r_a).norm();
            }

            // Total forces on a:
            // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;
            ////////////////////////////////
            totalForceOnA = viscoelaticforceOnA + gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;
            ////////////////////////////////

            // Total torque a and b:
            torqueA = r_a.cross(slideForceOnA + rollForceA);
            torqueB = r_b.cross(-slideForceOnA + rollForceA); // original code

            // omp_set_lock(&writelock);
            // #pragma omp critical
            // {
            vec3 aaccA = (1/moi[A])*torqueA;
            vec3 aaccB = (1/moi[B])*torqueB;

            aaccsq[A*num_particles+B].x = aaccA.x;
            aaccsq[A*num_particles+B].y = aaccA.y;
            aaccsq[A*num_particles+B].z = aaccA.z;
            aaccsq[B*num_particles+A].x = aaccB.x;
            aaccsq[B*num_particles+A].y = aaccB.y;
            aaccsq[B*num_particles+A].z = aaccB.z;


            if (writeStep) {
                // No factor of 1/2. Includes both spheres:
                // PE += -G * m[A] * m[B] * grav_scale / dist + 0.5 * k * overlap * overlap;
                // PE += -G * m[A] * m[B] / dist + 0.5 * k * overlap * overlap;

                // Van Der Waals + elastic:
                const double diffRaRb = R[A] - R[B];
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * R[A] * R[B];
                const double denom_sum = z * z - (sumRaRb * sumRaRb);
                const double denom_diff = z * z - (diffRaRb * diffRaRb);
                const double U_vdw =
                    -Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + 
                    log(denom_sum / denom_diff));
                pe += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
            }
        } else  // Non-contact forces:
        {

            // No collision: Include gravity and vdw:
            // const vec3 gravForceOnA = (G * m[A] * m[B] * grav_scale / (dist * dist)) * (rVecab / dist);
            const vec3 gravForceOnA = {0.0,0.0,0.0};
            // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
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
            // ==========================================
            // Test new vdw force equation with less division
            const double d1 = h2 + twoRah + twoRbh;
            const double d2 = d1 + 4 * Ra * Rb;
            const double numer = 64*Ha*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
            const double denomrecip = 1/(6*d1*d1*d2*d2);
            const vec3 vdwForceOnA = (numer*denomrecip)*rVecab.normalized();
            // ==========================================
            // const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
            //                          ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
            //                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
            //                                            (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
            //                          rVecab.normalized();
            // const vec3 vdwForceOnA = {0.0,0.0,0.0};
            /////////////////////////////
            totalForceOnA = vdwForceOnA + gravForceOnA;
            // totalForceOnA = vdwForceOnA;
            // totalForceOnA = gravForceOnA;
            /////////////////////////////
            if (writeStep) {
                // PE += -G * m[A] * m[B] * grav_scale / dist; // Gravitational
                const double diffRaRb = R[A] - R[B];
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * R[A] * R[B];
                const double denom_sum = z * z - (sumRaRb * sumRaRb);
                const double denom_diff = z * z - (diffRaRb * diffRaRb);
                const double U_vdw =
                    -Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                pe += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
            }

            // todo this is part of push_apart. Not great like this.
            // For pushing apart overlappers:
            // vel[A] = { 0,0,0 };
            // vel[B] = { 0,0,0 };
        }


        vec3 accA = (1/m[A])*totalForceOnA; 
        vec3 accB = -1.0*(1/m[B])*totalForceOnA; 

        accsq[A*num_particles+B].x = accA.x;
        accsq[A*num_particles+B].y = accA.y;
        accsq[A*num_particles+B].z = accA.z;
        accsq[B*num_particles+A].x = accB.x;
        accsq[B*num_particles+A].y = accB.y;
        accsq[B*num_particles+A].z = accB.z;

        distances[e] = dist;


        // #pragma acc update host(pe)

    }

    // #pragma acc loop seq
    #pragma acc parallel loop gang num_gangs(108) present(this,num_particles,acc[0:num_particles],aacc[0:num_particles],accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles])
    for (int i = 0; i < num_particles; i++)
    {
        #pragma acc loop seq
        for (int j = 0; j < num_particles; j++)
        {
            acc[i].x += accsq[i*num_particles+j].x;
            acc[i].y += accsq[i*num_particles+j].y;
            acc[i].z += accsq[i*num_particles+j].z;
            aacc[i].x += aaccsq[i*num_particles+j].x;
            aacc[i].y += aaccsq[i*num_particles+j].y;
            aacc[i].z += aaccsq[i*num_particles+j].z;
        }
    // #pragma acc update self(acc[0:num_particles],aacc[0:num_particles]) //if(write_step)
    // #pragma acc update self(acc[i],aacc[i]) //if(write_step)
    }

    #pragma acc update host(pe,acc[0:num_particles],aacc[0:num_particles])
    // std::cout<<aaccsq[0].x<<','<<aaccsq[0].y<<','<<aaccsq[0].z<<std::endl;
    PE = pe;

    #pragma acc exit data delete(pe,writeStep)
    // std::cerr<<"PEpre: "<<PE<<std::endl;
    // std::cerr<<"acc: "<<acc[0].x<<','<<acc[0].y<<','<<acc[0].z<<std::endl;

    #ifdef MPI_ENABLE
        MPI_Allreduce(MPI_IN_PLACE,acc,num_particles*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,aacc,num_particles*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if (writeStep)
        {
            double local_PE = PE;
            // std::cout<<"PE in rank "<<world_rank<<" : "<<PE<<std::endl;
            PE = 0.0;
            MPI_Reduce(&local_PE,&PE,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            // if (world_rank == 0)
            // {
            //     std::cerr<<"PE: "<<PE<<std::endl;
            //     std::cerr<<"acc: "<<acc[0].x<<','<<acc[0].y<<','<<acc[0].z<<std::endl;
            // }
        }
    #endif

    double t1 = omp_get_wtime();
    update_time += t1-t0;


    #pragma acc parallel loop gang worker num_gangs(108) num_workers(256) present(this,acc[0:num_particles],aacc[0:num_particles],w[0:num_particles],vel[0:num_particles],velh[0:num_particles],wh[0:num_particles],num_particles,dt)
    for (int Ball = 0; Ball < num_particles; Ball++) {
        // Velocity for next step:
        vel[Ball] = velh[Ball] + .5 * acc[Ball] * dt;
        w[Ball] = wh[Ball] + .5 * aacc[Ball] * dt;
    }  // THIRD PASS END


    #pragma acc update host(w[0:num_particles],vel[0:num_particles],pos[0:num_particles]) if(writeStep && world_rank == 0)
    if (writeStep && world_rank == 0) 
    {
        // std::cerr<<"Writing "<<num_particles<<" balls"<<std::endl;
        ballBuffer << '\n';  // Prepares a new line for incoming data.
        for (int Ball = 0; Ball < num_particles;  Ball++)
        {
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
            // Send positions and rotations to buffer:
            // std::cerr<<"Write ball "<<Ball<<std::endl;
            KE += .5 * m[Ball] * vel[Ball].normsquared() +
                    .5 * moi[Ball] * w[Ball].normsquared();  // Now includes rotational kinetic energy.
            mom += m[Ball] * vel[Ball];
            ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
        }
    }

    // #pragma acc exit data delete(ke,Mom,Ang_mom)

    // t.end_event("CalcVelocityforNextStep");
}  // one Step end