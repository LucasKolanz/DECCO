// This file was originally written for SpaceLab/DECCO

// The Ball_group class is meant to encompass all the balls in an aggregate collision or formation simulation

// Authors: Job Guidos and Lucas Kolanz

#pragma once
#ifndef BALL_GROUP_HPP
#define BALL_GROUP_HPP

#include <openacc.h>
#include "../external/json/single_include/nlohmann/json.hpp"
#include "../utilities/vec3.hpp"
#include "../utilities/MPI_utilities.hpp"
#include "../data/DECCOData.hpp"


#include <vector>

#ifdef MPI_ENABLE
    #include <mpi.h>
#endif

#ifdef EXPERIMENTAL_FILESYSTEM
    #include <experimental/filesystem>
    namespace fs = std::experimental::filesystem;
#else
    #include <filesystem>
    namespace fs = std::filesystem;
#endif

// using std::numbers::pi;
using json = nlohmann::json;
extern const int bufferlines;
enum distributions {constant, logNorm};
enum simType {BPCA, BCCA, BAPA, collider, relax};

constexpr double Kb = 1.380649e-16; //in erg/K
constexpr double pi = 3.14159265358979311599796346854;


//This struct is meant to encompass all values necessary to carry out a simulation, besides the physical
//attributes, such as position, velocity, etc. Note: it is important to add variables to the "=" operator definition
//if more attributes are added to this struct and they are expected to carry over from one sim to the next such as 
//during a BPCA growth simulation.
struct Ball_group_attributes
{
    std::string project_path = "";
    std::string output_folder = "";
    std::string data_directory = "";
    std::string projectileName = "";
    std::string targetName = "";
    std::string output_prefix = "";
    std::string random_folder_template = "";

    double radiiFraction = -1;
    bool debug = false;
    bool write_all = false;
    bool mid_sim_restart = false;

    // std::string out_folder;
    int num_particles = 0;
    int num_particles_added = 0;
    int MAXOMPthreads = 1;
    int OMPthreads = 1;
    int MAXMPInodes = 1;
    int MPInodes = 1;
    int start_index = 0;
    int start_step = 1;
    int relax_index = 0;

    int skip=-1;  // Steps thrown away before recording a step to the buffer. 500*.04 is every 20 seconds in sim.
    unsigned long long steps=0;

    double dt=-1;
    double kin=-1;  // Spring constant
    double kout=-1;

    //Don't copy these during add_particle. These are set at the beginning (or during) of sim_looper
    int world_rank = -1;
    int world_size = -1;
    int num_pairs = -1;
    bool write_step = false;

    const std::string sim_meta_data_name = "sim_info";

    int seed = -1;
    int output_width = -1;
    distributions radiiDistribution;
    simType typeSim = BPCA; //Default to BPCA for now
    double lnSigma = 0.2; //sigma for log normal distribution 

    //Variable only has an effect for BCCA simulations (typeSim = BCCA)
    //if true, projectile is copy of target, if false, projectile is taken from a 
    //different, random simulation that has already made it to this point.
    bool symmetric = true; 

    // Useful values:
    double r_min = -1;
    double r_max = -1;
    double m_total = -1;
    double initial_radius = -1;
    double v_collapse = 0;
    double v_max = -1;
    double v_max_prev = HUGE_VAL;
    double soc = -1;

    ///none of these should be negative so verify they were set before using
    bool dynamicTime = false;
    double G = 6.67e-08;  // in dyn*cm^2*g^-2// Gravitational constant defaults to the actual value, but you can change it
    double density = -1.0;
    double u_s = -1.0;                // Coeff of sliding friction
    double u_r = -1.0;               // Coeff of rolling friction
    double sigma = -1.0;              // Poisson ratio for rolling friction.
    double Y = -1.0;               // Young's modulus in erg/cm3
    double cor = -1.0;                // Coeff of restitution
    double simTimeSeconds = -1.0;  // Seconds
    double timeResolution = -1.0;    // Seconds - This is duration between exported steps.
    double fourThirdsPiRho = -1.0;  // for fraction of smallest sphere radius.
    double scaleBalls = -1.0;                         // base radius of ball.
    double maxOverlap = -1.0 ;                           // of scaleBalls
    double KEfactor = -1.0;                              // Determines collision velocity based on KE/PE
    double v_custom = -1.0;  // Velocity cm/s
    double temp = -1.0;          //tempurature of simulation in Kelvin
    double eta = -1.0;  //eta = KE/PE (used for setting speed of projectile based on a desired energy ratio)
    double kConsts = -1.0;
    double impactParameter = -1.0;  // Impact angle radians
    double Ha = -1.0;         // Hamaker constant for vdw force
    double h_min = -1.0;  // 1e8 * std::numeric_limits<double>::epsilon(), // 2.22045e-10 (epsilon is 2.22045e-16)
    double cone = -1.0;  // Cone of particles ignored moving away from center of mass. Larger angle ignores more.


    // Simulation Structure
    int properties = -1;  // Number of columns in simData file per ball
    int genBalls = -1;
    int attempts = -1.0;  // How many times to try moving every ball touching another in generator.


    double spaceRange = -1.0;  // Rough minimum space required
    double spaceRangeIncrement = -1.0;
    double z0Rot = -1.0;  // Cluster one z axis rotation
    double y0Rot = -1.0;  // Cluster one y axis rotation
    double z1Rot = -1.0;  // Cluster two z axis rotation
    double y1Rot = -1.0;  // Cluster two y axis rotation
    double simTimeElapsed = -1.0;

    int N=-1; //Number of balls to grow (if BPCA or BAPA)
    int M=-1; //BAPA size (fragment size of intermediate growth steps)

    const time_t start = time(nullptr);  // For end of program analysis
    time_t startProgress = 0;                // For progress reporting (gets reset)
    time_t lastWrite = 0;                    // For write control (gets reset)

    /////////////////////////////////
    const double h_min_physical = 2.1e-8; //prolly should make this a parameter/calculation
    const double max_mu = 0.5; // Make another parameter
    bool mu_scale = false;
    /////////////////////////////////
    // data_type 0 for hdf5 
    // data_type 1 for csv 
    int data_type = 0;
    std::string filetype = "h5";
    int num_writes = 0;


    // Overload the assignment operator
    Ball_group_attributes& operator=(const Ball_group_attributes& other) 
    {
        if (this != &other) // Protect against self-assignment 
        {  
            project_path = other.project_path;
            output_folder = other.output_folder;
            data_directory = other.data_directory;
            projectileName = other.projectileName;
            targetName = other.targetName;
            output_prefix = other.output_prefix;
            random_folder_template = other.random_folder_template;

            radiiFraction = other.radiiFraction;
            debug = other.debug;
            write_all = other.write_all;
            mid_sim_restart = other.mid_sim_restart;

            num_particles = other.num_particles;
            num_particles_added = other.num_particles_added;
            MAXOMPthreads = other.MAXOMPthreads;
            OMPthreads = other.OMPthreads;
            MAXMPInodes = other.MAXMPInodes;
            MPInodes = other.MPInodes;
            start_index = other.start_index;
            start_step = other.start_step;
            relax_index = other.relax_index;

            skip = other.skip;
            steps = other.steps;

            dt = other.dt;
            kin = other.kin;
            kout = other.kout;


            seed = other.seed;
            output_width = other.output_width;
            radiiDistribution = other.radiiDistribution;
            typeSim = other.typeSim;
            symmetric = other.symmetric;
            lnSigma = other.lnSigma;

            r_min = other.r_min;
            r_max = other.r_max;
            m_total = other.m_total;
            initial_radius = other.initial_radius;
            v_collapse = other.v_collapse;
            v_max = other.v_max;
            v_max_prev = other.v_max_prev;
            soc = other.soc;
            N = other.N;
            M = other.M;

            dynamicTime = other.dynamicTime;
            G = other.G;
            density = other.density;
            u_s = other.u_s;
            u_r = other.u_r;
            sigma = other.sigma;
            Y = other.Y;
            cor = other.cor;
            simTimeSeconds = other.simTimeSeconds;
            timeResolution = other.timeResolution;
            fourThirdsPiRho = other.fourThirdsPiRho;
            scaleBalls = other.scaleBalls;
            maxOverlap = other.maxOverlap;
            KEfactor = other.KEfactor;
            v_custom = other.v_custom;
            temp = other.temp;
            eta = other.eta;
            kConsts = other.kConsts;
            impactParameter = other.impactParameter;
            Ha = other.Ha;
            h_min = other.h_min;
            cone = other.cone;

            properties = other.properties;
            genBalls = other.genBalls;
            attempts = other.attempts;

            spaceRange = other.spaceRange;
            spaceRangeIncrement = other.spaceRangeIncrement;
            z0Rot = other.z0Rot;
            y0Rot = other.y0Rot;
            z1Rot = other.z1Rot;
            y1Rot = other.y1Rot;
            simTimeElapsed = other.simTimeElapsed;

            // start is const, no need to copy
            startProgress = other.startProgress;
            lastWrite = other.lastWrite;
            simTimeElapsed = other.simTimeElapsed;

            // h_min_physical and max_mu are const, no need to copy
            mu_scale = other.mu_scale;
            dynamicTime = other.dynamicTime;

            data_type = other.data_type;
            filetype = other.filetype;
            num_writes = other.num_writes;
        }
        return *this;
    }

};


/// @brief Facilitates the concept of a group of balls with physical properties.
class Ball_group
{
public:
    // enum distributions {constant, logNorm};
    // enum simType {BPCA, BCCA, collider, relax};
    

    Ball_group_attributes attrs;

    /////////////////////////////////
    // bool mu_scale = false;
    /////////////////////////////////


    //////////////////////////////////
    //The following attributes are the physical characteristics of all balls contained in Ball_group
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
    #ifdef GPU_ENABLE
        vec3* accsq = nullptr;
        vec3* aaccsq = nullptr;
    #endif
    double* R = nullptr;    ///< Radius
    double* m = nullptr;    ///< Mass
    double* moi = nullptr;  ///< Moment of inertia
    //////////////////////////////////

    //The DECCOData class takes care of reading and writing data in whatever format is specified in the input file 
    DECCOData* data = nullptr;
    
    //Buffers for energy and ball (pos,vel,etc.) data to avoid writing to disk too much
    std::vector<double> energyBuffer;
    std::vector<double> ballBuffer;

    //Constructors
    Ball_group() = default;
    Ball_group(std::string& path,bool test); //Testing constructor
    explicit Ball_group(const int nBalls);
    // explicit Ball_group(const std::string& path, const std::string& filename, int start_file_index);
    explicit Ball_group(std::string& path,const int index=-1);
    // explicit Ball_group(const std::string& path,const std::string& projectileName,const std::string& targetName,const double& customVel);
    Ball_group(const Ball_group& rhs);
    Ball_group& operator=(const Ball_group& rhs);

    //Functions directly related to changing the state of a Ball_group
    void kick(const vec3& vec) const;
    void move(const vec3& vec) const;
    void offset(const double& rad1, const double& rad2, const double& impactParam) const;
    void pushApart() const;
    void updatePE();
    void pos_and_vel_for_collision(Ball_group &projectile,Ball_group &target);
    void pos_and_vel_for_collision(Ball_group &projectile);
    void overwrite_v_custom(Ball_group &projectile,Ball_group &target);
    void overwrite_v_custom(Ball_group &projectile);
    void to_origin() const;
    void rotAll(const char axis, const double angle) const;
    // vec3 random_offset(const double3x3 local_coords,vec3 projectile_pos,vec3 projectile_vel,const double projectile_rad);
    vec3 random_offset(Ball_group &projectile, Ball_group &target);
    void comSpinner(const double& spinX, const double& spinY, const double& spinZ) const;
    void sim_one_step(int step);
    void sim_looper(unsigned long long start_step);
    
    //Functions which calculate/set values for Ball_group
    inline double calc_VDW_force_mag(const double Ra, const double Rb, const double h);
    // void calc_mu_scale_factor();
    void calibrate_dt(int const Step, const double& customSpeed);
    void calc_v_collapse();
    [[nodiscard]] double getVelMax();
    void calc_helpfuls(const bool includeRadius=true);
    double get_soc();    
    vec3 calc_momentum(const std::string& of) const;
    [[nodiscard]] double getRadius(const vec3& center) const;
    [[nodiscard]] vec3 getCOM() const;
    [[nodiscard]] vec3 getVCOM() const;
    void zeroVel() const;
    // void zeroSaveVals();
    void zeroAngVel() const;
    double calc_mass(const double& radius, const double& density);
    double calc_moi(const double& radius, const double& mass);
    double calc_max_bolt_velocity(double temp, double mass);
    double calc_eta_velocity(const double eta, Ball_group &projectile, Ball_group &target);
    double calc_group_noncontact_PE(Ball_group &projectile,Ball_group &target);
    double calc_noncontact_PE();
    [[nodiscard]] double getRmin();
    [[nodiscard]] double getRmax();
    [[nodiscard]] double getMmin() const;
    [[nodiscard]] double getMmax() const;
    double getMass();
    void set_seed_from_input(const std::string location);

    //Initializers
    Ball_group BPCA_projectile_init();
    Ball_group BCCA_projectile_init(const bool symmetric);
    Ball_group BAPA_projectile_init();
    void init_data(int counter);
    void relaxInit(const std::string path);
    void aggregationInit(const std::string path,const int index=-1);
    void colliderInit(const std::string path);
    void sim_init_write(int counter=0);
    void parse_input_file(std::string location);
    Ball_group spawn_particles(const int count);
    Ball_group add_projectile(const simType);
    void merge_ball_group(const Ball_group& src,const bool includeRadius=true);
    void freeMemory() const;
    std::string find_whole_file_name(std::string path,const int index=-1);
    int check_restart(std::string folder);
    #ifdef HDF5_ENABLE
        void loadDatafromH5(std::string path, std::string file);
    #endif
    std::string get_data_info();
    void parse_meta_data(std::string metadata);
    std::string find_file_name(std::string path,int index);
    int get_num_threads();
    std::string data_type_from_input(const std::string location);





    
    bool isAggregation();
    void allocate_group(const int nBalls);
    void init_conditions();
    void parseSimData(std::string line);
    void loadConsts(const std::string& path, const std::string& filename);
    [[nodiscard]] static std::string getLastLine(const std::string& path, const std::string& filename);
    // void simDataWrite(std::string& outFilename);
    void threeSizeSphere(const int nBalls);
    void generate_ball_field(const int nBalls);
    void loadSim(const std::string& path, const std::string& filename);
    void distSizeSphere(const int nBalls);
    void oneSizeSphere(const int nBalls);
    void placeBalls(const int nBalls);
    void updateDTK(const double& velocity);
    void simInit_cond_and_center(bool add_prefix);
    void sim_continue(const std::string& path);
    void sim_init_two_cluster(const std::string& path,const std::string& projectileName,const std::string& targetName);
    void verify_projectile(const std::string projectile_folder, const int index, const double max_wait_time);
private:

};


    
bool is_touching(Ball_group &projectile,Ball_group &target);
void moveApart(const vec3 &projectile_direction,Ball_group &projectile,Ball_group &target);

#endif