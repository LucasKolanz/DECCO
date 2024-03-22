#pragma once
#include "../default_files/dust_const.hpp"
// #include "dust_const.hpp"
#include "../external/json/single_include/nlohmann/json.hpp"
#include "../utilities/vec3.hpp"
#include "../utilities/linalg.hpp"
#include "../utilities/Utils.hpp"
#include "../utilities/MPI_utilities.hpp"
#include "../data/DECCOData.hpp"
#include "../timing/timing.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <limits.h>
#include <cstring>
#include <typeinfo>
#include <memory>
#include <random>
#include <omp.h>

#ifdef MPI_ENABLE
    #include <mpi.h>
#endif


// using std::numbers::pi;
using json = nlohmann::json;
namespace fs = std::filesystem;
extern const int bufferlines;


struct Ball_group_attributes
{


    std::string project_path = "";
    std::string output_folder = "";
    std::string data_directory = "";
    std::string projectileName = "";
    std::string targetName = "";
    std::string output_prefix = "";

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
    enum distributions {constant, logNorm};
    distributions radiiDistribution;
    enum simType {BPCA, collider, relax};
    simType typeSim = BPCA; //Default to BPCA for now
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

    int N=-1; //Number of balls to grow (if BPCA)

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
    Ball_group_attributes attrs;

    /////////////////////////////////
    // bool mu_scale = false;
    /////////////////////////////////

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
    // double* u_scale = nullptr; ///ADD TO COPY CONSTRUCTOR, ETC


    DECCOData* data = nullptr;
    // std::unique_ptr<DECCOData> data;
    
    std::vector<double> energyBuffer;
    std::vector<double> ballBuffer;

    Ball_group() = default;

    explicit Ball_group(const int nBalls);
    // explicit Ball_group(const std::string& path, const std::string& filename, int start_file_index);
    explicit Ball_group(std::string& path);
    explicit Ball_group(const std::string& path,const std::string& projectileName,const std::string& targetName,const double& customVel);
    Ball_group(const Ball_group& rhs);
    Ball_group& operator=(const Ball_group& rhs);
    void parse_input_file(std::string location);
    inline double calc_VDW_force_mag(const double Ra, const double Rb, const double h);
    // void calc_mu_scale_factor();
    // void zeroSaveVals();
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
    void sim_init_write(int counter);
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
    void freeMemory() const;
    std::string find_restart_file_name(std::string path);
    int check_restart(std::string folder);
    #ifdef HDF5_ENABLE
        void loadDatafromH5(std::string path, std::string file);
    #endif
    void init_data(int counter);
    std::string get_data_info();
    void parse_meta_data(std::string metadata);
    void relax_init(std::string path);
    void BPCA_init(std::string path);
    std::string find_file_name(std::string path,int index);
    int get_num_threads();


    void sim_one_step();
    #ifdef GPU_ENABLE
        void sim_one_step_GPU();
    #endif
    void sim_looper(unsigned long long start_step);


    
private:

    void allocate_group(const int nBalls);
    void init_conditions();
    [[nodiscard]] double getRmin();
    [[nodiscard]] double getRmax();
    [[nodiscard]] double getMassMax() const;
    void parseSimData(std::string line);
    void loadConsts(const std::string& path, const std::string& filename);
    [[nodiscard]] static std::string getLastLine(const std::string& path, const std::string& filename);
    // void simDataWrite(std::string& outFilename);
    [[nodiscard]] double getMass();
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
};

/// @brief For creating a new ballGroup of size nBalls
/// @param nBalls Number of balls to allocate.
Ball_group::Ball_group(const int nBalls)
{
    allocate_group(nBalls);
    for (size_t i = 0; i < nBalls; i++) {
        R[i] = 1;
        m[i] = 1;
        moi[i] = calc_moi(R[i], m[i]);
    }
}

/// @brief For generating a new ballGroup of any simTypes
/// @param path is a path to the job folder
Ball_group::Ball_group(std::string& path)
{
    parse_input_file(path);

    if (attrs.typeSim == attrs.BPCA)
    {
        BPCA_init(path);
    }
    else if (attrs.typeSim == attrs.collider)
    {
        MPIsafe_print(std::cerr,"COLLIDER NOT IMPLIMENTED. NOW EXITING . . .\n");
        MPIsafe_exit(-1);
    }
    else if (attrs.typeSim == attrs.relax)
    {
        if (attrs.relax_index > 0)
        {
            relax_init(path);
        }
        else
        {
            std::string message("ERROR: simType is relax but relax_index is ("+std::to_string(attrs.relax_index)+") < 0\n");
            MPIsafe_print(std::cerr,message);
        }
    }
        
}

// Initializes relax job (only new, not for restart)
void Ball_group::relax_init(std::string path)
{
    std::string filename = find_file_name(path,attrs.relax_index); 

    loadSim(path, filename.substr(filename.find_last_of('/')+1,filename.size()));

    zeroVel();
    // zeroAngVel();

    calc_v_collapse(); 
    if (attrs.dt < 0)
    {
        calibrate_dt(0, attrs.v_custom);
    }
    simInit_cond_and_center(false);

}


// Initializes BPCA job for restart or new job
void Ball_group::BPCA_init(std::string path)
{
    int restart = check_restart(path);

    // If the simulation is complete exit now. Otherwise, the call to 
    //find_restart_file_name will possibly delete one of the data files 
    if (restart == 2)
    {
        MPIsafe_print(std::cerr,"Simulation already complete. Now exiting. . .\n");
        MPIsafe_exit(0);
    }

    std::string filename = find_restart_file_name(path); 
    bool just_restart = false;

    if (filename != "")
    {
        if (filename.substr(filename.size()-4,filename.size()) == ".csv")
        {
            size_t _pos = filename.find_first_of("_");
            // int file_index = stoi(filename.substr(0,filename.find_first_of("_")));
            if (filename[_pos+1] == 'R')
            {
                just_restart = true;
            }
        }
    }

    

    if (!just_restart && restart==1)
    {
        MPIsafe_print(std::cerr,std::string("Loading sim "+path+filename+'\n'));
        loadSim(path, filename);
        calc_v_collapse(); 
        // getMass();
        if (attrs.dt < 0)
            calibrate_dt(0, attrs.v_custom);
        simInit_cond_and_center(false);
    }
    else if (restart == 0 || just_restart)
    {

        generate_ball_field(attrs.genBalls);
        // Hack - Override and creation just 2 balls position and velocity.
        if (attrs.genBalls > 0 && attrs.genBalls <= 2)
        {
            pos[0] = {0, R[0]+1.01e-6, 0};
            vel[0] = {0, 0, 0};
            if (attrs.genBalls > 1)
            {
                pos[1] = {0, -(R[1]+1.01e-6), 0};
                vel[1] = {0, 0, 0};
        
            }
        }
        else
        {
            MPIsafe_print(std::cerr,"ERROR: genBalls > 2 not yet implimented (right)?\n");
        }

        // if (mu_scale)
        // {
        //     calc_mu_scale_factor();
        // }
        // std::cerr<<initial_radius<<std::endl;

        attrs.m_total = getMass();
        calc_v_collapse();
        // std::cerr<<"INIT VCUSTOM "<<v_custom<<std::endl;
        calibrate_dt(0, attrs.v_custom);
        simInit_cond_and_center(true);
        
        
    }
    else
    {
        std::string message("ERROR: restart code '"+std::to_string(restart)+"' not recognized.\n");
        MPIsafe_print(std::cerr,message);

        MPIsafe_exit(-1);
    }
}

// /// @brief For continuing a sim.
// /// @param fullpath is the filename and path excluding the suffix _simData.csv, _constants.csv, etc.
// /// @param customVel To condition for specific vMax.
// Ball_group::Ball_group(const std::string& path, const std::string& filename, int start_file_index=0)
// {
//     parse_input_file(std::string(path));
//     sim_continue(path, filename,start_file_index);
//     calc_v_collapse();
//     calibrate_dt(0, v_custom);
//     simInit_cond_and_center(false);
// }

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
    parse_input_file(std::string(path));
    // std::cerr<<path<<std::endl;
    sim_init_two_cluster(path, projectileName, targetName);
    calc_v_collapse();
    if (customVel > 0){calibrate_dt(0, customVel);}
    else {calibrate_dt(0, attrs.v_custom);}
    simInit_cond_and_center(true);
}

Ball_group& Ball_group::operator=(const Ball_group& rhs)
{

    attrs = rhs.attrs;

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


    data = rhs.data;

   

    return *this;
}

Ball_group::Ball_group(const Ball_group& rhs)
{


    attrs = rhs.attrs;

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


    data = rhs.data;


}

void Ball_group::init_data(int counter = 0)
{


    if (data != nullptr)
    {
        delete data;
        data = nullptr; 
    }

    std::string type = "";
    
    if (attrs.typeSim == attrs.relax)
    {
        type = "RELAX";
    }

    std::string sav_file;
    if (attrs.data_type == 0) //h5
    {
        sav_file = attrs.output_folder+std::to_string(counter)+"_"+type+"data."+attrs.filetype;
    }
    else if (attrs.data_type == 1) //csv
    {
        sav_file = attrs.output_folder+std::to_string(counter)+"_"+type+".csv";
    }
    else
    {
        std::string message("ERROR: data_type '"+attrs.filetype+"' not supported.\n"); 
        MPIsafe_print(std::cerr,message);
        MPIsafe_exit(EXIT_FAILURE);
    }
    data = new DECCOData(sav_file,\
                        attrs.num_particles,attrs.steps/attrs.skip+1,attrs.steps);
    
}


//Parses input.json file that is in the same folder the executable is in
void Ball_group::parse_input_file(std::string location)
{

    if (location == "")
    {
        try {
            std::filesystem::path currentPath = std::filesystem::current_path();
            location = currentPath.string() + "/";
        } catch (const std::filesystem::filesystem_error& e) {
            MPIsafe_print(std::cerr,std::string("Error getting current directory: " + std::string(e.what()) + '\n'));
            exit(-1);
        }
    }
    // std::string s_location(location);
    std::string json_file = location + "input.json";
    MPIsafe_print(std::cerr,std::string("Parsing input file: "+json_file+'\n'));
    std::ifstream ifs(json_file);
    json inputs = json::parse(ifs);
    attrs.output_folder = inputs["output_folder"];
    attrs.data_directory = inputs["data_directory"];

    if (inputs["simType"] == "BPCA")
    {
        attrs.typeSim = attrs.BPCA;
    }
    else if (inputs["simType"] == "collider")
    {
        attrs.typeSim = attrs.collider;
    }
    else if (inputs["simType"] == "relax")
    {
        attrs.typeSim = attrs.relax;
        attrs.relax_index = inputs["relaxIndex"];
    }

    if (inputs["dataFormat"] == "h5" || inputs["dataFormat"] == "hdf5")
    {
        attrs.data_type = 0;
        attrs.filetype = "h5";
    }
    else if (inputs["dataFormat"] == "csv")
    {
        attrs.data_type = 1;
        attrs.filetype = "csv";
    }


    if (inputs["seed"] == std::string("default"))
    {
        attrs.seed = static_cast<unsigned int>(time(nullptr));
    }
    else
    {
        attrs.seed = static_cast<unsigned int>(inputs["seed"]);
    }

    if (getRank() == 0)
    {
        std::ofstream seedFile;
        seedFile.open(attrs.output_folder+"seedFile.txt",std::ios::app);
        seedFile<<attrs.seed<<std::endl;
        seedFile.close();
    }
    
    
    MPIsafe_print(std::cerr,std::string("Writing seed '"+std::to_string(attrs.seed)+"' to seedFile.txt\n"));
    
    random_generator.seed(attrs.seed);//This was in the else but it should be outside so random_generator is always seeded the same as srand (right?)
    srand(attrs.seed);

    std::string temporary_distribution = inputs["radiiDistribution"];
    std::transform(temporary_distribution.begin(), temporary_distribution.end(), temporary_distribution.begin(), ::tolower);
    if (temporary_distribution == "lognormal" || temporary_distribution == "lognorm")
    {
        attrs.radiiDistribution = attrs.logNorm;
    }
    else
    {
        attrs.radiiDistribution = attrs.constant;
    }
    attrs.N = inputs["N"];
    attrs.dynamicTime = inputs["dynamicTime"];
    attrs.G = inputs["G"];
    attrs.density = inputs["density"];
    attrs.u_s = inputs["u_s"];
    attrs.u_r = inputs["u_r"];
    attrs.sigma = inputs["sigma"];
    attrs.Y = inputs["Y"];
    attrs.cor = inputs["cor"];
    attrs.simTimeSeconds = inputs["simTimeSeconds"];
    attrs.timeResolution = inputs["timeResolution"];
    attrs.fourThirdsPiRho = 4. / 3. * pi * attrs.density;
    attrs.scaleBalls = inputs["scaleBalls"];
    attrs.maxOverlap = inputs["maxOverlap"];
    attrs.KEfactor = inputs["KEfactor"];

    attrs.MAXMPInodes = inputs["MPInodes"];
    attrs.MPInodes = attrs.MAXMPInodes;
    attrs.MAXOMPthreads = inputs["OMPthreads"];
    attrs.OMPthreads = attrs.MAXOMPthreads;
    
    if (inputs["v_custom"] == std::string("default"))
    {
        attrs.v_custom = 0.36301555459799423;
    }
    else
    {
        attrs.v_custom = inputs["v_custom"];
    }
    attrs.temp = inputs["temp"]; // this will modify v_custom in oneSizeSphere
    double temp_kConst = inputs["kConsts"];
    attrs.kConsts = temp_kConst * (attrs.fourThirdsPiRho / (attrs.maxOverlap * attrs.maxOverlap));
    attrs.impactParameter = inputs["impactParameter"];
    attrs.Ha = inputs["Ha"];
    double temp_h_min = inputs["h_min"];
    attrs.h_min = temp_h_min * attrs.scaleBalls;
    if (inputs["cone"] == std::string("default"))
    {
        attrs.cone = pi/2;
    }
    else
    {
        attrs.cone = inputs["cone"];
    }
    attrs.properties = inputs["properties"];
    attrs.genBalls = inputs["genBalls"];
    attrs.attempts = inputs["attempts"];
    attrs.skip = inputs["skip"];
    attrs.steps = inputs["steps"];
    attrs.dt = inputs["dt"];
    attrs.kin = inputs["kin"];
    attrs.kout = inputs["kout"];
    if (inputs["spaceRange"] == std::string("default"))
    {
        attrs.spaceRange = 4 * std::pow(
                        (1. / .74 * attrs.scaleBalls * attrs.scaleBalls * attrs.scaleBalls * attrs.genBalls),
                        1. / 3.); 
    }
    else
    {
        attrs.spaceRange = inputs["spaceRange"];
    }
    if (inputs["spaceRangeIncrement"] == std::string("default"))
    {
        attrs.spaceRangeIncrement = attrs.scaleBalls * 3;
    }
    else
    {
        attrs.spaceRangeIncrement = inputs["spaceRangeIncrement"];
    }
    attrs.z0Rot = inputs["z0Rot"];
    attrs.y0Rot = inputs["y0Rot"];
    attrs.z1Rot = inputs["z1Rot"];
    attrs.y1Rot = inputs["y1Rot"];
    attrs.simTimeElapsed = inputs["simTimeElapsed"];

    attrs.projectileName = inputs["projectileName"];
    attrs.targetName = inputs["targetName"];
    attrs.output_prefix = inputs["output_prefix"];
    if (attrs.output_prefix == std::string("default"))
    {
        attrs.output_prefix = "";
    }

    attrs.radiiFraction = inputs["radiiFraction"];

    attrs.output_width = attrs.num_particles;
}

// @brief calculates the vdw force
inline double Ball_group::calc_VDW_force_mag(const double Ra,const double Rb,const double h)
{
    const double h2 = h * h;
    // constexpr double h2 = h * h;
    const double twoRah = 2 * Ra * h;
    const double twoRbh = 2 * Rb * h;
    return attrs.Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
             ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                               (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                               (h2 + twoRah + twoRbh + 4 * Ra * Rb)));

}

// // @breif calculates the mu scaling factor for all pairs of particle sizes
// void Ball_group::calc_mu_scale_factor()
// {
//     int e;
//     for (int A = 1; A < attrs.num_particles; ++A)
//     {
//         for (int B = 0; B < A; ++B)
//         {   
//             e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
//             u_scale[e] = calc_VDW_force_mag(R[A],R[B],h_min_physical)/
//                                 calc_VDW_force_mag(R[A],R[B],h_min);  
//         }
//     }
// }

////////////////////////////////////
// void Ball_group::zeroSaveVals()
// {
//     int size = attrs.num_particles*attrs.num_particles;
//     for (int i = 0; i < size; ++i)
//     {
//         vdwForce[i] = {0,0,0};
//         // elasticForce[i] = {0,0,0};
//         // slideForce[i] = {0,0,0};
//         // rollForce[i] = {0,0,0};
//         // torqueForce[i] = {0,0,0};
//     }
//     // for (int i = 0; i < num_particles; ++i)
//     // {
//     //     // if (i < num_particles)
//     //     // {
//     //     // distB3[i] = 0.0;
//     //     // }
//     //     // slidDir[i] = {0,0,0};
//     //     // rollDir[i] = {0,0,0};
//     //     // inout[i] = 0.0;
//     //     // slidB3[i] = {0,0,0};
//     //     // rollB3[i] = {0,0,0};
//     //     // // slidFric[i] = {0,0,0};
//     //     // rollFric[i] = {0,0,0};
//     // }
// }
////////////////////////////////////

void Ball_group::calibrate_dt(int const Step, const double& customSpeed = -1.)
{
    const double dtOld = attrs.dt;


    std::string message = "";
    if (customSpeed > 0.) {
        updateDTK(customSpeed);
        message += "CUSTOM SPEED: " + std::to_string(customSpeed) + '\n';
    } else {
        // std::cerr << vCollapse << " <- vCollapse | Lazz Calc -> " << M_PI * M_PI * G * pow(density, 4.
        // / 3.) * pow(mTotal, 2. / 3.) * rMax;

        attrs.v_max = getVelMax();

        message += '\n';

        // Take whichever velocity is greatest:
        message += std::to_string(attrs.v_collapse) + " = vCollapse | vMax = " + std::to_string(attrs.v_max);
        if (attrs.v_max < attrs.v_collapse) { attrs.v_max = attrs.v_collapse; }

        if (attrs.v_max < attrs.v_max_prev) {
            updateDTK(attrs.v_max);
            attrs.v_max_prev = attrs.v_max;
            message += "\nk: " + std::to_string(attrs.kin) + "\tdt: " + std::to_string(attrs.dt) + '\n';
        }
    }
    MPIsafe_print(std::cerr,message);

    message = "";

    if (Step == 0 or dtOld < 0) {
        attrs.steps = static_cast<unsigned long long>(attrs.simTimeSeconds / attrs.dt) + 1;
        // std::cout<<simTimeSeconds / dt - steps*1.0<<std::endl;
        // if (simTimeSeconds / dt == steps) //There is one too few writes in the sim if this is true
        // {
        //     std::cout<<"IT HAPPENED, numparts: "<<num_particles<<std::endl;
        //     steps += 1;
        // }
        if (attrs.steps < 0)
        {
            message += "ERROR: STEPS IS NEGATIVE.\n";
            message += "simTimeSeconds/dt = " + std::to_string(attrs.simTimeSeconds / attrs.dt)+'\n';
            message += "casted simTimeSeconds/dt (steps) = " + std::to_string(static_cast<int>(attrs.simTimeSeconds / attrs.dt))+'\n';
            message += "Exiting program now.\n";
            MPIsafe_print(std::cerr,message);
            exit(-1);
        }

        message += "\tInitial Steps: " + std::to_string(attrs.steps) + '\n';
    } else {
        attrs.steps = static_cast<unsigned long long>(dtOld / attrs.dt) * (attrs.steps - Step) + Step;
        if (attrs.steps < 0)
        {
            message += "ERROR: STEPS IS NEGATIVE.\n";
            message += "dtOld/dt = " + std::to_string(dtOld / attrs.dt) + '\n';
            message += "(steps - Step) = " + std::to_string(attrs.steps - Step) + '\n';
            message += "Step = " + std::to_string(Step) + '\n';
            message += "Final steps = " + std::to_string(static_cast<unsigned long long>(dtOld / attrs.dt) * (attrs.steps - Step) + Step) + '\n';
            message += "Exiting program now.'\n'";
            MPIsafe_print(std::cerr,message);
            exit(-1);
        }
        message += "\tSteps: " + std::to_string(attrs.steps);
    }
    MPIsafe_print(std::cerr,message);

    message = "";

    if (attrs.timeResolution / attrs.dt > 1.) {
        attrs.skip = static_cast<int>(floor(attrs.timeResolution / attrs.dt));
        message += "\tSkip: " + std::to_string(attrs.skip) + '\n';
    } else {
        message += "Desired time resolution is lower than dt. Setting to 1 second per skip.\n";
        attrs.skip = static_cast<int>(floor(1. / attrs.dt));
    }
    MPIsafe_print(std::cerr,message);
}

// todo - make bigger balls favor the middle, or, smaller balls favor the outside.
/// @brief Push balls apart until no overlaps
void Ball_group::pushApart() const
{
    MPIsafe_print(std::cerr,std::string("Separating spheres - Current max overlap:\n"));
    /// Using acc array as storage for accumulated position change.
    int* counter = new int[attrs.num_particles];
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        acc[Ball] = {0, 0, 0};
        counter[Ball] = 0;
    }

    double overlapMax = -1;
    const double pseudoDT = attrs.r_min * .1;
    int step = 0;

    while (true) {
        // if (step % 10 == 0)
        //{
        //  simDataWrite("pushApart_");
        //}

        for (int A = 0; A < attrs.num_particles; A++) {
            for (int B = A + 1; B < attrs.num_particles; B++) {
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

        for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
            if (counter[Ball] > 0) {
                pos[Ball] += acc[Ball].normalized() * pseudoDT;
                acc[Ball] = {0, 0, 0};
                counter[Ball] = 0;
            }
        }

        if (overlapMax > 0) {
            MPIsafe_print(std::cerr,std::string(std::to_string(overlapMax) + "                        \r"));//Why is there a \r here? Keeping until I know
        } else {
            MPIsafe_print(std::cerr,"\nSuccess!\n");
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
    while (position < attrs.initial_radius) {
        // todo - include vdw!!!
        attrs.v_collapse += attrs.G * attrs.m_total / (attrs.initial_radius * attrs.initial_radius) * 0.1;
        position += attrs.v_collapse * 0.1;
    }
    attrs.v_collapse = fabs(attrs.v_collapse);
}

/// get max velocity
[[nodiscard]] double Ball_group::getVelMax()
{
    attrs.v_max = 0;

    // todo - make this a manual set true or false to use soc so we know if it is being used or not.
    if (attrs.soc > 0) {
        int counter = 0;
        for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
            if (vel[Ball].norm() > attrs.v_max) 
            { 
                attrs.v_max = vel[Ball].norm();
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

        MPIsafe_print(std::cerr,'(' + std::to_string(counter) + " spheres ignored"+ ") ");
    } else {
        for (int Ball = 0; Ball < attrs.num_particles; Ball++) {

            if (vel[Ball].norm() > attrs.v_max) 
            { 
                attrs.v_max = vel[Ball].norm();
            }
        }

        // Is vMax for some reason unreasonably small? Don't proceed. Probably a finished sim.
        // This shouldn't apply to extremely destructive collisions because it is possible that no
        // particles are considered, so it will keep pausing.
        if (attrs.v_max < 1e-10) {
            MPIsafe_print(std::cerr,"\nMax velocity in system is less than 1e-10.\n");
            system("pause");
        }
    }

    return attrs.v_max;
}

double Ball_group::get_soc()
{
    return attrs.soc;
}

void Ball_group::calc_helpfuls()
{
    attrs.r_min = getRmin();
    attrs.r_max = getRmax();
    attrs.m_total = getMass();
    attrs.initial_radius = get_radius(getCOM());
    attrs.soc = 4 * attrs.r_max + attrs.initial_radius;
    // soc = -1;
}   

// Kick ballGroup (give the whole thing a velocity)
void Ball_group::kick(const vec3& vec) const
{
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) { vel[Ball] += vec; }
}


vec3 Ball_group::calc_momentum(const std::string& of = "") const
{
    vec3 pTotal = {0, 0, 0};
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) { pTotal += m[Ball] * vel[Ball]; }
    // fprintf(stderr, "%s Momentum Check: %.2e, %.2e, %.2e\n", of.c_str(), pTotal.x, pTotal.y, pTotal.z);
    return pTotal;
}

// offset cluster
void Ball_group::offset(const double& rad1, const double& rad2, const double& impactParam) const
{
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        pos[Ball].x += (rad1 + rad2) * cos(impactParam);
        pos[Ball].y += (rad1 + rad2) * sin(impactParam);
    }
}

/// Approximate the radius of the ballGroup.
[[nodiscard]] double Ball_group::get_radius(const vec3& center) const
{
    double radius = 0;
    if (attrs.num_particles > 1) {
        for (size_t i = 0; i < attrs.num_particles; i++) {
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

    if (attrs.num_particles > 1)  // Code below only necessary for effects between balls.
    {
        for (int A = 1; A < attrs.num_particles; A++) {
            for (int B = 0; B < A; B++) {
                const double sumRaRb = R[A] + R[B];
                const double dist = (pos[A] - pos[B]).norm();
                const double overlap = sumRaRb - dist;

                // Check for collision between Ball and otherBall.
                if (overlap > 0) {
                    PE +=
                        -attrs.G * m[A] * m[B] / dist + attrs.kin * ((sumRaRb - dist) * .5) * ((sumRaRb - dist) * .5);
                } else {
                    PE += -attrs.G * m[A] * m[B] / dist;
                }
            }
        }
    } else  // For the case of just one ball:
    {
        PE = 0;
    }
}

std::string Ball_group::get_data_info()
{
    std::ostringstream out_stream;
    out_stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    
    out_stream<<"steps:"<<attrs.steps<<",skip:"<<attrs.skip;
    out_stream<<",kin:"<<attrs.kin<<",kout:"<<attrs.kout<<",dt:"<<attrs.dt;

    return out_stream.str();
}

void Ball_group::sim_init_write(int counter=0)
{
    std::cerr<<"Sim init write for index: "<<counter<<std::endl;
    init_data(counter);

    // if (counter > 0) { filename.insert(0, std::to_string(counter) + '_'); }

    std::vector<double> constData(data->getWidth("constants")*attrs.num_particles);
    // Write constant data:
    int pt = 0;
    int jump = data->getSingleWidth("constants");
    for (int i = 0; i < attrs.num_particles; i++) 
    {
        constData[pt] = R[i];
        constData[pt+1] = m[i];
        constData[pt+2] = moi[i];
        pt += jump;
    }

    data->Write(constData,"constants"); //THIS MUST COME BEFORE THE NEXT WRITE METADATA otherwise the h5 file wont be initiated
    if (attrs.data_type == 0) //This meta write is for restarting jobs. Only necessary for hdf5
    {
        data->WriteMeta(get_data_info(),attrs.sim_meta_data_name,"constants");
    }


    energyBuffer = std::vector<double> (data->getWidth("energy"));
    energyBuffer[0] = attrs.simTimeElapsed;
    energyBuffer[1] = PE;
    energyBuffer[2] = KE;
    energyBuffer[3] = PE+KE;
    energyBuffer[4] = mom.norm();
    energyBuffer[5] = ang_mom.norm();
    data->Write(energyBuffer,"energy");

    // Reinitialize energies for next step:
    KE = 0;
    PE = 0;
    mom = {0, 0, 0};
    ang_mom = {0, 0, 0};

    // Send position and rotation to buffer:
    ballBuffer = std::vector<double> (data->getWidth("simData"));
    pt = 0;
    jump = data->getSingleWidth("simData");
    for (int i = 0; i < attrs.num_particles; i++) 
    {
        ballBuffer[pt] = pos[i].x;
        ballBuffer[pt+1] = pos[i].y;
        ballBuffer[pt+2] = pos[i].z;
        ballBuffer[pt+3] = w[i].x;
        ballBuffer[pt+4] = w[i].y;
        ballBuffer[pt+5] = w[i].z;
        ballBuffer[pt+6] = w[i].norm();
        ballBuffer[pt+7] = vel[i].x;
        ballBuffer[pt+8] = vel[i].y;
        ballBuffer[pt+9] = vel[i].z;
        ballBuffer[pt+10] = 0;
        pt += jump;
    }
    data->Write(ballBuffer,"simData",1);

    //Initialize ballBuffer and energyBuffer to the size they should be for actual sim
    energyBuffer.clear();
    ballBuffer.clear();

    energyBuffer = std::vector<double> (data->getWidth("energy")*bufferlines);
    ballBuffer = std::vector<double> (data->getWidth("simData")*bufferlines);

    //initialize num_writes
    attrs.num_writes = 0;

    std::cerr << "\nSimulating " << attrs.steps * attrs.dt / 60 / 60 << " hours.\n";
    std::cerr << "Total mass: " << attrs.m_total << '\n';
    std::cerr << "\n===============================================================\n";

}


[[nodiscard]] vec3 Ball_group::getCOM() const
{
    if (attrs.m_total > 0) {
        vec3 comNumerator;
        for (int Ball = 0; Ball < attrs.num_particles; Ball++) { comNumerator += m[Ball] * pos[Ball]; }
        vec3 com = comNumerator / attrs.m_total;
        return com;
    } else {
        std::cerr << "Mass of cluster is zero.\n";
        exit(EXIT_FAILURE);
    }
}

void Ball_group::zeroVel() const
{
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) { vel[Ball] = {0, 0, 0}; }
}

void Ball_group::zeroAngVel() const
{
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) { w[Ball] = {0, 0, 0}; }
}

void Ball_group::to_origin() const
{
    const vec3 com = getCOM();

    for (int Ball = 0; Ball < attrs.num_particles; Ball++) { pos[Ball] -= com; }
}

// Set velocity of all balls such that the cluster spins:
void Ball_group::comSpinner(const double& spinX, const double& spinY, const double& spinZ) const
{
    const vec3 comRot = {spinX, spinY, spinZ};  // Rotation axis and magnitude
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        vel[Ball] += comRot.cross(pos[Ball] - getCOM());
        w[Ball] += comRot;
    }
}

void Ball_group::rotAll(const char axis, const double angle) const
{
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
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
Ball_group Ball_group::spawn_particles(const int count)
{
    // Load file data:
    std::cerr << "Add Particle\n";

    // Random particle to origin
    Ball_group projectile(count);
    // Particle random position at twice radius of target:
    // We want the farthest from origin since we are offsetting form origin. Not com.
    const auto cluster_radius = 3;

    const vec3 projectile_direction = rand_vec3(1).normalized();
    projectile.pos[0] = projectile_direction * (cluster_radius + attrs.scaleBalls * 4);
    projectile.w[0] = {0, 0, 0};
    // Velocity toward origin:
    projectile.vel[0] = -attrs.v_custom * projectile_direction;
    projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
    projectile.m[0] = attrs.density * 4. / 3. * pi * std::pow(projectile.R[0], 3);
    projectile.moi[0] = calc_moi(projectile.R[0], projectile.m[0]);

    const double3x3 local_coords = local_coordinates(to_double3(projectile_direction));
    // to_vec3(local_coords.y).print();
    // to_vec3(local_coords.z).print();
    // projectile.pos[0].print();
    for (int i = 1; i < projectile.attrs.num_particles - 3; i++) {
        const auto rand_y = rand_between(-cluster_radius, cluster_radius);
        const auto rand_z = rand_between(-cluster_radius, cluster_radius);
        // projectile.pos[i] = projectile.pos[0] + perpendicular_shift(local_coords, rand_y, rand_z);
        projectile.pos[i] = projectile.pos[0] + perpendicular_shift(local_coords, rand_y, rand_z);
        // std::cout << rand_y << '\t' << to_vec3(local_coords.y * rand_y) <<'\t'<< projectile.pos[i] <<
        // '\n';
    }
    projectile.pos[projectile.attrs.num_particles - 3] = projectile_direction * 2;
    projectile.pos[projectile.attrs.num_particles - 2] = projectile_direction * 4;
    projectile.pos[projectile.attrs.num_particles - 1] = projectile_direction * 6;

    Ball_group new_group{projectile.attrs.num_particles + attrs.num_particles};

    new_group.merge_ball_group(*this);
    new_group.merge_ball_group(projectile);

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
        for (size_t i = 0; i < attrs.num_particles; i++) {
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
    // projectile.radiiDistribution = radiiDistribution;
    // projectile.radiiFraction = radiiFraction;
    // // projectile.data = data;
    // //carry over folders
    // projectile.project_path = project_path;
    // projectile.output_folder = output_folder;
    // projectile.data_directory = data_directory;
    // projectile.projectileName = projectileName;
    // projectile.targetName = targetName;
    // projectile.output_prefix = output_prefix;

    // projectile.skip = skip;
    // projectile.steps = steps;

    // projectile.dt=dt;
    // projectile.kin=kin;  // Spring constant
    // projectile.kout=kout;
    // Particle random position at twice radius of target:
    // We want the farthest from origin since we are offsetting form origin. Not com.
    const auto cluster_radius = get_radius(vec3(0, 0, 0));

    const vec3 projectile_direction = rand_unit_vec3();
    projectile.pos[0] = projectile_direction * (cluster_radius + attrs.scaleBalls * 4);
    if (attrs.radiiDistribution == attrs.constant)
    {
        // std::cout<<"radiiFraction: "<<radiiFraction<<std::endl;
        projectile.R[0] = attrs.scaleBalls;  //MAKE BOTH VERSIONS SAME
        // projectile.R[0] = scaleBalls/radiiFraction;  //limit of 1.4// rand_between(1,3)*1e-5;
        // std::cout<<"(constant) Particle added with radius of "<<projectile.R[0]<<std::endl;
    }
    else
    {
        projectile.R[0] = lognorm_dist(attrs.scaleBalls*std::exp(-5*std::pow(attrs.lnSigma,2)/2),attrs.lnSigma);
        // std::cout<<"(lognorm) Particle added with radius of "<<projectile.R[0]<<std::endl;
    }
    projectile.w[0] = {0, 0, 0};
    projectile.m[0] = attrs.density * 4. / 3. * pi * std::pow(projectile.R[0], 3);
    // Velocity toward origin:
    if (attrs.temp > 0)
    {
        double a = std::sqrt(Kb*attrs.temp/projectile.m[0]);
        attrs.v_custom = max_bolt_dist(a); 

        std::string message("v_custom set to "+std::to_string(attrs.v_custom)+ "cm/s based on a temp of "+
                std::to_string(attrs.temp)+" degrees K.\n"); 
    }
    projectile.vel[0] = -attrs.v_custom * projectile_direction;

    
    // projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
    projectile.moi[0] = calc_moi(projectile.R[0], projectile.m[0]);

  

    const double3x3 local_coords = local_coordinates(to_double3(projectile_direction));
    
    projectile.pos[0] = dust_agglomeration_offset(local_coords,projectile.pos[0],projectile.vel[0],projectile.R[0]);


    
    return projectile;
}

// Uses previous O as target and adds one particle to hit it:
Ball_group Ball_group::add_projectile()
{
    // Load file data:
    MPIsafe_print(std::cerr,"Add Particle\n");

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

    std::ostringstream oss;
    oss << "\nTarget Velocity: " << std::scientific << vel[0].norm()
        << "\nProjectile Velocity: " << projectile.vel[0].norm() << "\n\n";
    MPIsafe_print(std::cerr,oss.str());

    projectile.calc_momentum("Projectile");
    calc_momentum("Target");
    Ball_group new_group{projectile.attrs.num_particles + attrs.num_particles};

    int new_num_particles = projectile.attrs.num_particles + attrs.num_particles;

    new_group.merge_ball_group(*this);
    new_group.merge_ball_group(projectile);
    new_group.attrs = attrs;
    //The next line is important because the previous line overwrites the value of num_particles set in the Ball_group constructor
    new_group.attrs.num_particles = new_num_particles;

    // Hack - if v_custom is less than 1 there are problems if dt is calibrated to this
    //        if v_custom is greater than 1 you need to calibrate dt to that v_custom
    if (attrs.v_custom < 1)
    {
        new_group.calibrate_dt(0, 1);
    }
    else
    {
        new_group.calibrate_dt(0, attrs.v_custom);
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
    std::memcpy(
        &distances[attrs.num_particles_added], src.distances, sizeof(src.distances[0]) * src.attrs.num_particles);
    std::memcpy(&pos[attrs.num_particles_added], src.pos, sizeof(src.pos[0]) * src.attrs.num_particles);
    std::memcpy(&vel[attrs.num_particles_added], src.vel, sizeof(src.vel[0]) * src.attrs.num_particles);
    std::memcpy(&velh[attrs.num_particles_added], src.velh, sizeof(src.velh[0]) * src.attrs.num_particles);
    std::memcpy(&acc[attrs.num_particles_added], src.acc, sizeof(src.acc[0]) * src.attrs.num_particles);
    std::memcpy(&w[attrs.num_particles_added], src.w, sizeof(src.w[0]) * src.attrs.num_particles);
    std::memcpy(&wh[attrs.num_particles_added], src.wh, sizeof(src.wh[0]) * src.attrs.num_particles);
    std::memcpy(&aacc[attrs.num_particles_added], src.aacc, sizeof(src.aacc[0]) * src.attrs.num_particles);
    std::memcpy(&R[attrs.num_particles_added], src.R, sizeof(src.R[0]) * src.attrs.num_particles);
    std::memcpy(&m[attrs.num_particles_added], src.m, sizeof(src.m[0]) * src.attrs.num_particles);
    std::memcpy(&moi[attrs.num_particles_added], src.moi, sizeof(src.moi[0]) * src.attrs.num_particles);
    

    // Keep track of now loaded ball set to start next set after it:
    attrs.num_particles_added += src.attrs.num_particles;

    // num_particles_added += src.num_particles;
    // radiiDistribution = src.radiiDistribution;
    // radiiFraction = src.radiiFraction;

    // //carry over folders
    // project_path = src.project_path;
    // output_folder = src.output_folder;
    // data_directory = src.data_directory;
    // projectileName = src.projectileName;
    // targetName = src.targetName;
    // output_prefix = src.output_prefix;

    // skip = src.skip;
    // steps = src.steps;

    // dt=src.dt;
    // kin=src.kin;  // Spring constant
    // kout=src.kout;
    // data = src.data;

    calc_helpfuls();
}

/// Allocate balls
void Ball_group::allocate_group(const int nBalls)
{
    attrs.num_particles = nBalls;

    try {
        distances = new double[(attrs.num_particles * attrs.num_particles / 2) - (attrs.num_particles / 2)];

        pos = new vec3[attrs.num_particles];
        vel = new vec3[attrs.num_particles];
        velh = new vec3[attrs.num_particles];
        acc = new vec3[attrs.num_particles];
        w = new vec3[attrs.num_particles];
        wh = new vec3[attrs.num_particles];
        aacc = new vec3[attrs.num_particles];
        R = new double[attrs.num_particles];
        m = new double[attrs.num_particles];
        moi = new double[attrs.num_particles];

        
    } catch (const std::exception& e) {
        std::cerr << "Failed trying to allocate group. " << e.what() << '\n';
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
    delete[] w;
    delete[] wh;
    delete[] aacc;
    delete[] R;
    delete[] m;
    delete[] moi;
    #ifdef GPU_ENABLE
        delete[] aaccsq;
        delete[] accsq;
    #endif
    // delete data;
    
}


// Initialize accelerations and energy calculations:
void Ball_group::init_conditions()
{
    // SECOND PASS - Check for collisions, apply forces and torques:
    for (int A = 1; A < attrs.num_particles; A++)  // cuda
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
                k = attrs.kin;
                // Apply coefficient of restitution to balls leaving collision.
                // if (dist >= oldDist) {
                //     k = kout;
                // } else {
                //     k = kin;
                // }

                // Cohesion (in contact) h must always be h_min:
                // constexpr double h = h_min;
                const double h = attrs.h_min;
                const double Ra = R[A];
                const double Rb = R[B];
                const double h2 = h * h;
                // constexpr double h2 = h * h;
                const double twoRah = 2 * Ra * h;
                const double twoRbh = 2 * Rb * h;
                const vec3 vdwForceOnA =
                    attrs.Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                    ((h + Ra + Rb) /
                     ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                      (h2 + twoRah + twoRbh + 4 * Ra * Rb) * (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                    rVecab.normalized();

                // Elastic force:
                const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

                // Gravity force:
                const vec3 gravForceOnA = (attrs.G * m[A] * m[B] / (dist * dist)) * (rVecab / dist);

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
                    slideForceOnA = attrs.u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                }

                // Compute rolling friction force:
                const double w_diff_mag = w_diff.norm();
                if (w_diff_mag > 1e-13)  // Divide by zero protection.
                {
                    rollForceA =
                        -attrs.u_r * elastic_force_A_mag * (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
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
                    -attrs.Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                PE += U_vdw + 0.5 * k * overlap * overlap;

            } else  // Non-contact forces:
            {
                // No collision: Include gravity and vdw:
                // const vec3 gravForceOnA = (G * m[A] * m[B] / (dist * dist)) * (rVecab / dist);

                // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic
                // cancellation:
                double h = std::fabs(overlap);
                if (h < attrs.h_min)  // If h is closer to 0 (almost touching), use hmin.
                {
                    h = attrs.h_min;
                }
                const double Ra = R[A];
                const double Rb = R[B];
                const double h2 = h * h;
                const double twoRah = 2 * Ra * h;
                const double twoRbh = 2 * Rb * h;
                const vec3 vdwForceOnA =
                    attrs.Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
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
                    -attrs.Ha / 6 *
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
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        KE += .5 * m[Ball] * vel[Ball].dot(vel[Ball]) + .5 * moi[Ball] * w[Ball].dot(w[Ball]);
        mom += m[Ball] * vel[Ball];
        ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
    }
}

[[nodiscard]] double Ball_group::getRmin()
{
    attrs.r_min = R[0];
    for (int Ball = 1; Ball < attrs.num_particles; Ball++) {
        if (R[Ball] < attrs.r_min) { attrs.r_min = R[Ball]; }
    }
    return attrs.r_min;
}

[[nodiscard]] double Ball_group::getRmax()
{
    attrs.r_max = R[0];
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        if (R[Ball] > attrs.r_max) { attrs.r_max = R[Ball]; }
    }
    return attrs.r_max;
}


[[nodiscard]] double Ball_group::getMassMax() const
{
    double mMax = m[0];
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        if (m[Ball] > mMax) { mMax = m[Ball]; }
    }
    return mMax;
}


void Ball_group::parseSimData(std::string line)
{
    std::string lineElement;
    // Get number of balls in file
    int count = 54 / attrs.properties;
    attrs.num_particles = static_cast<int>((static_cast<int>(std::count(line.begin(), line.end(), ',')) + 1)/11);
    if (attrs.num_particles > 0)
    {
        count = attrs.num_particles;
    }
    // int count = std::count(line.begin(), line.end(), ',') / properties + 1;
    allocate_group(count);
    std::stringstream chosenLine(line);  // This is the last line of the read file, containing all data
                                         // for all balls at last time step
    // Get position and angular velocity data:
    for (int A = 0; A < attrs.num_particles; A++) {
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
        for (int i = 0; i < attrs.properties - 10; i++)  // We used 10 elements. This skips the rest.
        {
            std::getline(chosenLine, lineElement, ',');
        }
    }
}

/// Get previous sim constants by filename.
void Ball_group::loadConsts(const std::string& path, const std::string& filename)
{
    // Get radius, mass, moi:
    std::string constantsFilename = path + filename + "constants.csv";
    if (auto ConstStream = std::ifstream(constantsFilename, std::ifstream::in)) {
        std::string line, lineElement;
        for (int A = 0; A < attrs.num_particles; A++) {
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
        MPIsafe_print(std::cerr,"Could not open constants file: " + constantsFilename + " ... Exiting program.\n");
        MPIsafe_exit(EXIT_FAILURE);
    }
    // getMass();
}


//This used to be [[nodiscard]] static std::string ... but wont compile outside the actual class definition
/// Get last line of previous simData by filename.
[[nodiscard]] std::string Ball_group::getLastLine(const std::string& path, const std::string& filename)
{
    std::string simDataFilepath = path + filename + "simData.csv";

    if (auto simDataStream = std::ifstream(simDataFilepath, std::ifstream::in)) {
        MPIsafe_print(std::cerr,"\nParsing last line of data.\n");

        simDataStream.seekg(-1, std::ios_base::end);  // go to 
         // spot before the EOF

        bool keepLooping = true;
        bool first_run = true;
        while (keepLooping) {
            char ch = ' ';
            simDataStream.get(ch);  // Get current byte's data

            if (static_cast<int>(simDataStream.tellg()) <=
                1) {                     // If the data was at or before the 0th byte
                simDataStream.seekg(0);  // The first line is the last line
                keepLooping = false;     // So stop there
            } else if (ch == '\n' && not first_run) {     // If the data was a newline
                keepLooping = false;     // Stop at the current position (if we arent on the first character).
            } else {                     // If the data was neither a newline nor at the 0 byte
                simDataStream.seekg(-2, std::ios_base::cur);  // Move to the front of that data, then to
                                                              // the front of the data before it
            }
            first_run = false;
        }
        std::string line;
        std::getline(simDataStream, line);  // Read the current line
        return line;
    } else {

        std::string message("Could not open simData file: "+simDataFilepath+"... Exiting program.\n");
        MPIsafe_print(std::cerr,message);
        MPIsafe_exit(EXIT_FAILURE);
        return "ERROR"; //This shouldn't return but not returning anything is giving a warning
    }
}

// void Ball_group::simDataWrite(std::string& outFilename)
// {
//     // todo - for some reason I need checkForFile instead of just using ballWrite. Need to work out why.
//     // Check if file name already exists. If not, initialize
//     std::ifstream checkForFile;
//     checkForFile.open(output_folder + outFilename + "simData.csv", std::ifstream::in);
//     if (checkForFile.is_open() == false) {
//         sim_init_write(outFilename);
//     } else {
//         ballBuffer << '\n';  // Prepares a new line for incoming data.

//         for (int Ball = 0; Ball < num_particles; Ball++) {
//             // Send positions and rotations to buffer:
//             if (Ball == 0) {
//                 ballBuffer << pos[Ball][0] << ',' << pos[Ball][1] << ',' << pos[Ball][2] << ','
//                            << w[Ball][0] << ',' << w[Ball][1] << ',' << w[Ball][2] << ','
//                            << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y << ','
//                            << vel[Ball].z << ',' << 0;
//             } else {
//                 ballBuffer << ',' << pos[Ball][0] << ',' << pos[Ball][1] << ',' << pos[Ball][2] << ','
//                            << w[Ball][0] << ',' << w[Ball][1] << ',' << w[Ball][2] << ','
//                            << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y << ','
//                            << vel[Ball].z << ',' << 0;
//             }
//         }

//         // Write simData to file and clear buffer.
//         std::ofstream ballWrite;
//         ballWrite.open(output_folder + outFilename + "simData.csv", std::ofstream::app);
//         ballWrite << ballBuffer.rdbuf();  // Barf buffer to file.
//         ballBuffer.str("");               // Resets the stream for that balls to blank.
//         ballWrite.close();
//     }
//     checkForFile.close();
// }


[[nodiscard]] double Ball_group::getMass()
{
    attrs.m_total = 0;
    {
        for (int Ball = 0; Ball < attrs.num_particles; Ball++) { attrs.m_total += m[Ball]; }
    }
    return attrs.m_total;
}

void Ball_group::threeSizeSphere(const int nBalls)
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

        R[Ball] = 3. * attrs.scaleBalls;
        m[Ball] = attrs.density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(attrs.spaceRange);
    }

    for (int Ball = larges; Ball < (larges + mediums); Ball++) {
        R[Ball] = 2. * attrs.scaleBalls;  // std::pow(1. / (double)nBalls, 1. / 3.) * 2. * scaleBalls;
        m[Ball] = attrs.density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(attrs.spaceRange);
    }
    for (int Ball = (larges + mediums); Ball < nBalls; Ball++) {
        R[Ball] = 1. * attrs.scaleBalls;  // std::pow(1. / (double)nBalls, 1. / 3.) * 1. * scaleBalls;
        m[Ball] = attrs.density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(attrs.spaceRange);
    }

    attrs.m_total = 0;
    for (int i = 0; i < nBalls; i++)
    {
        attrs.m_total += m[i];
        std::cerr<<"Ball "<<i<<"\tmass is "<<m[i]<<"\t"<<"radius is "<<R[i]<<std::endl;
    }

    std::cerr << "Smalls: " << smalls << " Mediums: " << mediums << " Larges: " << larges << '\n';

    // Generate non-overlapping spherical particle field:
    int collisionDetected = 0;
    int oldCollisions = nBalls;

    for (int failed = 0; failed < attrs.attempts; failed++) {
        for (int A = 0; A < nBalls; A++) {
            for (int B = A + 1; B < nBalls; B++) {
                // Check for Ball overlap.
                const double dist = (pos[A] - pos[B]).norm();
                const double sumRaRb = R[A] + R[B];
                const double overlap = dist - sumRaRb;
                if (overlap < 0) {
                    collisionDetected += 1;
                    // Move the other ball:
                    pos[B] = rand_vec3(attrs.spaceRange);
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
        if (failed == attrs.attempts - 1 ||
            collisionDetected >
                static_cast<int>(
                    1.5 *
                    static_cast<double>(
                        nBalls)))  // Added the second part to speed up spatial constraint increase when
                                   // there are clearly too many collisions for the space to be feasible.
        {
            std::cerr << "Failed " << attrs.spaceRange << ". Increasing range " << attrs.spaceRangeIncrement
                      << "cm^3.\n";
            attrs.spaceRange += attrs.spaceRangeIncrement;
            failed = 0;
            for (int Ball = 0; Ball < nBalls; Ball++) {
                pos[Ball] = rand_vec3(
                    attrs.spaceRange);  // Each time we fail and increase range, redistribute all balls randomly
                                  // so we don't end up with big balls near mid and small balls outside.
            }
        }
        collisionDetected = 0;
    }

    std::cerr << "Final spacerange: " << attrs.spaceRange << '\n';
    std::cerr << "m_total: " << attrs.m_total << '\n';
    std::cerr << "Initial Radius: " << get_radius(getCOM()) << '\n';
    std::cerr << "Mass: " << getMass() << '\n';
}

void Ball_group::generate_ball_field(const int nBalls)
{
    MPIsafe_print(std::cerr,"CLUSTER FORMATION (with "+std::to_string(nBalls)+" balls)\n");

    allocate_group(nBalls);

    // Create new random number set.
        //This should be d
         // in parse_input_file
    // const int seedSave = static_cast<int>(time(nullptr));
    // srand(seedSave);
    if (attrs.radiiDistribution == attrs.constant)
    {
        oneSizeSphere(nBalls);
    }
    else
    {
        distSizeSphere(nBalls);
    }
    
    calc_helpfuls();
    // threeSizeSphere(nBalls);

    attrs.output_prefix = std::to_string(nBalls) + "_R" + scientific(get_radius(getCOM())) + "_v" +
                    scientific(attrs.v_custom) + "_cor" + rounder(sqrtf(attrs.cor), 4) + "_mu" + rounder(attrs.u_s, 3) +
                    "_rho" + rounder(attrs.density, 4);
}

/// Make ballGroup from file data.
void Ball_group::loadSim(const std::string& path, const std::string& filename)
{
    std::string file = filename;
    //file we are loading is csv file
    size_t _pos;
    int file_index;

    if (file.substr(file.size()-4,file.size()) == ".csv")
    {
        //decrease index by 1 so we have most recent finished sim
        _pos = file.find_first_of("_");
        size_t _lastpos = file.find_last_of("_");
        
        file_index = stoi(file.substr(0,_pos));

        file = std::to_string(file_index) + file.substr(_pos,_lastpos-(_pos-1));
        attrs.start_index = file_index+1;//shouldnt be file_index-1 because that is just the one we read, we will write to the next index

        parseSimData(getLastLine(path, file));
        loadConsts(path, file);
    }
    else if (file.substr(file.size()-3,file.size()) == ".h5")
    {
        #ifdef HDF5_ENABLE
            _pos = file.find_first_of("_");
            file_index = stoi(file.substr(0,_pos));
            loadDatafromH5(path,file);
        #else
            MPIsafe_print(std::cerr,"ERROR: HDF5 not enabled. Please recompile with -DHDF5_ENABLE and try again.\n");
            MPIsafe_exit(EXIT_FAILURE);
        #endif
    }
    else
    {
        MPIsafe_print(std::cerr,"ERROR: filename in loadSim is of unknown type.\n");
        MPIsafe_exit(EXIT_FAILURE);
    }


    calc_helpfuls();

    std::string message("Balls: " + std::to_string(attrs.num_particles) + '\n' + 
                        "Mass: " + dToSci(attrs.m_total) + '\n' +
                        "Approximate radius: " + dToSci(attrs.initial_radius) + " cm.\n");
    MPIsafe_print(std::cerr,message);
}

void Ball_group::parse_meta_data(std::string metadata)
{
    std::string subdata,data_t,intstr;
    size_t comma_pos,colon_pos;
    bool run = true;    

    while (run)
    {
        comma_pos = metadata.find_first_of(",");      
        subdata = metadata.substr(0,comma_pos);
        colon_pos = subdata.find_first_of(":");

        data_t = subdata.substr(0,colon_pos);
        intstr = subdata.substr(colon_pos+1,subdata.length());

        if (data_t == "steps")
        {
            attrs.steps = stoi(intstr);
        }
        else if (data_t == "skip")
        {
            attrs.skip = stoi(intstr);
        }
        else if (data_t == "kin")
        {
            std::istringstream in_stream(intstr);
            double retrieved_double;
            in_stream >> retrieved_double;
            attrs.kin = retrieved_double;
        }
        else if (data_t == "kout")
        {
            std::istringstream in_stream(intstr);
            double retrieved_double;
            in_stream >> retrieved_double;
            attrs.kout = retrieved_double;
        }
        else if (data_t == "dt")
        {
            std::istringstream in_stream(intstr);
            double retrieved_double;
            in_stream >> retrieved_double;
            attrs.dt = retrieved_double;
        }
        else
        {
            std::cerr<<"DECCO ERROR: sim metadata '"<<data_t<<"' doesn't exist."<<std::endl;
            exit(EXIT_FAILURE);
        }


        metadata = metadata.substr(comma_pos+1,metadata.length());

        if (comma_pos == std::string::npos)
        {
            run = false;
        }
    }

}




#ifdef HDF5_ENABLE
void Ball_group::loadDatafromH5(std::string path,std::string file)
{
    //read metadata to determine steps and skip variables
    std::string meta = HDF5Handler::readMetadataFromDataset("constants",path+file,attrs.sim_meta_data_name);
    size_t _pos = file.find_first_of("_");
    int file_index = stoi(file.substr(0,_pos));
    bool has_meta = true;
    //If this error happens then then we cannot restart from midway through a sim.
    //This is because the metadata containing the info needed was missing for somereason
  
    if (meta == ERR_RET_MET)  
    {
        has_meta = false;
        //If the highest sim is not finished, we need to load up the previous one and delete the partially completed sim
        if (!HDF5Handler::sim_finished(path,file))
        {
            std::string rmfile = file;

            #ifdef MPI_ENABLE
                MPI_Barrier(MPI_COMM_WORLD);
                
                int status;
                int send_result;
                //If multiple nodes, we don't want to delete until everyone has loaded
                if (getRank() == 0)
                {
                    status = remove(rmfile.c_str());
                    if (getSize() > 1)
                    {
                        for (int i = 1; i < getSize(); i++)
                        {
                            send_result = MPI_Send(&status, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                            if (send_result != MPI_SUCCESS)
                            {
                                std::cerr<<"ERROR: MPI_Send to node "<<i<<" errored with code "<<send_result<<std::endl;   
                                MPIsafe_exit(-1);
                            }
                        }

                    }
                }
                else
                {
                    MPI_Status mpistat;
                    MPI_Recv(&status, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &mpistat);
                    //verify Recv worked
                    if (mpistat.MPI_ERROR != MPI_SUCCESS)
                    {
                        std::cerr<<"ERROR: MPI_Recv for node "<<getRank()<<" errored with code "<<mpistat.MPI_ERROR<<std::endl;   
                        MPIsafe_exit(-1);
                    }
                }
            #else
                int status = remove(rmfile.c_str());
            #endif
            

            if (status != 0)
            {
                std::string message("File: '"+rmfile+"' could not be removed, now exiting with failure.\n");
                MPIsafe_print(std::cerr,message);
                MPIsafe_exit(EXIT_FAILURE);
            }
            file_index--;
            file = std::to_string(file_index) + file.substr(_pos,file.size());

        }
    }

    //This needs to be here because its used in the following function
    attrs.start_index = file_index;
    
    allocate_group(HDF5Handler::get_num_particles(path,file));
    
    //Load constants because this can be done without an initialized instance of DECCOData
    HDF5Handler::loadConsts(path,file,R,m,moi);

    int writes;
    if (attrs.typeSim != attrs.relax && has_meta) //Relax jobs should not read in the metadata for dt, steps, etc. That is for restarting jobs.
    {
        parse_meta_data(meta);

        //Now we have all info we need to initialze an instance of DECCOData.
        //However, data_written_so_far needs to be determined and set since this is a restart.
        //All this happens in the next two functions. 
        init_data(attrs.start_index);
        //writes is 0 if there is no writes so far (I don't think this should happen but if it does, more stuff needs to happen).
        //writes is >0 then that is how many writes there have been.
        //writes is -1 if there are writes and the sim is already finished. 
        writes = data->setWrittenSoFar(path,file);
    }
    else
    {

        // init_data(attrs.start_index);
        writes = -1;
    }
    // if (writes == 0)//This should really never happen. If it did then there is an empty h5 file
    // {
    //     std::cerr<<"not implimented"<<std::endl;
    //     exit(-1);
    // }
    if(writes > 0) //Works
    {
        //This cannot be done without an instance of DECCOData, that is why these are different than loadConsts
        data -> loadSimData(path,file,pos,w,vel);


        //initiate buffers since we won't call sim_init_write on a restart
        energyBuffer = std::vector<double> (data->getWidth("energy")*bufferlines);
        ballBuffer = std::vector<double> (data->getWidth("simData")*bufferlines);
        
        MPIsafe_print(std::cerr,"mid_sim_restart\n");
        attrs.mid_sim_restart = true;
        attrs.start_step = attrs.skip*(writes-1)+1;
        attrs.start_index++;
    }
    else if(writes == -1) //Works
    {
        if (attrs.typeSim != attrs.relax && has_meta)
        {
            data -> loadSimData(path,file,pos,w,vel);
        }
        else
        {
            HDF5Handler::loadh5SimData(path,file,pos,w,vel);
        }
        attrs.start_index++;
    }
    else
    {
        MPIsafe_print(std::cerr,"ERROR: in setWrittenSoFar() output of value '"+std::to_string(writes)+"'.\n");
        MPIsafe_exit(EXIT_FAILURE);
    }
     
}
#endif

void Ball_group::distSizeSphere(const int nBalls)
{
    for (int Ball = 0; Ball < nBalls; Ball++) {
        R[Ball] = lognorm_dist(attrs.scaleBalls*std::exp(-5*std::pow(attrs.lnSigma,2)/2),attrs.lnSigma);
        m[Ball] = attrs.density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(attrs.spaceRange);
    }

    attrs.m_total = getMass();

    placeBalls(nBalls);
}

void Ball_group::oneSizeSphere(const int nBalls)
{
    for (int Ball = 0; Ball < nBalls; Ball++) {
        R[Ball] = attrs.scaleBalls;
        m[Ball] = attrs.density * 4. / 3. * 3.14159 * std::pow(R[Ball], 3);
        moi[Ball] = .4 * m[Ball] * R[Ball] * R[Ball];
        w[Ball] = {0, 0, 0};
        pos[Ball] = rand_vec3(attrs.spaceRange);
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

    attrs.m_total = getMass();

    placeBalls(nBalls);
}

void Ball_group::placeBalls(const int nBalls)
{
    // Generate non-overlapping spherical particle field:
    int collisionDetected = 0;
    int oldCollisions = nBalls;

    if (nBalls == 1)
    {
        pos[0] = {0,1e-5,0};
    }

    for (int failed = 0; failed < attrs.attempts; failed++) {
        for (int A = 0; A < nBalls; A++) {
            for (int B = A + 1; B < nBalls; B++) {
                // Check for Ball overlap.
                const double dist = (pos[A] - pos[B]).norm();
                const double sumRaRb = R[A] + R[B];
                const double overlap = dist - sumRaRb;
                if (overlap < 0) {
                    collisionDetected += 1;
                    // Move the other ball:
                    pos[B] = rand_vec3(attrs.spaceRange);
                }
            }
        }
        if (collisionDetected < oldCollisions) {
            oldCollisions = collisionDetected;
            MPIsafe_print(std::cerr,"Collisions: "+std::to_string(collisionDetected)+'\n');
        }
        if (collisionDetected == 0) {
            MPIsafe_print(std::cerr,"Success!\n");
            break;
        }
        if (failed == attrs.attempts - 1 ||
            collisionDetected >
                static_cast<int>(
                    1.5 *
                    static_cast<double>(
                        nBalls)))  // Added the second part to speed up spatial constraint increase when
                                   // there are clearly too many collisions for the space to be feasible.
        {

            std::string message("Failed "+std::to_string(attrs.spaceRange)+". Increasing range "+std::to_string(attrs.spaceRangeIncrement)+"cm^3.\n");
            MPIsafe_print(std::cerr,message);
            attrs.spaceRange += attrs.spaceRangeIncrement;
            failed = 0;
            for (int Ball = 0; Ball < nBalls; Ball++) {
                pos[Ball] = rand_vec3(
                    attrs.spaceRange);  // Each time we fail and increase range, redistribute all balls randomly
                                  // so we don't end up with big balls near mid and small balls outside.
            }
        }
        collisionDetected = 0;
    }

    std::string message("Final spacerange: " + std::to_string(attrs.spaceRange)+'\n' +
                        "Initial Radius: "+std::to_string(get_radius(getCOM()))+'\n' +
                        "Mass: "+std::to_string(attrs.m_total)+'\n');
    MPIsafe_print(std::cerr,message);
}




void Ball_group::updateDTK(const double& velocity)
{
    calc_helpfuls();
    attrs.kin = attrs.kConsts * attrs.r_max * velocity * velocity;
    attrs.kout = attrs.cor * attrs.kin;
    const double h2 = attrs.h_min * attrs.h_min;
    const double four_R_min = 4 * attrs.r_min * attrs.h_min;
    const double vdw_force_max = attrs.Ha / 6 * 64 * attrs.r_min * attrs.r_min * attrs.r_min * attrs.r_min * attrs.r_min * attrs.r_min *
                                 ((attrs.h_min + attrs.r_min + attrs.r_min) / ((h2 + four_R_min) * (h2 + four_R_min) *
                                                             (h2 + four_R_min + 4 * attrs.r_min * attrs.r_min) *
                                                             (h2 + four_R_min + 4 * attrs.r_min * attrs.r_min)));
    // todo is it rmin*rmin or rmin*rmax
    const double elastic_force_max = attrs.kin * attrs.maxOverlap * attrs.r_min;
    const double regime = (vdw_force_max > elastic_force_max) ? vdw_force_max : elastic_force_max;
    const double regime_adjust = regime / (attrs.maxOverlap * attrs.r_min);

    // dt = .02 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min);
    attrs.dt = .01 * sqrt((attrs.fourThirdsPiRho / regime_adjust) * attrs.r_min * attrs.r_min * attrs.r_min); //NORMAL ONE
    // dt = .005 * sqrt((fourThirdsPiRho / regime_adjust) * r_min * r_min * r_min);
    std::stringstream message;
    message << "==================" << '\n';
    message << "dt set to: " << attrs.dt << '\n';
    message << "kin set to: " << attrs.kin << '\n';
    message << "kout set to: " << attrs.kout << '\n';
    message << "h_min set to: " << attrs.h_min << '\n';
    message << "Ha set to: " << attrs.Ha << '\n';
    message << "u_s set to: " << attrs.u_s << '\n';
    message << "u_r set to: " << attrs.u_r << '\n';
    if (vdw_force_max > elastic_force_max)
    {
        message << "In the vdw regime.\n";
    }
    else
    {
        message << "In the elastic regime.\n";
    }
    message << "==================" << std::endl;
    MPIsafe_print(std::cerr,message.str());
}


void Ball_group::simInit_cond_and_center(bool add_prefix)
{
    std::string message("==================\ndt: "
                        + dToSci(attrs.dt) + '\n'
                        + "k : " + std::to_string(attrs.kin) + '\n'
                        + "Skip: " + std::to_string(attrs.skip) + '\n'
                        + "Steps: " + std::to_string(attrs.steps) + '\n'
                        + "==================\n");

    MPIsafe_print(std::cerr,message);

    if (attrs.num_particles > 1)
    {
        to_origin();
    }

    calc_momentum("After Zeroing");  // Is total mom zero like it should be?

    // Compute physics between all balls. Distances, collision forces, energy totals, total mass:
    init_conditions();

    // Name the file based on info above:
    if (add_prefix)
    {   
        attrs.output_prefix += "_k" + scientific(attrs.kin) + "_Ha" + scientific(attrs.Ha) + "_dt" + scientific(attrs.dt) + "_";
    }
}





// Set's up a two cluster collision.
void Ball_group::sim_init_two_cluster(
    const std::string& path,
    const std::string& projectileName,
    const std::string& targetName)
{
    // Load file data:
    std::string message("TWO CLUSTER SIM\nFile 1: " + projectileName + "\tFile 2: " + targetName + '\n');
    MPIsafe_print(std::cerr,message);

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

    attrs.num_particles = projectile.attrs.num_particles + target.attrs.num_particles;
    
    MPIsafe_print(std::cerr,"Total number of particles in sim: "+std::to_string(attrs.num_particles) + '\n');

    // DO YOU WANT TO STOP EVERYTHING?
    // projectile.zeroAngVel();
    // projectile.zeroVel();
    // target.zeroAngVel();
    // target.zeroVel();


    // Calc info to determined cluster positioning and collisions velocity:
    projectile.updateGPE();
    target.updateGPE();

    projectile.offset(
        projectile.attrs.initial_radius, target.attrs.initial_radius + target.getRmax() * 2, attrs.impactParameter);

    //      const double PEsys = projectile.PE + target.PE + (-G * projectile.mTotal * target.mTotal /
    //(projectile.getCOM() - target.getCOM()).norm());

    // Collision velocity calculation:
    const double mSmall = projectile.attrs.m_total;
    const double mBig = target.attrs.m_total;
    //      const double mTot = mBig + mSmall;
    // const double vSmall = -sqrt(2 * KEfactor * fabs(PEsys) * (mBig / (mSmall * mTot))); // Negative
    // because small offsets right.
    const double vSmall = -attrs.v_custom;                // DART probe override.
    const double vBig = -(mSmall / mBig) * vSmall;  // Negative to oppose projectile.
    // const double vBig = 0; // Dymorphous override.

    if (std::isnan(vSmall) || std::isnan(vBig)) {
        MPIsafe_print(std::cerr,"A VELOCITY WAS NAN!!!!!!!!!!!!!!!!!!!!!!\n\n");
        MPIsafe_exit(EXIT_FAILURE);
    }

    projectile.kick(vec3(vSmall, 0, 0));
    target.kick(vec3(vBig, 0, 0));

    std::ostringstream oss;
    oss << "\nTarget Velocity: " << std::scientific << vBig
        << "\nProjectile Velocity: " << vSmall << "\n\n";
    MPIsafe_print(std::cerr,oss.str());
    //This is jobs line. Keeping it in case this new printing fails
    // fprintf(message, "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n", vBig, vSmall);

    projectile.calc_momentum("Projectile");
    target.calc_momentum("Target");

    allocate_group(projectile.attrs.num_particles + target.attrs.num_particles);

    merge_ball_group(target);
    merge_ball_group(projectile);  // projectile second so smallest ball at end and largest ball at front
                                   // for dt/k calcs.

    attrs.output_prefix = projectileName + targetName + "T" + rounder(attrs.KEfactor, 4) + "_vBig" +
                    scientific(vBig) + "_vSmall" + scientific(vSmall) + "_IP" +
                    rounder(attrs.impactParameter * 180 / 3.14159, 2) + "_rho" + rounder(attrs.density, 4);
}

// @brief checks if this is new job or restart.
// @returns 0 if this is starting from scratch
// @returns 1 if this is a restart
// @returns 2 if this job is already finished
int Ball_group::check_restart(std::string folder) 
{
    std::string file;
    // int tot_count = 0;
    // int file_count = 0;
    int largest_file_index = -1;
    int file_index=0;
    for (const auto & entry : fs::directory_iterator(folder))
    {
        file = entry.path();
        size_t pos = file.find_last_of("/");
        file = file.erase(0,pos+1);
        
        if (file.substr(0,file.size()-4) == "timing")
        {
            return 2;
        }

        //Is the data in csv format?
        if (file.substr(file.size()-4,file.size()) == ".csv")
        {
            // file_count++;
            size_t _pos = file.find_first_of("_");
            size_t _secpos = file.substr(_pos+1,file.size()).find_first_of("_");
            _secpos += _pos+1; //add 1 to account for _pos+1 in substr above
            file_index = stoi(file.substr(0,file.find_first_of("_")));
            if (file[_pos+1] == 'R')
            {
                file_index = 0;
            }

            if (file_index > largest_file_index)
            {
                largest_file_index = file_index;
            }
        }
        else if (file.substr(file.size()-3,file.size()) == ".h5")
        {
            size_t _pos = file.find_first_of("_");
            file_index = stoi(file.substr(0,file.find_first_of("_")));
            if (file_index > largest_file_index)
            {
                largest_file_index = file_index;
            }
        }
    }
    
    if (largest_file_index > -1)
    {
        return 1;
    }
    else
    {
        return 0;
    }

    
}

//with a known slope and intercept, givin N, the number of particles, what is the 
//optimum number of threads. The function then chooses the power of 2 that is closest
//to this optimum
int Ball_group::get_num_threads()
{
    int N = attrs.num_particles;
    // //This is from speed tests on COSINE
    // double slope = ;
    // double intercept = ;

    // double interpolatedValue = slope * n + intercept; // Linear interpolation
    // return std::min(closestPowerOf2(interpolatedValue),attrs.MAXOMPthreads);        // Find the closest power of 2

    //I could only test up to 16 threads so far. Not enough data for linear interp
    

    int threads;
    // if (N < 0)
    // {
    //     std::cerr<<"ERROR: negative number of particles."<<std::endl;
    //     exit(-1);
    // }
    // else if (N < 80)
    // {
    //     threads = 1;
    // }
    // else if (N < 100)
    // {
    //     threads = 2;
    // }
    // else
    // {
    //     threads = 16;
    // }

    // if (threads > attrs.MAXOMPthreads)
    // {
        threads = attrs.MAXOMPthreads;
    // }
    return threads;
}


std::string Ball_group::find_file_name(std::string path,int index)
{
    std::string file;
    const std::string simDatacsv = "simData.csv";
    const std::string datah5 = "data.h5";
    int file_index;

    for (const auto & entry : fs::directory_iterator(path))
    {

        file = entry.path();
        size_t slash_pos = file.find_last_of("/");
        file = file.erase(0,slash_pos+1);
        size_t _pos = file.find_first_of("_");

        if (_pos != std::string::npos) // only go in here if we found a data file
        {
            //Is the data in csv format? (but first verify the call to substr wont fail)
            if (file.size() >= simDatacsv.size() && file.substr(file.size()-simDatacsv.size(),file.size()) == simDatacsv)
            {
                int num_ = std::count(file.begin(),file.end(),'_');
                file_index = stoi(file.substr(0,file.find_first_of("_")));
                if (num_ > 1) // old name convention
                {
                    if (index == 0)
                    {
                        if (file[_pos+1] == 'R')
                        {
                            return path+file;
                        }
                    }
                    else if (index > 0)
                    {
                        if (file_index == index)
                        {
                            return path+file;
                        }
                    }
                }
                else if (num_ == 1) // new name convention ONLY TESTED WITH THIS CASE
                {
                    if (file_index == index)
                    {
                        return path+file;
                    }
                }
                else
                {
                    MPIsafe_print(std::cerr,"ERROR: filename convention is not recognized for file '"+file+"'\n");
                    exit(-1);
                }

            }
            else if (file.size() >= datah5.size() && file.substr(file.size()-datah5.size(),file.size()) == datah5)
            {
                file_index = stoi(file.substr(0,file.find_first_of("_")));
                if (file_index == index)
                {
                    return path+file;
                }
            }
        }
    }
    
    MPIsafe_print(std::cerr,"ERROR: file at path '"+path+"' with index '"+std::to_string(index)+"' not found. Now exiting . . .\n");
    exit(-1);
}


std::string Ball_group::find_restart_file_name(std::string path)
{
    std::string file;
    std::string largest_file_name;
    std::string second_largest_file_name;
    std::string simDatacsv = "simData.csv";
    std::string datah5 = "data.h5";

    int largest_file_index = -1;
    int second_largest_file_index = -1;
    int file_index=0;
    bool csv = false;
    for (const auto & entry : fs::directory_iterator(path))
    {
        file = entry.path();
        size_t pos = file.find_last_of("/");
        file = file.erase(0,pos+1);

        //Is the data in csv format?
        if (file.size() > simDatacsv.size() && file.substr(file.size()-simDatacsv.size(),file.size()) == simDatacsv)
        {
            // file_count++;
            size_t _pos = file.find_first_of("_");
            size_t _secpos = file.substr(_pos+1,file.size()).find_first_of("_");
            _secpos += _pos+1; //add 1 to account for _pos+1 in substr above
            file_index = stoi(file.substr(0,file.find_first_of("_")));
            if (file[_pos+1] == 'R')
            {
                file_index = 0;
            }
            if (file_index > largest_file_index)
            {
                second_largest_file_index = largest_file_index;
                second_largest_file_name = largest_file_name;
                largest_file_index = file_index;
                largest_file_name = file;
                csv = true;
            }
            else if (file_index > second_largest_file_index)
            {
                second_largest_file_name = file;
                second_largest_file_index = file_index;
                csv = true;
            }



        }
        else if (file.size() > datah5.size() && file.substr(file.size()-datah5.size(),file.size()) == datah5)
        {
            size_t _pos = file.find_first_of("_");
            file_index = stoi(file.substr(0,file.find_first_of("_")));
            if (file_index > largest_file_index)
            {
                largest_file_index = file_index;
                largest_file_name = file;
            }
        }
    }

    if (csv)
    {
        if (getRank() == 0)
        {
            std::string file1 = path + largest_file_name;
            std::string file2 = path + largest_file_name.substr(0,largest_file_name.size()-simDatacsv.size()) + "constants.csv";
            std::string file3 = path + largest_file_name.substr(0,largest_file_name.size()-simDatacsv.size()) + "energy.csv";

            std::string message("Removing the following files: \n"
                                +'\t'+file1+'\n'
                                +'\t'+file2+'\n'
                                +'\t'+file3+'\n');
            MPIsafe_print(std::cerr,message);

            int status1 = remove(file1.c_str());
            int status2 = remove(file2.c_str());
            int status3 = remove(file3.c_str());

            if (status1 != 0)
            {
                MPIsafe_print(std::cerr,"File: '"+file1+"' could not be removed, now exiting with failure.\n");
                MPIsafe_exit(EXIT_FAILURE);
            }
            else if (status2 != 0)
            {
                MPIsafe_print(std::cerr,"File: '"+file2+"' could not be removed, now exiting with failure.\n");
                MPIsafe_exit(EXIT_FAILURE);
            }
            else if (status3 != 0)
            {
                MPIsafe_print(std::cerr,"File: '"+file3+"' could not be removed, now exiting with failure.\n");
                MPIsafe_exit(EXIT_FAILURE);
            }
        }
        MPIsafe_barrier();
        largest_file_name = second_largest_file_name;
    }
    return largest_file_name;
}


void Ball_group::sim_one_step()
{
    int world_rank = getRank();
    int world_size = getSize();
    /// FIRST PASS - Update Kinematic Parameters:
    // t.start_event("UpdateKinPar");
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        // Update velocity half step:
        velh[Ball] = vel[Ball] + .5 * acc[Ball] * attrs.dt;

        // Update angular velocity half step:
        wh[Ball] = w[Ball] + .5 * aacc[Ball] * attrs.dt;

        // Update position:
        pos[Ball] += velh[Ball] * attrs.dt;

        // Reinitialize acceleration to be recalculated:
        acc[Ball] = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        aacc[Ball] = {0, 0, 0};
    }
    // t.end_event("UpdateKinPar");

    double Ha = attrs.Ha;
    double u_r = attrs.u_r;
    double u_s = attrs.u_s;
    double kin = attrs.kin;
    double kout = attrs.kout;
    double h_min = attrs.h_min;
    double dt = attrs.dt;
    int num_parts = attrs.num_particles;
    int threads = attrs.OMPthreads;
    bool write_step = attrs.write_step;

    
    long long A;
    long long B;
    long long pc;
    long long lllen = attrs.num_particles;
    double t0 = omp_get_wtime();
    #pragma omp declare reduction(vec3_sum : vec3 : omp_out += omp_in)
    #pragma omp parallel for num_threads(threads)\
            reduction(vec3_sum:acc[:num_parts],aacc[:num_parts]) reduction(+:PE) \
            shared(world_rank,world_size,Ha,write_step,lllen,R,pos,vel,m,w,\
                u_r,u_s,moi,kin,kout,distances,h_min,dt)\
            default(none) private(A,B,pc) 
    for (pc = world_rank + 1; pc <= (((lllen*lllen)-lllen)/2); pc += world_size)
    {
        long double pd = (long double)pc;
        pd = (sqrt(pd*8.0L+1.0L)+1.0L)*0.5L;
        pd -= 0.00001L;
        A = (long long)pd;
        B = (long long)((long double)pc-(long double)A*((long double)A-1.0L)*.5L-1.0L);

 
        const double sumRaRb = R[A] + R[B];
        const vec3 rVecab = pos[B] - pos[A];  // Vector from a to b.
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
        double oldDist = distances[e];
        /////////////////////////////
        // double inoutT;
        /////////////////////////////
        // Check for collision between Ball and otherBall.
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

            // const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
            //                              ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
            //                              rVecab.normalized();

            // ==========================================
            // Test new vdw force equation with less division
            const double d1 = h2 + twoRah + twoRbh;
            const double d2 = d1 + 4 * Ra * Rb;
            const double numer = 64*Ha*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
            const double denomrecip = 1/(6*d1*d1*d2*d2);
            const vec3 vdwForceOnA = (numer*denomrecip)*rVecab.normalized();
            // ==========================================

            // Elastic force:
            // vec3 elasticForceOnA{0, 0, 0};
            // if (std::fabs(overlap) > 1e-6)
            // {
            //     elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
            // }
            const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
            ///////////////////////////////
            // elasticForce[A] += elasticForceOnA;
            // elasticForce[B] -= elasticForceOnA;
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
            // if (rel_vel_mag > 1e-20)  // Divide by zero protection.
            // if (rel_vel_mag > 1e-8)  // Divide by zero protection.
            ////////////////////////////////////////// CALC THIS AT INITIALIZATION for all combos os Ra,Rb
            // const double u_scale = calc_VDW_force_mag(Ra,Rb,h_min_physical)/
            //                         vdwForceOnA.norm();         //Friction coefficient scale factor
            //////////////////////////////////////////
            if (rel_vel_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
            {
                // slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                // In the frame of A, B applies force in the direction of B's velocity.
                ///////////////////////////////////
                // if (mu_scale)
                // {
                //     if (u_scale[e]*u_s > max_mu)
                //     {
                //         slideForceOnA = max_mu * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                //     }
                //     else
                //     {
                //         slideForceOnA = u_scale[e] * u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                //     }
                // }
                // else
                // {
                    slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                // }
                ///////////////////////////////////
            }
            //////////////////////////////////////
            // slideForce[A] += slideForceOnA;
            // slideForce[B] -= slideForceOnA;
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
                // if (mu_scale)
                // {
                //     if (u_scale[e]*u_r > max_mu)
                //     {
                //         rollForceA = 
                //             -max_mu * elastic_force_A_mag * (w_diff).cross(r_a) / 
                //             (w_diff).cross(r_a).norm();
                //     }
                //     else
                //     {
                //         rollForceA = 
                //             -u_scale[e] * u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
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



            aacc[A] += torqueA / moi[A];
            aacc[B] += torqueB / moi[B];

            if (write_step) {
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
                PE += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
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

            // const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
            //                              ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
            //                              rVecab.normalized();
            // ==========================================
            // Test new vdw force equation with less division
            const double d1 = h2 + twoRah + twoRbh;
            const double d2 = d1 + 4 * Ra * Rb;
            const double numer = 64*Ha*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
            const double denomrecip = 1/(6*d1*d1*d2*d2);
            const vec3 vdwForceOnA = (numer*denomrecip)*rVecab.normalized();
            // ==========================================
           
            /////////////////////////////
            totalForceOnA = vdwForceOnA + gravForceOnA;
            // totalForceOnA = vdwForceOnA;
            // totalForceOnA = gravForceOnA;
            /////////////////////////////
            if (write_step) {
                // PE += -G * m[A] * m[B] * grav_scale / dist; // Gravitational

                const double diffRaRb = R[A] - R[B];
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * R[A] * R[B];
                const double denom_sum = z * z - (sumRaRb * sumRaRb);
                const double denom_diff = z * z - (diffRaRb * diffRaRb);
                const double U_vdw =
                    -Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                PE += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
            }

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

    #ifdef MPI_ENABLE
        MPI_Allreduce(MPI_IN_PLACE,acc,attrs.num_particles*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,aacc,attrs.num_particles*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        double local_PE = PE;
        PE = 0.0;
        MPI_Reduce(&local_PE,&PE,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    #endif

    // t.end_event("CalcForces/loopApplicablepairs");

    // if (write_step) {
    //     ballBuffer << '\n';  // Prepares a new line for incoming data.
    //     // std::cerr<<"Writing "<<num_particles<<" balls"<<std::endl;
    // }

    // THIRD PASS - Calculate velocity for next step:
    // t.start_event("CalcVelocityforNextStep");
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) 
    {
        // Velocity for next step:
        vel[Ball] = velh[Ball] + .5 * acc[Ball] * attrs.dt;
        w[Ball] = wh[Ball] + .5 * aacc[Ball] * attrs.dt;

        /////////////////////////////////
        // if (true) {
        /////////////////////////////////
        if (write_step && world_rank == 0) 
        {
            // Send positions and rotations to buffer:
            int start = data->getWidth("simData")*attrs.num_writes+Ball*data->getSingleWidth("simData");
            ballBuffer[start] = pos[Ball][0];
            ballBuffer[start+1] = pos[Ball][1];
            ballBuffer[start+2] = pos[Ball][2];
            ballBuffer[start+3] = w[Ball][0];
            ballBuffer[start+4] = w[Ball][1];
            ballBuffer[start+5] = w[Ball][2];
            ballBuffer[start+6] = w[Ball].norm();
            ballBuffer[start+7] = vel[Ball][0];
            ballBuffer[start+8] = vel[Ball][1];
            ballBuffer[start+9] = vel[Ball][2];
            ballBuffer[start+10] = 0;

            KE += .5 * m[Ball] * vel[Ball].normsquared() +
                    .5 * moi[Ball] * w[Ball].normsquared();  // Now includes rotational kinetic energy.
            mom += m[Ball] * vel[Ball];
            ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
        }
    }  // THIRD PASS END
    if (write_step && world_rank == 0)
    {
        attrs.num_writes ++;
    }
    // t.end_event("CalcVelocityforNextStep");
}  // one Step end

#ifdef GPU_ENABLE
void Ball_group::sim_one_step_GPU()
{
    
    /// FIRST PASS - Update Kinematic Parameters:
    // t.start_event("UpdateKinPar");
    #pragma acc parallel loop gang worker present(this,velh[0:num_particles],vel[0:num_particles],\
        acc[0:num_particles],dt,wh[0:num_particles],w[0:num_particles],aacc[0:num_particles],\
        pos[0:num_particles],attrs.num_particles)
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        // Update velocity half step:
        velh[Ball] = vel[Ball] + .5 * acc[Ball] * attrs.dt;

        // Update angular velocity half step:
        wh[Ball] = w[Ball] + .5 * aacc[Ball] * attrs.dt;

        // Update position:
        pos[Ball] += velh[Ball] * attrs.dt;

        // Reinitialize acceleration to be recalculated:
        acc[Ball] = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        aacc[Ball] = {0, 0, 0};
    }
    // t.end_event("UpdateKinPar");

    // int threads = attrs.OMPthreads;

    // #pragma acc update device(accsq[0:num_particles*num_particles], aaccsq[0:num_particles*num_particles])

    #pragma acc parallel loop gang worker num_gangs(108) \
        present(attrs.num_particles,accsq[0:num_particles*num_particles],\
            aaccsq[0:num_particles*num_particles])
    for (int i = 0; i < attrs.num_particles*attrs.num_particles; ++i)
    {
        accsq[i] = {0.0,0.0,0.0};
        aaccsq[i] = {0.0,0.0,0.0};
    }


    double pe = 0.0;
    #pragma acc enter data copyin(pe)
    // #pragma acc enter data copyin(writeS/tep)

    double t0 = omp_get_wtime();

    #pragma acc parallel loop gang worker num_gangs(108) num_workers(256) reduction(+:pe) \
        present(pe,accsq[0:num_particles*num_particles],\
            aaccsq[0:num_particles*num_particles],m[0:num_particles],\
            moi[0:num_particles],w[0:num_particles],vel[0:num_particles],\
            pos[0:num_particles],R[0:num_particles],distances[0:num_pairs],\
            attrs.num_pairs,attrs.num_particles,attrs.Ha,attrs.kin,attrs.kout,attrs.h_min,\
            attrs.u_s,attrs.u_r,attrs.write_step,attrs.world_rank,attrs.world_size)
    for (int pc = attrs.world_rank+1; pc <= attrs.num_pairs; pc += attrs.world_size)
    {

        double pd = (double)pc;
        pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
        pd -= 0.00001;
        int A = (int)pd;
        int B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);

 
        const double sumRaRb = R[A] + R[B];
        const vec3 rVecab = pos[B] - pos[A];  // Vector from a to b.
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
        double oldDist = distances[e];
        /////////////////////////////
        // double inoutT;
        /////////////////////////////
        // Check for collision between Ball and otherBall.
        if (overlap > 0) {



            double k;
            if (dist >= oldDist) {
                k = attrs.kout;
            } else {
                k = attrs.kin;
            }

            // Cohesion (in contact) h must always be h_min:
            // constexpr double h = h_min;
            const double h = attrs.h_min;
            const double Ra = R[A];
            const double Rb = R[B];
            const double h2 = h * h;
            // constexpr double h2 = h * h;
            const double twoRah = 2 * Ra * h;
            const double twoRbh = 2 * Rb * h;

            // const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
            //                              ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
            //                              rVecab.normalized();

            // ==========================================
            // Test new vdw force equation with less division
            const double d1 = h2 + twoRah + twoRbh;
            const double d2 = d1 + 4 * Ra * Rb;
            const double numer = 64*attrs.Ha*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
            const double denomrecip = 1/(6*d1*d1*d2*d2);
            const vec3 vdwForceOnA = (numer*denomrecip)*rVecab.normalized();
            // ==========================================

            // Elastic force:
            // vec3 elasticForceOnA{0, 0, 0};
            // if (std::fabs(overlap) > 1e-6)
            // {
            //     elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
            // }
            const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);
            ///////////////////////////////
            // elasticForce[A] += elasticForceOnA;
            // elasticForce[B] -= elasticForceOnA;
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
            // if (rel_vel_mag > 1e-20)  // Divide by zero protection.
            // if (rel_vel_mag > 1e-8)  // Divide by zero protection.
            ////////////////////////////////////////// CALC THIS AT INITIALIZATION for all combos os Ra,Rb
            // const double u_scale = calc_VDW_force_mag(Ra,Rb,h_min_physical)/
            //                         vdwForceOnA.norm();         //Friction coefficient scale factor
            //////////////////////////////////////////
            if (rel_vel_mag > 1e-13)  // NORMAL ONE Divide by zero protection.
            {
                // slideForceOnA = u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                // In the frame of A, B applies force in the direction of B's velocity.
                ///////////////////////////////////
                // if (mu_scale)
                // {
                //     if (u_scale[e]*u_s > max_mu)
                //     {
                //         slideForceOnA = max_mu * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                //     }
                //     else
                //     {
                //         slideForceOnA = u_scale[e] * u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                //     }
                // }
                // else
                // {
                    slideForceOnA = attrs.u_s * elastic_force_A_mag * (frame_A_vel_B / rel_vel_mag);
                // }
                ///////////////////////////////////
            }
            //////////////////////////////////////
            // slideForce[A] += slideForceOnA;
            // slideForce[B] -= slideForceOnA;
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
                // if (mu_scale)
                // {
                //     if (u_scale[e]*u_r > max_mu)
                //     {
                //         rollForceA = 
                //             -max_mu * elastic_force_A_mag * (w_diff).cross(r_a) / 
                //             (w_diff).cross(r_a).norm();
                //     }
                //     else
                //     {
                //         rollForceA = 
                //             -u_scale[e] * u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
                //             (w_diff).cross(r_a).norm();
                //     }
                // }
                // else
                // {
                    rollForceA = 
                        -attrs.u_r * elastic_force_A_mag * (w_diff).cross(r_a) / 
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

            vec3 aaccA = (1/moi[A])*torqueA;
            vec3 aaccB = (1/moi[B])*torqueB;

            aaccsq[A*attrs.num_particles+B].x = aaccA.x;
            aaccsq[A*attrs.num_particles+B].y = aaccA.y;
            aaccsq[A*attrs.num_particles+B].z = aaccA.z;
            aaccsq[B*attrs.num_particles+A].x = aaccB.x;
            aaccsq[B*attrs.num_particles+A].y = aaccB.y;
            aaccsq[B*attrs.num_particles+A].z = aaccB.z;

            // aacc[A] += torqueA / moi[A];
            // aacc[B] += torqueB / moi[B];

            if (attrs.write_step) {
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
                    -attrs.Ha / 6 *
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
            if (h < attrs.h_min)  // If h is closer to 0 (almost touching), use hmin.
            {
                h = attrs.h_min;
            }
            const double Ra = R[A];
            const double Rb = R[B];
            const double h2 = h * h;
            const double twoRah = 2 * Ra * h;
            const double twoRbh = 2 * Rb * h;

            // const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
            //                              ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
            //                                                (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
            //                              rVecab.normalized();
            // ==========================================
            // Test new vdw force equation with less division
            const double d1 = h2 + twoRah + twoRbh;
            const double d2 = d1 + 4 * Ra * Rb;
            const double numer = 64*attrs.Ha*Ra*Ra*Ra*Rb*Rb*Rb*(h+Ra+Rb);
            const double denomrecip = 1/(6*d1*d1*d2*d2);
            const vec3 vdwForceOnA = (numer*denomrecip)*rVecab.normalized();
            // ==========================================
           
            /////////////////////////////
            totalForceOnA = vdwForceOnA + gravForceOnA;
            // totalForceOnA = vdwForceOnA;
            // totalForceOnA = gravForceOnA;
            /////////////////////////////
            if (attrs.write_step) {
                // PE += -G * m[A] * m[B] * grav_scale / dist; // Gravitational

                const double diffRaRb = R[A] - R[B];
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * R[A] * R[B];
                const double denom_sum = z * z - (sumRaRb * sumRaRb);
                const double denom_diff = z * z - (diffRaRb * diffRaRb);
                const double U_vdw =
                    -attrs.Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                pe += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
            }

            // todo this is part of push_apart. Not great like this.
            // For pushing apart overlappers:
            // vel[A] = { 0,0,0 };
            // vel[B] = { 0,0,0 };
        }

        // Newton's equal and opposite forces applied to acceleration of each ball:
        vec3 accA = (1/m[A])*totalForceOnA; 
        vec3 accB = -1.0*(1/m[B])*totalForceOnA; 

        accsq[A*attrs.num_particles+B].x = accA.x;
        accsq[A*attrs.num_particles+B].y = accA.y;
        accsq[A*attrs.num_particles+B].z = accA.z;
        accsq[B*attrs.num_particles+A].x = accB.x;
        accsq[B*attrs.num_particles+A].y = accB.y;
        accsq[B*attrs.num_particles+A].z = accB.z;


        // So last distance can be known for COR:
        distances[e] = dist;

    }

    // #pragma acc loop seq
    #pragma acc parallel loop gang num_gangs(108) \
        present(attrs.num_particles,acc[0:num_particles],aacc[0:num_particles],\
            accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles])
    for (int i = 0; i < attrs.num_particles; i++)
    {
        #pragma acc loop seq
        for (int j = 0; j < attrs.num_particles; j++)
        {
            acc[i].x += accsq[i*attrs.num_particles+j].x;
            acc[i].y += accsq[i*attrs.num_particles+j].y;
            acc[i].z += accsq[i*attrs.num_particles+j].z;
            aacc[i].x += aaccsq[i*attrs.num_particles+j].x;
            aacc[i].y += aaccsq[i*attrs.num_particles+j].y;
            aacc[i].z += aaccsq[i*attrs.num_particles+j].z;
        }
    // #pragma acc update self(acc[0:num_particles],aacc[0:num_particles]) //if(write_step)
    // #pragma acc update self(acc[i],aacc[i]) //if(write_step)
    }

    #pragma acc update host(acc[0:num_particles],aacc[0:num_particles])
    // std::cout<<aaccsq[0].x<<','<<aaccsq[0].y<<','<<aaccsq[0].z<<std::endl;
    
    if (attrs.write_step)
    {
        #pragma acc update host(pe)
        PE = pe;
    }

    #ifdef MPI_ENABLE
        MPI_Allreduce(MPI_IN_PLACE,acc,attrs.num_particles*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,aacc,attrs.num_particles*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        double local_PE = PE;
        PE = 0.0;
        MPI_Reduce(&local_PE,&PE,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
        #pragma acc update device(acc[0:num_particles],aacc[0:num_particles])
    #endif


    // t.end_event("CalcForces/loopApplicablepairs");
    #pragma acc parallel loop gang worker num_gangs(108) num_workers(256) \
        present(acc[0:num_particles],aacc[0:num_particles],w[0:num_particles],\
            vel[0:num_particles],velh[0:num_particles],wh[0:num_particles],attrs.num_particles,attrs.dt)
    for (int Ball = 0; Ball < attrs.num_particles; Ball++) {
        // Velocity for next step:
        vel[Ball] = velh[Ball] + .5 * acc[Ball] * attrs.dt;
        w[Ball] = wh[Ball] + .5 * aacc[Ball] * attrs.dt;
    }  // THIRD PASS END

    // THIRD PASS - Calculate velocity for next step:
    // t.start_event("CalcVelocityforNextStep");
    if (attrs.write_step && attrs.world_rank == 0) 
    {
        #pragma acc update host(w[0:num_particles],vel[0:num_particles],pos[0:num_particles])// if(attrs.write_step && attrs.world_rank == 0)
       
        for (int Ball = 0; Ball < attrs.num_particles; Ball++) 
        {
            // Send positions and rotations to buffer:
            int start = data->getWidth("simData")*attrs.num_writes+Ball*data->getSingleWidth("simData");
            ballBuffer[start] = pos[Ball][0];
            ballBuffer[start+1] = pos[Ball][1];
            ballBuffer[start+2] = pos[Ball][2];
            ballBuffer[start+3] = w[Ball][0];
            ballBuffer[start+4] = w[Ball][1];
            ballBuffer[start+5] = w[Ball][2];
            ballBuffer[start+6] = w[Ball].norm();
            ballBuffer[start+7] = vel[Ball][0];
            ballBuffer[start+8] = vel[Ball][1];
            ballBuffer[start+9] = vel[Ball][2];
            ballBuffer[start+10] = 0;

            KE += .5 * m[Ball] * vel[Ball].normsquared() +
                    .5 * moi[Ball] * w[Ball].normsquared();  // Now includes rotational kinetic energy.
            mom += m[Ball] * vel[Ball];
            ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
        }
    }  // THIRD PASS END
    if (attrs.write_step && attrs.world_rank == 0)
    {
        attrs.num_writes ++;
    }
    // t.end_event("CalcVelocityforNextStep");
}  // one Step end
#endif

void
Ball_group::sim_looper(unsigned long long start_step=1)
{

    attrs.world_rank = getRank();
    attrs.world_size = getSize();
    attrs.num_pairs = static_cast<int>(attrs.num_particles*(attrs.num_particles-1)/2);

    attrs.num_writes = 0;
    unsigned long long Step;
    // attrs.writeStep = false;

    if (attrs.world_rank == 0)
    {   
        attrs.startProgress = time(nullptr);
    }

    std::string message(
        "Beginning simulation...\nstart step:" +
        std::to_string(start_step)+'\n' +
        "Stepping through "+std::to_string(attrs.steps)+" steps.\n" + 
        "Simulating "+dToSci(attrs.simTimeSeconds)+" seconds per sim.\n" + 
        "Writing out every "+dToSci(attrs.timeResolution)+" seconds.\n" +
        "For a total of "+dToSci(attrs.simTimeSeconds/attrs.timeResolution)+" timesteps saved per sim.\n");
    MPIsafe_print(std::cerr,message);

    //Set the number of threads to be appropriate
    #ifndef GPU_ENABLE
        attrs.OMPthreads = get_num_threads();
    #else
        #pragma acc enter data create(accsq[0:attrs.num_particles*attrs.num_particles],aaccsq[0:attrs.num_particles*attrs.num_particles])

        // #pragma acc enter data copyin(this) 
        #pragma acc enter data copyin(moi[0:num_particles],m[0:num_particles],\
            w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],\
            distances[0:num_pairs]) 
        #pragma acc enter data copyin(accsq[0:num_particles*num_particles],\
            aaccsq[0:num_particles*num_particles],acc[0:num_particles],aacc[0:num_particles],\
            velh[0:num_particles],wh[0:num_particles]) 
        #pragma acc enter data copyin(attrs.dt,attrs.num_pairs,attrs.num_particles,attrs.Ha,attrs.kin,attrs.kout,attrs.h_min,\
            attrs.u_s,attrs.u_r,attrs.world_rank,attrs.world_size,attrs.write_step)
    #endif

    for (Step = start_step; Step < attrs.steps; Step++)  // Steps start at 1 for non-restart because the 0 step is initial conditions.
    {
        // simTimeElapsed += dt; //New code #1
        // Check if this is a write step:
        if (Step % attrs.skip == 0) {
            // if (world_rank == 0)
            // {
            //     t.start_event("writeProgressReport");
            // }
            attrs.write_step = true;
            #ifdef GPU_ENABLE
                #pragma acc update device(attrs.write_step)
            #endif
            // std::cerr<<"Write step "<<Step<<std::endl;

            /////////////////////// Original code #1
            attrs.simTimeElapsed += attrs.dt * attrs.skip;
            ///////////////////////

            if (attrs.world_rank == 0)
            {
                // Progress reporting:
                float eta = ((time(nullptr) - attrs.startProgress) / static_cast<float>(attrs.skip) *
                             static_cast<float>(attrs.steps - Step)) /
                            3600.f;  // Hours.
                float real = (time(nullptr) - attrs.start) / 3600.f;
                float simmed = static_cast<float>(attrs.simTimeElapsed / 3600.f);
                float progress = (static_cast<float>(Step) / static_cast<float>(attrs.steps) * 100.f);
                fprintf(
                    stderr,
                    "%llu\t%2.0f%%\tETA: %5.2lf\tReal: %5.2f\tSim: %5.2f hrs\tR/S: %5.2f\n",
                    Step,
                    progress,
                    eta,
                    real,
                    simmed,
                    real / simmed);
                // fprintf(stdout, "%u\t%2.0f%%\tETA: %5.2lf\tReal: %5.2f\tSim: %5.2f hrs\tR/S: %5.2f\n", Step,
                // progress, eta, real, simmed, real / simmed);
                fflush(stdout);
                attrs.startProgress = time(nullptr);
                // t.end_event("writeProgressReport");
            }
        } else {
            attrs.write_step = attrs.debug;
        }

        // Physics integration step:
        #ifndef GPU_ENABLE
            sim_one_step();
        #else
            sim_one_step_GPU();
        #endif

        if (attrs.write_step) {
            // t.start_event("writeStep");
            // Write energy to stream:
            ////////////////////////////////////
            //TURN THIS ON FOR REAL RUNS!!!
            // energyBuffer = std::vector<double> (data->getWidth("energy"));
            // std::cerr<<"start,num_writes: "<<start<<','<<num_writes<<std::endl;
            if (attrs.world_rank == 0)
            {    
                int start = data->getWidth("energy")*(attrs.num_writes-1);
                energyBuffer[start] = attrs.simTimeElapsed;
                energyBuffer[start+1] = PE;
                energyBuffer[start+2] = KE;
                energyBuffer[start+3] = PE+KE;
                energyBuffer[start+4] = mom.norm();
                energyBuffer[start+5] = ang_mom.norm();

                if (Step / attrs.skip % 10 == 0) 
                {

                    std::cerr << "vMax = " << getVelMax() << " Steps recorded: " << Step / attrs.skip << '\n';
                    std::cerr << "Data Write to "<<data->getFileName()<<"\n";
                    
                    data->Write(ballBuffer,"simData",bufferlines);

                    ballBuffer.clear();
                    ballBuffer = std::vector<double>(data->getWidth("simData")*bufferlines);
                    data->Write(energyBuffer,"energy");
                    energyBuffer.clear();
                    energyBuffer = std::vector<double>(data->getWidth("energy")*bufferlines);

                    attrs.num_writes = 0;

                }  // Data export end
                
                attrs.lastWrite = time(nullptr);
            }
            
            // Reinitialize energies for next step:
            KE = 0;
            PE = 0;
            mom = {0, 0, 0};
            ang_mom = {0, 0, 0};

            if (attrs.dynamicTime) { calibrate_dt(Step, false); }
            // t.end_event("writeStep");
        }  // writestep end
    }

    #ifdef GPU_ENABLE
        #pragma acc exit data delete(accsq[0:num_particles*num_particles],\
            aaccsq[0:num_particles*num_particles],acc[0:num_particles],aacc[0:num_particles])
        #pragma acc exit data delete(m[0:num_particles],w[0:num_particles],vel[0:num_particles],\
            pos[0:num_particles],R[0:num_particles],distances[0:num_pairs])
        #pragma acc exit data delete(attrs.dt,attrs.num_pairs,attrs.num_particles,attrs.Ha,\
            attrs.kin,attrs.kout,attrs.h_min,attrs.u_s,attrs.u_r,attrs.world_rank,attrs.world_size,attrs.write_step)
        // #pragma acc exit data delete(this)
    #endif

    if (attrs.world_rank == 0)
    {
        const time_t end = time(nullptr);

        std::cerr << "Simulation complete! \n"
                  << attrs.num_particles << " Particles and " << Step << '/' << attrs.steps << " Steps.\n"
                  << "Simulated time: " << attrs.steps * attrs.dt << " seconds\n"
                  << "Computation time: " << end - attrs.start << " seconds\n";
        std::cerr << "\n===============================================================\n";
    }


}  // end simLooper
