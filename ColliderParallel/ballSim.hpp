#include "../external/json/single_include/nlohmann/json.hpp"
#include "../Utils.hpp"
#include "../vec3.hpp"
#include "../timing/timing.hpp"
// #include <tbb>
#include <mutex>
#include <execution>
#include <random>


std::mutex g_mutex;

using json = nlohmann::json;
using std::numbers::pi;

struct P;
struct P_pair;
struct P_group;
constexpr double Kb = 1.380649e-16; //in erg/K

// template <typename T>
// ostream& operator<<(ostream& os, const vector<T>& v)
// {
//     os << "[";
//     for (int i = 0; i < v.size(); ++i) {
//         os << v[i];
//         if (i != v.size() - 1)
//             os << ", ";
//     }
//     os << "]\n";
//     return os;
// }

//Every ball is represented with struct P
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

    std::vector<int> id = std::vector<int> (3);
    int unique_id;

    double
        R = 0,
        m = 0,
        moi = 0;

    bool operator==(const P& p) const
    {
        return (unique_id == p.unique_id);
    }
    bool operator!=(const P& p) const
    {
        return !(*this == p);
    }
    // P() = default;
    // P(vec3 p, vec3 v, vec3 ww, 
    //     double rad, double mass, double moia)
    //     :pos(p),vel{v},w(ww),R(rad),m(mass),w(ww){}
};


//Every ball will feel forces from every other ball (with gravity)
//So this struct is used to find all pairs in order to loop over them
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


//This struct defines the whole group of balls in the simulation
//including important values for the simulation
struct P_group
{
    enum distributions {constant, logNorm};
    double lnSigma = 0.2; //sigma for log normal distribution 
    bool parallel = true;

    std::vector<double> distances; ///temperary for testing why this code isnt the same as the original
    
    int n; // Number of particles in the group
    int mapElemets = -1;
    distributions radiiDistribution;
    std::vector<P> p_group; // Group of particles
    std::vector<P_pair> pairs;
    std::vector<int> applicable_pairs;
    std::vector<int> used_ids;
    // std::vector<std::vector<int>> gridIDs;
    std::unordered_map<std::string,std::vector<P>> IDToGrid;
    // double U; // Potential energy
    // double T; // Kinetic Energy
    vec3 mom = {0.0,0.0,0.0};
    vec3 ang_mom = {0.0,0.0,0.0};
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
    std::string output_prefix = "dataFile";
    std::string filenum;

    double u_s = -1.0, u_r = -1.0;
    double cor = -1.0;
    double simTimeSeconds = -1.0;
    double timeResolution = -1.0;
    double fourThirdsPiRho = -1.0;
    double scaleBalls = -1.0; //ball radius
    double maxOverlap = -1.0;
    double v_custom = -10.0;
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
    int ball_counter = 0;

    bool writeStep = false;
    bool dynamicTime;

    time_t start = time(nullptr);        // For end of program analysis
    time_t startProgress; // For progress reporting (gets reset)
    time_t lastWrite;     // For write control (gets reset)

    timey t;

    void parse_input_file(char const* location);
    void makeSpheres(const int nBalls, distributions radiiDistribution);
    // void oneSizeSphere(const int nBalls);
    // void distSizeSphere(const int nBalls);
    void placeBalls(const int nBalls);
    [[nodiscard]] vec3 getCOM() const;
    [[nodiscard]] double get_radius(const vec3& center) const;
    void calibrate_dt(int const Step, const double& customSpeed = -1.);
    void updateDTK(const double& velocity);
    void calc_v_collapse();
    [[nodiscard]] double getVelMax();
    [[nodiscard]] double getRmax();
    [[nodiscard]] double getRmin();
    void calc_helpfuls();
    void simInit_cond_and_center();
    vec3 calc_momentum(const std::string& of = "") const;
    void to_origin();
    void init_conditions();
    void make_pairs();
    void compute_acceleration(P_pair& p_pair);
    void compute_distance(P_pair& p_pair);
    void compute_energy(P& P);
    void compute_velocity(P& P);
    void update_kinematics(P& P);
    void generate_ball_field(const int nBalls);
    double getMass();
    void sim_looper();
    void sim_one_step(const bool write_step);
    void sim_continue(const std::string& path, int start_file_index=0);
    void write_to_buffer();
    void zeroAngVel();
    void zeroVel();
    void loadSim(const std::string& path, const std::string& filename);
    // [[nodiscard]] static std::string getLastLine(const std::string& path, 
    [[nodiscard]] std::string getLastLine(const std::string& path, 
                    const std::string& filename);
    void parseSimData(const std::string& line,
                    const std::string& path,const std::string& filename);
    void loadConsts(const std::string& path, const std::string& filename);
    void sim_init_write(std::string filename, int counter = 0);
    P dust_agglomeration_particle_init();
    inline vec3 calc_momentum(const P& partile) const;
    vec3 calc_momentum(const P_group& group) const;
    inline double getMass(const P& p);
    double getMass(const P_group& pg);
    inline void kick(const vec3& vec, P& p);
    void kick(const vec3& vec, P_group& pg);
    void kick(const vec3& vec);
    void add_projectile();
    vec3 dust_agglomeration_offset(
            const double3x3 local_coords, vec3 projectile_pos,
            vec3 projectile_vel, const double projectile_rad);
    inline double calc_moi(const double& radius, const double& mass);
    void findGroup(P &p); //TODO impliment grid 
    void findGroups();
    void mapGroups();
    inline std::string getKey(std::vector<int> v);
    inline std::string getKey(int x, int y, int z);
    std::vector<P> getBalls(P ball);
    void findPairs(P &p);
    void getRelaventPairs();
    void printPairs();
    int get_unique_id();
    void write_timing();

    P_group() = default;
    P_group(int n) : n(n) {}
    P_group(std::vector<P> p_group) : p_group(p_group) {}
    P_group(const bool generate, const char* path);
    // P_group(const bool generate, const char* path, timey* time);
    P_group(const std::string& path, int start_file_index=0);
    // P_group(const std::string& path,timey* time, int start_file_index=0);

    void sim_one_step_single_core(const bool write_step);///temperary for testing why this code isnt the same as the original
    void compute_acceleration_speed_testing(P_pair& p_pair);
    void init_conditions_single_core();


    std::stringstream ballBuffer;
    std::stringstream energyBuffer;
private:

    
};


/// @brief For starting a sim.
/// @param generate is a bool (Is this necessary???)
/// @param path is the path to where the sim will save.
/// @param time is the instance of the timing class used for this sim
P_group::P_group(const bool generate, const char* path)
{
    // t = time;
    parse_input_file(path);//should be first in constructor
    energyBuffer.precision(12);  // Need more precision on momentum.
    
    generate_ball_field(genBalls);

    p_group[0].pos = {0, 1.101e-5, 0};
    p_group[0].vel = {0, 0, 0};
    if (genBalls > 1)
    {
        p_group[1].pos = {0, -1.101e-5, 0};
        p_group[1].vel = {0, 0, 0};
    }
    mTotal = getMass();
    calc_v_collapse();
    std::cerr<<"INIT VCUSTOM "<<v_custom<<std::endl;
    calibrate_dt(0, v_custom);
    simInit_cond_and_center();



    // std::cerr<<p_group[0].pos<<std::endl;
    // std::cerr<<p_group[0].vel<<std::endl;

    // std::cerr<<p_group[1].pos<<std::endl;
    // std::cerr<<p_group[1].vel<<std::endl;
}


/// @brief For continuing a sim.
/// @param path is the filename and path excluding the suffix _simData.csv, _constants.csv, etc.
/// @param time is the instance of the timing class used for this sim
/// @param start_file_index tells the sim where to start if you are restarting
P_group::P_group(const std::string& path, int start_file_index)
{
    // t = time;
    parse_input_file(path.c_str());
    energyBuffer.precision(12);  // Need more precision on momentum.
    num_particles = genBalls + start_file_index;
    sim_continue(path,start_file_index-1);
    calc_v_collapse();
    calibrate_dt(0, v_custom);
    simInit_cond_and_center();
}


/// @brief Make ballGroup from file data.
/// @param path is the full path to the folder you want the sim in
/// @param filename excluding the suffix _simData.csv, _constants.csv, etc. 
void P_group::loadSim(const std::string& path, const std::string& filename)
{

    parseSimData(getLastLine(path, filename),path,filename);

    calc_helpfuls();

    std::cerr << "Balls: " << num_particles << '\n';
    std::cerr << "Mass: " << mTotal << '\n';
    std::cerr << "Approximate radius: " << initialRadius << " cm.\n";
}


/// @brief returns the last line of a given data file
/// @param path is the full path to the folder you want the sim in
/// @param filename excluding the suffix _simData.csv, _constants.csv, etc. 
[[nodiscard]] std::string P_group::getLastLine(const std::string& path, const std::string& filename)
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
        std::cerr << "Could not open simData file: " << simDataFilepath << "... Exiting program."
                  << '\n';
        exit(EXIT_FAILURE);
    }
    return "-1";
}


/// @brief get and set previous sim constants by filename.
/// @param path is the full path to the folder you want the sim in
/// @param filename excluding the suffix _simData.csv, _constants.csv, etc. 
void P_group::loadConsts(const std::string& path, const std::string& filename)
{
    // Get radius, mass, moi:
    std::string constantsFilename = path + filename + "constants.csv";
    if (auto ConstStream = std::ifstream(constantsFilename, std::ifstream::in)) {
        std::string line, lineElement;
        for (int A = 0; A < num_particles; A++) {
            std::getline(ConstStream, line);  // Ball line.
            std::stringstream chosenLine(line);
            std::getline(chosenLine, lineElement, ',');  // Radius.
            p_group[A].R = std::stod(lineElement);
            std::getline(chosenLine, lineElement, ',');  // Mass.
            p_group[A].m = std::stod(lineElement);
            std::getline(chosenLine, lineElement, ',');  // Moment of inertia.
            p_group[A].moi = std::stod(lineElement);
        }
    } else {
        std::cerr << "Could not open constants file: " << constantsFilename << "... Existing program."
                  << '\n';
        exit(EXIT_FAILURE);
    }
}

int P_group::get_unique_id()
{
    ball_counter++;
    return ball_counter-1;
}


/// @brief makes new P_group from a line of a previously existing sim file
/// @param line is the line of data you want parsed
/// @param path is the full path to the folder you want the sim in
/// @param filename excluding the suffix _simData.csv, _constants.csv, etc. 
void P_group::parseSimData(const std::string& line,const std::string& path,const std::string& filename)
{
    std::string lineElement;

    std::stringstream chosenLine(line);  // This is the last line of the read file, containing all data
                                         // for all balls at last time step
    // Get position and angular velocity data:
    for (int A = 0; A < num_particles; A++) {
        P new_p;
        new_p.unique_id = get_unique_id();
        for (int i = 0; i < 3; i++)  // Position
        {
            std::getline(chosenLine, lineElement, ',');
            new_p.pos[i] = std::stod(lineElement);
        }
        for (int i = 0; i < 3; i++)  // Angular Velocity
        {
            std::getline(chosenLine, lineElement, ',');
            new_p.w[i] = std::stod(lineElement);
        }
        std::getline(chosenLine, lineElement, ',');  // Angular velocity magnitude skipped
        for (int i = 0; i < 3; i++)                  // velocity
        {
            std::getline(chosenLine, lineElement, ',');
            new_p.vel[i] = std::stod(lineElement);
        }
        for (int i = 0; i < properties - 10; i++)  // We used 10 elements. This skips the rest.
        {
            std::getline(chosenLine, lineElement, ',');
        }
        p_group.push_back(new_p);
    }
    loadConsts(path, filename);
}


/// @brief returns the mass of a given P
inline double P_group::getMass(const P& p)
{
    return p.m;        
}


/// @brief sets the whole P_group's velocity to zero
void P_group::zeroVel()
{
    auto lambda = [=](P &b){ b.vel = {0,0,0}; };
    
    if (parallel)
    {
        std::for_each(std::execution::par_unseq, p_group.begin(), p_group.end(),
                lambda);
    }
    else
    {
        std::for_each(std::execution::unseq, p_group.begin(), p_group.end(),
                lambda);
    }
}


/// @brief Kick ballGroup (give the whole thing a velocity) given a velocity vector
void P_group::kick(const vec3& vec)
{
    for (int Ball = 0; Ball < num_particles; Ball++) 
    { 
        p_group[Ball].vel = p_group[Ball].vel + vec; 
    }
}


/// @brief Kick a given ballGroup (give the whole thing a velocity) given a velocity vector
void P_group::kick(const vec3& vec, P_group& pg) 
{
    for (int Ball = 0; Ball < num_particles; Ball++) 
    { 
        pg.p_group[Ball].vel = pg.p_group[Ball].vel + vec; 
    }
}


/// @brief Kick a given ball (P) given a velocity vector
inline void P_group::kick(const vec3& vec, P& p) 
{
    p.vel = p.vel + vec;
}


/// @brief zeros the angular velocity for the p_group
void P_group::zeroAngVel()
{
    auto lambda = [&](P &b){ b.w = {0,0,0}; };
    if (parallel)
    {
        std::for_each(std::execution::par_unseq, p_group.begin(), p_group.end(),
                lambda);
    }
    else
    {
        std::for_each(std::execution::unseq, p_group.begin(), p_group.end(),
                lambda);   
    }

    // for (int Ball = 0; Ball < num_particles; Ball++) { p_group[Ball].w = {0, 0, 0}; }
}


/// @brief moves the center of mass of the P_group to the origin
void P_group::to_origin()
{
    const vec3 com = getCOM();

    auto lambda = [=](P &p){ p.pos = p.pos - com; };
    if (parallel)
    {
        std::for_each(std::execution::par_unseq, p_group.begin(), p_group.end(),
                lambda);
    }
    else
    {
        std::for_each(std::execution::unseq, p_group.begin(), p_group.end(),
                lambda);
    }
    // for (int Ball = 0; Ball < num_particles; Ball++) { p_group[Ball].pos -= com; }
}


/// @brief returns the total mass of the P_group
double P_group::getMass()
{
    auto lambda = [&](double sum, const P &b){return sum + b.m; };
    return std::accumulate(p_group.begin(), p_group.end(), 0.0, lambda);        
}


/// @brief returns the total mass of a given P_group
double P_group::getMass(const P_group& pg)
{
    auto lambda = [&](double sum, const P &b){return sum + b.m; };
    return std::accumulate(pg.p_group.begin(), pg.p_group.end(), 0.0, lambda);        
}


/// @brief returns the total momentum of the P_group
/// @param of does nothing for the function but is useful to keep track of what you want the momentum of
//TODO: make this parallel
vec3 P_group::calc_momentum(const std::string& of) const
{
    // vec3 init = p_group[0].m * p_group[0].vel;
    // auto lambda = [&](vec3& sum, const P &b){return sum + b.m * b.vel; };
    // return std::accumulate(p_group.begin(), p_group.end(), init, lambda);  

    vec3 pTotal = {0, 0, 0};
    for (int Ball = 0; Ball < num_particles; Ball++) 
    { 
        pTotal += p_group[Ball].m * p_group[Ball].vel; 
    }
    // fprintf(stderr, "%s Momentum Check: %.2e, %.2e, %.2e\n", of.c_str(), pTotal.x, pTotal.y, pTotal.z);
    return pTotal;
}


/// @brief returns the total momentum of a given P_group
vec3 P_group::calc_momentum(const P_group& group) const
{
    vec3 pTotal = {0, 0, 0};
    for (int Ball = 0; Ball < num_particles; Ball++)
    { 
        pTotal += group.p_group[Ball].m * group.p_group[Ball].vel;
    }
    // fprintf(stderr, "%s Momentum Check: %.2e, %.2e, %.2e\n", of.c_str(), pTotal.x, pTotal.y, pTotal.z);
    return pTotal;
}


/// @brief returns the momentum of a single given particle (ball)
inline vec3 P_group::calc_momentum(const P& partile) const
{
    return partile.m*partile.vel;
}

void P_group::init_conditions_single_core()
{
    make_pairs();
    // SECOND PASS - Check for collisions, apply forces and torques:
    for (int A = 1; A < num_particles; A++)  // cuda
    {
        // DONT DO ANYTHING HERE. A STARTS AT 1.
        for (int B = 0; B < A; B++) {
            const double sumRaRb = p_group[A].R + p_group[B].R;
            const vec3 rVecab = p_group[B].pos - p_group[A].pos;  // Vector from a to b.
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
                const double Ra = p_group[A].R;
                const double Rb = p_group[B].R;
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
                const vec3 gravForceOnA = (G * p_group[A].m * p_group[B].m / (dist * dist)) * (rVecab / dist);

                // Sliding and Rolling Friction:
                vec3 slideForceOnA{0, 0, 0};
                vec3 rollForceA{0, 0, 0};
                vec3 torqueA{0, 0, 0};
                vec3 torqueB{0, 0, 0};

                // Shared terms:
                const double elastic_force_A_mag = elasticForceOnA.norm();
                const vec3 r_a = rVecab * p_group[A].R / sumRaRb;  // Center to contact point
                const vec3 r_b = rVecba * p_group[B].R / sumRaRb;
                const vec3 w_diff = p_group[A].w - p_group[B].w;

                // Sliding friction terms:
                const vec3 d_vel = p_group[B].vel - p_group[A].vel;
                const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
                                           p_group[A].w.cross(r_a) - p_group[B].w.cross(r_a);

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

                p_group[A].aacc += torqueA / p_group[A].moi;
                p_group[B].aacc += torqueB / p_group[B].moi;


                // No factor of 1/2. Includes both spheres:
                // PE += -G * m[A] * m[B] / dist + 0.5 * k * overlap * overlap;

                // Van Der Waals + elastic:
                const double diffRaRb = p_group[A].R - p_group[B].R;
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * p_group[A].R * p_group[B].R;
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
                const double Ra = p_group[A].R;
                const double Rb = p_group[B].R;
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

                const double diffRaRb = p_group[A].R - p_group[B].R;
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * p_group[A].R * p_group[B].R;
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
            p_group[A].acc += totalForceOnA / p_group[A].m;
            p_group[B].acc -= totalForceOnA / p_group[B].m;

            // So last distance can be known for COR:
            distances[e] = dist;
        }
        // DONT DO ANYTHING HERE. A STARTS AT 1.
    }

    // Calc energy:
    for (int Ball = 0; Ball < num_particles; Ball++) {
        KE += .5 * p_group[Ball].m * p_group[Ball].vel.dot(p_group[Ball].vel) + .5 * p_group[Ball].moi * p_group[Ball].w.dot(p_group[Ball].w);
        mom += p_group[Ball].m * p_group[Ball].vel;
        ang_mom += p_group[Ball].m * p_group[Ball].pos.cross(p_group[Ball].vel) + p_group[Ball].moi * p_group[Ball].w;
    }
}


/// @brief sets initial conditions for a simulation
void P_group::init_conditions()
{

    t.start_event("init_conditions");
    make_pairs();
    //calc init accelerations
    writeStep = true;
    if (parallel)
    {
        std::for_each(std::execution::par, pairs.begin(), pairs.end(),
                    std::bind_front(&P_group::compute_acceleration, this));
    }
    else
    {
        std::for_each(std::execution::seq, pairs.begin(), pairs.end(),
                    std::bind_front(&P_group::compute_acceleration, this));        
                    // std::bind_front(&P_group::compute_acceleration_speed_testing, this));        
    }
    writeStep = false;
    //calc init energy:
    for (int i = 0; i < num_particles; ++i)
    {
        compute_energy(p_group[i]);
    }
    // std::for_each(std::execution::par_unseq, p_group.begin(), p_group.end(),
                // std::bind_front(&P_group::compute_energy, this));
    t.end_event("init_conditions");
}


//@brief returns new position of particle after it is given random offset
//@param local_coords is plane perpendicular to direction of projectile
//@param projectile_pos is projectile's position before offset is applied
//@param projectile_vel is projectile's velocity
//@param projectile_rad is projectile's radius
vec3 P_group::dust_agglomeration_offset(
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
            if (line_sphere_intersect(test_pos, projectile_vel, p_group[i].pos, p_group[i].R + projectile_rad)) {
                new_position = test_pos;
                intersect = true;
                break;
            }
        }
    } while (!intersect);
    return new_position;
}


/// @brief returns the moment of inertia of a solid sphere given radius and mass 
inline double P_group::calc_moi(const double& radius, const double& mass) 
    { return .4 * mass * radius * radius; }


/// @brief returns the projectile for a BPCA style sim 
P P_group::dust_agglomeration_particle_init()
{
    P projectile;
    projectile.unique_id = get_unique_id();
    const auto cluster_radius = get_radius(vec3(0, 0, 0));
    const vec3 projectile_direction = rand_unit_vec3();
    projectile.pos = projectile_direction * (cluster_radius + scaleBalls * 4);
    if (radiiDistribution == constant)
    {
        projectile.R = scaleBalls;  // rand_between(1,3)*1e-5;
    }
    else
    {
        projectile.R = lognorm_dist(scaleBalls*std::exp(-5*std::pow(lnSigma,2)/2),lnSigma);
    }
    projectile.w = {0, 0, 0};
    projectile.m = density * 4. / 3. * pi * std::pow(projectile.R, 3);
    // Velocity toward origin:
    if (temp > 0)
    {
        double a = std::sqrt(Kb*temp/projectile.m);
        v_custom = max_bolt_dist(a); 
        std::cerr<<"v_custom set to "<<v_custom<< "cm/s based on a temp of "
                <<temp<<" degrees K."<<std::endl; 
    }
    projectile.vel = -v_custom * projectile_direction;
    // projectile.R[0] = 1e-5;  // rand_between(1,3)*1e-5;
    projectile.moi = calc_moi(projectile.R, projectile.m);

    const double3x3 local_coords = local_coordinates(to_double3(projectile_direction));

    projectile.pos = dust_agglomeration_offset(local_coords,projectile.pos,projectile.vel,projectile.R);
    return projectile;
}


/// @brief Uses P_group as target and adds one particle to hit it:
void P_group::add_projectile()
{
    // Load file data:
    std::cerr << "Add Particle\n";

    auto projectile = dust_agglomeration_particle_init();

    
    // Collision velocity calculation:
    const vec3 p_target{calc_momentum("p_target")};
    const vec3 p_projectile{calc_momentum(projectile)};
    const vec3 p_total{p_target + p_projectile};
    const double m_target{getMass()};
    const double m_projectile{getMass(projectile)};
    const double m_total{m_target + m_projectile};
    const vec3 v_com = p_total / m_total;

    // Negate total system momentum:
    kick(-v_com,projectile);
    kick(-v_com);

    fprintf(
        stderr,
        "\nTarget Velocity: %.2e\nProjectile Velocity: %.2e\n",
        p_group[0].vel.norm(),
        projectile.vel.norm());

    std::cerr << '\n';

    p_group.push_back(projectile);
    num_particles = num_particles + 1;

    calc_helpfuls();

    // Hack - if v_custom is less than 1 there are problems if dt is calibrated to this
    //        if v_custom is greater than 1 you need to calibrate dt to that v_custom
    if (v_custom < 1)
    {
        calibrate_dt(0, 1);
    }
    else
    {
        calibrate_dt(0, v_custom);
    }

    distances = std::vector<double> ((int(num_particles*num_particles-num_particles)/2));

    // init_conditions();
    init_conditions_single_core();
    
    to_origin();

    if (true)
    {
        for (int i = 0; i < num_particles; i++)
        {
            std::cerr<<"=================="<<std::endl;
            std::cerr<<p_group[i].pos<<std::endl;
            std::cerr<<p_group[i].vel<<std::endl;
            std::cerr<<"=================="<<std::endl;
        }
    }

}


/// @brief generates a random field of balls given the number of balls
void P_group::generate_ball_field(const int nBalls)
{
    std::cerr << "CLUSTER FORMATION\n";
    makeSpheres(nBalls, radiiDistribution);
    distances = std::vector<double> ((int(nBalls*nBalls-nBalls)/2));

    // if (radiiDistribution == constant)
    // {
    //     oneSizeSphere(nBalls);
    // }
    // else
    // {
    //     distSizeSphere(nBalls);
    // }
    calc_helpfuls();
}


/// Not sure if this will be useful in the future. Leaving for now
// void P_group::write_to_buffer()
// {
//     for (size_t i = 0; i < num_particles; i++)
//     {
//         P& cp = p_group[i]; // Current particle

//         // Send positions and rotations to buffer:
//         if (i == 0)
//         {
//             ballBuffer
//                 << cp.pos[0] << ','
//                 << cp.pos[1] << ','
//                 << cp.pos[2] << ','
//                 << cp.w[0] << ','
//                 << cp.w[1] << ','
//                 << cp.w[2] << ','
//                 << cp.w.norm() << ','
//                 << cp.vel.x << ','
//                 << cp.vel.y << ','
//                 << cp.vel.z << ','
//                 << 0;
//         }
//         else
//         {
//             ballBuffer
//                 << ',' << cp.pos[0] << ','
//                 << cp.pos[1] << ','
//                 << cp.pos[2] << ','
//                 << cp.w[0] << ','
//                 << cp.w[1] << ','
//                 << cp.w[2] << ','
//                 << cp.w.norm() << ','
//                 << cp.vel.x << ','
//                 << cp.vel.y << ','
//                 << cp.vel.z << ','
//                 << 0;
//         }

//         T += .5 * cp.m * cp.vel.normsquared() + .5 * cp.moi * cp.w.normsquared(); // Now includes rotational kinetic energy.
//         mom += cp.m * cp.vel;
//         ang_mom += cp.m * cp.pos.cross(cp.vel) + cp.moi * cp.w;
//     }
// }


/// @brief steps kinematic variables one time step given a P
void P_group::update_kinematics(P& P)
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


/// @brief steps velocites one time step given a P
void P_group::compute_velocity(P& P)
{
    // Velocity for next step:
    P.vel = P.velh + .5 * P.acc * dt;
    P.w = P.wh + .5 * P.aacc * dt;
}

/// @brief computes the energy associated with a given P and adds it to the total energy
void P_group::compute_energy(P& P)
{
    KE += .5 * P.m * P.vel.dot(P.vel) + .5 * P.moi * P.w.dot(P.w);
    mom += P.m * P.vel;
    ang_mom += P.m * P.pos.cross(P.vel) + P.moi * P.w;
}

void P_group::compute_distance(P_pair& p_pair)
{
    p_pair.dist = (p_pair.A->pos - p_pair.B->pos).norm();
}

/// @brief Computes the acceleration of a given P_pair based on Van der waals and collision forces
void P_group::compute_acceleration(P_pair& p_pair)
{
    if (true)
    // if (p_pair.dist <= 6*scaleBalls)
    {
        const double sumRaRb = p_pair.A->R + p_pair.B->R;
        const vec3 rVecab = p_pair.B->pos - p_pair.A->pos;  // Vector from a to b.
        const vec3 rVecba = -rVecab;
        const double dist = (rVecab).norm();

        // Check for collision between Ball and otherBall:
        double overlap = sumRaRb - dist;

        vec3 totalForceOnA{0, 0, 0};

        double oldDist = p_pair.dist;

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
            const double Ra = p_pair.A->R;
            const double Rb = p_pair.B->R;
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
            const vec3 r_a = rVecab * p_pair.A->R / sumRaRb;  // Center to contact point
            const vec3 r_b = rVecba * p_pair.B->R / sumRaRb;
            const vec3 w_diff = p_pair.A->w - p_pair.B->w;

            // Sliding friction terms:
            const vec3 d_vel = p_pair.B->vel - p_pair.A->vel;
            const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
                                       p_pair.A->w.cross(r_a) - p_pair.B->w.cross(r_a);

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
                    -u_r * elastic_force_A_mag * 
                    (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
            }

            // Total forces on a:
            //Took out gravity force
            totalForceOnA = elasticForceOnA + slideForceOnA + vdwForceOnA;
            // totalForceOnA = gravForceOnA + elasticForceOnA + slideForceOnA + vdwForceOnA;

            // Total torque a and b:
            torqueA = r_a.cross(slideForceOnA + rollForceA);
            torqueB = r_b.cross(-slideForceOnA + rollForceA);

            if (parallel)
            {
                {
                    const std::lock_guard<std::mutex> lock(g_mutex);
                    p_pair.A->aacc += torqueA / p_pair.A->moi;
                }
                {
                    const std::lock_guard<std::mutex> lock(g_mutex);
                    p_pair.B->aacc += torqueB / p_pair.B->moi;
                }
            }
            else
            {
                p_pair.B->aacc += torqueB / p_pair.B->moi;
                p_pair.A->aacc += torqueA / p_pair.A->moi;
            }

            if (writeStep) {
                // No factor of 1/2. Includes both spheres:
                // O.PE += -G * O.m[A] * O.m[B] / dist + 0.5 * k * overlap * overlap;

                // Van Der Waals + elastic:
                const double diffRaRb = p_pair.A->R - p_pair.B->R;
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * p_pair.A->R * p_pair.B->R;
                const double denom_sum = z * z - (sumRaRb * sumRaRb);
                const double denom_diff = z * z - (diffRaRb * diffRaRb);
                const double U_vdw =
                    -Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                PE += U_vdw + 0.5 * k * overlap * overlap;
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
            const double Ra = p_pair.A->R;
            const double Rb = p_pair.B->R;
            const double h2 = h * h;
            const double twoRah = 2 * Ra * h;
            const double twoRbh = 2 * Rb * h;
            const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                     ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                       (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                       (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                     rVecab.normalized();

            totalForceOnA = vdwForceOnA;  // +gravForceOnA;
            if (writeStep) {
                // O.PE += -G * O.m[A] * O.m[B] / dist; // Gravitational

                const double diffRaRb = p_pair.A->R - p_pair.B->R;
                const double z = sumRaRb + h;
                const double two_RaRb = 2 * p_pair.A->R * p_pair.B->R;
                const double denom_sum = z * z - (sumRaRb * sumRaRb);
                const double denom_diff = z * z - (diffRaRb * diffRaRb);
                const double U_vdw =
                    -Ha / 6 *
                    (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                PE += U_vdw;  // Van Der Waals
            }

            // todo this is part of push_apart. Not great like this.
            // For pushing apart overlappers:
            // O.vel[A] = { 0,0,0 };
            // O.vel[B] = { 0,0,0 };
        }

        // Newton's equal and opposite forces applied to acceleration of each ball:
        if (parallel)
        {
            {
                const std::lock_guard<std::mutex> lock(g_mutex);
                p_pair.A->acc += totalForceOnA / p_pair.A->m;
            }
            {
                const std::lock_guard<std::mutex> lock(g_mutex);
                p_pair.B->acc -= totalForceOnA / p_pair.B->m;
            }
        }
        else
        {
            p_pair.A->acc += totalForceOnA / p_pair.A->m;
            p_pair.B->acc -= totalForceOnA / p_pair.B->m;   
        }

        // So last distance can be known for COR:
        p_pair.dist = dist;
    }
} 


/// This version of compute_acceleration should be more efficient but has a bug
void P_group::compute_acceleration_speed_testing(P_pair& p_pair)
{
    const double Ra = p_pair.A->R;
    const double Rb = p_pair.B->R;
    const vec3 rVecab = p_pair.B -> pos - p_pair.A -> pos;
    const vec3 rVecba = -rVecab;
    const double m_a = p_pair.A->m;
    const double m_b = p_pair.B->m;
    const double sumRaRb = Ra + Rb;
    // vec3 rVec = p_pair.B->pos - p_pair.A->pos; // Start with rVec from a to b.
    const double dist = (rVecab).norm();

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
        rVecab.normalized();
    

    vec3 totalForce = vdwForce;


    // Check for collision between Ball and otherBall:
    double overlap = sumRaRb - dist;

    double oldDist = p_pair.dist;
    

    if (writeStep)
    {
        const double diffRaRb = Ra - Rb;
        const double z = sumRaRb + h;
        const double two_RaRb = 2 * Ra * Rb;
        const double denom_sum = z * z - (sumRaRb * sumRaRb);
        const double denom_diff = z * z - (diffRaRb * diffRaRb);
        const double U_vdw =
            -Ha / 6 *
            (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
        PE = PE + U_vdw;  // Van Der Waals potential
    }

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

        if (writeStep)
        {
            PE = PE + 0.5 * k * overlap * overlap;
        }

        // Elastic a:
        vec3 elasticForce = -k * overlap * .5 * (rVecab / dist);
        const double elastic_force_A_mag = elasticForce.norm();
        // Friction a:
        vec3 dVel = p_pair.B->vel - p_pair.A->vel;
        // vec3 frictionForce = { 0, 0, 0 };
        const vec3 r_a = rVecab * Ra / sumRaRb;  // Center to contact point
        const vec3 r_b = rVecba * Rb / sumRaRb;
        const vec3 w_diff = p_pair.A -> w - p_pair.B -> w;
        const double w_diff_mag = w_diff.norm();
        const vec3 relativeVelOfA = dVel - dVel.dot(rVecab) * (rVecab / (dist * dist)) - 
                                    p_pair.A->w.cross(r_a) - 
                                    p_pair.B->w.cross(r_b);
        double relativeVelMag = relativeVelOfA.norm();
        
        vec3 slideForceOnA, rollForceA;
        if (relativeVelMag > 1e-10) // When relative velocity is very low, dividing its vector components by its magnitude below is unstable.
        {
            slideForceOnA = u_s * (elastic_force_A_mag) *
                            (relativeVelOfA / relativeVelMag);
            // frictionForce = mu * (elasticForce.norm() + vdwForce.norm()) *
            //              (relativeVelOfA / relativeVelMag);
        }

        if (w_diff_mag > 1e-13)  // Divide by zero protection.
        {
            rollForceA = -u_r * elastic_force_A_mag * 
                        (w_diff).cross(r_a) / (w_diff).cross(r_a).norm();
        }

        // Torque a:
        const vec3 aTorque = r_a.cross(slideForceOnA + rollForceA);
        const vec3 bTorque = r_b.cross(-slideForceOnA + rollForceA);

        // Gravity on a:
        // const vec3 gravForceOnA = (G * p_pair.A->m * p_pair.B->m / (dist * dist)) * (rVec / dist);

        // Total forces on a:
        totalForce = totalForce + elasticForce + slideForceOnA;

        {
            const std::lock_guard<std::mutex> lock(g_mutex);
            p_pair.A->aacc += aTorque / p_pair.A->moi;
        }
        {
            const std::lock_guard<std::mutex> lock(g_mutex);
            p_pair.B->aacc += bTorque / p_pair.B->moi;
        }

    }

    // Newton's equal and opposite forces applied to acceleration of each ball:
    {
        const std::lock_guard<std::mutex> lock(g_mutex);
        p_pair.A->acc += totalForce / m_a;
    }
    {
        const std::lock_guard<std::mutex> lock(g_mutex);
        p_pair.B->acc -= totalForce / m_b;
    }
}

void
P_group::sim_one_step_single_core(const bool write_step)
{
    /// FIRST PASS - Update Kinematic Parameters:
    t.start_event("UpdateKinPar");
    for (int Ball = 0; Ball < num_particles; Ball++) {
        // Update velocity half step:
        p_group[Ball].velh = p_group[Ball].vel + .5 * p_group[Ball].acc * dt;

        // Update angular velocity half step:
        p_group[Ball].wh = p_group[Ball].w + .5 * p_group[Ball].aacc * dt;

        // Update position:
        p_group[Ball].pos += p_group[Ball].velh * dt;

        // Reinitialize acceleration to be recalculated:
        p_group[Ball].acc = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        p_group[Ball].aacc = {0, 0, 0};
    }
    t.end_event("UpdateKinPar");

    /// SECOND PASS - Check for collisions, apply forces and torques:
    t.start_event("CalcForces/loopApplicablepairs");
    for (int A = 1; A < num_particles; A++)  
    {
        /// DONT DO ANYTHING HERE. A STARTS AT 1.
        for (int B = 0; B < A; B++) {
            const double sumRaRb = p_group[A].R + p_group[B].R;
            const vec3 rVecab = p_group[B].pos - p_group[A].pos;  // Vector from a to b.
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

                // if (!contact && A == O.num_particles-1)
                // {
                //     // std::cout<<"CONTACT MADE"<<std::endl;
                //     contact = true;
                //     contactBuffer<<A<<','<<simTimeElapsed<<'\n';
                // }

                double k;
                if (dist >= oldDist) {
                    k = kout;
                } else {
                    k = kin;
                }

                // Cohesion (in contact) h must always be h_min:
                // constexpr double h = h_min;
                const double h = h_min;
                const double Ra = p_group[A].R;
                const double Rb = p_group[B].R;
                const double h2 = h * h;
                // constexpr double h2 = h * h;
                const double twoRah = 2 * Ra * h;
                const double twoRbh = 2 * Rb * h;
                const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                         ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                         rVecab.normalized();
                

        
                const vec3 elasticForceOnA = -k * overlap * .5 * (rVecab / dist);

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
                const vec3 r_a = rVecab * p_group[A].R / sumRaRb;  // Center to contact point
                const vec3 r_b = rVecba * p_group[B].R / sumRaRb;
                const vec3 w_diff = p_group[A].w - p_group[B].w;

                // Sliding friction terms:
                const vec3 d_vel = p_group[B].vel - p_group[A].vel;
                const vec3 frame_A_vel_B = d_vel - d_vel.dot(rVecab) * (rVecab / (dist * dist)) -
                                           p_group[A].w.cross(r_a) - p_group[B].w.cross(r_a);

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


                p_group[A].aacc += torqueA / p_group[A].moi;
                p_group[B].aacc += torqueB / p_group[B].moi;

                if (write_step) {
                    // No factor of 1/2. Includes both spheres:
                    // O.PE += -G * O.m[A] * O.m[B] * grav_scale / dist + 0.5 * k * overlap * overlap;
                    // O.PE += -G * O.m[A] * O.m[B] / dist + 0.5 * k * overlap * overlap;

                    // Van Der Waals + elastic:
                    const double diffRaRb = p_group[A].R - p_group[B].R;
                    const double z = sumRaRb + h;
                    const double two_RaRb = 2 * p_group[A].R * p_group[B].R;
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
                // const vec3 gravForceOnA = (G * O.m[A] * O.m[B] * grav_scale / (dist * dist)) * (rVecab / dist);
                const vec3 gravForceOnA = {0.0,0.0,0.0};
                // Cohesion (non-contact) h must be positive or h + Ra + Rb becomes catastrophic cancellation:
                double h = std::fabs(overlap);
                if (h < h_min)  // If h is closer to 0 (almost touching), use hmin.
                {
                    h = h_min;
                }
                const double Ra = p_group[A].R;
                const double Rb = p_group[B].R;
                const double h2 = h * h;
                const double twoRah = 2 * Ra * h;
                const double twoRbh = 2 * Rb * h;
                const vec3 vdwForceOnA = Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb *
                                         ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb) *
                                                           (h2 + twoRah + twoRbh + 4 * Ra * Rb))) *
                                         rVecab.normalized();

                totalForceOnA = vdwForceOnA + gravForceOnA;

                if (write_step) {
                    // O.PE += -G * O.m[A] * O.m[B] * grav_scale / dist; // Gravitational

                    const double diffRaRb = p_group[A].R - p_group[B].R;
                    const double z = sumRaRb + h;
                    const double two_RaRb = 2 * p_group[A].R * p_group[B].R;
                    const double denom_sum = z * z - (sumRaRb * sumRaRb);
                    const double denom_diff = z * z - (diffRaRb * diffRaRb);
                    const double U_vdw =
                        -Ha / 6 *
                        (two_RaRb / denom_sum + two_RaRb / denom_diff + log(denom_sum / denom_diff));
                    PE += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
                }

                // todo this is part of push_apart. Not great like this.
                // For pushing apart overlappers:
                // O.vel[A] = { 0,0,0 };
                // O.vel[B] = { 0,0,0 };
            }

            // Newton's equal and opposite forces applied to acceleration of each ball:
            p_group[A].acc += totalForceOnA / p_group[A].m;
            p_group[B].acc -= totalForceOnA / p_group[B].m;

            // So last distance can be known for COR:
            distances[e] = dist;

        }
        // DONT DO ANYTHING HERE. A STARTS AT 1.
    }

   

    t.end_event("CalcForces/loopApplicablepairs");

    if (write_step) {
        ballBuffer << '\n';  // Prepares a new line for incoming data.
        // std::cerr<<"Writing "<<O.num_particles<<" balls"<<std::endl;
    }

    // THIRD PASS - Calculate velocity for next step:
    t.start_event("CalcVelocityforNextStep");
    for (int Ball = 0; Ball < num_particles; Ball++) {
        // Velocity for next step:
        p_group[Ball].vel = p_group[Ball].velh + .5 * p_group[Ball].acc * dt;
        p_group[Ball].w = p_group[Ball].wh + .5 * p_group[Ball].aacc * dt;

        /////////////////////////////////
        // if (true) {
        /////////////////////////////////
        if (write_step) {
            // Send positions and rotations to buffer:
            // std::cerr<<"Write ball "<<Ball<<std::endl;
            if (Ball == 0) {
                ballBuffer << p_group[Ball].pos[0] << ',' << p_group[Ball].pos[1] << ',' << p_group[Ball].pos[2] << ','
                           << p_group[Ball].w[0] << ',' << p_group[Ball].w[1] << ',' << p_group[Ball].w[2] << ','
                           << p_group[Ball].w.norm() << ',' << p_group[Ball].vel.x << ',' << p_group[Ball].vel.y << ','
                           << p_group[Ball].vel.z << ',' << 0;
            } else {
                ballBuffer << ',' << p_group[Ball].pos[0] << ',' << p_group[Ball].pos[1] << ',' << p_group[Ball].pos[2] << ','
                           << p_group[Ball].w[0] << ',' << p_group[Ball].w[1] << ',' << p_group[Ball].w[2] << ','
                           << p_group[Ball].w.norm() << ',' << p_group[Ball].vel.x << ',' << p_group[Ball].vel.y << ','
                           << p_group[Ball].vel.z << ',' << 0;
            }

            KE += .5 * p_group[Ball].m * p_group[Ball].vel.normsquared() +
                    .5 * p_group[Ball].moi * p_group[Ball].w.normsquared();  // Now includes rotational kinetic energy.
            mom += p_group[Ball].m * p_group[Ball].vel;
            ang_mom += p_group[Ball].m * p_group[Ball].pos.cross(p_group[Ball].vel) + p_group[Ball].moi * p_group[Ball].w;
        }
    }  // THIRD PASS END
    t.end_event("CalcVelocityforNextStep");
}  // one Step end



/// @brief steps the sim forward by one time step and writes out data if bool write_step is true
void P_group::sim_one_step(const bool write_step)
{
    // make_pairs(); //This parallel block seems to be very slow
    //uncomment following line to use grid
    // getRelaventPairs();
    ///////////////////

    if (parallel)
    {
        std::for_each(std::execution::par_unseq, p_group.begin(), p_group.end(),
                    std::bind_front(&P_group::update_kinematics, this)); //0.07hr/step
        std::for_each(std::execution::par, pairs.begin(), pairs.end(),
                    std::bind_front(&P_group::compute_acceleration, this)); //0.38hr/step
        // std::for_each(std::execution::par, pairs.begin(), pairs.end(),
                    // std::bind_front(&P_group::compute_distance, this)); //0.38hr/step
        std::for_each(std::execution::par_unseq, p_group.begin(), p_group.end(),
                    std::bind_front(&P_group::compute_velocity, this)); //0.05hr/step
    }
    else
    {
        std::for_each(std::execution::unseq, p_group.begin(), p_group.end(),
                    std::bind_front(&P_group::update_kinematics, this));
        std::for_each(std::execution::seq, pairs.begin(), pairs.end(),
                    std::bind_front(&P_group::compute_acceleration, this));
                    // std::bind_front(&P_group::compute_acceleration_speed_testing, this));
        // std::for_each(std::execution::seq, pairs.begin(), pairs.end(),
                    // std::bind_front(&P_group::compute_distance, this));
        std::for_each(std::execution::unseq, p_group.begin(), p_group.end(),
                    std::bind_front(&P_group::compute_velocity, this));    
    }
    if (write_step)
    {
        t.start_event("WriteToBuffer");
        ballBuffer << '\n';  // Prepares a new line for incoming data.
        for (auto p : p_group)  
        {  
            if (p != *p_group.begin()) 
            {
                ballBuffer << ',' << p.pos[0] << ',' << p.pos[1] << ',' << p.pos[2] << ','
                           << p.w[0] << ',' << p.w[1] << ',' << p.w[2] << ','
                           << p.w.norm() << ',' << p.vel.x << ',' << p.vel.y << ','
                           << p.vel.z << ',' << 0;
            } 
            else 
            {
                ballBuffer << p.pos[0] << ',' << p.pos[1] << ',' << p.pos[2] << ','
                           << p.w[0] << ',' << p.w[1] << ',' << p.w[2] << ','
                           << p.w.norm() << ',' << p.vel.x << ',' << p.vel.y << ','
                           << p.vel.z << ',' << 0;
            }

            KE += .5 * p.m * p.vel.normsquared() +
                    .5 * p.moi * p.w.normsquared();  // Now includes rotational kinetic energy.
            mom += p.m * p.vel;
            ang_mom += p.m * p.pos.cross(p.vel) + p.moi * p.w;
        }
        t.end_event("WriteToBuffer");
    }
}


/// @brief loops through a whole simulation keeping track of timing
void P_group::sim_looper()
{
    std::cerr << "Beginning simulation...\n";

    startProgress = time(nullptr);

    for (int Step = 1; Step < steps; Step++)  // Steps start at 1 because the 0 step is initial conditions.
    {
        // Check if this is a write step:
        if (Step % skip == 0) {
            t.start_event("writeProgressReport");
            writeStep = true;

            simTimeElapsed += dt * skip;

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
            t.end_event("writeProgressReport");
        } else {
            writeStep = false;
        }

        // Physics integration step:
        t.start_event("OneStep");
        sim_one_step_single_core(writeStep);
        // sim_one_step(writeStep);
        t.end_event("OneStep");
        // if (writeStep)
        // {
        //     t.print_events();
        //     exit(0);
        // }


        if (writeStep) {
            // Write energy to stream:
            energyBuffer << '\n'
                         << simTimeElapsed << ',' << PE << ',' << KE << ',' << PE + KE << ','
                         << mom.norm() << ','
                         << ang_mom.norm();  // the two zeros are bound and unbound mass

            // Reinitialize energies for next step:
            KE = 0;
            PE = 0;
            mom = {0, 0, 0};
            ang_mom = {0, 0, 0};
            // unboundMass = 0;
            // boundMass = massTotal;

            // Data Export. Exports every 10 writeSteps (10 new lines of data) and also if the last write was
            // a long time ago.
            // if (time(nullptr) - lastWrite > 1800 || Step / skip % 10 == 0) {
            if (Step / skip % 10 == 0) {
                // Report vMax:
                std::cerr<<"Step: "<<Step<<std::endl;
                std::cerr<<"skip: "<<skip<<std::endl;
                std::cerr << "vMax = " << getVelMax() << " Steps recorded: " << Step / skip << '\n';
                std::cerr << "Data Write to "<<s_location<<"\n";


                // Write simData to file and clear buffer.
                t.start_event("WriteToFile");
                std::ofstream ballWrite;

                ballWrite.open(s_location + filenum + output_prefix + "simData.csv", std::ofstream::app);
                ballWrite << ballBuffer.rdbuf();  // Barf buffer to file.
                ballBuffer.str("");               // Empty the stream for next filling.
                ballWrite.close();

                // Write Energy data to file and clear buffer.
                std::ofstream energyWrite;
                energyWrite.open(s_location + filenum + output_prefix + "energy.csv", std::ofstream::app);
                energyWrite << energyBuffer.rdbuf();
                energyBuffer.str("");  // Empty the stream for next filling.
                energyWrite.close();
                t.end_event("WriteToFile");

                lastWrite = time(nullptr);
            }  // Data export end


            if (dynamicTime) { calibrate_dt(Step, false); }
        }  // writestep end
    }
    // t.print_events();
    const time_t end = time(nullptr);

    std::cerr << "Simulation complete!\n"
              << num_particles << " Particles and " << steps << " Steps.\n"
              << "Simulated time: " << steps * dt << " seconds\n"
              << "Computation time: " << end - start << " seconds\n";
    std::cerr << "\n===============================================================\n";
}  // end simLooper

void P_group::makeSpheres(const int nBalls, distributions radiiDistribution)
{
    for (int Ball = 0; Ball < nBalls; ++Ball) 
    {
        P new_p;
        new_p.unique_id = get_unique_id();
        if (radiiDistribution == constant)
        {
            new_p.R = scaleBalls;
        }
        else
        {
            new_p.R = lognorm_dist(scaleBalls*std::exp(-5*std::pow(lnSigma,2)/2),lnSigma);
        }
        new_p.m = density * 4. / 3. * 3.14159 * std::pow(new_p.R, 3);
        new_p.moi = .4 * new_p.m * new_p.R * new_p.R;
        new_p.w = {0,0,0};
        new_p.vel = {0,0,0}; //?????
        new_p.acc = {0,0,0}; //?????
        new_p.aacc = {0,0,0}; //?????
        new_p.pos = rand_vec3(spaceRange);
        p_group.push_back(new_p);
    }

    mTotal = getMass();

    placeBalls(nBalls);
}

// void P_group::distSizeSphere(const int nBalls)
// {
//     for (int Ball = 0; Ball < nBalls; ++Ball) 
//     {
//         P new_p;
//         new_p.unique_id = get_unique_id();
//         new_p.R = lognorm_dist(scaleBalls*std::exp(-5*std::pow(lnSigma,2)/2),lnSigma);
//         new_p.m = density * 4. / 3. * 3.14159 * std::pow(new_p.R, 3);
//         new_p.moi = .4 * new_p.m * new_p.R * new_p.R;
//         new_p.w = {0,0,0};
//         new_p.vel = {0,0,0}; //?????
//         new_p.pos = rand_vec3(spaceRange);
//         p_group.push_back(new_p);
//     }

//     mTotal = getMass();

//     placeBalls(nBalls);
// }

// ///@brief makes a new P_group, all with the same sized radius, given the number of balls in the group
// void P_group::oneSizeSphere(const int nBalls)
// {

//     for (int Ball = 0; Ball < nBalls; ++Ball) 
//     {
//         P new_p;
//         new_p.unique_id = get_unique_id();
//         new_p.R = scaleBalls;
//         new_p.m = density * 4. / 3. * 3.14159 * std::pow(new_p.R, 3);
//         new_p.moi = .4 * new_p.m * new_p.R * new_p.R;
//         new_p.w = {0,0,0};
//         new_p.vel = {0,0,0}; //?????
//         new_p.pos = rand_vec3(spaceRange);
//         p_group.push_back(new_p);
//     }

//     mTotal = getMass();

//     placeBalls(nBalls);
    
// }


void P_group::placeBalls(const int nBalls)
{
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


/// @brief parses an input file, setting relavent values for the sim
/// @param location is the path to the input file
void P_group::parse_input_file(char const* location)
{
    std::string loc(location);
    s_location = loc;
    std::string json_file = loc + "input.json";
    std::ifstream ifs(json_file);
    json inputs = json::parse(ifs);

    // if (inputs["parallel"] == "true")
    // {
    //     parallel = true;
    //     safe = std::execution::par;
    //     danger = std::execution::par_unseq;
    // }
    // else
    // {
    //     parallel = false;
    //     safe = std::execution::seq;
    //     danger = std::execution::unsequenced_policy;
    // }

    if (inputs["seed"] == std::string("default"))
    {
        seed = static_cast<int>(time(nullptr));
    }
    else
    {
        seed = static_cast<int>(inputs["seed"]);
        random_generator.seed(seed);
    }
    srand(seed);  // srand(seedSave);

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
    if (inputs["v_custom"] == std::string("default"))
    {
        v_custom = 0.36301555459799423;
    }
    else
    {
        v_custom = inputs["v_custom"];
    }
    temp = inputs["temp"]; // this will modify v_custom in dust_agglomeration_particle_init
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
    num_particles = genBalls;
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

/// @brief returns approximate radius of the P_group given a center.
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


/// @brief returns center of mass vector of the P_group
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


/// @brief calculates the vecocity a P_group will collapse at 
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


/// @brief determins how many time steps are needed for a given simulation to end based on simulation parameters
/// @param Step is what step the simulation is currently on. Useful so sim knows how many steps should be left
/// @param customSpeed will update time step parameters if you want to have a custom speed
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


/// @brief updates the time step length and K (the constant of restitution)
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

/// @brief returns max velocity of a P from the P_group
[[nodiscard]] double P_group::getVelMax()
{
    vMax = 0;

    // todo - make this a manual set true or false to use soc so we know if it is being used or not.
    // if (false) {
    std::cerr<<"SOC: "<<soc<<std::endl;
    if (soc > 0) {
        std::cerr<<"WE NOOOOT HERERER"<<std::endl;
        int counter = 0;
        for (int Ball = 0; Ball < num_particles; Ball++) {
            if (p_group[Ball].vel.norm() > vMax) 
            { 
                vMax = p_group[Ball].vel.norm(); 
                std::cerr<<"V_MAX for ball "<<Ball<<" = "<<vMax<<std::endl; 
                
            }
            /////////////////SECTION COMMENTED FOR ACCURACY TESTS
            // Only consider balls moving toward com and within 4x initial radius around it.
            // const vec3 fromCOM = p_group[Ball].pos - getCOM();
            // if (acos(p_group[Ball].vel.normalized().dot(fromCOM.normalized())) > cone && fromCOM.norm() < soc) {
            //     if (p_group[Ball].vel.norm() > vMax) { vMax = p_group[Ball].vel.norm(); }
            // } else {
            //     counter++;
            // }
        }
        std::cerr << '(' << counter << " spheres ignored"
                  << ") ";
    } else {
        std::cerr<<"WE IN HERERER"<<std::endl;
        for (int Ball = 0; Ball < num_particles; Ball++) {
            if (p_group[Ball].vel.norm() > vMax) 
            { 
                vMax = p_group[Ball].vel.norm(); 
                std::cerr<<"V_MAX for ball "<<Ball<<" = "<<vMax<<std::endl; 
                
            }
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

//TODO make next two functions parallel
/// @brief returns the smallest radius of any P from P_group
[[nodiscard]] double P_group::getRmin()
{
    rMin = p_group[0].R;
    for (int Ball = 1; Ball < num_particles; Ball++) {
        if (p_group[Ball].R < rMin) { rMin = p_group[Ball].R; }
    }
    return rMin;
}

/// @brief returns the largest radius of any P from P_group
[[nodiscard]] double P_group::getRmax()
{
    rMax = p_group[0].R;
    for (int Ball = 0; Ball < num_particles; Ball++) {
        if (p_group[Ball].R > rMax) { rMax = p_group[Ball].R; }
    }
    return rMax;
}


/// @brief calculates helpful parameters for a simulation. 
///         Can be called at the beginning of a sim or when stuff changes (after a time step is completed)
void P_group::calc_helpfuls()
{
    rMin = getRmin();
    rMax = getRmax();
    mTotal = getMass();
    initialRadius = get_radius(getCOM());
    soc = 4 * rMax + initialRadius;
}


void P_group::simInit_cond_and_center()
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
    init_conditions_single_core();
    // init_conditions();
}


void P_group::getRelaventPairs()
{
    applicable_pairs.clear();
    used_ids.clear();
    findGroups();
    mapGroups();
    for (auto pt : p_group)
    {
        findPairs(pt);
    }
    if (applicable_pairs.size() == 0)
    {
        std::iota (std::begin(applicable_pairs), std::end(applicable_pairs), 0); // Fill with 0, 1, ..., 99.;
    }
    // printPairs();
}

void P_group::printPairs()
{
    std::cout<<"Start printPairs\n";
    for (P_pair pair : pairs)
    {
        std::cout<<"("<<pair.A->unique_id<<", "<<pair.B->unique_id<<")\n";
    } 
    std::cout<<"End printPairs"<<std::endl;
}

/// @brief makes an unordered_map where the key is the grid space and it maps to all P in that grid space
void P_group::mapGroups()
{
    for (auto pt : p_group)
    {
        std::string key = getKey(pt.id);
        if (IDToGrid.find(key) != IDToGrid.end()) // key present
        {
            IDToGrid[key].push_back(pt);
        }
        else//key not present
        {
            std::vector<P> points;
            IDToGrid.insert({key, points});
        }
    }
    // printMap();
    return;
}

inline std::string P_group::getKey(std::vector<int> v)
{
    return std::to_string(v[0]) + std::to_string(v[1]) + std::to_string(v[2]);
}

inline std::string P_group::getKey(int x, int y, int z)
{
    return std::to_string(x) + std::to_string(y) + std::to_string(z);
}


/// @brief returns a vector of all the balls relavent (close enough) to a given ball
std::vector<P> P_group::getBalls(P ball)
{

    int currx, curry, currz;
    currx = ball.id[0];
    curry = ball.id[1]; 
    currz = ball.id[2]; 

    std::vector<P> points;
    for (int x = -1; x < 2; ++x)
    {
        for (int y = -1; y < 2; ++y)
        {
            for (int z = -1; z < 2; ++z)
            {
                std::string key = getKey(currx+x,curry+y,currz+z);
                if (IDToGrid.find(key) != IDToGrid.end())
                {
                    // for (auto ind_it = begin(IDToGrid[key]); ind_it != end(IDToGrid[key]); ++ind_it)
                    for (auto pt : IDToGrid[key])
                    {
                        // if (ball_indicies.find(*ind_it) == ball_indicies.end())
                        if(std::find(std::execution::par_unseq, points.begin(), points.end(), pt) == points.end())
                        {
                            points.push_back(pt);
                        }
                    }
                }       
            }
        }
    }       
    return points;
}


/// @brief sets a ball's (P) grid id 
void P_group::findGroup(P &p)
{
    p.id[0] = floor(p.pos[0]/gridSize);
    p.id[1] = floor(p.pos[1]/gridSize);
    p.id[2] = floor(p.pos[2]/gridSize);
}


/// @brief sets all P's groups of the P_group
void P_group::findGroups()
{
    std::for_each(std::execution::par_unseq,p_group.begin(),p_group.end(),
                std::bind_front(&P_group::findGroup, this));
    // #pragma omp parallel
    // {
    //     #pragma omp parallel for default(none) shared(gridIDs)
    //     for (int i = 0; i < numBalls; ++i)
    //     {   
    //         gridIDs[i] = std::vector<int>(3);
    //         for (int j = 0; j < 3; j++)
    //         {
    //             gridIDs[i][j] = floor(pos[i][j]/gridSize);
    //         }
    //     }
    // }
    // return;
}


/// @brief makes relavent particle pairs based on how close balls are to the given P
void P_group::findPairs(P &p)
{
    std::vector<P> points = getBalls(p);
    for (auto point : points)
    {
        if (p.unique_id != point.unique_id)
        {
            if (std::find(used_ids.begin(), used_ids.end(), point.unique_id) == used_ids.end())
            {
                double dist = (point.pos - p.pos).norm();
                P_pair pair{ &p_group[p.unique_id], &p_group[point.unique_id], dist };
                pairs.push_back(pair); 
            }
        }
    }    
    used_ids.push_back(p.unique_id);
}


/// @brief makes ALL possible pairs of P's 
void P_group::make_pairs() 
{
    // pairs.clear(); // All particle pairs

    int n = p_group.size();
    int n_pairs = n * (n - 1) / 2;
    std::vector<P_pair> p_pairs(n_pairs); // All particle pairs
    for (size_t i = 0; i < n_pairs; i++)
    {
        // Pair Combinations [A,B] [B,C] [C,D]... [A,C] [B,D] [C,E]... ...
        int A = i % n;
        int stride = 1 + i / n; // Stride increases by 1 after each full set of pairs
        int B = (A + stride) % n;

        // Create particle* pair
        double dist = (p_group[B].pos - p_group[A].pos).norm();
        p_pairs[i] = { &p_group[A], &p_group[B], dist };
        // p_pairs[i] = { &p_group[A], &p_group[B] };
    }
    pairs = p_pairs;
}


/// @brief continues a sim given the path to the sim and what index you want it to restart at
void P_group::sim_continue(const std::string& path, int start_file_index)
{
    // Load file data:
    // num_particles = 3 + start_file_index;
    if (start_file_index == 0)
    {
        std::cerr << "Continuing Sim...\nFile: " << path << start_file_index << "_*" << '\n';
        loadSim(path,"");
    }
    else
    {
        // output_prefix = std::to_string(start_file_index);
        std::cerr << "Continuing Sim...\nFile: " << start_file_index << '_' << '\n';
        loadSim(path, std::to_string(start_file_index) + '_');
    }

    std::cerr << '\n';
    calc_momentum("O");
}


//TODO remove parameter filename
/// @brief initializes file streams and writes initial data to them
/// @param filename doesn't do anything anymore
/// (optional) @param counter is the file index the sim should start at for continuing a simulation.  
void P_group::sim_init_write(std::string filename, int counter)
{
    std::ifstream checkForFile;
    if (counter != 0)
    {
        filenum = std::to_string(counter) + '_';
    }
    else
    {
        filenum = "";
    }

    checkForFile.open(s_location + filenum + output_prefix + "simData.csv", std::ifstream::in);
    
    // Add a counter to the file name until it isn't overwriting anything:
    while (checkForFile.is_open()) {
        counter++;
        checkForFile.close();
        checkForFile.open(s_location + std::to_string(counter) + '_' + output_prefix + "_simData.csv", std::ifstream::in);
    }

    // Complete file names:
    std::string simDataFilename = s_location + filenum + output_prefix + "simData.csv";
    std::string energyFilename = s_location + filenum + output_prefix + "energy.csv";
    std::string constantsFilename = s_location + filenum + output_prefix + "constants.csv";
    // std::cout<<s_location + filenum + output_prefix + "simData.csv"<<std::endl;

    std::cerr << "New file tag: " << filenum + output_prefix;

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
        constWrite << p_group[Ball].R << ',' << p_group[Ball].m << ',' << 
                    p_group[Ball].moi << '\n';
    }

    // Write energy data to buffer:
    energyBuffer << '\n'
                 << simTimeElapsed << ',' << PE << ',' << KE << ','
                 << PE + KE << ',' << mom.norm() << ','
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
    ballBuffer << p_group[0].pos.x << ',' << p_group[0].pos.y << ',' << p_group[0].pos.z << ',' 
               << p_group[0].w.x << ',' << p_group[0].w.y << ','
               << p_group[0].w.z << ',' << p_group[0].w.norm() << ',' 
               << p_group[0].vel.x << ',' << p_group[0].vel.y << ',' << p_group[0].vel.z
               << ',' << 0;  // bound[0];
    for (int Ball = 1; Ball < num_particles; Ball++) {
        ballBuffer << ',' << p_group[Ball].pos.x
                   << ','  // Needs comma start so the last bound doesn't have a dangling comma.
                   << p_group[Ball].pos.y << ',' << p_group[Ball].pos.z << ',' 
                   << p_group[Ball].w.x << ',' << p_group[Ball].w.y << ',' << p_group[Ball].w.z 
                   << ',' << p_group[Ball].w.norm() << ',' 
                   << p_group[Ball].vel.x << ',' << p_group[Ball].vel.y
                   << ',' << p_group[Ball].vel.z << ',' << 0;  // bound[Ball];
    }
    // Write position and rotation data to file:
    ballWrite << ballBuffer.rdbuf();
    ballBuffer.str("");  // Resets the stream buffer to blank.

    // Close Streams for user viewing:
    energyWrite.close();
    ballWrite.close();
    constWrite.close();

    std::cerr << "\nSimulating " << steps * dt / 60 / 60 << " hours.\n";
    std::cerr << "Total mass: " << mTotal << '\n';
    std::cerr << "\n===============================================================\n";
}

void P_group::write_timing()
{
    t.save_events(s_location+"ballSim_timing.txt");
}