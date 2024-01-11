
#include <omp.h>
#include <openacc.h>
#include <fstream>
#include "vec3.hpp"
#include "../SpaceLab/external/json/single_include/nlohmann/json.hpp"
#include "Utils.hpp"
#include "linalg.hpp"

// double PE = 0.0;
const double Ha = 4.7e-12;
double u_r = 1e-5;
double u_s = 0.1;
double kin = 3.7;
double kout = 1.5;
double h_min = 0.1;
double dt = 1e-5;

using json = nlohmann::json;


class lass
{
public:
	int particles = 400;
	double *R = new double[particles];
	double *m = new double[particles];
	double *moi = new double[particles];
	vec3 *acc = new vec3[particles];
	vec3 *aacc = new vec3[particles];
	vec3 *pos = new vec3[particles];
	vec3 *vel = new vec3[particles];
	vec3 *w = new vec3[particles];
	int num_pairs = (((particles*particles)-particles)/2);
	double *distances = new double[num_pairs];

	int num_particles = particles;
	bool write_step;

	int accum = 0;

    // double *PE = new double[1];
	double PE = 0.0;
    int world_rank,world_size;

	void init();
	void tofu();
	void parse_input_file(char const* location);
    void loop_one_step(bool write_step);
};

void lass::init()
{
	for (int i = 0; i < particles; ++i)
	{
		acc[i] = {0.0,0.0,0.0};
		aacc[i] = {0.0,0.0,0.0};
		pos[i] = {static_cast<double>(i/particles),static_cast<double>((i+1)/particles),static_cast<double>((i+2)/particles)};
		vel[i] = {static_cast<double>(i/particles),static_cast<double>((i-1)/particles),static_cast<double>((i-2)/particles)};
		w[i] = {static_cast<double>(i/particles),static_cast<double>((i*2)/particles),static_cast<double>((i*3)/particles)};
		R[i] = 1e-5;
		m[i] = 7.07e-10;
		moi[i] = (2/5)*m[i]*R[i]*R[i];
	}

	for (int i = 0; i < num_pairs; i++)
	{
		distances[i] = i*1.5;
	}
    PE = 0.0;
	// PE[0] = 0.0;
}


void lass::tofu()
{
	int outerLoop = 2000;

	// int pc;

	// int A,B;
	
	
	world_rank = 0;
	world_size = 1;

    

	
    double t0 = omp_get_wtime();
    // #pragma acc data reduction(+:PE) copyin(this) copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copy(acc[0:num_particles],aacc[0:num_particles]) copyin(A,B,pc,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    #pragma acc enter data copyin(this) 
    #pragma acc enter data copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) 
    #pragma acc enter data copyin(acc[0:num_particles],aacc[0:num_particles]) 
    #pragma acc enter data copyin(num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r)
    // #pragma acc enter data create(A,B,pc,world_rank,world_size)
    // {                                                   
	for (int k = 0; k < outerLoop; k++)
	{
        // printf("HREERREERERE");
        if(k % 5 == 0)
        {
            write_step = true;
        }
        else
        {
            write_step = true;
        }

        // #pragma update device(write_step)
	    
        loop_one_step(write_step);
        // #pragma acc update host(pe)
        // PE = pe;
        // std::cerr<<"PE: "<<PE<<std::endl;
        // std::cerr<<"pe: "<<pe<<std::endl;
    }

    // }


    // #pragma acc update device()
	for (int i = 0; i < 1; i++)
	{
		// acc[i] = {0.0,0.0,0.0};
		std::cout<<acc[i]<<std::endl;
		std::cout<<aacc[i]<<std::endl;
	}
	
    #pragma acc exit data delete(acc[0:num_particles],aacc[0:num_particles])
	#pragma acc exit data delete(m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs])
	#pragma acc exit data delete(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
	#pragma acc exit data delete(this)
	// delete[] acc;
	// free(dacc);
	
	
	double t1 = omp_get_wtime();
	std::cerr<<"GPU took "<<t1-t0<<" seconds"<<std::endl;
    std::cerr<<"PE REDUCTION: "<<PE<<std::endl;
	// std::cerr<<"PE REDUCTION: "<<PE[0]<<std::endl;
	// std::cerr<<"TEST"<<std::endl;
}
// #pragma omp end declare target
// #pragma omp end declare target

void lass::loop_one_step(bool writeStep)
{
    
    double pe = 0.0;
    // vec3 *ACC = new vec3[particles];
    // for (int i = 0; i < particles; i++)
    // {
    //     ACC[i] = acc[i];
    // }
    // vec3 *AACC = new vec3[particles];
    // #pragma acc declare copyin(pe)
    #pragma acc enter data copyin(pe)
    #pragma acc enter data copyin(writeStep)
    // std::cerr<<"CONFIRMED"<<std::endl;
    // #pragma acc parallel loop gang worker  num_gangs(100) num_workers(60) private(pc,A,B) present_or_copy(PE,acc,aacc,m,w,vel,pos,R,distances) present_or_copy(this,pc,A,B,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)  
    // #pragma acc kernels private(pc,A,B)// reduction(+:PE[0:1]) //present(this,acc,aacc,PE,m,w,vel,pos,R,distances,num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)   
    // #pragma acc parallel loop gang reduction(+:PE[0:1]) private(pc,A,B) present_or_copy(PE[0:1],acc[0:num_particles],aacc[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) present_or_copy(this,pc,A,B,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)  
    // #pragma acc parallel loop gang worker private(pc,A,B) reduction(+:PE[0:1]) copyin(this) copy(acc[0:num_particles],aacc[0:num_particles]) copyin(m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copyin(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    // #pragma acc parallel loop gang worker reduction(+:PE[0:1]) num_workers(50) num_gangs(50) private(pc,A,B) copyin(this) copy(PE[0:1],acc[0:num_particles],aacc[0:num_particles]) copyin(m[0:num_particles],moi[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copyin(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    // #pragma acc parallel loop gang worker reduction(+:PE) num_workers(100) num_gangs(50) private(pc,A,B) copyin(this) copy(PE,acc[0:num_particles],aacc[0:num_particles]) copyin(m[0:num_particles],moi[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copyin(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    // { 
        // #pragma acc loop gang worker 
    #pragma acc parallel loop num_gangs(108) num_workers(256) reduction(+:pe) present(this,acc,aacc,m,moi,w,vel,pos,R,distances,num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,writeStep,world_rank,world_size)
    for (int pc = world_rank + 1; pc <= num_pairs; pc += world_size)
    {
        // #pragma acc atomic
        if (writeStep)
        {
            pe += 1;
        }
        // long double pd = (long double)pc;
        // pd = (sqrt(pd*8.0L+1.0L)+1.0L)*0.5L;
        // pd -= 0.00001L;
        // A = (long long)pd;
        // B = (long long)((long double)pc-(long double)A*((long double)A-1.0L)*.5L-1.0L);
        double pd = (double)pc;
        pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
        pd -= 0.00001;
        int A = (int)pd;
        int B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);

        const double sumRaRb =  R[A] + R[B];
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
        // int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
        int e = pc-1;
        double oldDist = distances[e];

        // Check for collision between Ball and otherBall.
        // if (overlap > 0) {
        if (true) {

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

            // for (int i = 0; i < 3; i ++)
            // {
            //     #pragma omp atomic
            //         aacc[A][i] += aaccA[i];
            //     #pragma omp atomic
            //         aacc[B][i] += aaccB[i];

            // }
            // #pragma acc atomic
            //     aacc[A].x += aaccA.x;
            // #pragma acc atomic
            //     aacc[A].y += aaccA.y;
            // #pragma acc atomic
            //     aacc[A].z += aaccA.z;
            // #pragma acc atomic
            //     aacc[B].x += aaccB.x;
            // #pragma acc atomic
            //     aacc[B].y += aaccB.y;
            // #pragma acc atomic
            //     aacc[B].z += aaccB.z;

            #pragma acc atomic
                aacc[A].x += 1;
            #pragma acc atomic
                aacc[A].y += 1;
            #pragma acc atomic
                aacc[A].z += 1;
            #pragma acc atomic
                aacc[B].x += 1;
            #pragma acc atomic
                aacc[B].y += 1;
            #pragma acc atomic
                aacc[B].z += 1;
            // }
            // omp_unset_lock(&writelock);


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
                // #pragma omp critical
                // #pragma acc atomic
                // PE += 1; ///TURN ON FOR REAL SIM
                // std::cerr<<"PE: "<<PE<<std::endl;
                // printf("HERE");
                // PE[0] += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
                // pe += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
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
                // #pragma omp critical
                // printf("HERE1");
                // PE += 1.0;  // Van Der Waals TURN ON FOR REAL SIM
                // pe += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
                // PE[0] += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
            }

            // todo this is part of push_apart. Not great like this.
            // For pushing apart overlappers:
            // vel[A] = { 0,0,0 };
            // vel[B] = { 0,0,0 };
        }


        vec3 accA = (1/m[A])*totalForceOnA; 
        vec3 accB = (1/m[B])*totalForceOnA; 
        // #pragma acc atomic
        //     acc[A].x += accA.x;
        // #pragma acc atomic
        //     acc[A].y += accA.y;
        // #pragma acc atomic
        //     acc[A].z += accA.z;
        // #pragma acc atomic
        //     acc[B].x -= accB.x;
        // #pragma acc atomic
        //     acc[B].y -= accB.y;
        // #pragma acc atomic
        //     acc[B].z -= accB.z;

        // #pragma acc atomic
        //     ACC[A].x += 1;
        // #pragma acc atomic
        //     ACC[A].y += 1;
        // #pragma acc atomic
        //     ACC[A].z += 1;
        // #pragma acc atomic
        //     ACC[B].x += 1;
        // #pragma acc atomic
        //     ACC[B].y += 1;
        // #pragma acc atomic
        //     ACC[B].z += 1;

        #pragma acc atomic
            acc[A].x += 1;
        #pragma acc atomic
            acc[A].y += 1;
        #pragma acc atomic
            acc[A].z += 1;
        #pragma acc atomic
            acc[B].x += 1;
        #pragma acc atomic
            acc[B].y += 1;
        #pragma acc atomic
            acc[B].z += 1;

        distances[e] = dist;

        // #pragma acc update device(PE)
    // #pragma acc update self(acc[0:particles],aacc[0:particles],PE) //if(write_step)
    // #pragma acc update self(PE) //if(write_step)
    // #pragma acc update host(ACC[0:particles])
        #pragma acc update host(pe) //WHY DOES THIS LINE MAKE EVERYTHING WORK
    }
        // #pragma acc update self(acc[0:particles],aacc[0:particles],PE) //if(write_step)
            // #pragma acc serial present(pe)
            // {
    PE = pe;
    // acc = ACC;
    // std::cerr<<ACC[0].x<<std::endl;
    // delete[] acc;
    #pragma acc exit data delete(pe,writeStep)
    // std::cerr<<"PE: "<<PE<<std::endl;
    // std::cerr<<"pe: "<<pe<<std::endl;
            // }

            // #pragma acc data ex 

}

int main()
{

	// MPI_Init(NULL, NULL);
    // int world_rank, world_size;

    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	
	lass clas;
	

	

    std::cerr<<"default device: "<<omp_get_default_device()<<std::endl;
	clas.init();

    std::cerr<<"num devices   : "<<omp_get_num_devices()<<std::endl;
	clas.tofu();






	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Finalize();

}

void lass::parse_input_file(char const* location)
{
    std::string s_location(location);
    std::string json_file = s_location + "input.json";
    std::ifstream ifs(json_file);
    // std::cerr<<json_file<<std::endl;
    //// CANNOT USE json::parse() IF YOU RDBUF TOO
    // std::cerr<<ifs.rdbuf()<<std::endl;
    json inputs = json::parse(ifs);

    // int world_rank = getRank();
    
    // if (world_rank == 0)
    // {
    //     if (inputs["seed"] == std::string("default"))
    //     {
    //         seed = static_cast<int>(time(nullptr));
    //     }
    //     else
    //     {
    //         seed = static_cast<int>(inputs["seed"]);
    //     }

    // }
    // MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);

    // srand(seed);
    // random_generator.seed(seed);
    // std::cout<<"SEED: "<<seed<<" for rank: "<<world_rank<<std::endl;

    // OMPthreads = inputs["OMPthreads"];
    // total_balls_to_add = inputs["N"];
    // std::string temp_typeSim = inputs["simType"];
    // if (temp_typeSim == "BPCA")
    // {
    //     typeSim = BPCA;
    // }
    // else
    // {
    //     typeSim = collider;
    // }
    // std::string temp_distribution = inputs["radiiDistribution"];
    // if (temp_distribution == "logNormal")
    // {
    //     radiiDistribution = logNorm;
    // }
    // else
    // {
    //     radiiDistribution = constant;
    // }
    // dynamicTime = inputs["dynamicTime"];
    // G = inputs["G"];
    // density = inputs["density"];
    u_s = inputs["u_s"];
    u_r = inputs["u_r"];
    // sigma = inputs["sigma"];
    // Y = inputs["Y"];
    // cor = inputs["cor"];
    // simTimeSeconds = inputs["simTimeSeconds"];
    // timeResolution = inputs["timeResolution"];
    // fourThirdsPiRho = 4. / 3. * pi * density;
    // scaleBalls = inputs["scaleBalls"];
    // maxOverlap = inputs["maxOverlap"];
    // KEfactor = inputs["KEfactor"];
    // if (inputs["v_custom"] == std::string("default"))
    // {
    //     v_custom = 0.36301555459799423;
    // }
    // else
    // {
    //     v_custom = inputs["v_custom"];
    // }
    // temp = inputs["temp"]; // this will modify v_custom in oneSizeSphere
    // double temp_kConst = inputs["kConsts"];
    // kConsts = temp_kConst * (fourThirdsPiRho / (maxOverlap * maxOverlap));
    // impactParameter = inputs["impactParameter"];
    // Ha = inputs["Ha"];
    // double temp_h_min = inputs["h_min"];
    // h_min = temp_h_min * scaleBalls;
    // if (inputs["cone"] == std::string("default"))
    // {
    //     cone = pi/2;
    // }
    // else
    // {
    //     cone = inputs["cone"];
    // }
    // properties = inputs["properties"];
    // genBalls = inputs["genBalls"];
    // attempts = inputs["attempts"];
    // skip = inputs["skip"];
    // steps = inputs["steps"];
    // dt = inputs["dt"];
    // kin = inputs["kin"];
    // kout = inputs["kout"];
    // if (inputs["spaceRange"] == std::string("default"))
    // {
    //     spaceRange = 4 * std::pow(
    //                     (1. / .74 * scaleBalls * scaleBalls * scaleBalls * genBalls),
    //                     1. / 3.); 
    // }
    // else
    // {
    //     spaceRange = inputs["spaceRange"];
    // }
    // if (inputs["spaceRangeIncrement"] == std::string("default"))
    // {
    //     spaceRangeIncrement = scaleBalls * 3;
    // }
    // else
    // {
    //     spaceRangeIncrement = inputs["spaceRangeIncrement"];
    // }
    // z0Rot = inputs["z0Rot"];
    // y0Rot = inputs["y0Rot"];
    // z1Rot = inputs["z1Rot"];
    // y1Rot = inputs["y1Rot"];
    // simTimeElapsed = inputs["simTimeElapsed"];
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
    // projectileName = inputs["projectileName"];
    // targetName = inputs["targetName"];
    // output_prefix = inputs["output_prefix"];
    // if (output_prefix == std::string("default"))
    // {
    //     output_prefix = "";
    // }

    // radiiFraction = inputs["radiiFraction"];

    // output_width = num_particles;
}
