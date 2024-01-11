
#include <omp.h>
#include <openacc.h>
#include <fstream>
#include "vec3.hpp"
// #include "../SpaceLab/external/json/single_include/nlohmann/json.hpp"
// #include "Utils.hpp"
// #include "linalg.hpp"

// double PE = 0.0;
const double Ha = 4.7e-12;
double u_r = 1e-5;
double u_s = 0.1;
double kin = 3.7;
double kout = 1.5;
double h_min = 0.1;
double dt = 1e-5;

// using json = nlohmann::json;


class lass
{
public:
	int num_particles = 400;
	double *R = new double[num_particles];
	double *m = new double[num_particles];
	double *moi = new double[num_particles];
    vec3 *acc = new vec3[num_particles];
	vec3 *aacc = new vec3[num_particles];
	vec3 *pos = new vec3[num_particles];
	vec3 *vel = new vec3[num_particles];
	vec3 *w = new vec3[num_particles];
	int num_pairs = (((num_particles*num_particles)-num_particles)/2);
	double *distances = new double[num_pairs];

	bool write_step;

    vec3 *accsq = new vec3[num_particles*num_particles];
	vec3 *aaccsq = new vec3[num_particles*num_particles];

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
	for (int i = 0; i < num_particles; ++i)
	{
		acc[i] = {0.0,0.0,0.0};
		aacc[i] = {0.0,0.0,0.0};
		pos[i] = {static_cast<double>(i/num_particles),static_cast<double>((i+1)/num_particles),static_cast<double>((i+2)/num_particles)};
		vel[i] = {static_cast<double>(i/num_particles),static_cast<double>((i-1)/num_particles),static_cast<double>((i-2)/num_particles)};
		w[i] = {static_cast<double>(i/num_particles),static_cast<double>((i*2)/num_particles),static_cast<double>((i*3)/num_particles)};
		R[i] = 1e-5;
		m[i] = 7.07e-10;
		moi[i] = (2/5)*m[i]*R[i]*R[i];
        for (int j = 0; j < num_particles; j++)
        {
            accsq[i*num_particles+j] = {0.0,0.0,0.0};
            aaccsq[i*num_particles+j] = {0.0,0.0,0.0};
        }
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
	int outerLoop = 1000;

	// int pc;

	// int A,B;
	
	
	world_rank = 0;
	world_size = 1;

    

	
    double t0 = omp_get_wtime();
    // #pragma acc data reduction(+:PE) copyin(this) copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copy(acc[0:num_particles],aacc[0:num_particles]) copyin(A,B,pc,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    #pragma acc enter data copyin(this) 
    #pragma acc enter data copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) 
    #pragma acc enter data copyin(acc[0:num_particles],aacc[0:num_particles],accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles]) 
    #pragma acc enter data copyin(num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,world_rank,world_size)
    // #pragma acc enter data create(A,B,pc,world_rank,world_size)
    // {                                                   
	for (int k = 0; k < outerLoop; k++)
	{
        // printf("HREERREERERE");
        if(k % 5 == 0)
        {
            write_step = false;
        }
        else
        {
            write_step = true;
        }

        // #pragma update device(write_step)
	    init();
        loop_one_step(write_step);
        // #pragma acc update host(pe)
        // PE = pe;
        // std::cerr<<"PE: "<<PE<<std::endl;
        // std::cerr<<"pe: "<<pe<<std::endl;
    }

    // }


    // #pragma acc update host(acc[0:num_particles],aacc[0:num_particles])
	for (int i = 0; i < 1; i++)
	{
		// acc[i] = {0.0,0.0,0.0};
		std::cout<<acc[i]<<std::endl;
		std::cout<<aacc[i]<<std::endl;
	}
	
    #pragma acc exit data delete(acc[0:num_particles],aacc[0:num_particles],accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles])
	#pragma acc exit data delete(m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs])
	#pragma acc exit data delete(num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,world_rank,world_size)
	#pragma acc exit data delete(this)
	// delete[] acc;
	// free(dacc);

    delete[] acc;
    delete[] aacc;
    delete[] accsq;
    delete[] aaccsq;
    delete[] vel;
    delete[] pos;
    delete[] w;
    delete[] m;
    delete[] moi;
    delete[] R;
    delete[] distances;
	
	
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

    #pragma acc enter data copyin(pe)
    #pragma acc enter data copyin(writeStep)
    
    // #pragma acc parallel loop num_gangs(108) num_workers(256) reduction(+:pe) present(pe,this,accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],num_particles,writeStep)
    #pragma acc parallel loop gang worker num_gangs(108) num_workers(256) reduction(+:pe) present(pe,this,accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],m[0:num_particles],moi[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs],num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,writeStep,world_rank,world_size)
    for (int pc = 1; pc <= num_pairs; pc ++)
    {
        if (writeStep)
        {
            pe += 1;
        }
        
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



            aaccsq[A*num_particles+B].x = 1;
            aaccsq[A*num_particles+B].y = 1;
            aaccsq[A*num_particles+B].z = 1;
            aaccsq[B*num_particles+A].x = -1;
            aaccsq[B*num_particles+A].y = -1;
            aaccsq[B*num_particles+A].z = -1;


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

        accsq[A*num_particles+B].x = 1;
        accsq[A*num_particles+B].y = 1;
        accsq[A*num_particles+B].z = 1;
        accsq[B*num_particles+A].x = -1;
        accsq[B*num_particles+A].y = -1;
        accsq[B*num_particles+A].z = -1;

        // distances[e] = dist;


        // #pragma acc update host(pe)

    }

    #pragma acc parallel loop gang worker present(this,acc[0:num_particles],aacc[0:num_particles],accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles])
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

    PE = pe;

    #pragma acc exit data delete(pe,writeStep)

}

int main()
{

	
	lass clas;
	

    std::cerr<<"default device: "<<omp_get_default_device()<<std::endl;
	clas.init();

    std::cerr<<"num devices   : "<<omp_get_num_devices()<<std::endl;
	clas.tofu();

}


