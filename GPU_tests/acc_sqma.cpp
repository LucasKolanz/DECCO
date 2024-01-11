
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
	vec3 *velh = new vec3[num_particles];
    vec3 *w = new vec3[num_particles];
	vec3 *wh = new vec3[num_particles];
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
        velh[i] = {0.0,0.0,0.0};
		wh[i] = {0.0,0.0,0.0};
		pos[i] = {static_cast<double>(i),static_cast<double>((i+1)),static_cast<double>((i+2))};
		vel[i] = {0,0,0};
		w[i] = {0,0,0};
		R[i] = 1e-5;
		m[i] = 7.07e-10;
		moi[i] = 1;//(2/5)*m[i]*R[i]*R[i];
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
	int outerLoop = 10;

	// int pc;

	// int A,B;
	
	
	world_rank = 0;
	world_size = 1;

    

	
    double t0 = omp_get_wtime();
    // #pragma acc data reduction(+:PE) copyin(this) copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copy(acc[0:num_particles],aacc[0:num_particles]) copyin(A,B,pc,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    #pragma acc enter data copyin(this) 
    #pragma acc enter data copyin(moi[0:num_particles],m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) 
    #pragma acc enter data copyin(accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],acc[0:num_particles],aacc[0:num_particles],velh[0:num_particles],wh[0:num_particles]) 
    #pragma acc enter data copyin(dt,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,world_rank,world_size)
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
	    // init();
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
	
    #pragma acc exit data delete(accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],acc[0:num_particles],aacc[0:num_particles])
    #pragma acc exit data delete(m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],velh[0:num_particles],wh[0:num_particles],distances[0:num_pairs])
    #pragma acc exit data delete(dt,num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,world_rank,world_size)
    #pragma acc exit data delete(this)


    delete[] acc;
    delete[] aacc;
    delete[] accsq;
    delete[] aaccsq;
    delete[] vel;
    delete[] velh;
    delete[] pos;
    delete[] w;
    delete[] wh;
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
    
    /// FIRST PASS - Update Kinematic Parameters:
    #pragma acc parallel loop gang present(this,velh[0:num_particles],vel[0:num_particles],acc[0:num_particles],dt,wh[0:num_particles],w[0:num_particles],aacc[0:num_particles],pos[0:num_particles],num_particles)
    for (int Ball = 0; Ball < num_particles; Ball++) {
        // Update velocity half step:
        velh[Ball] = vel[Ball] + .5 * acc[Ball] * dt;

        // Update angular velocity half step:
        wh[Ball] = w[Ball] + .5 * aacc[Ball] * dt;

        // Update position:
        pos[Ball] += velh[Ball] * dt;

        // Reinitialize acceleration to be recalculated:
        acc[Ball] = {0, 0, 0};

        // Reinitialize angular acceleration to be recalculated:
        aacc[Ball] = {0, 0, 0};
    }

    double pe = 0.0;
    #pragma acc enter data copyin(pe)
    #pragma acc enter data copyin(writeStep)

    double t0 = omp_get_wtime();
    // std::cerr<<"IN simonestep"<<std::endl;

    // #pragma acc parallel loop gang worker private(pc,A,B) reduction(+:PE[0:1]) copyin(this) copy(PE[0:1],acc[0:num_particles],aacc[0:num_particles]) copyin(m[0:num_particles],moi[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs]) copyin(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size)
    #pragma acc parallel loop gang worker num_gangs(108) num_workers(256) reduction(+:pe) present(pe,this,accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles],m[0:num_particles],moi[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs],num_pairs,num_particles,Ha,k_in,k_out,h_min,u_s,u_r,writeStep,world_rank,world_size)
    for (int pc = 1; pc <= num_pairs; pc ++)
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



            // aaccsq[0].x = torqueA.x;
            // aaccsq[0].y = torqueA.y;
            // aaccsq[0].z = torqueA.z;
            // if (std::isnan(aaccA.x))
            // {
            //     aaccA.x = 10;
            // }
            // if (std::isnan(aaccA.y))
            // {
            //     aaccA.y = 10;
            // }
            // if (std::isnan(aaccA.z))
            // {
            //     aaccA.z = 10;
            // }
            // if (std::isnan(aaccB.x))
            // {
            //     aaccB.x = 10;
            // }
            // if (std::isnan(aaccB.y))
            // {
            //     aaccB.y = 10;
            // }
            // if (std::isnan(aaccB.z))
            // {
            //     aaccB.z = 10;
            // }
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
                // #pragma omp critical
                // #pragma acc atomic
                // PE += 1; ///TURN ON FOR REAL SIM
                // std::cerr<<"PE: "<<PE<<std::endl;
                // printf("HERE");
                // PE[0] += U_vdw + 0.5 * k * overlap * overlap; ///TURN ON FOR REAL SIM
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
                // #pragma omp critical
                // printf("HERE1");
                // PE += 1.0;  // Van Der Waals TURN ON FOR REAL SIM
                pe += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
                // PE[0] += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
            }

            // todo this is part of push_apart. Not great like this.
            // For pushing apart overlappers:
            // vel[A] = { 0,0,0 };
            // vel[B] = { 0,0,0 };
        }


        vec3 accA = (1/m[A])*totalForceOnA; 
        vec3 accB = (1/m[B])*totalForceOnA; 

        // accsq[0].x = totalForceOnA.x;
        // accsq[0].y = totalForceOnA.y;
        // accsq[0].z = totalForceOnA.z;
        accsq[A*num_particles+B].x = accA.x;
        accsq[A*num_particles+B].y = accA.y;
        accsq[A*num_particles+B].z = accA.z;
        accsq[B*num_particles+A].x = accA.x;
        accsq[B*num_particles+A].y = accA.y;
        accsq[B*num_particles+A].z = accA.z;

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

    #pragma acc update host(pe,acc[0:num_particles],aacc[0:num_particles],accsq[0:num_particles*num_particles],aaccsq[0:num_particles*num_particles])
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
            if (world_rank == 0)
            {
                std::cerr<<"PE: "<<PE<<std::endl;
            }
        }
    #endif

    double t1 = omp_get_wtime();


    
    // double ke = 0.0;
    // vec3 Mom,Ang_mom;
    // Mom = {0.0,0.0,0.0};
    // Ang_mom = {0.0,0.0,0.0};
    // #pragma acc enter data copyin(ke,Mom,Ang_mom)
    #pragma acc parallel loop num_gangs(108) num_workers(256) present(this,acc[0:num_particles],aacc[0:num_particles],w[0:num_particles],vel[0:num_particles],velh[0:num_particles],wh[0:num_particles],num_particles,dt)
    for (int Ball = 0; Ball < num_particles; Ball++) {
        // Velocity for next step:
        vel[Ball] = velh[Ball] + .5 * acc[Ball] * dt;
        w[Ball] = wh[Ball] + .5 * aacc[Ball] * dt;

        /////////////////////////////////
        // if (true) {
        /////////////////////////////////
        
        // #pragma acc update host(ke,Mom,Ang_mom)
    }  // THIRD PASS END


    // if (write_step && world_rank == 0) {
    //     // std::cerr<<"Writing "<<num_particles<<" balls"<<std::endl;
        
    // }
    #pragma acc update host(w[0:num_particles],vel[0:num_particles],pos[0:num_particles]) if(writeStep && world_rank == 0)
    // if (writeStep && world_rank == 0) 
    // {
    //     ballBuffer << '\n';  // Prepares a new line for incoming data.
    //     for (int Ball = 0; Ball < num_particles;  Ball++)
    //     {
    //         if (Ball == 0) {
    //             ballBuffer << pos[Ball][0] << ',' << pos[Ball][1] << ',' << pos[Ball][2] << ','
    //                        << w[Ball][0] << ',' << w[Ball][1] << ',' << w[Ball][2] << ','
    //                        << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y << ','
    //                        << vel[Ball].z << ',' << 0;
    //         } else {
    //             ballBuffer << ',' << pos[Ball][0] << ',' << pos[Ball][1] << ',' << pos[Ball][2] << ','
    //                        << w[Ball][0] << ',' << w[Ball][1] << ',' << w[Ball][2] << ','
    //                        << w[Ball].norm() << ',' << vel[Ball].x << ',' << vel[Ball].y << ','
    //                        << vel[Ball].z << ',' << 0;
    //         }
    //         // Send positions and rotations to buffer:
    //         // std::cerr<<"Write ball "<<Ball<<std::endl;
    //         KE += .5 * m[Ball] * vel[Ball].normsquared() +
    //                 .5 * moi[Ball] * w[Ball].normsquared();  // Now includes rotational kinetic energy.
    //         mom += m[Ball] * vel[Ball];
    //         ang_mom += m[Ball] * pos[Ball].cross(vel[Ball]) + moi[Ball] * w[Ball];
    //     }
    // }

}

int main()
{

	
	lass clas;
	

    std::cerr<<"default device: "<<omp_get_default_device()<<std::endl;
	clas.init();

    std::cerr<<"num devices   : "<<omp_get_num_devices()<<std::endl;
	clas.tofu();

}


