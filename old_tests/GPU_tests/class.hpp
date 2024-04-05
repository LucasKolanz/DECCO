
#pragma omp declare target
#include "vec3.hpp"
#pragma omp end declare target


double PE = 0.0;
constexpr double Ha = 4.7e-12;
double u_r = 1e-5;
double u_s = 0.1;
double kin = 3.7;
double kout = 1.5;
double h_min = 0.1;
double dt = 1e-5;
#pragma omp declare target
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
	int lllen = particles;//static_cast<long long>(particles);
	int num_pairs = (((lllen*lllen)-lllen)/2);
	double *distances = new double[num_pairs];
	std::string writeFileName = "timing.csv";


	void init();
	void tofu();
};
#pragma omp end declare target

#pragma omp declare target
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
}
#pragma omp end declare target


// #pragma omp declare target
void lass::tofu()
{
	int accum = 0;
	std::cerr<<"accum: "<<accum<<std::endl;
	#pragma omp target
	for (int l = 0; l < 10000; l++)
	{
		accum += l;
	}

	std::cerr<<"accum: "<<accum<<std::endl;

	int outerLoop = 50;

	int pc;
	bool write_step;
	int world_rank = 0;
	int world_size = 1;
	int A,B;
    double t0 = omp_get_wtime();
	for (int k = 0; k < outerLoop; k++)
	{

		if (k%5==0)
		{
			write_step = true;
		}
		else
		{
			write_step = false;
		}
		
		#pragma omp target //defaultmap(none) map(to:vel[0:particles],pos[0:particles],m[0:particles],w[0:particles],Ha,world_rank,world_size,lllen,u_r,u_s,kin,kout,num_pairs,write_step,A,B,R[0:particles],moi[0:particles],h_min,dt) map(tofrom:PE,acc[0:particles],aacc[0:particles],distances[0:num_pairs]) 
	    {
        	#pragma omp parallel for //reduction(+:PE) default(none) private(A,B,pc) shared(k,num_pairs,acc,aacc,world_rank,world_size,Ha,write_step,lllen,R,pos,vel,m,w,u_r,u_s,moi,kin,kout,distances,h_min,dt)
		    for (pc = world_rank+1; pc <= num_pairs; pc+=world_size)
		    {
		    	// std::cout<<omp_get_thread_num()<<std::endl;
		    	// int threadNum = omp_get_thread_num();
		    	double pd = (double)pc;
		    	pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
		    	pd -= 0.00001;
		    	// i = (long long)pd;
		    	A = (int)pd;
	            B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);

	            const double sumRaRb = R[A] + R[B];
	            const vec3 rVecab = pos[B] - pos[A];  // Vector from a to b.
	            const vec3 rVecba = -rVecab;
	            const double dist = (rVecab).norm();




	            // Check for collision between Ball and otherBall:
	            double overlap = sumRaRb - dist;

	            vec3 totalForceOnA{0, 0, 0};

	            // Distance array element: 1,0    2,0    2,1    3,0    3,1    3,2 ...
	            // int e = static_cast<unsigned>(A * (A - 1) * .5) + B;  // a^2-a is always even, so this works.
	            int e = pc-1;
	            double oldDist = distances[e];

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

	                for (int i = 0; i < 3; i ++)
	                {
	                    #pragma omp atomic
	                        aacc[A][i] += aaccA[i];
	                    #pragma omp atomic
	                        aacc[B][i] += aaccB[i];

	                }



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
	                    PE += U_vdw;  // Van Der Waals TURN ON FOR REAL SIM
	                }

	                // todo this is part of push_apart. Not great like this.
	                // For pushing apart overlappers:
	                // vel[A] = { 0,0,0 };
	                // vel[B] = { 0,0,0 };
	            }


	            vec3 accA = (1/m[A])*totalForceOnA; 
	            vec3 accB = (1/m[B])*totalForceOnA; 
	            #pragma omp atomic
	                acc[A][0] += accA[0];
	            #pragma omp atomic
	                acc[A][1] += accA[1];
	            #pragma omp atomic
	                acc[A][2] += accA[2];
	            #pragma omp atomic
	                acc[B][0] -= accB[0];
	            #pragma omp atomic
	                acc[B][1] -= accB[1];
	            #pragma omp atomic
	                acc[B][2] -= accB[2];


	            distances[e] = dist;

	        }
		}


		for (int i = 0; i < particles; i++)
		{
			acc[i] = {0.0,0.0,0.0};
			aacc[i] = {0.0,0.0,0.0};
		}
		

	}
	double t1 = omp_get_wtime();
	std::cout<<"GPU, MPI, and OMP took "<<t1-t0<<" seconds"<<std::endl;
}
// #pragma omp end declare target
// #pragma omp end declare target