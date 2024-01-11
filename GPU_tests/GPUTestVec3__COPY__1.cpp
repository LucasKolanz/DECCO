// #include <iostream>
// #include <fstream>
#include <omp.h>
// #include <mpi.h>
// #include <cmath>
// #include <random>
// #include <string>
// #include <iostream>
// #include <sstream>
// #include <cassert>
// #include <string>
// #include <random>
// #include <sstream>
// #include <iomanip>
// #include <cstdlib>      // To resolve std::abs ambiguity on clang
// #include <cstdint>      // For implementing namespace linalg::aliases
// #include <array>        // For std::array
// #include <iosfwd>       // For forward definitions of std::ostream
// #include <type_traits>  // For std::enable_if, std::is_same, std::declval
// #include <functional>   // For std::hash declaration
// #include "../SpaceLab/vec3.hpp"
// #include "include.hpp"
// #pragma omp declare target
#include "vec3.hpp"
// #pragma omp end declare target
// #include "class.hpp"


// double PE = 0.0;
// constexpr double Ha = 4.7e-12;
// double u_r = 1e-5;
// double u_s = 0.1;
// double kin = 3.7;
// double kout = 1.5;
// double h_min = 0.1;
// double dt = 1e-5;

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
	std::cerr<<"HERERREERREER"<<std::endl;
	int accum = 0;
	std::cerr<<"accum: "<<accum<<std::endl;
	// #pragma omp target
	// for (int l = 0; l < 10000; l++)
	// {
	// 	accum += l;
	// }

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
		
		// #pragma omp target defaultmap(none) map(to:vel[0:particles],pos[0:particles],m[0:particles],w[0:particles],Ha,world_rank,world_size,lllen,u_r,u_s,kin,kout,num_pairs,write_step,A,B,pc,R[0:particles],moi[0:particles],h_min,dt) map(tofrom:PE,acc[0:particles],aacc[0:particles],distances[0:num_pairs]) 
		#pragma omp target map(tofrom:acc[0:particles]) map(to:A,B,pc)
	    {
        	#pragma omp parallel for default(none) private(A,B,pc) shared(num_pairs,acc,write_step,world_rank,world_size)
		    for (pc = world_rank+1; pc <= num_pairs; pc+=world_size)
		    {
		    	// // std::cout<<omp_get_thread_num()<<std::endl;
		    	// // int threadNum = omp_get_thread_num();
		    	// double pd = (double)pc;
		    	// pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
		    	// pd -= 0.00001;
		    	// // i = (long long)pd;
		    	// A = (int)pd;
	            // B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);

	            

				#pragma omp atomic
	                acc[A][0] += 1;
	            #pragma omp atomic
	                acc[A][1] += 1;
	            #pragma omp atomic
	                acc[A][2] += 1;
	            #pragma omp atomic
	                acc[B][0] -= 1;
	            #pragma omp atomic
	                acc[B][1] -= 1;
	            #pragma omp atomic
	                acc[B][2] -= 1;


	            // distances[e] = dist;

	        }
		}


		for (int i = 0; i < particles; i++)
		{
			acc[i] = {0.0,0.0,0.0};
		}
		

	}
	double t1 = omp_get_wtime();
	std::cout<<"GPU, MPI, and OMP took "<<t1-t0<<" seconds"<<std::endl;
}
// #pragma omp end declare target
// #pragma omp end declare target

int main()
{

	// MPI_Init(NULL, NULL);
    // int world_rank, world_size;

    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	
	lass clas;
	

    std::cerr<<"default device: "<<omp_get_default_device()<<std::endl;
    std::cerr<<"num devices   : "<<omp_get_num_devices()<<std::endl;
	

	clas.init();

	clas.tofu();






	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Finalize();

}
