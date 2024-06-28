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
#include "vec3.hpp"
// #include "class.hpp"




// #pragma omp declare target
class lass
{
public:
	int particles = 400;

	vec3 *acc = new vec3[particles];
	double *aacc = new double[particles];

	int num_pairs = (((particles*particles)-particles)/2);


	void init();
	void tofu();
};
// #pragma omp end declare target

void lass::init()
{
	for (int i = 0; i < particles; ++i)
	{
		acc[i] = {0.0,0.0,0.0};
		aacc[i] = 0.0;
	}

}


// #pragma omp declare target
void lass::tofu()
{
	std::cerr<<"Start main loop execution"<<std::endl;

	int outerLoop = 50;

	int pc;

	int A,B;
    double t0 = omp_get_wtime();
	for (int k = 0; k < outerLoop; k++)
	{
		// #pragma omp target map(tofrom:aacc[0:particles]) map(to:A,B,pc,num_pairs)
		#pragma omp target map(tofrom:acc[0:particles],aacc[0:particles]) map(to:A,B,pc,num_pairs)
	    {
        	// #pragma omp parallel for default(none) private(A,B,pc) shared(num_pairs,aacc)
        	#pragma omp parallel for default(none) private(A,B,pc) shared(num_pairs,acc,aacc)
		    for (pc = 1; pc <= num_pairs; pc++)
		    {
		    	
		    	//////////////////////////////////////
		    	double pd = (double)pc;
		    	pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
		    	pd -= 0.00001;
		    	// i = (long long)pd;
		    	i = (int)pd;
		    	j = (int)((double)pc-(double)i*((double)i-1.0)*.5-1.0);
		    	//////////////////////////////////////
				#pragma omp atomic
	                acc[A].x += 1;
	            #pragma omp atomic
	                acc[A].y += 1;
	            #pragma omp atomic
	                acc[A].z += 1;
	            #pragma omp atomic
	                acc[B].x -= 1;
	            #pragma omp atomic
	                acc[B].y -= 1;
	            #pragma omp atomic
	                acc[B].z -= 1;

	            // #pragma omp atomic
	            // 	acc[A] += {1,1,1};

	            aacc[A] += 1;


	            // distances[e] = dist;

	        }
		}


		for (int i = 0; i < particles; i++)
		{
			acc[i] = {0.0,0.0,0.0};
			aacc[i] = 0.0;
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
