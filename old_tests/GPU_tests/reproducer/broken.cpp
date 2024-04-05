
#include <omp.h>

#include "vec3.hpp"




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

void lass::init()
{
	for (int i = 0; i < particles; ++i)
	{
		acc[i] = {0.0,0.0,0.0};
		aacc[i] = 0.0;
	}

}


void lass::tofu()
{
	int outerLoop = 50;

	int pc;

	int A,B;
    double t0 = omp_get_wtime();
	for (int k = 0; k < outerLoop; k++)
	{
		#pragma omp target map(tofrom:acc[0:particles],aacc[0:particles],this) map(to:A,B,pc,num_pairs)
	    {
        	#pragma omp parallel for default(none) private(A,B,pc) shared(num_pairs,acc,aacc)
		    for (pc = 1; pc <= num_pairs; pc++)
		    {
		    	
		    	//////////////////////////////////////
		    	double pd = (double)pc;
		    	pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
		    	pd -= 0.00001;
		    	
		    	A = (int)pd;
		    	B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);
		    	//////////////////////////////////////
		    	
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

	        }
		}


		for (int i = 0; i < particles; i++)
		{
			acc[i] = {0.0,0.0,0.0};
		}
		

	}
	double t1 = omp_get_wtime();
	std::cout<<"GPU took "<<t1-t0<<" seconds"<<std::endl;
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
