
#include <omp.h>
#include <openacc.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include "vec3.hpp"

class lass
{
public:
    int num_particles = 400;
	int num_pairs = (num_particles*num_particles-num_particles)*0.5;

    vec3 *acc = new vec3[num_particles];
    vec3 *accaccum = new vec3[num_particles];
    vec3 *accsq = new vec3[num_particles*num_particles];
	double PE = 0.0;

	void init();
	void tofu();
    void loop_one_step();
};

void lass::init()
{
	for (int i = 0; i < num_particles; ++i)
	{
        // acc[i] = 0.0;
        acc[i] = {0.0,0.0,0.0};
		accaccum[i] = {0.0,0.0,0.0};
        for (int j = 0; j < num_particles; j++)
        {
            accsq[i*num_particles+j] = {0.0,0.0,0.0};
        }
	}

    PE = 0.0;
}


void lass::tofu()
{
	int outerLoop = 2000;

    	
    double t0 = omp_get_wtime();
  
    #pragma acc enter data copyin(this,accaccum[0:num_particles],acc[0:num_particles],accsq[0:num_particles*num_particles])    
	for (int k = 0; k < outerLoop; k++)
	{
        init();
        loop_one_step();
    }

	for (int i = 0; i < 10; i++)
	{
        std::cout<<acc[i]<<std::endl;
	}
    std::cout<<std::endl;
    for (int i = 0; i < 10; i++)
    {
		std::cout<<accaccum[i]<<std::endl;
    }
	
    #pragma acc exit data delete(acc[0:num_particles],accaccum[0:num_particles],accsq[0:num_particles*num_particles])

    delete[] acc;
	
	double t1 = omp_get_wtime();
    std::cerr<<"PE REDUCTION: "<<PE<<std::endl;
	std::cerr<<"Parallel part and data movement took "<<t1-t0<<" seconds"<<std::endl;

}


void lass::loop_one_step()
{
    
    double pe = 0.0;
    #pragma acc enter data copyin(pe)

    #pragma acc parallel loop num_gangs(108) num_workers(256) reduction(+:pe) present(this,pe,accaccum[0:num_particles],acc[0:num_particles],accsq[0:num_particles*num_particles])
    for (int pc = 1; pc <= num_pairs; pc++)
    {
        pe += 1;
        double pd = (double)pc;
        pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
        pd -= 0.00001;
        int A = (int)pd;
        int B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);

        
        // acc[A] += 1;
        // acc[B] += 1;
        // if (pc % 2 == 0)
        {
            #pragma acc atomic
            acc[A].x += 1;
            #pragma acc atomic
            acc[A].y += 1;
            #pragma acc atomic
            acc[A].z += 1;
            #pragma acc atomic
            acc[B].x -= 1;
            #pragma acc atomic
            acc[B].y -= 1;
            #pragma acc atomic
            acc[B].z -= 1;
        }

        accsq[A*num_particles+B].x = 1;
        accsq[A*num_particles+B].y = 1;
        accsq[A*num_particles+B].z = 1;
        accsq[B*num_particles+A].x = -1;
        accsq[B*num_particles+A].y = -1;
        accsq[B*num_particles+A].z = -1;

        // #pragma acc update host(pe) //Why does it work here

    }

    #pragma acc parallel loop gang worker present(this,acc[0:num_particles],accsq[0:num_particles*num_particles])
    for (int i = 0; i < num_particles; i++)
    {
        #pragma acc loop seq
        for (int j = 0; j < num_particles; j++)
        {
            accaccum[i].x += accsq[i*num_particles+j].x;
            accaccum[i].y += accsq[i*num_particles+j].y;
            accaccum[i].z += accsq[i*num_particles+j].z;
        }
    }
    
    #pragma acc update host(pe,accaccum[0:num_particles],acc[0:num_particles],accsq[0:num_particles*num_particles]) //but doesnt work here?

    PE = pe;

    // std::cout<<"PE: "<<PE<<std::endl;
    // std::cout<<"acc[0]: "<<acc[0].x<<','<<acc[0].y<<','<<acc[0].z<<std::endl;

    #pragma acc exit data delete(pe)

}

int main()
{
	lass clas;
	
    std::cerr<<"default device: "<<omp_get_default_device()<<std::endl;
	clas.init();

    std::cerr<<"num devices   : "<<omp_get_num_devices()<<std::endl;
	clas.tofu();
}


