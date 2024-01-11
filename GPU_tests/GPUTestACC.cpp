#include <iostream>
#include <fstream>
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include "../SpaceLab/vec3.hpp"

#ifndef PSEUDOPARTICLES
	#define PSEUDOPARTICLES 1000
#endif

#ifndef THREADS
	#define THREADS 4
#endif


// bool arrComp(double *x,double *y,double *z,vec3 *a2)
// {
// 	bool enter = false;
// 	for (int i = 0; i < PSEUDOPARTICLES; i++)
// 	{
// 		// std::cerr<<i<<std::endl;
// 		// std::cerr<<a2[0][0]<<std::endl;
// 		// std::cerr<<a1[0][0]<<std::endl;
// 		enter = true;
// 		if (x[i] != a2[i][0] || y[i] != a2[i][1] || z[i] != a2[i][2])
// 		{
// 			return false;
// 		}
// 	}
// 	return enter;
// }

bool arrComp(vec3 *a1,vec3 *a2)
{
	bool enter = false;
	for (int i = 0; i < PSEUDOPARTICLES; i++)
	{
		// std::cerr<<i<<std::endl;
		// std::cerr<<a2[0][0]<<std::endl;
		// std::cerr<<a1[0][0]<<std::endl;
		enter = true;
		if (a1[i].x != a2[i].x || a1[i].y != a2[i].y || a1[i].z != a2[i].z)
		{
			return false;
		}
	}
	return enter;
}



int main()
{
	// MPI_Init(NULL, NULL);
    int world_rank, world_size;
    world_rank = 0;
    world_size = 1;
    // double accum1 = 0;
    // double accum2 = 0;
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	int outerLoop = 10;

	vec3 *arr_comp = new vec3[PSEUDOPARTICLES];
	vec3* arr = new vec3[PSEUDOPARTICLES];
	// double *arrx = new double[PSEUDOPARTICLES];
	// double *arry = new double[PSEUDOPARTICLES];
	// double *arrz = new double[PSEUDOPARTICLES];
	// vec3 *arr2 = new vec3[PSEUDOPARTICLES];
	// vec3 *arr3 = new vec3[PSEUDOPARTICLES];
	int lllen = PSEUDOPARTICLES;//static_cast<long long>(PSEUDOPARTICLES);
	int num_pairs = (((lllen*lllen)-lllen)/2);
	std::string writeFileName = "timing.csv";
	
	// MPI Initialization
    
    // int prob_size = ceil(num_pairs/world_size);
    // int local_pairs = prob_size;
    // if (world_rank == world_size-1)
    // {
    // 	local_pairs = (world_rank*prob_size) - num_pairs;
    // }
    // vec3 *local_arr = new vec3[prob_size]
	

    if (world_rank == 0)
    {
		std::cout<<"PSEUDOPARTICLES:\t"<<PSEUDOPARTICLES<<std::endl;
		std::cout<<"THREADS:\t\t"<<THREADS<<std::endl;
		std::cout<<"NODES:\t\t\t"<<world_size<<std::endl;
		std::cout<<"total iterations:\t"<<(((lllen*lllen)-lllen)/2)<<std::endl;
	}
	for (int i = 0; i < PSEUDOPARTICLES; ++i)
	{
		arr[i] = {0,0,0};
		arr_comp[i] = {0,0,0};
	}


	// MPI_Bcast(arr1,PSEUDOPARTICLES*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	// MPI_Bcast(arr2,PSEUDOPARTICLES*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	// MPI_Scatter(&arr1,PSEUDOPARTICLES*3,MPI_DOUBLE,&local_arr,prob_size*3,MPI_DOUBLE,0,MPI_COMM_WORLD);

    //Start regular openmp 
    double PE = 0.0;
    double t0 = omp_get_wtime();
    int A;
    int B;
	int pc;
	// int start = (world_rank*prob_size)+1;
	// int stop =  prob_size+start-1;
	// if (world_rank == world_size-1)
	// {
	// 	// std::cout<<"HERE: "<<prob_size<<" "<<world_rank<<" "<<std::endl;
	// 	stop = num_pairs;
	// }
	// std::cout<<"rank: "<<world_rank<<"\t\tstart: "<<start<<"\tstop: "<<stop<<std::endl;
	// omp_set_num_threads(THREADS);
	// vec3* reduceMe = new vec3[THREADS*PSEUDOPARTICLES];
	// for (int l = 0; l < THREADS; l++)
	// {
	// 	for (int m = 0; m < PSEUDOPARTICLES; m++)
	// 	{
	// 		reduceMe[l*PSEUDOPARTICLES+m] = {0,0,0};
	// 		// reduceMe[l*PSEUDOPARTICLES+m].x = 0;
	// 		// reduceMe[l*PSEUDOPARTICLES+m].y = 0;
	// 		// reduceMe[l*PSEUDOPARTICLES+m].z = 0;
	// 	}
	// }

	for (int k = 0; k < outerLoop; k++)
	{
		
		// for(int i = 0; i < PSEUDOPARTICLES; ++i) {
		// 	arr1[i] = {0,0,0};
		// }
		// #pragma omp declare reduction(vec3_sum : vec3 : omp_out += omp_in)
		// #pragma omp parallel for default(none) reduction(+:accum1) reduction(vec3_sum:arr1[:PSEUDOPARTICLES]) num_threads(THREADS) shared(num_pairs,world_size,std::cerr,world_rank) private(prob_size,i,j,pc)
		    // for (pc = world_rank+1; pc <= num_pairs; pc+=world_size)
	    // #pragma omp parallel for default(none) reduction(+:arr1[:PSEUDOPARTICLES]) shared(world_rank,world_size,num_pairs) private(i,j,pc)
	    // #pragma omp target map(tofrom:arr1[0:PSEUDOPARTICLES])
	    // std::cout<<"BEFORE TARGET"<<std::endl;
		// #pragma acc data pcopy(this)
	    // #pragma acc data pcopyin(num_pairs,num_particles,A,B,pc,Ha,k_in,k_out,h_min,u_s,u_r,write_step,world_rank,world_size) 
	    // #pragma acc data pcopyin(m[0:num_particles],w[0:num_particles],vel[0:num_particles],pos[0:num_particles],R[0:num_particles],distances[0:num_pairs])
	    #pragma acc data pcopy(arr[0:PSEUDOPARTICLES],PE)
	    {

			// PE += 1;
			#pragma acc parallel loop gang worker private(pc,A,B) reduction(+:PE) 
	    	for (pc = 1; pc <= num_pairs; pc++)
		    {
		    	// std::cout<<omp_get_thread_num()<<std::endl;
		    	// int threadNum = omp_get_thread_num();
		    	double pd = (double)pc;
		    	pd = (sqrt(pd*8.0+1.0)+1.0)*0.5;
		    	pd -= 0.00001;
		    	// i = (long long)pd;
		    	A = (int)pd;
		    	B = (int)((double)pc-(double)A*((double)A-1.0)*.5-1.0);
		    	// j = (long long)((long double)pc-(long double)i*((long double)i-1.0L)*.5L-1.0L);

		    	// if (world_rank == 0)
		    	// {
		    	// 	std::cerr<<"i: "<<i<<"\tj: "<<j<<std::endl;
		    	// }
		    	// accum1 += pc*1.0;
		        // arr1[i] += {1,1,1};
		        // arr1[j] += {1,1,1};

		    	// if (pc == 1)
		    	// {
			    // 	// int tid = omp_get_default_device();
			    // 	int nthreads = omp_get_num_devices();
			    // 	// printf("def device %d out of %d devices\n",tid,nthreads);
			    // 	printf("%d devices\n",nthreads);
			    	// printf("total threads: %d\n",omp_get_num_threads());
		    	// }
		    	// getchar();
 				PE += 1.0;
 				#pragma acc atomic
		        arr[A].x += 1.0;
 				#pragma acc atomic
		        arr[B].x += 1.0;

 				#pragma acc atomic
		        arr[A].y += 1.0;
 				#pragma acc atomic
		        arr[B].y += 1.0;

 				#pragma acc atomic
		        arr[A].z += 1.0;
 				#pragma acc atomic
		        arr[B].z += 1.0;

		        // reduceMe[threadNum*PSEUDOPARTICLES+j] += {1,1,1};
		        // reduceMe[threadNum*PSEUDOPARTICLES+i] += {1,1,1};
		        // reduceMe[threadNum*PSEUDOPARTICLES+j].x += 1;
		        // reduceMe[threadNum*PSEUDOPARTICLES+j].y += 1;
		        // reduceMe[threadNum*PSEUDOPARTICLES+j].z += 1;
		        // reduceMe[threadNum*PSEUDOPARTICLES+i].x += 1;
		        // reduceMe[threadNum*PSEUDOPARTICLES+i].y += 1;
		        // reduceMe[threadNum*PSEUDOPARTICLES+i].z += 1;
		    	// std::cerr<<arr1[i]<<std::endl;
		    }
		}
	    // getchar();
		// MPI_Allreduce(arr1,arr2,PSEUDOPARTICLES*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
		
		// #pragma omp parallel for //shared(arr1) default(none)
		// for (int a = 0; a < PSEUDOPARTICLES; a++)
		// {
		// 	vec3 sum = {0.0,0.0,0.0};
		// 	// vec3 sum;
		// 	// sum.x = 0;
		// 	// sum.y = 0;
		// 	// sum.z = 0;

		// 	for (int b = 0; b < THREADS; b++)
		// 	{
		// 		sum += reduceMe[b*PSEUDOPARTICLES+a];
		// 		reduceMe[b*PSEUDOPARTICLES+a] = {0.0,0.0,0.0};
		// 		// sum.x += reduceMe[b*PSEUDOPARTICLES+a].x;
		// 		// sum.y += reduceMe[b*PSEUDOPARTICLES+a].y;
		// 		// sum.z += reduceMe[b*PSEUDOPARTICLES+a].z;
		// 		// reduceMe[b*PSEUDOPARTICLES+a].x = 0.0;
		// 		// reduceMe[b*PSEUDOPARTICLES+a].y = 0.0;
		// 		// reduceMe[b*PSEUDOPARTICLES+a].z = 0.0;
		// 	}
		// 	arr1[a] += sum;
		// 	// arr1[a].x += sum.x;
		// 	// arr1[a].y += sum.y;
		// 	// arr1[a].z += sum.z;
		// 	// arr1[i] = {0,0,0};
		// 	// arr3[i] += arr2[i];
		// 	// arr2[i] = {0,0,0};
		// }
		
	    // MPI_Barrier(MPI_COMM_WORLD);
	}

	// for (int i = 0; i < PSEUDOPARTICLES; i++)
	// {
	// 	arr1[i] = arr3[i];
	// }

	// std::cout<<"world_size: "<<world_size<<std::endl;

	// if (world_rank == 0)
	// {
	// 	for (int i = 1; i < world_size; i++)
	// 	{
	// 		MPI_Recv(&arr2,PSEUDOPARTICLES*3,MPI_DOUBLE,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	// 		std::cout<<"Recieved node "<<i<<std::endl;
	// 		if (arrComp(arr1,arr2))
	// 		{
	// 			std::cout<<"Same arrays for rank 0 and "<<i<<std::endl;
	// 		}
	// 		else
	// 		{
	// 			std::cout<<"DIFFERENT arrays for rank 0 and "<<i<<std::endl;
	// 		}
	// 	}
	// }
	// else
	// {
	// 	MPI_Send(&arr1,PSEUDOPARTICLES*3,MPI_DOUBLE,0,world_rank,MPI_COMM_WORLD);
	// }

	
	// std::cerr<<"MPI_Finalize "<<MPI_Finalize()<<std::endl;

	// std::cout<<"HERE1"<<std::endl;
	double PE1 = 0.0;
	if (world_rank == 0)
	{
		//Start regular single core run
		double t1 = omp_get_wtime();
		for (int k = 0; k < outerLoop; k++)
		{
			for (int c = 1; c < PSEUDOPARTICLES; c++)
			{
				for (int d = 0; d < c; ++d)
				{
					PE1 += 1;
					arr_comp[d] += {1,1,1};
		        	arr_comp[c] += {1,1,1};
		        	// arr_comp[d].x += 1;
		        	// arr_comp[d].y += 1;
		        	// arr_comp[d].z += 1;
		        	// arr_comp[c].x += 1;
		        	// arr_comp[c].y += 1;
		        	// arr_comp[c].z += 1;
				}
			}
		}	
		double t2 = omp_get_wtime();

		bool same1 = arrComp(arr,arr_comp);
		
		std::cout<<"same = "<<same1<<std::endl;
		std::cout<<"GPU, MPI, and OMP took "<<t1-t0<<" seconds"<<std::endl;
		std::cout<<"regular adding up took "<<t2-t1<<" seconds"<<std::endl;

		if (same1)
		{
			// std::fstream write;
			// write.open(writeFileName,std::ios::app);
			// // write<<"threads"<<','<<"particles"<<','<<"PSe time"<<','<<"regular omp time"<<','<<"Man red time (NO reduce)"<<','<<"Man red time (reduce)"<<','<<"reg time"<<std::endl;
			// write<<THREADS<<','<<PSEUDOPARTICLES<<','<<t1-t0<<','<<t2-t1<<','<<manRedTotTime<<','<<t3-t2<<','<<t4-t3<<std::endl;
			std::cerr<<"Arrays are the same"<<std::endl;
		}
		else
		{
			std::cerr<<"ARRAYS NOT THE SAME"<<std::endl;
		}
		int count = 0;
		for (int e = 0; e < PSEUDOPARTICLES; ++e)
		{
			// if (arr1[e].x != arr_comp[e].x || arr1[e].y != arr_comp[e].y || arr1[e].z != arr_comp[e].z)
			if (arr[e][0] != arr_comp[e][0] || arr[e][1] != arr_comp[e][1] || arr[e][2] != arr_comp[e][2])
			{
				count++;
				std::cout<<"element "<<e<<", array: "<<arr[e]<<", arrcomp: "<<arr_comp[e]<<std::endl;
				if (count > 100)
				{
					break;
				}
			}
		}
		std::cout<<count<<" elements are different."<<std::endl;
	}
	std::cout<<"PE, PE1: "<<PE<<", "<<PE1<<std::endl;
	// MPI_Barrier(MPI_COMM_WORLD);
	// MPI_Finalize();

}
