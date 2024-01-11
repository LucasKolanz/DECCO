#include <iostream>
#include <omp.h>

int main(){
    int particles = 2416;

    double *acc, *aacc, *vel, *pos, *m, *distances; 
    double *dacc, *daacc, *dvel, *dpos, *dm, *ddistances; 

    acc   = (double*)malloc(sizeof(double) * particles*3);
    aacc   = (double*)malloc(sizeof(double) * particles*3);
    vel = (double*)malloc(sizeof(double) * particles*3);
    pos = (double*)malloc(sizeof(double) * particles*3);
    m = (double*)malloc(sizeof(double) * particles);
    distances = (double*)malloc(sizeof(double) * int((particles*(particles-1))*0.5));

    cudaMalloc((void **) &dacc, sizeof(double) * particles*3);
    cudaMalloc((void **) &daacc, sizeof(double) * particles*3);
    cudaMalloc((void **) &dvel, sizeof(double) * particles*3);
    cudaMalloc((void **) &dpos, sizeof(double) * particles*3);
    cudaMalloc((void **) &dm, sizeof(double) * particles);
    cudaMalloc((void **) &ddistances, sizeof(double) * int((particles*(particles-1))*0.5));

    for (int i = 0; i < int((particles*(particles-1))*0.5); i++)
    {
        if (i < particles*3)
        {
            acc[i] = i*0.5;
            aacc[i] = i*0.5;
            vel[i] = i*0.5;
            pos[i] = i*0.5;
        }
        if (i < particles)
        {
            m[i] = i*1.5;
        }
        distances[i] = i*2.5;
    }


    double time0 = omp_get_wtime(); 
    cudaMemcpy(dacc,acc,sizeof(double)*particles*3,cudaMemcpyHostToDevice);
    cudaMemcpy(daacc,aacc,sizeof(double)*particles*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dvel,vel,sizeof(double)*particles*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dpos,pos,sizeof(double)*particles*3,cudaMemcpyHostToDevice);
    cudaMemcpy(dm,m,sizeof(double)*particles,cudaMemcpyHostToDevice);
    cudaMemcpy(ddistances,distances,sizeof(double)*int((particles*(particles-1))*0.5),cudaMemcpyHostToDevice);

    cudaMemcpy(acc,dacc,sizeof(double)*particles*3,cudaMemcpyDeviceToHost);
    cudaMemcpy(aacc,daacc,sizeof(double)*particles*3,cudaMemcpyDeviceToHost);

    double time1 = omp_get_wtime();

    cudaFree(dacc);
    cudaFree(daacc);
    cudaFree(dpos);
    cudaFree(dvel);
    cudaFree(dm);
    cudaFree(ddistances);

    delete [] acc;
    delete [] aacc;
    delete [] vel;
    delete [] pos;
    delete [] m;
    delete [] distances;


    std::cout<<"Equivalent copying took "<<(time1-time0)*994<<" seconds"<<std::endl;


}