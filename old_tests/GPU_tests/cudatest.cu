#include <iostream>
#include <omp.h>
#include <assert.h>

inline
cudaError_t checkCuda(cudaError_t result)
{
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
  return result;
}

void noCuda()
{
    int size = 10;

    double* matrix = new double[size*size];
    double* vector = new double[size];
    double* output = new double[size];

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i*size+j] = 1.0;
        }
        vector[i] = 2.0;
        output[i] = 0.0;
    }    

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            output[i] += matrix[i*size+j]*vector[i];
        }
    }

    for (int i = 0; i < size; i++)
    {
        std::cout<<output[i]<<std::endl;
    }
}


__global__ void naivMVM(double* matrix,double* vector,double* output,int len)
{
    int row = blockIdx.y*blockDim.y+threadIdx.y;
    // int col = blockIdx.x*blockDim.x+threadIdx.x;
    double sum = 0;
    for (int i = 0; i < len; i++)
    {
        sum += matrix[row*len+i]*vector[i];
    }
    output[row] = sum;
}

void naivCuda()
{
    int block_size = 32;
    int size = block_size*31;

    double *h_output,*h_matrix,*h_vector;
    h_matrix = (double*)malloc(sizeof(double)*size*size);
    h_vector = (double*)malloc(sizeof(double)*size);
    h_output = (double*)malloc(sizeof(double)*size);

    double *d_output,*d_matrix,*d_vector;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            h_matrix[i*size+j] = 1.0;
        }
        h_vector[i] = 2.0;
        h_output[i] = 0.0;
    }    


    dim3 dimGrid(size/block_size, size/block_size, 1);
    dim3 dimBlock(block_size, block_size, 1);

    checkCuda(cudaMalloc(&d_matrix,sizeof(double)*size*size));
    checkCuda(cudaMalloc(&d_vector,sizeof(double)*size));
    checkCuda(cudaMalloc(&d_output,sizeof(double)*size));
    cudaDeviceSynchronize();
    
    checkCuda(cudaMemcpy(d_matrix, h_matrix, sizeof(double)*size*size, cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(d_vector, h_vector, sizeof(double)*size, cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();

    naivMVM<<<dimGrid,dimBlock>>>(d_matrix,d_vector,d_output,size);

    cudaDeviceSynchronize();
    checkCuda(cudaMemcpy(h_output, d_output, sizeof(double)*size, cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();
    

    for (int i = 0; i < size; i++)
    {
        std::cout<<h_output[i]<<',';
    }
    std::cout<<std::endl;

    delete[] h_matrix;
    delete[] h_vector;
    delete[] h_output;

    cudaFree(&d_output);
    cudaFree(&d_matrix);
    cudaFree(&d_vector);
}

#define TILE

__global__ void goodMVM(double* matrix,double* vector,double* output,int len)
{
    const int tileDim = 64;
    const int blockRows = 8;
    __shared__ double mat_tile[tileDim][tileDim];
    __shared__ double vec_tile[blockRows];

    int row = blockIdx.y*tileDim+threadIdx.y;
    int col = blockIdx.x*tileDim+threadIdx.x;
    int width = gridDim.x*tileDim;

        // vec_tile[0] = 3;
    for (int i = 0; i < tileDim; i+=blockRows)
    {

        // mat_tile[i][threadIdx.x] = matrix[(row+i)*width+col];
        // mat_tile[0][0] = 3;
        mat_tile[threadIdx.y+i][threadIdx.x] = matrix[(row+i)*width+col];
    }
    for (int i = 0; i < blockRows; i++)
    {
        vec_tile[i] = vector[threadIdx.y+i];
    }

    __syncthreads();

    double sum = 0;
    for (int i = 0; i < tileDim; i++)
    {
        sum += mat_tile[row][i]*vec_tile[i];
    }
    atomicAdd((output+row),sum);
}

void goodCuda()
{
    int block_size = 32;
    int size = block_size*32;

    double *h_output,*h_matrix,*h_vector;
    h_matrix = (double*)malloc(sizeof(double)*size*size);
    h_vector = (double*)malloc(sizeof(double)*size);
    h_output = (double*)malloc(sizeof(double)*size);

    double *d_output,*d_matrix,*d_vector;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            h_matrix[i*size+j] = 1.0;
        }
        h_vector[i] = 2.0;
        h_output[i] = 0.0;
    }    


    dim3 dimGrid(size/block_size, size/block_size, 1);
    dim3 dimBlock(block_size, block_size, 1);

    checkCuda(cudaMalloc(&d_matrix,sizeof(double)*size*size));
    checkCuda(cudaMalloc(&d_vector,sizeof(double)*size));
    checkCuda(cudaMalloc(&d_output,sizeof(double)*size));
    cudaDeviceSynchronize();
    
    checkCuda(cudaMemcpy(d_matrix, h_matrix, sizeof(double)*size*size, cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(d_vector, h_vector, sizeof(double)*size, cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();

    goodMVM<<<dimGrid,dimBlock>>>(d_matrix,d_vector,d_output,size);

    cudaDeviceSynchronize();
    checkCuda(cudaMemcpy(h_output, d_output, sizeof(double)*size, cudaMemcpyDeviceToHost));
    cudaDeviceSynchronize();
    

    for (int i = 0; i < size; i++)
    {
        std::cout<<h_output[i]<<',';
    }
    std::cout<<std::endl;

    delete[] h_matrix;
    delete[] h_vector;
    delete[] h_output;

    cudaFree(&d_output);
    cudaFree(&d_matrix);
    cudaFree(&d_vector);
}


int main()
{
    // naivCuda();
    goodCuda();
}