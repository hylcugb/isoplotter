#define CUDA
#include "isoplotter.h"

#include <cuda.h>
#include <iostream>

#define xcuda(stmt) {                                                   \
        cudaError_t err = stmt;                                         \
        if (err != cudaSuccess) {                                       \
            std::cerr << __FILE__ << ":" << __LINE__ << ": Failed to run " << #stmt << ". Reason: " << cudaGetErrorString(err) << std::endl; \
            exit(1);                                                    \
        }                                                               \
    }

static const uint Threads_Per_Block = 512;

__device__ __host__ double entropy(double gcmean, double atmean) {
    return -(gcmean * log2(gcmean) + atmean * log2(atmean));
}

__device__ __host__ double entropy(double gcmean) {
    return entropy(gcmean, 1 - gcmean);
}

__global__ void find_best_split_kernel(gc_sum_t gc_sum,
                                       double segment_entropy,
                                       double *Djs,
                                       uint64_t *is) {
    double best_Djs = 0.0;
    uint64_t best_i = 0;
    uint64_t n = gc_sum.length();

    for(uint64_t i = 2 + (blockDim.x * blockIdx.x) + threadIdx.x;
        i < n;
        i += gridDim.x * blockDim.x) {

        uint64_t len_left = i;
        uint64_t len_right = n - i;
        double entropy_left = entropy( gc_sum.get(i-1) / len_left );
        double entropy_right = entropy( gc_sum.get_reverse(i) / len_right );
        double weighted_entropy_left = double(len_left) / n * entropy_left;
        double weighted_entropy_right = double(len_right) / n * entropy_right;

        double candidateDjs = segment_entropy - (weighted_entropy_left + weighted_entropy_right);
        if(best_i == 0 || candidateDjs > best_Djs) {
            best_Djs = candidateDjs;
            best_i = i;
        }
    }

    __shared__ double Djs_shared[Threads_Per_Block];
    __shared__ uint64_t is_shared[Threads_Per_Block];

    int tid = threadIdx.x;
    Djs_shared[tid] = best_Djs;
    is_shared[tid] = best_i;

    for(int stride = Threads_Per_Block / 2; stride > 0; stride >>= 1) {
        __syncthreads();

        if(tid < stride) {
            if(Djs_shared[tid + stride] > Djs_shared[tid]) {
                Djs_shared[tid] = Djs_shared[tid + stride];
                is_shared[tid] = is_shared[tid + stride];
            }
        }
    }

    if(tid == 0) {
        Djs[blockIdx.x] = Djs_shared[0];
        is[blockIdx.x] = is_shared[0];
    }

}

uint64_t find_best_split_cuda(segment_t segment) {
    uint64_t n = segment.len() - 2;
    uint64_t nblocks = (n - 1) / Threads_Per_Block + 1;
    if(nblocks > 2048)
        nblocks = 2048;

    double *Djs = new double[nblocks];
    double *gpu_Djs;
    uint64_t *is = new uint64_t[nblocks];
    uint64_t *gpu_is;

    xcuda( cudaMalloc((void**)&gpu_Djs, nblocks*sizeof(double)) );
    xcuda( cudaMalloc((void**)&gpu_is, nblocks*sizeof(uint64_t)) );

    find_best_split_kernel<<<nblocks, Threads_Per_Block>>>(segment.gc_sum.gpu_copy(),
                                                           segment.entropy,
                                                           gpu_Djs,
                                                           gpu_is);
	xcuda( cudaPeekAtLastError() );
	xcuda( cudaThreadSynchronize() );

    xcuda( cudaMemcpy(Djs, gpu_Djs, nblocks * sizeof(double), cudaMemcpyDeviceToHost) );    
    xcuda( cudaMemcpy(is, gpu_is, nblocks * sizeof(uint64_t), cudaMemcpyDeviceToHost) );    

    double best_Djs = Djs[0];
    uint64_t best_i = is[0];

    for(uint64_t i = 1; i < nblocks; i++) {
        double d = Djs[i];
        if(d > best_Djs) {
            best_Djs = d;
            best_i = is[i];
        }
    }

    delete [] Djs;
    delete [] is;
    xcuda( cudaFree(gpu_Djs) );
    xcuda( cudaFree(gpu_is) );

    return best_i;
}

