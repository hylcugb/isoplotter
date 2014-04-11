#include "gc_sum.h"

#include <assert.h>
#include <cuda.h>
#include <iostream>

#define xcuda(stmt) {                                                   \
        cudaError_t err = stmt;                                         \
        if (err != cudaSuccess) {                                       \
            std::cerr << __FILE__ << ":" << __LINE__ << ": Failed to run " << #stmt << ". Reason: " << cudaGetErrorString(err) << std::endl; \
            exit(1);                                                    \
        }                                                               \
    }

gc_sum_t create_gc_sum(double *gc_win_sum, uint64_t nwins, bool gpu) {
    gc_sum_t result;
    result.cumsum = gc_win_sum + 1;
    result.begin = 0;
    result.end = nwins;
    result.sum_begin = 0;
    result.sum_end = gc_win_sum[nwins];

    if(gpu) {
        result.gpu_alloc();
    } else {
        result.gpu_cumsum = NULL;
    }

    return result;
}

void dispose_gc_sum(gc_sum_t gc_sum) {
    delete [] (gc_sum.cumsum - 1);
    if(gc_sum.gpu_cumsum) {
        gc_sum.gpu_dispose();
    }
}

void gc_sum_t::gpu_alloc() {
    size_t sizeof_ = sizeof(double) * (end - begin + 1);
    xcuda( cudaMalloc((void**)&gpu_cumsum, sizeof_) );
    xcuda( cudaMemcpy(gpu_cumsum, cumsum - 1, sizeof_, cudaMemcpyHostToDevice) );
}

void gc_sum_t::gpu_dispose() {
    xcuda( cudaFree(gpu_cumsum) );
}

gc_sum_t gc_sum_t::gpu_copy() {
    gc_sum_t result = *this;
    result.cumsum = gpu_cumsum + 1;
    return result;
}

__device__ __host__ void gc_sum_t::split(uint64_t midpoint, gc_sum_t &left, gc_sum_t &right) const {
    left.cumsum = this->cumsum;
    left.gpu_cumsum = this->gpu_cumsum;
    left.begin = this->begin;
    left.end = this->begin + midpoint;
    left.sum_begin = this->sum_begin;
    left.sum_end = left.cumsum[left.end - 1];

    right.cumsum = this->cumsum;
    right.gpu_cumsum = this->gpu_cumsum;
    right.begin = left.end;
    right.end = this->end;
    right.sum_begin = left.sum_begin + left.get(midpoint - 1);
    right.sum_end = right.cumsum[right.end - 1];
}

double gc_sum_t::range(uint64_t begin, uint64_t end) const {
    return cumsum[this->begin + end - 1] - cumsum[this->begin + begin - 1];
}
