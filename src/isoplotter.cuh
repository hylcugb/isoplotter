#pragma once

#include <stdint.h>

#ifndef __device__
#define __device__
#define __host__
#endif

__device__ __host__ double entropy(double gcmean, double atmean);
__device__ __host__ double entropy(double gcmean);

struct segment_t { // todo: make a "result" segment that has gc_mean... make this a purely
                   //       intermediate/computational structure.
    uint64_t start;
    uint64_t end;
    double entropy;
    double gc_mean;
    double gc_stddev;
    gc_sum_t gc_sum;
    gc_sum_t gc_sum2;

    __device__ __host__ uint64_t len() const {
        return end - start;
    }

    void split(uint64_t midpoint, segment_t &left, segment_t &right) const;
};

uint64_t find_best_split_cuda(segment_t segment);
