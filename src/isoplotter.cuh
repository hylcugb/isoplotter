#pragma once

#include <stdint.h>

#ifndef __device__
#define __device__
#define __host__
#endif

__device__ __host__ double entropy(double gcmean, double atmean);
__device__ __host__ double entropy(double gcmean);

struct segment_t {
    uint64_t start;
    uint64_t end;
    double entropy;
    gc_sum_t gc_sum;
    gc_sum_t gc_sum2;

    uint64_t len() const {
        return end - start;
    }

    __device__ __host__ void split(uint64_t midpoint, segment_t &left, segment_t &right) const;
};
