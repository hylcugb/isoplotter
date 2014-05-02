#pragma once

#include <stdint.h>

#ifndef __device__
#define __device__
#define __host__
#endif

class gc_sum_t {
public:
    __device__ __host__ void split(uint64_t midpoint, gc_sum_t &left, gc_sum_t &right) const;

    void gpu_alloc();
    void gpu_dispose();
    gc_sum_t gpu_copy();

    double range(uint64_t begin, uint64_t end) const;

    __device__ __host__ inline double get(uint64_t offset) const {
        return cumsum[begin + offset] - sum_begin;
    }
    __device__ __host__ inline double get_reverse(uint64_t offset) const {
        return sum_end - cumsum[begin + offset - 1];
    }
    __device__ __host__ inline uint64_t length() const {
        return end - begin;
    }

    double *cumsum;
    double *gpu_cumsum;
    uint64_t begin;
    uint64_t end;
    double sum_begin;
    double sum_end;
};

gc_sum_t create_gc_sum(double *gc_win_sum, uint64_t nwins, bool gpu);
void dispose_gc_sum(gc_sum_t gc_sum);

