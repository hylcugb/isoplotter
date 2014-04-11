#pragma once

#include "gc_sum.h"

#include <list>
#include <stdint.h>

typedef enum {
    A = 0x0, C = 0x1, G = 0x2, T = 0x3, __nbase = 4
} base_t;

typedef struct {
    uint64_t n[__nbase];
} base_count_t;

struct n_segment_t {
    uint64_t start;
    uint64_t end;

    uint64_t len() {
        return end - start;
    }
};

struct segment_t {
    uint64_t start;
    uint64_t end;
    double entropy;
    gc_sum_t gc_sum;
    gc_sum_t gc_sum2;

    uint64_t len() const {
        return end - start;
    }

    void split(uint64_t midpoint, segment_t &left, segment_t &right) const;
};

std::list<segment_t> find_isochores(char *seq, size_t seqlen, size_t winlen = 32, size_t mindomainlen = 3008, size_t min_n_domain_len = 50000);
