#pragma once

#include <list>
#include <stdint.h>

typedef enum {
    A = 0x0, C = 0x1, G = 0x2, T = 0x3, __nbase = 4
} base_t;

typedef struct {
    uint64_t n[__nbase];
} base_count_t;

struct segment_t {
    uint64_t start;
    uint64_t end;
    double entropy;

    uint64_t len() {
        return end - start;
    }
};

std::list<segment_t> find_isochores(char *seq, size_t seqlen, size_t winlen = 32, size_t mindomainlen = 3008);
