#pragma once

#include "gc_sum.h"
#include "isoplotter.cuh"

#include <list>

std::list<segment_t> find_isochores(char *seq,
                                    size_t seqlen,
                                    size_t winlen = 32,
                                    size_t mindomainlen = 3008,
                                    size_t min_n_domain_len = 50000);
