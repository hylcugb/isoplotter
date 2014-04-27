#include <math.h>

#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>

#include "gc_sum.h"
#include "isoplotter.cuh"
#include "isoplotter.h"
#include "util.h"

using namespace std;

struct n_segment_t {
    uint64_t start;
    uint64_t end;

    uint64_t len() {
        return end - start;
    }
};

void segment_t::split(uint64_t midpoint, segment_t &left, segment_t &right) const {
    assert( (midpoint > 0) && (midpoint < len()) );

    left.start = this->start;
    left.end = this->start + midpoint;
    left.entropy = ::entropy( gc_sum.get(midpoint-1) / left.len() );

    right.start = left.end;
    right.end = this->end;
    right.entropy = ::entropy( gc_sum.get_reverse(midpoint) / right.len() );

    gc_sum.split(midpoint, left.gc_sum, right.gc_sum);
    gc_sum2.split(midpoint, left.gc_sum2, right.gc_sum2);
}

double mean(gc_sum_t gc_sum) {
    uint64_t n = gc_sum.length();
    double sum = gc_sum.get(n - 1);

    return sum / n;
}

double stddev(gc_sum_t gc_sum, gc_sum_t gc_sum2, bool sample = false) {
    uint64_t n = gc_sum.length();
    double sum = gc_sum.get(n - 1);
    double sum2 = gc_sum2.get(n - 1);

    double divisor = sample ? n - 1 : n;
    
    return sqrt( (sum2 - sum*sum/n) / divisor );
}

double dynamic_threshold(segment_t segment, size_t winlen) {
    return double(-0.975*log(segment.len()*winlen) + (0.699*log(stddev(segment.gc_sum, segment.gc_sum2)) + 1.1866));
}

segment_t create_segment(size_t start, size_t end, double entropy) {
    return {start, end, entropy};
}


uint64_t find_best_split_cpu(segment_t segment) {
    double Djs = 0.0f;
    size_t i_best = 2;

    for(size_t i = 2; i < segment.len(); i++) {
        uint64_t len_left = i;
        uint64_t len_right = segment.len() - i;
        double entropy_left = entropy( segment.gc_sum.get(i-1) / len_left );
        double entropy_right = entropy( segment.gc_sum.get_reverse(i) / len_right );
        double weighted_entropy_left = double(len_left) / segment.len() * entropy_left;
        double weighted_entropy_right = double(len_right) / segment.len() * entropy_right;

        double candidateDjs = segment.entropy - (weighted_entropy_left + weighted_entropy_right);
        if(i == 2 || candidateDjs > Djs) {
              Djs = candidateDjs;
            i_best = i;
        }
    }

    return i_best;
}

pair<segment_t, segment_t> divide_segment(segment_t segment, double &Djs) {
    uint64_t i;
    if(segment.len() > 2000) {
        i = find_best_split_cuda(segment);
    } else {
        i = find_best_split_cpu(segment);
    }

    pair<segment_t, segment_t> result;
    segment.split(i, result.first, result.second);        

    pair<double, double> weightedEntropy = {
        (double(result.first.len())/segment.len()) * result.first.entropy,
        (double(result.second.len())/segment.len()) * result.second.entropy
    };

    Djs = segment.entropy - (weightedEntropy.first + weightedEntropy.second);

    return result;
}


/*
pair<segment_t, segment_t> divide_segment(segment_t segment, double &Djs) {
    Djs = 0.0f;
    pair<segment_t, segment_t> result;

    for(size_t i = 2; i < segment.len(); i++) {
        pair<segment_t, segment_t> candidate;
        segment.split(i, candidate.first, candidate.second);        

        pair<double, double> weightedEntropy = {
            (double(candidate.first.len())/segment.len()) * candidate.first.entropy,
            (double(candidate.second.len())/segment.len()) * candidate.second.entropy
        };

        double candidateDjs = segment.entropy - (weightedEntropy.first + weightedEntropy.second);
        if(i == 2 || candidateDjs > Djs) {
            result = candidate;
            Djs = candidateDjs;
        }
    }

    return result;
}
*/

list<segment_t> merge(list<segment_t> segments,
                      list<n_segment_t> n_segments,
                      uint64_t min_n_domain_len) {
    list<segment_t> result;
    
    for(auto n_segment: n_segments) {
        bool small_n_segment = n_segment.len() <= min_n_domain_len;
        bool merged = false;

        for(auto it_segment = segments.begin(); it_segment != segments.end(); ) {
            segment_t &segment = *it_segment;

            if(segment.start >= n_segment.start) {
                // Segment is after N island, so simply shift the segment forward
                // by size of N island.
                segment.start += n_segment.len();
                segment.end += n_segment.len();

                if( small_n_segment && (n_segment.end == segment.start) ) {
                    // It's a small N island, so just prepend it to this segment.
                    segment.start = n_segment.start;
                    merged = true;
                }

                ++it_segment;
            } else if(segment.end <= n_segment.start) {
                // Segment is before N island, so it requires no further processing.

                if( small_n_segment && (n_segment.start == segment.end) ) {
                    // It's a small N island, so just append it to this segment.
                    segment.end = n_segment.end;
                    merged = true;
                }

                result.push_back(segment);
                auto it_erase = it_segment++;
                segments.erase(it_erase);
            } else {
                // N island is in middle of segment.
                        
                if( small_n_segment ) {
                    // It's a small N island, so just insert it into this segment.
                    segment.end += n_segment.len();
                    merged = true;
                } else {
                    segment_t segment_before_n = {segment.start, n_segment.start, segment.entropy};
                    result.push_back(segment_before_n);

                    uint64_t len_after_n = segment.end - segment_before_n.end;
                    segment.start = n_segment.end;
                    segment.end = segment.start + len_after_n;
                }
                ++it_segment;
            }
        }

        if(!merged) {
            result.push_back( {n_segment.start, n_segment.end, 0.0, 0.0, 0.0} );
        }
    }

    // There will be some remaining if the sequence didn't end with an N island.
    for(auto &segment: segments) {
        result.push_back(segment);
    }

    if(result.size() > 1) {
        // Since we drop the last 1 or 2 windows initially, there's going to be
        // a gap if the last segment is an N segment. Consistent with Matlab version,
        // we just insert a tiny segment with the same entropy.
        auto it_last = --result.end();
        auto it_prev = --(--result.end());

        if(it_prev->end != it_last->start) {
            segment_t filler;
            filler.start = it_prev->end;
            filler.end = it_last->start;
            filler.entropy = it_prev->entropy;
            result.insert(it_last, filler);
        }
    }

    return result;
}

void create_win_gc(char *seq, size_t seqlen, size_t winlen, double **out_gc_, double **out_gc2, size_t &out_nwins, list<n_segment_t> &out_n_segments) {
    size_t win_bases_count = 0;
    double gc_mean_accum = 0.0;
    double gc_mean_accum2 = 0.0;
    uint64_t gc_count_win = 0;
    
    out_nwins = 0;
    *out_gc_ = new double[seqlen / winlen + 1];
    (*out_gc_)[0] = 0.0;

    *out_gc2 = new double[seqlen / winlen + 1];
    (*out_gc2)[0] = 0.0;
    
    size_t i = 0;
    while(i < seqlen) {
        char base = seq[i];

        if(base == 'N') {
            n_segment_t n_segment;
            n_segment.start = i++;
            for(; (i < seqlen) && (seq[i] == 'N'); i++) {
            }
            n_segment.end = i;
            out_n_segments.push_back(n_segment);
        } else {
            if( (base == 'G') || (base == 'C') ) {
                gc_count_win++;
            }
            if(++win_bases_count == winlen) {
                double gc_mean = double(gc_count_win) / winlen;
                gc_mean_accum += gc_mean;
                gc_mean_accum2 += gc_mean * gc_mean;
                (*out_gc_)[out_nwins + 1] = gc_mean_accum;
                (*out_gc2)[out_nwins + 1] = gc_mean_accum2;
                out_nwins++;
                win_bases_count = 0;
                gc_count_win = 0;
            }

            i++;
        }
    }

    // Note: we ignore any leftover bases beyond a multiple of window length.
    // Also, we ignore the final window in order to recreate results of Matlab version.
    out_nwins--;
}

list<segment_t> find_isochores(char *seq, size_t seqlen, size_t winlen, size_t mindomainlen, size_t min_n_domain_len) {
    mindomainlen /= winlen;

    double *gc_;
    double *gc2;
    size_t nwins;
    list<n_segment_t> n_segments;

    create_win_gc(seq, seqlen, winlen, &gc_, &gc2, nwins, n_segments);

    gc_sum_t sum = create_gc_sum(gc_, nwins, true);
    gc_sum_t sum2 = create_gc_sum(gc2, nwins, false);

    list<segment_t> segments = { create_segment(0, nwins, entropy( sum.get(nwins-1) / nwins )) };
    segments.front().gc_sum = sum;
    segments.front().gc_sum2 = sum2;

    for(list<segment_t>::iterator it = segments.begin(); it != segments.end(); ) {
        segment_t segment = *it;
        if(segment.len() <= (mindomainlen*2)) {
            // No point in trying to divide it.
            ++it;
            continue;
        }

        double Djs;
        pair<segment_t, segment_t> subsegments = divide_segment(segment, Djs);

        double dt = dynamic_threshold(segment, winlen);
        double log_Djs = log(Djs);

        //cout << "(" << (segment.start+1) << "," << segment.end << ") -> " << subsegments.first.len() << "; dt=" << dt << ", log(Djs)=" << log_Djs << endl;
	    
        if( (subsegments.first.len() <= mindomainlen)
            || (subsegments.second.len() <= mindomainlen) ) {
            // One of the segments is smaller than allowed, so don't accept
            // the division and move on to the next segment.
            ++it;
        } else if( log_Djs <= dt ) {
            // Divergence not high enough, so reject division and move on
            // to next segment.
            ++it;
        } else {
            // Division is acceptable.
            list<segment_t>::iterator it_insert = it;

            // Overwrite current segment in list with left partition.
            *it_insert = subsegments.first;
            
            // Add right partition to list.
            segments.insert(++it_insert, subsegments.second);
        }
    }

    for(auto &seg: segments) {
        seg.start = seg.start * winlen;
        seg.end = seg.end * winlen;
    }
    segments.back().end = nwins * winlen;

    for(auto &seg: segments) {
        seg.gc_mean = mean(seg.gc_sum);
        seg.gc_stddev = stddev(seg.gc_sum, seg.gc_sum2);
    }
    dispose_gc_sum(sum);
    dispose_gc_sum(sum2);

    list<segment_t> result = merge(segments, n_segments, min_n_domain_len);

    // Sanity check
    assert(result.front().start == 0);
    for(auto it = result.begin(); next(it) != result.end(); it++) {
        assert(it->end == next(it)->start);
    }
    // The end is a bit unpredictable with N island injection. Just make sure it's sane.
    // If it ends with an N island, it will have the original seqlen.
    // Otherwise, it will be truncated to a window border.
    assert( abs(seqlen - result.back().end) < (winlen * 3) );

    // Round down to window border.
    result.back().end = (result.back().end / winlen) * winlen;

    return result;
}
