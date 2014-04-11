#include <math.h>

#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>

#include "gc_sum.h"
#include "isoplotter.h"
#include "util.h"

using namespace std;

#define WINLEN 32 //tmp

double gcmean(uint64_t gcsum, size_t winlen, uint64_t nwins) {
    return double(gcsum) / (winlen * nwins);
}

double entropy(double gcmean, double atmean) {
    return -(gcmean * log2(gcmean) + atmean * log2(atmean));
}

double entropy(double gcmean) {
    return entropy(gcmean, 1 - gcmean);
}

void segment_t::split(uint64_t midpoint, segment_t &left, segment_t &right) const {
    assert( (midpoint > 0) && (midpoint < len()) );

    left.start = this->start;
    left.end = this->start + midpoint;
    left.entropy = ::entropy( ::gcmean( gc_sum.get(midpoint-1), WINLEN, left.len()) );

    right.start = left.end;
    right.end = this->end;
    right.entropy = ::entropy( ::gcmean( gc_sum.get_reverse(midpoint), WINLEN, right.len() ) );

    gc_sum.split(midpoint, left.gc_sum, right.gc_sum);
    gc_sum2.split(midpoint, left.gc_sum2, right.gc_sum2);
}

double stddev(gc_sum_t gc_sum, gc_sum_t gc_sum2) {
    uint64_t n = gc_sum.length();
    double sum = gc_sum.get(n - 1) / double(n);
    double sum2 = gc_sum2.get(n - 1) / double(n * n);
    
    return sqrt( (sum2 - sum*sum/n) / n );
}

double dynamic_threshold(segment_t segment, size_t winlen) {
    return double(-0.975*log(segment.len()*winlen) + (0.699*log(stddev(segment.gc_sum, segment.gc_sum2)) + 1.1866));
}

segment_t create_segment(size_t start, size_t end, double entropy) {
    return {start, end, entropy};
}

pair<segment_t, segment_t> divide_segment(segment_t segment, double &Djs) {
    Djs = 0.0f;
    pair<segment_t, segment_t> result;

    for(size_t i = 0; i < segment.len() - 1; i++) {
        size_t midpoint = segment.start + i + 1;

        pair<segment_t, segment_t> candidate;
        segment.split(i+1, candidate.first, candidate.second);        

        pair<double, double> weightedEntropy = {
            (double(candidate.first.len())/segment.len()) * candidate.first.entropy,
            (double(candidate.second.len())/segment.len()) * candidate.second.entropy
        };

        double candidateDjs = segment.entropy - (weightedEntropy.first + weightedEntropy.second);
        if(i == 1 || candidateDjs > Djs) {
            result = candidate;
            Djs = candidateDjs;
        }
    }

    return result;
}

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
            result.push_back( {n_segment.start, n_segment.end, 0.0} );
        }
    }

    // There will be some remaining if the sequence didn't end with an N island.
    for(auto &segment: segments) {
        result.push_back(segment);
    }

    {
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

void create_win_gc(char *seq, size_t seqlen, size_t winlen, uint64_t **out_gc_, uint64_t **out_gc2, size_t &out_nwins, list<n_segment_t> &out_n_segments) {
    int win_bases_count = 0;
    uint64_t gc_count_ = 0;
    
    out_nwins = 0;
    *out_gc_ = new uint64_t[seqlen / winlen + 1];
    (*out_gc_)[0] = 0;

    *out_gc2 = new uint64_t[seqlen / winlen + 1];
    (*out_gc2)[0] = 0;
    
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
                gc_count_++;
            }
            if(++win_bases_count == winlen) {
                (*out_gc_)[out_nwins + 1] = gc_count_;
                uint64_t prev_win_count = gc_count_ - (*out_gc_)[out_nwins];
                (*out_gc2)[out_nwins + 1] = (*out_gc2)[out_nwins] + prev_win_count * prev_win_count;
                out_nwins++;
                win_bases_count = 0;
            }

            i++;
        }
    }

    // Note: we ignore any leftover bases beyond a multiple of window length.
    // Also, we ignore the final window in order to recreate results of Matlab version.
    out_nwins--;
}

list<segment_t> find_isochores(char *seq, size_t seqlen, size_t winlen, size_t mindomainlen, size_t min_n_domain_len) {
    assert(winlen==WINLEN); // tmp

    mindomainlen /= winlen;

    uint64_t *gc_;
    uint64_t *gc2;
    size_t nwins;
    list<n_segment_t> n_segments;

    create_win_gc(seq, seqlen, winlen, &gc_, &gc2, nwins, n_segments);

    gc_sum_t sum = create_gc_sum(gc_, nwins);
    gc_sum_t sum2 = create_gc_sum(gc2, nwins);

    list<segment_t> segments = { create_segment(0, nwins, entropy( gcmean(sum.get(nwins-1), winlen, nwins) )) };
    segments.front().gc_sum = sum;
    segments.front().gc_sum2 = sum2;

    for(list<segment_t>::iterator it = segments.begin(); it != segments.end(); ) {
        segment_t segment = *it;
        double Djs;
        pair<segment_t, segment_t> subsegments = divide_segment(segment, Djs);

        double dt = dynamic_threshold(segment, winlen);
        double log_Djs = log(Djs);

        //cout << "(" << (segment.start+1) << "," << segment.end << ") -> " << subsegments.first.len() << "; dt=" << dt << ", log(Djs)=" << log_Djs << endl;
	    
        if( (subsegments.first.len() <= mindomainlen)
            || (subsegments.second.len() <= mindomainlen) ) {
            //cout << "  REJECT" << endl;
            ++it;
        } else if( log_Djs <= dt ) {
            ++it;
            //cout << "  REJECT" << endl;
        } else {
            list<segment_t>::iterator it_insert = it;

            *it_insert = subsegments.first;
            segments.insert(++it_insert, subsegments.second);
        }
    }

    delete [] gc_;
    delete [] gc2;

    for(auto &seg: segments) {
        seg.start = seg.start * winlen;
        seg.end = seg.end * winlen;
    }
    segments.back().end = nwins * winlen;

    list<segment_t> result = merge(segments, n_segments, min_n_domain_len);

    // Sanity check
    assert(result.front().start == 0);
    for(auto it = result.begin(); next(it) != result.end(); it++) {
        assert(it->end == next(it)->start);
    }
    // If it ends with an N island, it will have the original seqlen.
    // Otherwise, it will be truncated to a window border.
    assert( (result.back().end == seqlen) || (result.back().end == ( (seqlen / winlen - 1) * winlen )) );

    // Round down to window border.
    result.back().end = (result.back().end / winlen) * winlen;

    return result;
}
