#include <math.h>

#include <assert.h>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>

#include "isoplotter.h"
#include "util.h"

using namespace std;

static base_count_t count(char *seq, size_t from, size_t to) {
    base_count_t result = {0,0,0,0};

    for(size_t i = from; i < to; i++) {
        switch(seq[i]) {
        case 'A': result.n[A]++; break;
        case 'C': result.n[C]++; break;
        case 'T': result.n[T]++; break;
        case 'G': result.n[G]++; break;
        }
    }

    return result;
}

static size_t base_count_len(base_count_t bc) {
    return bc.n[G] + bc.n[A] + bc.n[T] + bc.n[C];
}


static bool entropy(base_count_t bc, float *ret_e) {
    size_t n = base_count_len(bc);
    float e = 0;

    for(base_t base = (base_t)0; base < __nbase; base = (base_t)(base + 1)) {
        size_t count = bc.n[base];
        if(count == 0) {
            *ret_e = 0;
            return false;
        }

        float avg = (float)count / n;
        e += avg * log2(avg);
    }

    *ret_e = -e;

    return true;
}

float tgc(char *seq, size_t seqLen) {
    size_t n = seqLen;
    base_count_t nbases = count(seq, 0, n);
    size_t gc = nbases.n[G] + nbases.n[C];
    float e;

    dbf("nG=%lu nA=%lu nT=%lu nC=%lu", nbases.n[G], nbases.n[A], nbases.n[T], nbases.n[C]);

    if(entropy(nbases, &e)) {
        dbf("len=%lu gc=%lu", base_count_len(nbases), gc);
        dbf("entropy=%f", e);
        return e;
    } else {
        db("Failed finding e");
        return -1;
    }
}

template<typename T>
float stddevf(T *data, size_t n) {
    double sum = 0;
    float mean;

    for(size_t i = 0; i < n; i++) {
        sum += data[i];
    }
    mean = (float)(sum / n);

    sum = 0;
    for(size_t i = 0; i < n; i++) {
        float diff = data[i] - mean;
        sum += diff * diff;
    }

    return (float)sqrt(sum / n);
}

size_t win_count(char *seq, size_t seqlen, size_t winlen) {
    return (seqlen / winlen) + (seqlen % winlen ? 1 : 0);
}

double *win_gc(char *seq, size_t seqlen, size_t winlen) {
    size_t nwins = win_count(seq, seqlen, winlen);
    double *gc = new double[nwins];

    for(size_t iwin = 0; iwin < nwins; iwin++) {
        size_t win_from = iwin * winlen;
        size_t win_to = min(win_from + winlen, seqlen);
        base_count_t nbases = count(seq, win_from, win_to);
	
        gc[iwin] = double(nbases.n[G] + nbases.n[C]) / winlen;
    }

    return gc;
}

double *cumulative_gcmean(double *gcwins, size_t len, bool reverse) {
    double *result = new double[len];
    double accum = 0;

    if(!reverse) {
        for(size_t i = 0; i < len; i++) {
            accum += gcwins[i];
            result[i] = accum / (i+1);
            errif_(result[i] < 0.0f || result[i] > 1.0f);
        }
    } else {
        for(size_t i = 1; i <= len; i++) {
            size_t index = len - i;
            accum += gcwins[index];
            result[index] = accum / i;
            errif_(result[index] < 0.0f || result[index] > 1.0f);
        }
    }

    return result;
}


double entropy(double gcmean, double atmean) {
    return -(gcmean * log2(gcmean) + atmean * log2(atmean));
}

double entropy(double gcmean) {
    return entropy(gcmean, 1 - gcmean);
}

double entropy(double *gc, size_t len) {
    double gcsum = 0;
    double atsum = 0;

    for(size_t i = 0; i < len; i++) {
        gcsum += gc[i];
        atsum += 1 - gc[i];
    }

    double gcmean = gcsum / len;
    double atmean = atsum / len;	

    return entropy(gcmean, atmean);
}

template<typename T>
double stddev(T *data, size_t n) {
    double sum = 0;
    double mean;

    for(size_t i = 0; i < n; i++) {
        sum += data[i];
    }
    mean = (double)(sum / n);

    sum = 0;
    for(size_t i = 0; i < n; i++) {
        double diff = data[i] - mean;
        sum += diff * diff;
    }

    return (double)sqrt(sum / n);
}

double stddev(segment_t segment, double *gc) {
    return stddev(gc + segment.start, segment.len());
}

double dynamic_threshold(segment_t segment, double *gc, size_t winlen) {
    return double(-0.975*log(segment.len()*winlen) + (0.699*log(stddev(segment, gc)) + 1.1866));
}

void test_cumulative_gcmean(double *gc, size_t len) {
    len = min(len, 10000ul);

    double *cum = cumulative_gcmean(gc, len, false);
    for(size_t i = 0; i < len; i++) {
        double e = entropy(gc, i+1);
        double ecum = entropy(cum[i]);
        errif(!equals(e, ecum, 0.000001), "Mismatch at %ld, %f %f", i, e, ecum);
    }

    cum = cumulative_gcmean(gc, len, true);
    for(size_t i = len; i > 0; i--) {
        double e = entropy(gc + i - 1, len - i + 1);
        double ecum = entropy(cum[i - 1]);
        errif(!equals(e, ecum, 0.000001), "Reverse mismatch at %ld, %f %f", i, e, ecum);
    }
}

segment_t create_segment(size_t start, size_t end, double entropy) {
    return {start, end, entropy};
}

pair<segment_t, segment_t> divide_segment(double *gc, segment_t segment, double &Djs) {
    Djs = 0.0f;
    pair<segment_t, segment_t> result;
    size_t result_i = 0;
    double *Djs_vec = new double[segment.len()];

    pair<double *, double *> gcmean = {
        cumulative_gcmean(gc + segment.start, segment.len(), false),
        cumulative_gcmean(gc + segment.start, segment.len(), true)
    };

    for(size_t i = 0; i < segment.len() - 1; i++) {
        size_t midpoint = segment.start + i + 1;
        pair<segment_t, segment_t> candidate = {
            create_segment(segment.start, midpoint, entropy(gcmean.first[i])),
            create_segment(midpoint, segment.end, entropy(gcmean.second[i+1]))
        };
        pair<double, double> weightedEntropy = {
            (double(candidate.first.len())/segment.len()) * candidate.first.entropy,
            (double(candidate.second.len())/segment.len()) * candidate.second.entropy
        };

        double candidateDjs = segment.entropy - (weightedEntropy.first + weightedEntropy.second);
        Djs_vec[i] = candidateDjs;
        if(i == 1 || candidateDjs > Djs) {
            result = candidate;
            result_i = i;
            Djs = candidateDjs;
        }
    }

/*
  cout << setprecision(12);
  cout << "  Hs1=" << result.first.entropy << ", Hs2=" << result.second.entropy;
  cout << ", Djs(" << (result_i - 9) << ":" << (result_i + 11) << ")=[";
  for(size_t i = result_i - 10; i <= result_i + 10; i++) {
  cout << Djs_vec[i] << ",";
  }
  cout << "]" << endl;
*/
    delete [] Djs_vec;
    delete [] gcmean.first;
    delete [] gcmean.second;

    return result;
}

list<segment_t> merge(list<segment_t> segments,
                      list<n_segment_t> n_segments,
                      uint64_t min_n_domain_len) {
    list<segment_t> result;
    
    #if 0
    cout << "N" << endl;
    for(auto n_segment: n_segments) {
        cout << n_segment.start << "\t" << n_segment.end << endl;
    }

    cout << endl;
    cout << "non-N" << endl;
    for(auto segment: segments) {
        cout << segment.start << "\t" << segment.end << endl;
    }
    #endif

    for(auto n_segment: n_segments) {
        bool small_n_segment = n_segment.len() <= min_n_domain_len;

        //cout << "N: " << n_segment.start << "\t" << n_segment.end << endl;
        for(auto it_segment = segments.begin(); it_segment != segments.end(); ) {
            segment_t &segment = *it_segment;

            //cout << "  seg: " << segment.start << "\t" << segment.end << endl;

            if(segment.start >= n_segment.start) {
                // Segment is after N island, so simply shift the segment forward
                // by size of N island.
                segment.start += n_segment.len();
                segment.end += n_segment.len();

                if( small_n_segment && (n_segment.end == segment.start) ) {
                    segment.start = n_segment.start;
                }

                ++it_segment;
            } else if(segment.end <= n_segment.start) {
                // Segment is before N island, so it requires no further processing.

                if( small_n_segment && (n_segment.start == segment.end) ) {
                    segment.end = n_segment.end;
                }

                result.push_back(segment);
                auto it_erase = it_segment++;
                segments.erase(it_erase);
            } else {
                // N island is in middle of segment.
                //assert( (segment.start < n_segment.start) && (segment.end >= n_segment.start) );
                        
                if( small_n_segment ) {
                    segment.end += n_segment.len();
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

        if(!small_n_segment) {
            result.push_back( {n_segment.start, n_segment.end, 0.0} );
        }
    }

    return result;
}

void create_win_gc(char *seq, size_t seqlen, size_t winlen, double **out_gc, size_t &out_nwins, list<n_segment_t> &out_n_segments) {
    int gc_count = 0;
    int win_bases_count = 0;
    
    out_nwins = 0;
    *out_gc = new double[seqlen / winlen];
    
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
                gc_count++;
            }
            if(++win_bases_count == winlen) {
                (*out_gc)[out_nwins++] = double(gc_count) / winlen;
                gc_count = 0;
                win_bases_count = 0;
            }

            i++;
        }
    }

    // Note: we ignore any leftover bases for the final window
}

list<segment_t> find_isochores(char *seq, size_t seqlen, size_t winlen, size_t mindomainlen, size_t min_n_domain_len) {
    mindomainlen /= winlen;

    double *gc;
    size_t nwins;
    list<n_segment_t> n_segments;

    create_win_gc(seq, seqlen, winlen, &gc, nwins, n_segments);

    list<segment_t> segments = { create_segment(0, nwins, entropy(gc, nwins)) };

    for(list<segment_t>::iterator it = segments.begin(); it != segments.end(); ) {
        segment_t segment = *it;
        double Djs;
        pair<segment_t, segment_t> subsegments = divide_segment(gc, segment, Djs);

        double dt = dynamic_threshold(segment, gc, winlen);
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

    delete [] gc;

    for(auto &seg: segments) {
        seg.start = seg.start * winlen;
        seg.end = seg.end * winlen;
    }
    segments.back().end = nwins * winlen;

    list<segment_t> result = merge(segments, n_segments, min_n_domain_len);

    result.back().end = (result.back().end / winlen) * winlen;

    return result;
}
