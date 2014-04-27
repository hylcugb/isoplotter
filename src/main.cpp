#include "util.h"

#include <assert.h>

#include "isoplotter.h"
#include "pna.h"

int main(int argc, char* argv[])
{
    assert(argc == 2);

    const char *path = argv[1];

    seqio::PnaReader reader(path);

    for(int i = 0; i < (int)reader.getSequenceCount(); i++) {
    //int i = 13; {
        std::shared_ptr<seqio::PnaSequenceReader> seq = reader.openSequence(i);
        char *bases = new char[seq->size()];
        seq->read(bases, seq->size());

        std::list<segment_t> segments = find_isochores(bases,
                                                       seq->size(),
                                                       32,
                                                       3008,
                                                       50000);
        for(auto seg: segments) {
            dbf("%d\t%10lu\t%10lu\t%10lu\t%.3f\t%.4f", i+1, seg.start+1, seg.end, seg.len(), seg.gc_mean, seg.gc_stddev);
        }

        delete [] bases;
    }
}
