#include "util.h"

#include <assert.h>

#include "isoplotter.h"
#include "pna.h"

int main(int argc, char* argv[])
{
    assert(argc == 2);

    //const char *path = "test/seq/chr22.pna";
    const char *path = argv[1];

    seqio::PnaReader reader(path);
    std::shared_ptr<seqio::PnaSequenceReader> seq = reader.openSequence(0);
    char *bases = new char[seq->size()];
    seq->read(bases, seq->size());

    std::list<segment_t> segments = find_isochores(bases,
                                                   seq->size(),
                                                   32,
                                                   3008,
                                                   50000);
    //dbf("%10s  %10s  %010s  %10s", "Start", "End", "Len", "Entropy");
    for(auto seg: segments) {
        dbf("%10lu\t%10lu\t%10lu\t%10f", seg.start+1, seg.end, seg.len(), seg.entropy);
    }
}
