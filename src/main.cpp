#include "util.h"

#include "isoplotter.h"
#include "pna.h"

int main(int argc, char* argv[])
{
    const char *path = "/home/sean/tmp/no_ns.pna";

    seqio::PnaReader reader(path);
    std::shared_ptr<seqio::PnaSequenceReader> seq = reader.openSequence(0);
    char *bases = new char[seq->size()];
    seq->read(bases, seq->size());

    std::list<segment_t> segments = find_isochores(bases,
                                                   seq->size(),
                                                   32,
                                                   3008);
    //dbf("%10s  %10s  %010s  %10s", "Start", "End", "Len", "Entropy");
    for(auto seg: segments) {
        dbf("%10lu  %10lu  %10lu  %10f", seg.start+1, seg.end, seg.len(), seg.entropy);
    }
}
