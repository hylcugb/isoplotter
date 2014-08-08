#include "util.h"
#include "isoplotter.h"

#include <iostream>

#include <assert.h>
#include <seqio.h>

using namespace std;

int main(int argc, char* argv[])
{
    if(argc == 1) {
        cerr << "usage: " << argv[0] << " input_path" << endl;
        cerr << endl;
        cerr << "  example: ./isoplotter /genomes/human.pna > results.out" << endl;
        return 1;
    }

    assert(argc == 2);

    const char *sequence_path = argv[1];

    char *sequence_buffer = nullptr;
    uint64_t sequence_buffer_length;
    uint64_t sequence_length;

    seqio_sequence_iterator iterator;
    seqio_sequence_options sequence_options = SEQIO_DEFAULT_SEQUENCE_OPTIONS;
    sequence_options.base_transform = SEQIO_BASE_TRANSFORM_CAPS_GATCN;
    seqio_create_sequence_iterator(sequence_path,
                                   sequence_options,
                                   &iterator);

    seqio_sequence sequence;
    for(int i = 0;
        (SEQIO_SUCCESS == seqio_next_sequence(iterator, &sequence))
            && (sequence != nullptr);
        i++) {

        char const *name, *comment;
        seqio_const_dictionary metadata;
        seqio_get_metadata(sequence, &metadata);
        seqio_get_value(metadata, SEQIO_KEY_NAME, &name);
        seqio_get_value(metadata, SEQIO_KEY_COMMENT, &comment);

        seqio_read_all(sequence, &sequence_buffer, &sequence_buffer_length, &sequence_length);

        std::list<segment_t> segments = find_isochores(sequence_buffer,
                                                       sequence_length,
                                                       32,
                                                       3008,
                                                       50000);
        for(auto seg: segments) {
            dbf("%d\t%10lu\t%10lu\t%10lu\t%.3f\t%.4f\t%d", i+1, seg.start+1, seg.end, seg.len(), seg.gc_mean, seg.gc_stddev, seg.homo);
        }

        seqio_dispose_sequence(&sequence);
    }

    seqio_dispose_sequence_iterator(&iterator);
    seqio_dispose_buffer(&sequence_buffer);
}
