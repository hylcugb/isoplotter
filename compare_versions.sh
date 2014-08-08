#!/bin/bash

#
# This script will run both versions of isoplotter and compare their output.
#
# Warning! Version 1 can only process uncompressed FASTA. Also, please only
#          use one sequence per FASTA.
#

set -e


#
# Analysis configuration
#

# useful if you just want to run the comparison code.
skip_analysis=false


#
# Input files
#
if [ -z $1 ]; then
    # See if we're on Sean's machine.
    if [ -e /media/student/poly/genome/fa/chr1.fa ]; then
        fastas=$(echo /media/student/poly/genome/fa/chr{{1..22},X,Y}.fa)
    else
        echo "usage: $(basename $0) fasta_file..."
        exit 1
    fi
fi

#
# Setup filesystem structure for version1
#
outdir=/tmp/isoplotter-compare

function v1_indir() {
    local fasta=$1
    echo $outdir/$(basename $fasta)/input
}

function v1_outdir() {
    local fasta=$1
    echo $outdir/$(basename $fasta)/output
}

function bounds() {
    local fasta=$1
    local version=$2

    echo $outdir/bounds/$(basename $fasta).v$version
}

if ! $skip_analysis; then
    #
    # Perform analysis
    #
    rm -rf $outdir
    mkdir -p $outdir
    mkdir -p $outdir/bounds

    for fasta in $fastas; do
        mkdir -p $(v1_indir $fasta)
        mkdir -p $(v1_outdir $fasta)
        ln -s $fasta $(v1_indir $fasta)/
    done

    # version 1
    (
        for fasta in $fastas; do
            version1/isoplotter pipeline $(v1_indir $fasta) $(v1_outdir $fasta)
            cat $(v1_outdir $fasta)/IsoPlotter_ns_H.txt | 
            awk '{print $2 " " $3 " " $4}' > $(bounds $fasta 1)
        done
    )

    # version 2
    (
        for fasta in $fastas; do
            version2/isoplotter $fasta |
            awk '{print $2 " " $3 " " $4}' > $(bounds $fasta 2)
        done
    )
fi

#
# Compare results
#
total_ndiffs=0
total_nlines=0

function report_diffs() {
    local label=$1
    local ndiffs=$2
    local nlines=$3

    echo "$label: Found differences in $ndiffs out of $nlines" $(python -c "print '(%.6f%%)' % (float($ndiffs)/$nlines*100)")
}

for fasta in $fastas; do
    bounds1=$(bounds $fasta 1)
    bounds2=$(bounds $fasta 2)
    ndiffs=$(diff $bounds1 $bounds2 | grep ">" | wc -l)
    nlines=$(cat $fasta | wc -l)
    if (( $ndiffs != 0 )); then
        report_diffs $(basename $fasta) $ndiffs $nlines
    fi
    total_ndiffs=$((total_ndiffs + ndiffs))
    total_nlines=$((total_nlines + nlines))
done

report_diffs "TOTAL" $total_ndiffs $total_nlines
