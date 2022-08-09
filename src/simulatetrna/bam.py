import os
import pysam
import pandas as pd


def iterate_reads(inreads):
    last_read_name = None
    reads = ()

    for read in inreads:

        if read.is_unmapped:
            continue

        read_name = read.qname

        if read_name != last_read_name:

            if last_read_name is not None:
                yield reads

            reads = set((read,))

        else:
            reads.add(read)

        last_read_name = read_name

def keep_random_alignment(infile, outfile):

    inbam = pysam.Samfile(infile, 'r')
    unsorted_outfile = outfile + 'tmp.bam'
    outbam = pysam.AlignmentFile(unsorted_outfile, "wb", template=inbam)

    for reads in iterate_reads(inbam):
        outbam.write(reads.pop())

    outbam.close()

    pysam.sort(unsorted_outfile, "-o", outfile)
    os.unlink(unsorted_outfile)
    pysam.index(outfile)


def filter_sam(infile, outfile):
    inbam = pysam.Samfile(infile, 'r')
    outbam = pysam.AlignmentFile(outfile, "w", template=inbam)

    for reads in iterate_reads(inbam):

        top_as = max([read.get_tag('AS') for read in reads])

        updated_reads = list(reads)
        for read in reads:
            if read.get_tag('AS') < top_as:
                updated_reads.remove(read)

        for read in updated_reads:
            outbam.write(read)

    outbam.close()
