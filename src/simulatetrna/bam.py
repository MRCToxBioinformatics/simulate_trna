import pysam
import pandas as pd

def iterate_reads(inreads):
    last_read_name = None
    reads = ()

    for read in inreads:

        read_name = read.qname

        if read_name != last_read_name:

            if last_read_name is not None:
                yield reads

            reads = set((read,))

        else:
            reads.add(read)

        last_read_name = read_name

def filter_sam(infile, outfile):
    inbam = pysam.Samfile(infile, 'r')
    outbam = pysam.AlignmentFile(outfile, "w", template=inbam)
    inreads = inbam.fetch()


    for reads in iterate_reads(inbam):

        top_as = max([read.get_tag('AS') for read in reads])

        updated_reads = list(reads)
        for read in reads:
            if read.get_tag('AS') < top_as:
                updated_reads.remove(read)

        for read in updated_reads:
            outbam.write(read)


    outbam.close()
