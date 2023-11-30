import os
import pysam
import pandas as pd

def iterate_reads(inreads, allow_multimapping=True, remove_trx=False, remove_mt=False):
    '''
    Yield groups of reads with the same read name

    inreads: pysam.Samfile object

    allow_multimapping: Set True to include multimapped reads

    remove_trx: Set True to remove tRNAs with 'tRX' in name

    remove_mt: Set True to remove tRNAs with 'MTtRNA' in name
    '''

    last_read_name = None
    reads = ()

    for read in inreads:

        if read.is_unmapped:
            continue

        if remove_trx and "tRX" in read.reference_name:
            continue

        if remove_mt and 'MTtRNA' in read.reference_name:
            continue

        read_name = read.qname

        if read_name != last_read_name:

            if last_read_name is not None:
                if allow_multimapping or len(reads)==1:
                    yield reads

            reads = set((read,))

        else:
            reads.add(read)

        last_read_name = read_name

    if allow_multimapping or len(reads)==1:
        yield reads


def keep_random_alignment(infile, outfile):

    '''
    Take a samfile and retain a single random alignment for multimapped reads. Will also sort and index the output

    infile: filepath for output samfile

    outfile: filepath for output bamfile
    '''

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
    '''
    Take a samfile and for each read, retain only the alignments with an alignment score equal to the best alignment score for the read. 

    infile: filepath for output samfile

    outfile: filepath for output bamfile
    '''


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
