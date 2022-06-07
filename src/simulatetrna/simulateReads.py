import numpy as np

from numpy.random import choice

import gzip

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from scipy.stats import norm
from scipy.stats import binom

import simulatetrna.alignmentSummary as alignmentSummary
import pysam
from random import choices

from importlib import reload

from collections import defaultdict, Counter

def make_gt(fasta_infile, mu=100, sd=50):
    '''make a random ground truth for each feature in a fasta_infile
    Values below zero are replaced with zero
    '''
    records = list(SeqIO.parse(fasta_infile, "fasta"))

    counts = norm.rvs(mu, sd, size=len(records))
    counts[counts<0]=0

    return(dict(zip([record.name for record in records],
                    [round(c) for c in counts])))

def addSeqError(base, size=None):
    'Returns a random error for a given (DNA) base'

    if not base in ['a', 't', 'c', 'g',
                    'A', 'T', 'C', 'G']:
        raise ValueError('base (%s) should be a, t, c, g, A, T, C or G' % base)

    base2error = {'a':['t', 'c', 'g'],
                  't':['a', 'c', 'g'],
                  'c':['a', 't', 'g'],
                  'g':['a', 't', 'c'],
                  'A':['T', 'C', 'G'],
                  'T':['A', 'C', 'G'],
                  'C':['A', 'T', 'G'],
                  'G':['A', 'T', 'C']}

    error = choice(base2error[base], size=size)

    return(error)

addErrors = np.frompyfunc(addSeqError, 1, 1)

class simulated_read_sequences:
    '''Base class for simulated read sequences


    Attributes:
        seq: tRNA sequence
        seq_name: tRNA sequence name
        reads: numpy array, where rows=sequences
        read_names: names for simulated reads
        start: read start position relative to sequence
        end: read end position relative to sequence
        alignment_summary: alignmentSummary.clustalwtrnaAlignmentSummary object
        summary_level: alignment summary at anticodon or sequence level

    '''

    def __init__(self,
                 seq_name,
                 seq,
                 count,
                 alignment_summary=None,
                 summary_level=None):
        self.seq = seq
        self.seq_name = seq_name
        self.read_names = ['%s_%s' % (seq_name, ix) for ix in range(0, count)]
        self.reads = np.broadcast_to(np.array(seq), (count, len(np.array(seq))))
        self.start = (0,)*count
        self.end = (len(seq),)*count

        # so that we can assign new values to sim_sequences
        self.reads.setflags(write=1)

        if alignment_summary is not None:
            if not isinstance(alignment_summary, alignmentSummary.clustalwtrnaAlignmentSummary):
                raise ValueError('alignment_summary must be an alignmentSummary.clustalwtrnaAlignmentSummary object')
            if summary_level not in ['anticodon', 'sequence']:
                raise ValueError('summary_level must be anticodon or sequence')

        self.alignment_summary = alignment_summary
        self.summary_level = summary_level


    def getFeatureName(self):
        'When working at anticodon level, the feature is the anticodon, otherwise, seq_name'
        if self.summary_level == 'anticodon':
            feature_name = self.alignment_summary.trna_to_anticodon[self.seq_name]
        else:
            feature_name = self.seq_name

        return(feature_name)

    def addSeqError(self, error_rate):
        'Add random errors to simulate sequencing error'
        sim_errors = np.random.binomial(1, error_rate, len(self.reads.flatten()))
        sim_errors.shape = self.reads.shape

        mask = sim_errors > 0
        self.reads[mask] = addErrors(self.reads[mask])


    def addMutations(self, mutation_threshold):
        'Add mutations at frequencies from alignment_summary'

        assert self.alignment_summary is not None, (
            'Need to use the alignment_summary when initiating object if you want to add mutations')

        feature_name = self.getFeatureName()

        # if feature_name is not in contig_base_frequencies, no bases will be mutated
        if feature_name in self.alignment_summary.contig_base_frequencies.keys():

            mut_sig = defaultdict(int, self.alignment_summary.contig_base_frequencies[feature_name])

            # if at anticodon level, need to convert to tRNA seq level
            if self.summary_level == 'anticodon':
                mut_sig = {x : mut_sig[self.alignment_summary.trna_pos_to_alignment_pos[self.seq_name][x]]
                        for x in self.alignment_summary.trna_pos_to_alignment_pos[self.seq_name].keys()}

            for ix, ref in enumerate(self.seq):

                # If alignment_summary has been pickled, mut_sig[ix] is a dict, otherwise,
                # it's a Counter and will return values regardless if key present.
                # try-except is to catch no key for dict, if-elif-else is to
                # catch no key for Counter
                try:
                    base_mut_sig = mut_sig[ix][ref]

                    if base_mut_sig[ref]>0:
                        mutation_rate = 1 - (base_mut_sig[ref] / sum(base_mut_sig.values()))

                    elif sum(base_mut_sig.values()) == 0: # No data for position. For now, assert mutation rate is zero!
                        mutation_rate = 0

                    else: # Reference never seen, so mutation rate = 1
                        mutation_rate = 1

                except: # No data for position. For now, assert mutation rate is zero!
                    mutation_rate = 0


                if mutation_rate >= mutation_threshold:

                    mutated_bases = choices(list(base_mut_sig.keys()),
                                            weights=base_mut_sig.values(),
                                            k=self.reads.shape[0])

                    # Need to work out why line below is needed!!
                    self.reads = self.reads.copy()

                    self.reads[:,ix] = np.array(mutated_bases)


    def truncate(self):
        '''Truncate reads using frequencies in alignment_summary.

        Truncation only changes the self.start and self.end attributes, not the self.reads array'''

        assert self.alignment_summary is not None, (
            'Need to use the alignment_summary when initiating object if you want to truncate reads')

        feature_name = self.getFeatureName()

        feature_truncation_sig = self.alignment_summary.alignment_coordinates[feature_name]

        # if at anticodon level, need to convert to tRNA seq level
        if self.summary_level == 'anticodon':
            final_feature_truncation_sig = defaultdict(Counter)

            for end_pos, start_pos_counter in feature_truncation_sig.items():

                for start_pos, count in start_pos_counter.items():


                    updated_end_pos = self.alignment_summary.alignment_pos_to_trna_pos[self.seq_name][end_pos]
                    updated_start_pos = self.alignment_summary.alignment_pos_to_trna_pos[self.seq_name][start_pos]

                    if np.isnan(updated_end_pos) or np.isnan(updated_start_pos):
                        # alignment position does not exist in tRNA sequence
                        # (e.g gap relative to multiple alignment)
                        continue

                    final_feature_truncation_sig[updated_end_pos][updated_start_pos] = count

        else:
            final_feature_truncation_sig = feature_truncation_sig

        end_weights = {x: sum(final_feature_truncation_sig[x].values())
                       for x in final_feature_truncation_sig.keys()}

        # Otherwise, no truncation information for tRNA. For now, assert no truncations!!
        if sum(end_weights.values()) > 0:

            self.end = choices(list(end_weights.keys()),
                               weights=end_weights.values(),
                               k=self.reads.shape[0])

            self.start = [choices(list(final_feature_truncation_sig[end].keys()),
                                  weights=final_feature_truncation_sig[end].values())[0]
                          for end in self.end]


    def __str__(self):
        obj_summary = 'This is a summary of the object attributes\n'
        for attr, value in vars(self).items():
            if attr == 'read_names':
                continue
            obj_summary += ''.join(map(str, (attr, ': ', value, '\n')))
        return(obj_summary)

    def write(self, filehandle):
        for ix in range(0, self.reads.shape[0]):
            start = self.start[ix]
            end = self.end[ix]
            sequence = ''.join(self.reads[ix,start:end])

            sim_record = SeqRecord(
                    Seq(sequence),
                    letter_annotations={"phred_quality":[40] * len(sequence)},
                    id=self.read_names[ix],
                    name=self.read_names[ix],
                    description=self.read_names[ix])

            SeqIO.write(sim_record, filehandle, "fastq")


def simulate_reads(infile,
                   outfile,
                   ground_truth,
                   error_rate=1/100,
                   mutation_threshold=None,
                   truncate=False,
                   alignment_summary=None,
                   summary_level='anticodon'):

    records = list(SeqIO.parse(infile, "fasta"))

    # Sanity checks for correct imput when using alignment summary for mutations/truncation signatures
    if mutation_threshold is not None or truncate:
        assert alignment_summary is not None, 'Need to supply alignment_summary'

    if(outfile.endswith('.gz')):
        fastq_out = gzip.open(outfile, "wt")
    else:
        fastq_out = open(outfile, "w")

    for feature in records:

        feature_count = ground_truth[feature.name]

        sim_sequences = simulated_read_sequences(
            seq_name=feature.name,
            seq=feature.seq,
            count=feature_count,
            summary_level=summary_level,
            alignment_summary=alignment_summary)

        if error_rate > 0:
             sim_sequences.addSeqError(error_rate)

        if mutation_threshold is not None:
            sim_sequences.addMutations(mutation_threshold)

        if truncate:
            sim_sequences.truncate()

        sim_sequences.write(fastq_out)

    fastq_out.close()
