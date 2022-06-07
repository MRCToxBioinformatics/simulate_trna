"""
Summarising the tRNA alignments with respect to mutational signatures and truncations
"""

import pysam

from pysamstats import stat_variation

from numpy import nan

from tempfile import NamedTemporaryFile

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

from os import unlink

from collections import Counter
from collections import defaultdict

from re import sub

class trnaAlignmentSummary:
    '''Base class for alignment summary


    Attributes:
        trna_pos_to_alignment_pos: Map from trna seq position to aligned sequence position
        alignment_pos_to_trna_pos: Map from aligned sequence position to trna seq position
        alignment_coordinates: Counter for alignment coordinates (e.g to find sites of truncations)
        contig_base_frequencies: Counter for base frequencies
        alignments: pysam.AlignmentFile
        anticodons: Anticodons
        trna_to_anticodon: Map from tRNA name to anticodon group
        anticodon_to_trnas: Map from anticodon group to tRNAs
    '''

    def __init__(self):

        self.trna_pos_to_alignment_pos =  defaultdict(defaultdict)
        self.alignment_pos_to_trna_pos =  defaultdict(defaultdict)

        self.alignment_coordinates = defaultdict(lambda: defaultdict(Counter))

        self.contig_base_frequencies = defaultdict(
            lambda : defaultdict(
                lambda : defaultdict(
                    lambda: Counter())))

        self.alignments = None
        self.anticodons = None
        self.trna_records = None
        self.trna_to_anticodon = {}
        self.anticodon_to_trnas = defaultdict(set)

    def loadSam(self, samfile_path):
        self.alignments = pysam.Samfile(samfile_path, 'r')

    def loadTrnaFasta(self, trna_seq_fasta_filepath):
        self.trna_records = SeqIO.parse(trna_seq_fasta_filepath, "fasta")
        self.trna_records = {x.name:x for x in self.trna_records}

    def getAntiCodons(self):
        'Assumes naming convention. May need to update'
        if self.trna_records is None:
            raise ValueError('Need to first parse the tRNA fasta with: loadTrnaFasta')
        self.anticodons = set([sub('-\d+-\d+', '', x) for x in self.trna_records.keys()])
        for anticodon in self.anticodons:
            for trna in self.trna_records.keys():
                if(sub('-\d+-\d+', '', trna) == anticodon):
                    self.trna_to_anticodon[trna] = anticodon
                    self.anticodon_to_trnas[anticodon].add(trna)


    def getErrorProfile(self, trna_name):
        if self.contig_base_frequencies is None:
            raise ValueError('Need to first summarise the alignments with: build')

        if trna_name not in self.contig_base_frequencies.keys():
            raise ValueError('%s is not in the fasta file used to summarise the alignments' % trna_name)
        else:
            return(self.contig_base_frequencies[trna_name])

    def getAlignmentCoordinates(self, trna_name):
        if self.alignment_coordinates is None:
            raise ValueError('Need to first summarise the alignments with: build')

        if trna_name not in self.alignment_coordinates.keys():
            raise ValueError('%s is not in the fasta file used to summarise the alignments' % trna_name)
        else:
            return(self.alignment_coordinates[trna_name])


class clustalwtrnaAlignmentSummary(trnaAlignmentSummary):
    '''
    trnaAlignmentSummary using clustalw
    '''
    def build(self, samfile_path, trna_seq_fasta_filepath):

        self.loadSam(samfile_path)

        self.loadTrnaFasta(trna_seq_fasta_filepath)

        self.getAntiCodons()


        for anticodon in self.anticodons:
            tRNAs = self.anticodon_to_trnas[anticodon]
            #[y for (x,y) in self.trna_records.items() if anticodon in x]

            # Some anticodons have a single tRNA sequence, in which case, no need to cluster!
            if len(tRNAs) > 1:
                tmp_file = NamedTemporaryFile(delete=False)

                with open(tmp_file.name, 'w') as outf:
                    for tRNA in tRNAs:
                        SeqIO.write(self.trna_records[tRNA], outf, "fasta")

                tmp_file.close()
                clustaw_cline = ClustalwCommandline(infile=tmp_file.name, outfile=tmp_file.name, quiet=True)
                clustaw_cline()

                clustaw_alignments = AlignIO.read(tmp_file.name, format='clustal')

                name2alignment = {}

                for trna_seq_aligned in clustaw_alignments:
                    name2alignment[trna_seq_aligned.description] = trna_seq_aligned.seq
            else:
                tRNA = tRNAs.pop()
                name2alignment = {tRNA:self.trna_records[tRNA]}

            # Tom: All the lines below could be generic. Consider extracting to base class level methods
            # if other instances of the class are tested, e.g clustal Omega.
            # Get map from tRNA sequence position to multiple alignment position
            for contig in name2alignment.keys():
                trna_ix = 0
                for pos_ix, pos in enumerate(name2alignment[contig]):
                    if(pos!='-'):
                        self.trna_pos_to_alignment_pos[contig][trna_ix] = pos_ix
                        self.alignment_pos_to_trna_pos[contig][pos_ix] = trna_ix
                        trna_ix += 1
                    else:
                        self.alignment_pos_to_trna_pos[contig][pos_ix] = nan

            # Summarise alignment coordinates
            for contig in name2alignment.keys():
                for read in self.alignments.fetch(contig):
                    # reference_end points to one past the last aligned residue (https://pysam.readthedocs.io/en/latest/api.html)
                    trna_pos_end = read.reference_end-1
                    trna_pos_start = read.reference_start

                    malignment_pos_end = self.trna_pos_to_alignment_pos[contig][trna_pos_end]
                    malignment_pos_start = self.trna_pos_to_alignment_pos[contig][trna_pos_start]


                    self.alignment_coordinates[contig][trna_pos_end][trna_pos_start] += 1

                    self.alignment_coordinates[anticodon][trna_pos_end][trna_pos_start] += 1

            for contig in name2alignment.keys():

                pysamstats_records = list(stat_variation(
                    self.alignments, trna_seq_fasta_filepath, chrom=contig))

                for pysamstat_record in pysamstats_records:

                    ref = pysamstat_record['ref']

                    trna_pos = pysamstat_record['pos']
                    malignment_pos = self.trna_pos_to_alignment_pos[contig][trna_pos]

                    self.contig_base_frequencies[contig][trna_pos][ref]['A'] += pysamstat_record['A']
                    self.contig_base_frequencies[contig][trna_pos][ref]['T'] += pysamstat_record['T']
                    self.contig_base_frequencies[contig][trna_pos][ref]['C'] += pysamstat_record['C']
                    self.contig_base_frequencies[contig][trna_pos][ref]['G'] += pysamstat_record['G']

                    self.contig_base_frequencies[anticodon][malignment_pos][ref]['A'] += pysamstat_record['A']
                    self.contig_base_frequencies[anticodon][malignment_pos][ref]['T'] += pysamstat_record['T']
                    self.contig_base_frequencies[anticodon][malignment_pos][ref]['C'] += pysamstat_record['C']
                    self.contig_base_frequencies[anticodon][malignment_pos][ref]['G'] += pysamstat_record['G']
