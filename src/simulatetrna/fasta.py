from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import collections

def santitiseGtRNAdbFasta(infile, outfile):
    records = list(SeqIO.parse(infile, "fasta"))

    recode_seq_to_name = collections.defaultdict(list)

    for record in records:
        recode_seq_to_name[record.seq].append(record.name)

    with open(outfile, 'w') as outf:
        for trna_seq, names in recode_seq_to_name.items():
            name = names[0] # take first name for all tRNA with same seqeunce
            santitised_tRNA = SeqRecord(
                trna_seq.replace('U', 'T'),
                id=name,
                name=name,description=name)
            SeqIO.write(santitised_tRNA, outf, "fasta")
