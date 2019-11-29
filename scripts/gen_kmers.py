import sys
import re
from Bio import SeqIO

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[x] for x in seq[::-1]])

def gen_kmers(seqs, length, step):
    """
    Takes a list of seqs and constructs the set of all kmers of provided length
    at given step values
    """
    kmers = []

    for seq in seqs:
        new_kmers = [seq[x:x+length] for x in range(0, len(seq) - length + 1, step)]
        kmers += new_kmers
        kmers += [revcomp(kmer) for kmer in new_kmers]

    #Remove kmers that move into masked regions
    kmers = [kmer for kmer in kmers if not kmer.find('N') == '-1']

    return set(kmers)

if __name__ == "__main__":
    fasta = sys.argv[1]
    outfile = sys.argv[2]
    length = int(sys.argv[3])
    step = int(sys.argv[4])

    seqs = []
    for record in SeqIO.parse(fasta, 'fasta'):
        seqs.append(str(record.seq))

    kmers = gen_kmers(seqs, length, step)

    with open(outfile, 'w') as of:
        for kmer in kmers:
            of.write(f"{kmer}\n")
