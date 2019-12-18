import sys
from Bio import pairwise2 #Biopython

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[x] for x in seq[::-1]])

def inserts_to_pairs(in_fq, out_fq1, out_fq2, length = 125):
    """
    Takes a fastq containing simulated inserts, and converts them to paired reads based on the requested
    length

    Input:
        in_fq (str): Filename of a fastq file containing a list of reads.
        out_fq1 (str): Filename of output fastq that will hold the first read in each pair
        out_fq2 (str): Filename of output fastq that will hold the second read in each pair
        length (int): Paired read length. Must be shorter or equal to the shortest read in the input fastq.

    Output:
        Writes paired reads to output fastq files.
    """
    fq1 = open(out_fq1, 'w')
    fq2 = open(out_fq2, 'w')

    i = 0
    with open(in_fq, 'r') as fq:
        entry_id = ""
        seq = ""
        optional = ""
        quality = ""
        for line in fq:
            if i % 100 == 0:
                print(i)
            i += 1
            entry_id = line.strip().split()[0]
            seq = next(fq).strip()
            optional = next(fq).strip()
            quality = next(fq).strip()
            seq1 = seq[:length]
            seq2 = revcomp(seq[len(seq) - length:])
            qual1 = quality[:length]
            qual2 = quality[len(quality) - length: ][::-1]
            fq1.write(f"{entry_id}\n")
            fq1.write(f"{seq1}\n")
            fq1.write(f"{optional}\n")
            fq1.write(f"{qual1}\n")
            fq2.write(f"{entry_id}\n")
            fq2.write(f"{seq2}\n")
            fq2.write(f"{optional}\n")
            fq2.write(f"{qual2}\n")
            
    fq1.close()
    fq2.close()

if __name__ == "__main__":
    in_fq = sys.argv[1]
    out_fq1 = sys.argv[2]
    out_fq2 = sys.argv[3]
    length = int(sys.argv[4])

    inserts_to_pairs(in_fq, out_fq1, out_fq2, length)