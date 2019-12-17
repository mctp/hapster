import sys
from Bio import pairwise2 #Biopython

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[x] for x in seq[::-1]])

def enrich_targets(probes_file, in_fq, out_fq1, out_fq2, similarity = .82, length = 126):
    """
    Takes a fastq containing simulated inserts, and checks them for similarity
    against a probe to determine if they would have been pulled down during target
    enrichment. For fragments that pass, consider them to be inserts in paired sequencing
    and take paired reads from each end.

    Input:
        probes_file (str): Filename of file containing a list of probes for enrichment. Currently expects
                            files of the format provided by Agilent for their v4 probe set.
        in_fq (str): Filename of a fastq file containing a list of reads.
        out_fq1 (str): Filename of output fastq that will hold the first read in each pair
        out_fq2 (str): Filename of output fastq that will hold the second read in each pair
        similarity (double): Percent similarity of a probes sequence to a read to consider that read targeted.
                             Tolerance of up to 18% dissimilarity (so 82% similarity) assumed based on 
                             a guideline from Norman et.al., 2017, doi:10.1101/gr.213538.116.
        length (int): Paired read length. Must be shorter or equal to the shortest read in the input fastq.

    Output:
        Writes paired reads to output fastq files.
    """
    #Assumes format provided by Agilent v4 probe set
    probes = []
    with open(probes_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            probe_id = line[1]
            probe_seq = line[2]
            probes.append([probe_id, probe_seq])
    probe_length = len(probes[0][1])

    #Alignment scores are equal to 1 for every match, and -1 for every mismatch/indel
    #Thus an 82% similarity means 18% mismatches, and a score 36% below the max since
    #you both lose the +1 for the match and get a -1 penalty for the mismatch.
    #Max score is the length of the probe.
    threshold = probe_length - (int(probe_length * (1-similarity)) * 2)
    #threshold = 0

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

            #If probe alignment score exceeds threshold, consider it an insert and add seqs to our paired fastqs
            for probe in probes:
                #This alignment method might not be optimal
                #The intent is to create an alignment with integer scores so it is easy to count sequence similarity.
                #Should probably test to see if this leads to unacceptable false positives/false negatives
                if pairwise2.align.globalms(probe[1], seq, 1, -1, -1, -1, penalize_end_gaps = False, score_only = True) > threshold:
                    seq1 = seq[:length]
                    seq2 = revcomp(seq[len(seq) - length:])
                    qual1 = quality[:length]
                    qual2 = quality[len(quality) - length: ][::-1]
                    fq1.write(f"{entry_id + '_' + probe[0]}\n")
                    fq1.write(f"{seq1}\n")
                    fq1.write(f"{optional}\n")
                    fq1.write(f"{qual1}\n")
                    fq2.write(f"{entry_id + '_' + probe[0]}\n")
                    fq2.write(f"{seq2}\n")
                    fq2.write(f"{optional}\n")
                    fq2.write(f"{qual2}\n")
                    break
    
    fq1.close()
    fq2.close()

if __name__ == "__main__":
    probes_file = sys.argv[1]
    in_fq = sys.argv[2]
    out_fq1 = sys.argv[3]
    out_fq2 = sys.argv[4]
    similarity = float(sys.argv[5])
    length = int(sys.argv[6])

    enrich_targets(probes_file, in_fq, out_fq1, out_fq2, similarity, length)