import sys
from Bio import SeqIO
from Bio import Seq

if __name__ == "__main__":
    fasta_filename = sys.argv[1]
    gtf_filename = sys.argv[2]
    outfile = sys.argv[3]

    seqs = {record.id:record for record in SeqIO.parse(fasta_filename, 'fasta')}
    spliced_seqs = {allele:"" for allele in seqs.keys()}

    with open(gtf_filename, 'r') as f:
        for line in f:
            line = line.strip().split()
            allele = line[0]
            feature = line[2]
            start = int(line[3]) - 1
            end = int(line[4]) - 1
            if feature == "CDS":
                spliced_seqs[allele] += str(seqs[allele].seq[start:(end+1)])
    
    proteins = {allele:str(Seq.Seq(spliced_seqs[allele]).translate()) for allele in spliced_seqs.keys()}

    with open(outfile, 'w') as of:
        for allele, protein in proteins.items():
            of.write(f">{allele}\n")
            i = 0
            while i < len(protein):
                protein_segment = protein[i:i+100]
                of.write(f"{protein_segment}\n")
                i += 100