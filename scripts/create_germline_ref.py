import sys
import subprocess
import re

class vcf_entry:
    def __init__(self, line):
        line = line.split()
        self.chromosome = line[0]
        self.position = int(line[1])
        self.id = line[2]
        self.reference = line[3]
        self.alternate = line[4]
        self.quality = line[5]
        self.filter = line[6]
        self.info = {x[0]:x[1] for x in [x.split('=') for x in line[7].split(';')] if x != ["STR"]}
        self.format = line[8].split(':')
        self.normal = {k:v for k, v in zip(self.format, line[9].split(':'))}

    def __str__(self):
        info = ';'.join(['='.join(x) for x in self.info.items()])
        format = ':'.join(self.format)
        normal = ':'.join([x for x in self.normal.values()])
        return '\t'.join([self.chromosome, str(self.position), self.id, self.reference, self.alternate, self.quality, self.filter, info, format, normal])

class vcf:
    def __init__(self, filename):
        with open(filename, 'r') as f:
            self.meta_information = []
            self.data = []
            for line in f:
                line = line.strip()
                if line.startswith("##"):
                    self.meta_information.append(line[2:])
                elif line.startswith('#'):
                    self.header = line[1:].split()
                else:
                    self.data.append(vcf_entry(line))

    def __str__(self):
        meta_information = '\n'.join([f"##{x}" for x in self.meta_information])
        header = '#' + '\t'.join(self.header)
        data = '\n'.join([str(x) for x in self.data])

        return '\n'.join([meta_information, header, data])

class fasta:
    def __init__(self, filename):
        with open(filename, 'r') as f:
            self.seqs = {}
            first = True
            seq_id = ""
            seq = ""
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if not first:
                        self.seqs[seq_id] = seq
                        seq = ""
                    else:
                        first = False
                    seq_id = line[1:]
                else:
                    seq += line
            self.seqs[seq_id] = seq

if __name__ == "__main__":
    haplotype_ref_file = sys.argv[1]
    vcf_file = sys.argv[2]
    germline_ref_file = sys.argv[3]

    hap_ref = fasta(haplotype_ref_file)
    mutations = vcf(vcf_file)

    #Offset keeps track of how previous indels will affect the position in the new string
    offset = 0
    cur_pos = 0
    chrom = ""
    for mutation in mutations.data:
        if mutation.filter == "PASS":
            #Make sure offset and current pos is reset on new allele
            if chrom == "":
                chrom = mutation.chromosome
            if chrom != mutation.chromosome:
                chrom = mutation.chromosome
                offset = 0
                cur_pos = 0
            
            seq = hap_ref.seqs[mutation.chromosome]
            pos = mutation.position - 1
            ref = mutation.reference
            alt = mutation.alternate

            #If alt is not made up of only nucleotides, skip it
            if re.search("[^ACGT]", alt) != None:
                continue
            
            #If identified variant occurs at less than 70% frequency, skip it
            alt_depth = int(mutation.normal["AD"].split(',')[1])
            total_depth = int(mutation.normal["DP"])
            if total_depth == 0:
                continue
            if alt_depth/total_depth < .7:
                continue
            
            #If pos + offset is less than cur_pos, that means it was in a region that got deleted by a previous variant
            #So we skip it
            if pos + offset > cur_pos:
                #Replace reference base with germline call
                hap_ref.seqs[mutation.chromosome] = seq[:pos + offset] + alt + seq[pos + offset + len(ref):]
                cur_pos = pos + offset
                offset += len(alt) - len(ref)

    raw_filename = re.sub("germline", "raw_germline", germline_ref_file)
    with open(f"{raw_filename}", 'w') as of:
        for allele_id, seq in hap_ref.seqs.items():
            of.write(f">{allele_id}\n")
            of.write(f"{seq}\n")

    subprocess.call(["picard", "NormalizeFasta", f"I={raw_filename}", f"O={germline_ref_file}"])
    subprocess.call(["bwa", "index", germline_ref_file])
    subprocess.call(["samtools", "faidx", germline_ref_file])
    subprocess.call(["picard", "CreateSequenceDictionary", f"R={germline_ref_file}"])
    