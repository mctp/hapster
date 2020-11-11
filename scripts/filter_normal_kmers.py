import sys
import re
import subprocess
import copy

def revcomp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[x] for x in seq[::-1]])

class vcf_entry:
    def __init__(self, line, tumor_col, normal_col):
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
        #Tumor column won't exist if its a germline VCF
        if tumor_col != -1:
            self.tumor = {k:v for k, v in zip(self.format, line[tumor_col].split(':'))}
        self.normal = {k:v for k, v in zip(self.format, line[normal_col].split(':'))}

    def __str__(self):
        info = ';'.join(['='.join(x) for x in self.info.items()])
        format = ':'.join(self.format)
        normal = ':'.join([x for x in self.normal.values()])
        #Tumor column won't exist if its a germline VCF
        try:
            tumor = ':'.join([x for x in self.tumor.values()])
        except AttributeError:
            return '\t'.join([self.chromosome, str(self.position), self.id, self.reference, self.alternate, self.quality, self.filter, info, format, normal])
        return '\t'.join([self.chromosome, str(self.position), self.id, self.reference, self.alternate, self.quality, self.filter, info, format, tumor, normal])

class vcf:
    def __init__(self, filename):
        with open(filename, 'r') as f:
            self.meta_information = []
            self.data = []
            self.tumor_sample = ""
            self.tumor_col = -1
            self.normal_sample = ""
            self.normal_col = -1
            for line in f:
                line = line.strip()
                if line.startswith("##"):
                    self.meta_information.append(line[2:])
                    info = line[2:]
                    if info.startswith("normal_sample"):
                        self.normal_sample = info.split('=')[1]
                    if info.startswith("tumor_sample"):
                        self.tumor_sample = info.split('=')[1]
                elif line.startswith('#'):
                    info = line[1:].split()
                    if len(info) == 10:
                        self.normal_col = 9
                    elif len(info) == 11:
                        if info[9] == self.normal_sample:
                            self.normal_col = 9
                            self.tumor_col = 10
                        else:
                            self.normal_col = 10
                            self.tumor_col = 9
                    self.header = info
                else:
                    self.data.append(vcf_entry(line, self.tumor_col, self.normal_col))

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

def gen_kmers(chrom, pos, alt, ref, ref_seqs, k, mutate = True):
    kmers = []
    seq = ref_seqs.seqs[chrom]
    offset = 8

    #VCF pos is 1-indexed, so change to 0 indexed
    pos = pos-1
    print(pos)

    if mutate:
        #Generate mutated sequence
        #Takes into account indels
        mut_seq = ''.join([seq[:pos], alt, seq[pos + len(ref):]])
    else:
        mut_seq = seq
    
    #Gen kmers
    if pos - k + 1 >= 0:
        start = pos - k + len(ref) + offset
    else:
        start = 0
    
    if pos + k <= len(mut_seq):
        end = pos
    else:
        end = len(mut_seq) - k - offset

    kmers = [mut_seq[x:x+k] for x in range(start, end + 1)]
    kmer_revcomps = [revcomp(x) for x in kmers]
    
    return kmers + kmer_revcomps

if __name__ == "__main__":
    vcf_file = sys.argv[1]
    bam_file = sys.argv[2]
    ref_file = sys.argv[3]
    k = int(sys.argv[4])
    out_file = sys.argv[5]

    mutations = vcf(vcf_file)
    germ_ref = fasta(ref_file)

    #Keep track of just passing mutations, and their offset with respect to the somatic_ref
    passing_mutations = [[mutation, 0] for mutation in mutations.data if mutation.filter == "PASS"]

    #Make ref seq with all passing "somatic" mutations
    #Offset keeps track of how previous indels will affect the position in the new string
    offset = 0
    cur_pos = 0
    chrom = ""
    somatic_ref = copy.deepcopy(germ_ref)
    for i, mutation in enumerate(passing_mutations):
        #If alt is not made up of only nucleotides, skip it
        alt = mutation[0].alternate
        if re.search("[^ACGT]", alt) != None:
            continue

        #Make sure offset and current pos is reset on new allele
        if chrom == "":
            chrom = mutation[0].chromosome
        if chrom != mutation[0].chromosome:
            chrom = mutation[0].chromosome
            offset = 0
            cur_pos = 0
        
        #Update offset value with respect to somatic_ref
        passing_mutations[i][1] = offset

        seq = somatic_ref.seqs[mutation[0].chromosome]
        pos = mutation[0].position - 1
        ref = mutation[0].reference
        
        #If pos + offset is less than cur_pos, that means it was in a region that got deleted by a previous variant
        #So we skip it
        if pos + offset > cur_pos:
            #Replace reference base with germline call
            somatic_ref.seqs[mutation[0].chromosome] = seq[:pos + offset] + alt + seq[pos + offset + len(ref):]
            offset += len(alt) - len(ref)
            cur_pos = pos

    #Gen kmers and count occurences
    i = 0
    for mutation, offset in passing_mutations:
        alt = mutation.alternate
        if re.search("[^ACGT]", alt) != None:
            continue
        
        #Gen kmers for single mutation
        kmers = gen_kmers(mutation.chromosome, mutation.position, mutation.alternate, mutation.reference, germ_ref, k)

        #Gen kmers for mutation with respect to somatically mutated ref
        kmers += gen_kmers(mutation.chromosome, mutation.position + offset, mutation.alternate, mutation.reference, somatic_ref, k, mutate = False)
        
        #Remove duplicates
        #kmers = set(kmers)

        #Write kmers to temp file for grep
        temp_kmers = f'temp/temp_kmers_{mutations.header[9]}_{vcf_file[-20:-4]}_{i}.txt'
        with open(temp_kmers, 'w') as of:
            for kmer in kmers:
                of.write(f"{kmer}\n")
        i += 1

        #Check all kmers against sam file using grep
        #Grep errors out if it doesn't find a match,
        #so on failure that means no kmers were found in the normal
        view = subprocess.Popen(['samtools', 'view', bam_file], stdout = subprocess.PIPE)
        try:
            test = subprocess.check_output(['grep', '-f', temp_kmers], stdin = view.stdout)
        except:
            found = False
        else:
            found = True
            #Minus 1 for trailing newline
            print(mutation)
            test_p = test.decode("utf-8")
            print(f"{test_p}\n")
            n_found = len(test.decode("utf-8").split('\n')) - 1

        #Add kmer counts as part of the info field
        if found:
            mutation.info["KMER"] = str(n_found)
        else:
            mutation.info["KMER"] = '0'
    
    with open(out_file, 'w') as of:
        of.write(str(mutations))