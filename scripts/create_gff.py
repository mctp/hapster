import sys
import re

def create_gff(seqs_file, output_file, offset):
    alleles = []
    with open(seqs_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            line[1] = re.sub("[*.]", '', line[1]).split('|')
            alleles.append(line)

    of = open(output_file, 'w')
    of.write('##gff-version 3\n')
    #Creating GFF file for bcftools csq
    # The program looks for "CDS", "exon", "three_prime_UTR" and "five_prime_UTR" lines,
    # looks up their parent transcript (determined from the "Parent=transcript:" attribute),
    # the gene (determined from the transcript's "Parent=gene:" attribute), and the biotype
    # (the most interesting is "protein_coding").
    #
    # Attributes required for
    #   gene lines:
    #   - ID=gene:<gene_id>
    #   - biotype=<biotype>
    #   - Name=<gene_name>      [optional]
    #
    #   transcript lines:
    #   - ID=transcript:<transcript_id>
    #   - Parent=gene:<gene_id>
    #   - biotype=<biotype>
    #
    #   other lines (CDS, exon, five_prime_UTR, three_prime_UTR):
    #   - Parent=transcript:<transcript_id>
    #
    # Supported biotypes:
    #   - see the function gff_parse_biotype() in bcftools/csq.c 
    for j, entry in enumerate(alleles, start = 100000 + offset):

        seq_id, seq = entry
        seqname = seq_id
        seq_id = re.sub('[*:]', '', seq_id)
        source = 'create_gff.py'
        score = '.'
        strand = '+'

        #Gene entry - our whole region is the gene
        of.write(f'{seqname}\t{source}\t{"gene"}\t{1}\t{sum(len(segment) for segment in seq)}\t{score}\t{strand}\t{"."}\tID=gene:ENSG{j:011d};biotype=protein_coding;Name={seq_id}\n')
        
        #Transcript entry - our whole gene is the transcript
        of.write(f'{seqname}\t{source}\t{"transcript"}\t{1}\t{sum(len(segment) for segment in seq)}\t{score}\t{strand}\t{"."}\tID=transcript:ENST{j:011d};Parent=gene:ENSG{j:011d};biotype=protein_coding\n')
        

        #Generate entries for UTRs, CDSs, exons
        total_length = 0
        coding_length = 0
        for i, segment in enumerate(seq):
            #First segment (if it exists) is 5' UTR
            if i == 0 and len(segment) > 0:
                of.write(f'{seqname}\t{source}\t{"five_prime_UTR"}\t{1}\t{len(segment)}\t{score}\t{strand}\t{"."}\tParent=transcript:ENST{j:011d}\n')
            
            #Last segment (if it exists) is 3' UTR
            if i == len(seq) - 1 and len(segment) > 0:
                of.write(f'{seqname}\t{source}\t{"three_prime_UTR"}\t{total_length + 1}\t{total_length + len(segment)}\t{score}\t{strand}\t{"."}\tParent=transcript:ENST{j:011d}\n')

            #Our seqs start with a leading non-exon region, so 
            #exons are odd numbered segments
            #considering all exons to be full coding regions (CDS)
            if i % 2 == 1:
                start = total_length + 1
                end = total_length + len(segment)

                coding_mod = coding_length % 3
                phase = 0
                if coding_mod == 1:
                    phase = 2
                elif coding_mod == 2:
                    phase = 1

                of.write(f'{seqname}\t{source}\t{"exon"}\t{start}\t{end}\t{score}\t{strand}\t{"."}\tParent=transcript:ENST{j:011d}\n')
                of.write(f'{seqname}\t{source}\t{"CDS"}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tParent=transcript:ENST{j:011d}\n')
                coding_length += len(segment)
            
            total_length += len(segment)

    return alleles

def validate_alleles(input_file):
    alleles = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            line[1] = re.sub("[*.]", '', line[1]).split('|')
            alleles.append(line)
    
    for id, seq in alleles:
        total_length = 0
        coding_region = ""
        for i, segment in enumerate(seq):
            #Our seqs start with a leading non-exon region, so 
            #exons are odd numbered segments
            if i % 2 == 1:
                total_length += len(segment)
                coding_region += segment
        if not coding_region.startswith('ATG'):
            print(id)
            print(total_length)
            print(coding_region)
        if not coding_region.endswith('TGA'):
            print(total_length)
            print(coding_region)
        if not total_length % 3 == 0:
            print(total_length)
            print(coding_region)

    return alleles

if __name__ == "__main__":
    seqs_file = sys.argv[1]
    output_file = sys.argv[2]
    offset = int(sys.argv[3])

    create_gff(seqs_file, output_file, offset)