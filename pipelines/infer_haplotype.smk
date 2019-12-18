import os
import pysam
from Bio import SeqIO
from pathlib import Path

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

#inputs
patient = config['patient']
sample = config['sample']
nm = config['nm']
protocol = config['protocol']

# polytect global references
genes = config['genes']
fasta_files = {gene: str(PD / config['gene_prefix'] / "alts" / f"{gene}.fa") for gene in genes}
bams = {gene: str(PD / "results" / patient / "alignments" / sample / f"{sample}_{gene}_complete.bam") for gene in genes}
gene_cor_cutoffs = config['gene_cor_cutoffs']
likelihood_files = {gene: str(PD / config['gene_prefix'] / "matrices" / f"{gene}_likelihoods_{protocol}.csv") for gene in genes}

#Derived variables
allele_ids = {gene:[record.id for record in SeqIO.parse(fasta, 'fasta')] for gene, fasta in fasta_files.items()}

rule all:
    input:
        haplotypes = f"results/{patient}/{sample}_haplotype.csv"

rule count_reads:
    input:
        bamfile = lambda w: bams[w.gene]
    output:
        counts = temp(f"temp/{sample}/{sample}_{{gene}}_counts.txt")
    params:
        alleles = lambda w: allele_ids[w.gene]
    run:
        """
        Takes a BAM file that has been processed by bwa-postalt.js 
        and finds the number of read pairs that align with nm <= threshold
        to each allele in 'alleles'
        """
        region_counts = {allele:0 for allele in params.alleles}
        bam = pysam.AlignmentFile(input.bamfile, 'rb')
        pairs = {}
        n_pairs = 1
        cur_read = "none"

        #For each read, track all entries
        #Add NM score from each member of the pair to pairs
        #At the end, for each pair that aligned with total NM <= threshold, add the probability to the probability output
        for read in bam:
            if read.query_name == "none":
                cur_read = read.query_name
            if read.query_name != cur_read:
                for region, values in pairs.items():
                    #Values[0] is boolean for whether or not each member of pair mapped, so values[0] == 2 means both mapped
                    #Values[1] is NM scores for each member, so check that total sum is <= threshold
                    if sum(values[0]) == 2 and sum(values[1]) <= nm and region in region_counts:
                        region_counts[region] += 1
                pairs = {}
                n_pairs += 1
                cur_read = read.query_name
            if read.has_tag('NM'):
                region = bam.get_reference_name(read.tid)
                #SAM flags dont transfer properly after bwa-postalt so we cant use read flags to make sure reads aligned in pairs
                #Each entry in pairs is a list of lists, where the first list is boolean for whether or not
                #the primary or secondary member of a read pair mapped, and the second list is the NM scores
                #of those reads if mapped
                if region in pairs:
                    if read.is_read1:
                        pairs[region][0][0] = 1
                        pairs[region][1][0] = read.get_tag('NM')
                    else:
                        pairs[region][0][1] = 1
                        pairs[region][1][1] = read.get_tag('NM')
                else:
                    if read.is_read1:
                        pairs[region] = [[1, 0], [read.get_tag('NM'), 9999999]]
                    else:
                        pairs[region] = [[0, 1], [9999999, read.get_tag('NM')]]

        with open(output.counts, 'w') as of:
            of.write(f"allele,count\n")
            for allele, count in region_counts.items():
                of.write(f"{allele},{count}\n")

rule call_haplotype:
    input:
        likelihoods = lambda w: likelihood_files[w.gene],
        counts = rules.count_reads.output.counts
    output:
        haplotype = temp(f"temp/{sample}/{sample}_{{gene}}_haplotype.txt")
    params:
        cutoff = lambda w: gene_cor_cutoffs[w.gene]
    script:
        "../scripts/type_alleles.R"

rule consolidate_haplotypes:
    input:
        expand("temp/{sample}/{sample}_{gene}_haplotype.txt", sample = sample, gene = genes)
    output:
        haplotypes = f"results/{patient}/{sample}_haplotype.csv"
    run:
        with open(output.haplotypes, 'w') as of:
            of.write(f"{sample}")
            for g in genes:
                with open(f"temp/{sample}/{sample}_{g}_haplotype.txt", 'r') as infile:
                    for line in infile:
                        line = line.strip()
                        of.write(f",{line}")
        of.close()