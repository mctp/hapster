"""
Simulates reads to construct the A matrices for haplotyping
"""
import os
import re
import subprocess
import pysam
from Bio import SeqIO
from pathlib import Path

# hapster runtime
PD = Path(config['HAPSTER_DIR'])
NCORES = config['NCORES']

#inputs
protocol = config['protocol']
min_insert_length = config['min_insert_length']
max_insert_length = config['max_insert_length']
read_length = config['read_length']
n_reads = config['n_reads']
nm = config['nm']
capture_targets = config['capture_targets']
similarity = config['similarity']

# hapster global references
genes = config['genes']
fasta_files = {gene: str(PD / config['gene_prefix'] / "alts" / f"{gene}.fa") for gene in genes}
regions = str(PD / config['gene_prefix'] / "sim" / "regions.txt")
base_fasta = str(PD / config['gene_prefix'] / "sim" / "base.fa")
complete_reference = str(PD / config['gene_prefix'] / "fa" / "complete.fa")
complete_reference_alt = str(PD / config['gene_prefix'] / "fa" / "complete.fa.alt")
extraction_reference = str(PD / config['gene_prefix'] / "fa" / "extraction.fa")
extraction_reference_alt = str(PD / config['gene_prefix'] / "fa" / "extraction.fa.alt")

# derived variables
full_regions = ""
with open(regions, 'r') as f:
    for line in f:
        line = line.strip().split()
        gene = line[0].replace(':', '')
        region = ' '.join(line[1:])
        full_regions += region + ' '

allele_ids = {gene:[record.id for record in SeqIO.parse(fasta, 'fasta')] for gene, fasta in fasta_files.items()}
allele_filenames = {k:[f"{re.sub('[^0-9a-zA-Z]+', '_', x)}" for x in v] for k, v in allele_ids.items()}
faidx_id = {k:{f"{re.sub('[^0-9a-zA-Z]+', '_', x)}":x for x in v} for k, v in allele_ids.items()}

rule all:
    input:
        #expand("temp/sim/{allele}_counts.txt", allele = allele_filenames),
        matrices = [PD / config["gene_prefix"] / "matrices" / f"{gene}_likelihoods_{protocol}.csv" for gene in genes]

#Makes a fasta file for use for simming reads
#This is where you can add pseudogenes to your real genes if you want
#to include them in the simulation
rule make_fasta:
    input:
        b_fasta = base_fasta,
        alt_fasta = PD / config["gene_prefix"] / "alts" / f"{{gene}}.fa"
    output:
        sim_fasta = temp(f"temp/sim/{{gene}}/{{allele}}_sim.fa")
    params:
        faidx_id = lambda w: faidx_id[w.gene][w.allele]
    shell:
        """
        cp {input.b_fasta} {output.sim_fasta}
        samtools faidx {input.alt_fasta} {params.faidx_id} >> {output.sim_fasta}
        """

#This rule simulates inserts with size based on the experimental
#design of the experiment it is simulating
#bbmap writes to local files that interfere with each other if run
#in parallel, so I move each instance to its own dir
rule sim_reads:
    input:
        sim_fasta = rules.make_fasta.output.sim_fasta
    output:
        insert_fastq = temp(f"temp/sim/{{gene}}/{{allele}}_inserts.fq")
    threads: 1
    shell:
        """
        mkdir {wildcards.allele}
        cd {wildcards.allele}
        randomreads.sh \
            ref=../{input.sim_fasta} \
            reads={n_reads} \
            minlength={min_insert_length} \
            maxlength={max_insert_length} \
            out=../{output.insert_fastq}
        cd ../
        rm -r {wildcards.allele}
        """

rule sim_exome_capture:
    input:
        capture_targets = capture_targets,
        insert_fastq = rules.sim_reads.output.insert_fastq
    output:
        exome_fq1 = temp(f"temp/sim/{{gene}}/{{allele}}_reads_wes.fq1"),
        exome_fq2 = temp(f"temp/sim/{{gene}}/{{allele}}_reads_wes.fq2")
    shell:
        """
        python scripts/targeted_enrichment.py \
            {input.capture_targets} \
            {input.insert_fastq} \
            {output.exome_fq1} \
            {output.exome_fq2} \
            {similarity} \
            {read_length}
        """

rule inserts_to_reads:
    input:
        insert_fastq = rules.sim_reads.output.insert_fastq
    output:
        wgs_fq1 = temp(f"temp/sim/{{gene}}/{{allele}}_reads_wgs.fq1"),
        wgs_fq2 = temp(f"temp/sim/{{gene}}/{{allele}}_reads_wgs.fq2")
    shell:
        """
        python scripts/inserts_to_reads.py \
            {input.insert_fastq} \
            {output.wgs_fq1} \
            {output.wgs_fq2} \
            {read_length}
        """

rule remove_ignored:
    input:
        fq1 = f"temp/sim/{{gene}}/{{allele}}_reads_{protocol}.fq1",
        fq2 = f"temp/sim/{{gene}}/{{allele}}_reads_{protocol}.fq2",
        extraction_ref = extraction_reference,
        full_ref = complete_reference
    output:
        extracted_sorted_bam = temp(f"temp/sim/{{gene}}/{{allele}}_extracted_sorted.bam"),
        extracted_sorted_bai = temp(f"temp/sim/{{gene}}/{{allele}}_extracted_sorted.bam.bai"),
        no_ignored_bam = temp(f"temp/sim/{{gene}}/{{allele}}_extracted_no_ignored.bam"),
        out_bam = temp(f"temp/sim/{{gene}}/{{allele}}_postalt_no_ignored.bam"),
        fq1 = temp(f"temp/sim/{{gene}}/{{allele}}_1_2.fq"),
        fq2 = temp(f"temp/sim/{{gene}}/{{allele}}_2_2.fq"),
        singles = temp(f"temp/sim/{{gene}}/{{allele}}_singles_2.fq"),
        orphan1 = temp(f"temp/sim/{{gene}}/{{allele}}_orphans_1_2.fq"),
        orphan2 = temp(f"temp/sim/{{gene}}/{{allele}}_orphans_2_2.fq")
    shell:
        """
        bwa mem {input.extraction_ref} {input.fq1} {input.fq2} | \
            samtools sort > {output.extracted_sorted_bam}
        samtools index {output.extracted_sorted_bam}
        samtools view -hb {output.extracted_sorted_bam} {full_regions} > {output.no_ignored_bam}
        bamtofastq \
            collate=1 \
            filename={output.no_ignored_bam} \
            gz=0 \
            inputformat=bam \
            F={output.fq1} \
            F2={output.fq2} \
            S={output.singles} \
            O={output.orphan1} \
            O2={output.orphan2}
        bwa mem {input.full_ref} {output.fq1} {output.fq2} | \
            k8 ./bin/bwa-postalt.js {complete_reference_alt} | \
            samtools view -hb | \
            samtools sort -n -o {output.out_bam}
        """

rule count_reads:
    input:
        bamfile = rules.remove_ignored.output.out_bam
    output:
        counts = temp(f"temp/sim/{{gene}}/{{allele}}_counts.txt")
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
            of.write(f"{wildcards.allele}")
            for count in region_counts.values():
                of.write(f",{count/n_pairs:.5f}")

rule consolidate_counts:
    input:
        in_counts = lambda w: expand("temp/sim/{{gene}}/{allele}_counts.txt", allele = allele_filenames[w.gene])
    output:
        matrix_raw = temp(f"temp/matrices/{{gene}}_likelihoods.csv")
    params:
        allele_ids = lambda w: allele_ids[w.gene]
    run:
        with open(output.matrix_raw, 'w') as of:
            for name in params.allele_ids:
                of.write(f",{name}",)
            of.write('\n')
            for filename in input.in_counts:
                with open(filename, 'r') as f:
                    for line in f:
                        of.write(f"{line}\n")

rule make_symmetrical:
    input:
        matrix_raw = rules.consolidate_counts.output.matrix_raw
    output:
        matrix = PD / config["gene_prefix"] / "matrices" / f"{{gene}}_likelihoods_{protocol}.csv"
    script:
        "../scripts/make_symmetrical.R"
