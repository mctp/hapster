# Pipeline to tune the correlation cutoff parameters for a given experimental setup
# Takes a sample, and runs infer_haplotypes + call_germline mutations once for each
# correlation cutoff value from .70 to .99

import os
from pathlib import Path

cutoffs = [x/100 for x in range(50, 100, 1)]

configfilename = workflow.overwrite_configfiles[0]

# hapster runtime
PD = Path(config['HAPSTER_DIR'])
NCORES = config['NCORES']

#Inputs
protocol = config['protocol']
patient = config['patient']
sample = config['sample']
aligned_file = config['aligned_file']
nm = config['nm']
cram_reference = config['cram_reference']
extraction_regions = config['extraction_regions']

#Algorithm parameters
genes = config['genes']

rule all:
    input:
        hap_fa = expand(f"results/cutoffs/{patient}/" + "{cutoff}" + f"/{sample}_original.fa", cutoff = cutoffs),
        germline_vcf = expand(f"results/cutoffs/{patient}/" + "{cutoff}" + f"/{sample}_germline_filtered.vcf", cutoff = cutoffs),
        germline_counts = f"results/cutoffs/{patient}/germline_counts.csv"
    shell:
        """
        rm -r temp/bambamtemps
        """

rule extract_reads:
    input:
        alignment = aligned_file
    output:
        bam_for_kmer_filter = f"results/{patient}/alignments/{sample}/{sample}_extracted_with_blacklist.bam",
        fq1 = expand(f"results/{patient}/seqs/{sample}/{sample}_" + "{gene}_1.fq", gene = genes),
        fq2 = expand(f"results/{patient}/seqs/{sample}/{sample}_" + "{gene}_2.fq", gene = genes),
        bam = expand(f"results/{patient}/alignments/{sample}/{sample}_" + "{gene}_complete.bam", gene = genes),
        fq1_consolidated = f"results/{patient}/seqs/{sample}_1.fq",
        fq2_consolidated = f"results/{patient}/seqs/{sample}_2.fq"
    threads: NCORES
    shell:
        """
        snakemake --unlock \
            --configfile {configfilename} \
            --config HAPSTER_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     aligned_file={aligned_file} \
                     cram_reference={cram_reference} \
                     extraction_regions={extraction_regions} \
            --snakefile {PD}/pipelines/extract_reads.smk \
            --directory {PD} \
            --cores {NCORES} --notemp
        snakemake \
            --configfile {configfilename} \
            --config HAPSTER_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     aligned_file={aligned_file} \
                     cram_reference={cram_reference} \
                     extraction_regions={extraction_regions} \
            --snakefile {PD}/pipelines/extract_reads.smk \
            --directory {PD} \
            --cores {NCORES} --notemp
        """

rule infer_haplotype:
    input:
        bams = expand(f"results/{patient}/alignments/{sample}/{sample}_" + "{gene}_complete.bam", gene = genes)
    output:
        haplotypes_cutoff = f"results/cutoffs/{patient}/{{cutoff}}/{sample}_haplotype.csv"
    params:
       haplotypes = f"results/{patient}/{sample}_haplotype.csv"
    threads: NCORES
    shell:
        """
        snakemake --unlock \
            --configfile {configfilename} \
            --config HAPSTER_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     nm={nm} \
                     protocol={protocol} \
                     gene_cor_cutoffs={wildcards.cutoff} \
            --snakefile {PD}/pipelines/infer_haplotype.smk \
            --directory {PD} \
            --cores {NCORES}
        snakemake \
            --configfile {configfilename} \
            --config HAPSTER_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     nm={nm} \
                     protocol={protocol} \
                     gene_cor_cutoffs={wildcards.cutoff} \
            --snakefile {PD}/pipelines/infer_haplotype.smk \
            --directory {PD} \
            --cores {NCORES}
        mv {params.haplotypes} {output.haplotypes_cutoff}
        """

rule call_germline_mutations:
    input:
        fq1 = expand(f"results/{patient}/seqs/{sample}/{sample}_" + "{gene}_1.fq", gene = genes),
        fq2 = expand(f"results/{patient}/seqs/{sample}/{sample}_" + "{gene}_2.fq", gene = genes),
        haplotype = f"results/cutoffs/{patient}/{{cutoff}}/{sample}_haplotype.csv"
    output:
        hap_fa_cutoff = f"results/cutoffs/{patient}/{{cutoff}}/{sample}_original.fa",
        germline_vcf_cutoff = f"results/cutoffs/{patient}/{{cutoff}}/{sample}_germline_filtered.vcf"
    params:
        move_fa = lambda w: f"mv results/{patient}/refs/{sample}_original.fa results/cutoffs/{patient}/{w.cutoff}/{sample}_original.fa",
        move_vcf = lambda w: f"mv results/{patient}/calls/{sample}_germline_filtered.vcf results/cutoffs/{patient}/{w.cutoff}/{sample}_germline_filtered.vcf"
    threads: NCORES
    shell:
        """
        snakemake \
            --configfile {configfilename} \
            --config HAPSTER_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     haplotype={input.haplotype} \
            --snakefile {PD}/pipelines/call_germline_mutations.smk \
            --directory {PD} \
            --cores {NCORES}
        snakemake --unlock \
            --configfile {configfilename} \
            --config HAPSTER_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     haplotype={input.haplotype} \
            --snakefile {PD}/pipelines/call_germline_mutations.smk \
            --directory {PD} \
            --cores {NCORES}
        {params.move_fa}
        {params.move_vcf}
        """

rule count_mutations:
    input:
        germline_vcfs = expand(f"results/cutoffs/{patient}/" + "{cutoff}" + f"/{sample}_germline_filtered.vcf", cutoff = cutoffs)
    output:
        germline_counts = f"results/cutoffs/{patient}/germline_counts.csv"
    run:
        with open(output.germline_counts, 'w') as of:
            of.write(f"patient,cutoff,gene,mutations\n")
            for cutoff in cutoffs:
                with open(f"results/cutoffs/{patient}/{cutoff}/{sample}_germline_filtered.vcf", 'r') as f:
                    gene_counts = {gene:0 for gene in genes}
                    for line in f:
                        if (not line.startswith('#')) and ("PASS" in line):
                            line = line.strip().split('*')
                            gene = line[0]
                            gene_counts[gene] += 1
                    for gene in gene_counts.keys():
                        count = gene_counts[gene]
                        of.write(f"{patient},{cutoff},{gene},{count}\n")
