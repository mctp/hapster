# Pipeline to tune the correlation cutoff parameters for a given experimental setup
# Takes a sample, and runs infer_haplotypes + call_germline mutations once for each
# correlation cutoff value from .70 to .99

import os
from pathlib import Path

cutoffs = [x/100 for x in range(70, 100, 1)]

configfilename = workflow.overwrite_configfiles[0]

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

#Inputs
patient = config['patient']
sample = config['sample']
aligned_file = config['aligned_file']
nm = config['nm']
cram_reference = config['cram_reference']

#Algorithm parameters
genes = config['genes']

rule all:
    input:
        hap_fa = expand(f"results/cutoffs/{patient}/" + "{cutoff}/{gene}_original.fa", cutoff = cutoffs, gene = genes),
        germline_vcf = expand(f"results/cutoffs/{patient}/" + "{cutoff}/{gene}_germline_filtered.vcf", cutoff = cutoffs, gene = genes),
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
        bam = expand(f"results/{patient}/alignments/{sample}/{sample}_" + "{gene}_complete.bam", gene = genes)
    threads: NCORES
    shell:
        """
        snakemake --unlock \
            --configfile {configfilename} \
            --config POLYTECT_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     aligned_file={aligned_file} \
                     cram_reference={cram_reference} \
            --snakefile {PD}/pipelines/extract_reads.smk \
            --directory {PD} \
            --cores {NCORES}
        snakemake \
            --configfile {configfilename} \
            --config POLYTECT_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     aligned_file={aligned_file} \
                     cram_reference={cram_reference} \
            --snakefile {PD}/pipelines/extract_reads.smk \
            --directory {PD} \
            --cores {NCORES}
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
            --config POLYTECT_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     nm={nm} \
            --snakefile {PD}/pipelines/infer_haplotype.smk \
            --directory {PD} \
            --cores {NCORES}
        snakemake \
            --configfile {configfilename} \
            --config POLYTECT_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     nm={nm} \
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
        hap_fa_cutoff = expand(f"results/cutoffs/{patient}/" + "{{cutoff}}/{gene}_original.fa", gene = genes),
        germline_vcf_cutoff = expand(f"results/cutoffs/{patient}/" + "{{cutoff}}/{gene}_germline_filtered.vcf", gene = genes)
    params:
        move_command_fa = "\n".join([f"mv {x} {y}" for x, y in zip(expand(f"results/{patient}/refs/{sample}_" + "{gene}_original.fa", gene = genes), expand(f"results/cutoffs/{patient}/" + "{{cutoff}}/{gene}_original.fa", gene = genes))]),
        move_command_vcf = "\n".join([f"mv {x} {y}" for x, y in zip(expand(f"results/{patient}/calls/{sample}_" + "{gene}_germline_filtered.vcf", gene = genes), expand(f"results/cutoffs/{patient}/" + "{{cutoff}}/{gene}_germline_filtered.vcf", gene = genes))])
    threads: NCORES
    shell:
        """
        snakemake \
            --configfile {configfilename} \
            --config POLYTECT_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     haplotype={input.haplotype} \
            --snakefile {PD}/pipelines/call_germline_mutations.smk \
            --directory {PD} \
            --cores {NCORES}
        snakemake --unlock \
            --configfile {configfilename} \
            --config POLYTECT_DIR={PD} \
                     NCORES={NCORES} \
                     patient={patient} \
                     sample={sample} \
                     haplotype={input.haplotype} \
            --snakefile {PD}/pipelines/call_germline_mutations.smk \
            --directory {PD} \
            --cores {NCORES}
        {params.move_command_fa}
        {params.move_command_vcf}
        """

rule count_mutations:
    input:
        germline_vcfs = expand(f"results/cutoffs/{patient}/" + "{cutoff}/{gene}_germline_filtered.vcf", cutoff = cutoffs, gene = genes)
    output:
        germline_counts = f"results/cutoffs/{patient}/germline_counts.csv"
    run:
        with open(output.germline_counts, 'w') as of:
            of.write(f"patient,cutoff,gene,mutations\n")
            for cutoff in cutoffs:
                for gene in genes:
                    with open(f"results/cutoffs/{patient}/{cutoff}/{gene}_germline_filtered.vcf", 'r') as f:
                        count = 0
                        for line in f:
                            if (not line.startswith('#')) and ("PASS" in line):
                                count += 1
                        of.write(f"{patient},{cutoff},{gene},{count}\n")