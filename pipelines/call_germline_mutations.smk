import pysam
import re
import os
from pathlib import Path

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

# inputs
patient = config['patient']
sample = config['sample']
haplotype = config['haplotype']

# polytect global references
genes = config['genes']
fq1s = {gene: str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_1.fq") for gene in genes}
fq2s = {gene: str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_2.fq") for gene in genes}
fasta_files = {gene: str(PD / config['gene_prefix'] / "alts" / f"{gene}.fa") for gene in genes}

#Derived variables
haplotype_dict = {}
with open(haplotype, 'r') as f:
    alleles = next(f).strip().split(',')[1:]
    for allele in alleles:
        gene = allele.split('*')[0]
        if gene in haplotype_dict:
            haplotype_dict[gene].append(allele)
        else:
            haplotype_dict[gene] = [allele]

rule all:
    input:
        hap_fa = expand(f"results/{patient}/refs/{sample}_" + "{gene}_original.fa", gene = genes),
        germline_vcf = expand(f"results/{patient}/calls/{sample}_" + "{gene}_germline_filtered.vcf", gene = genes)

rule create_regions_list:
    input:
        haplotype=haplotype
    output:
        regions_list=temp(f'temp/{sample}/{sample}_regions.list')
    run:
        haplotype = []
        with open(input.haplotype, 'r') as f:
            line = next(f).strip().split(',')
            haplotype += line[1:]
        haplotype = set(haplotype)
        haplotype = '\n'.join(set(haplotype))
        of = open(output.regions_list, 'w')
        of.write(haplotype)

rule make_haplotype_ref:
    input:
        haplotype = haplotype
    output:
        hap_fa = temp(f"temp/{sample}/{sample}_{{gene}}.fa"),
        hap_fa_amb = temp(f"temp/{sample}/{sample}_{{gene}}.fa.amb"),
        hap_fa_ann = temp(f"temp/{sample}/{sample}_{{gene}}.fa.ann"),
        hap_fa_bwt = temp(f"temp/{sample}/{sample}_{{gene}}.fa.bwt"),
        hap_fa_pac = temp(f"temp/{sample}/{sample}_{{gene}}.fa.pac"),
        hap_fa_fai = temp(f"temp/{sample}/{sample}_{{gene}}.fa.fai"),
        hap_fa_sa = temp(f"temp/{sample}/{sample}_{{gene}}.fa.sa"),
        hap_dict = temp(f"temp/{sample}/{sample}_{{gene}}.dict"),
        hap_fa_saved = f"results/{patient}/refs/{sample}_{{gene}}_original.fa"
    params:
        allele_1 = lambda w: sorted(haplotype_dict[w.gene])[0],
        allele_2 = lambda w: sorted(haplotype_dict[w.gene])[1]
    shell:
        """
        samtools faidx {input.fa} {params.allele_1} >> {output.hap_fa}
        samtools faidx {input.fa} {params.allele_2} >> {output.hap_fa}
        bwa index {output.hap_fa}
        samtools faidx {output.hap_fa}
        picard CreateSequenceDictionary R={output.hap_fa}
        cp {output.hap_fa} {output.hap_fa_saved}
        """

rule realign_to_haplotype_ref:
    input:
        rules.make_haplotype_ref.output,
        hap_fa = rules.make_haplotype_ref.output.hap_fa,
        fq1 = lambda w: fq1s[w.gene],
        fq2 = lambda w: fq2s[w.gene]
    output:
        bam = temp(f"temp/{sample}/{sample}_{{gene}}_realigned.bam"),
        deduped_bam = temp(f"temp/{sample}/{sample}_{{gene}}_realigned_deduped.bam"),
        bam_bai = temp(f"temp/{sample}/{sample}_{{gene}}_realigned_deduped.bam.bai"),
        metrics = temp(f"temp/{sample}/{sample}_{{gene}}_realigned_deduped_metrics.txt")
    params:
        read_group = lambda w: f"@RG\\tSM:{sample}\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}"
    threads: NCORES
    shell:
        """
        bwa mem -t 4 {input.hap_fa} {input.fq1} {input.fq2} -R "{params.read_group}" | \
            samtools sort -@ 4 | \
            samtools view -@ 4 -hb > {output.bam}
        picard MarkDuplicates I={output.bam} O={output.deduped_bam} M={output.metrics}
        samtools index {output.deduped_bam}
        """

rule mapq_to_60:
    input:
        normal_bam=rules.realign_to_haplotype_ref.output.deduped_bam
    output:
        normal_bam_60=temp(f"temp/{sample}/{sample}_{{gene}}_realigned_deduped_60.bam")
    run:
        normal_bam = pysam.AlignmentFile(input.normal_bam, 'rb')
        normal_bam_60 = pysam.AlignmentFile(output.normal_bam_60, 'wb', template = normal_bam)
        for read in normal_bam:
            read.mapq = 60
            normal_bam_60.write(read)

rule index_60:
    input:
        normal_bam_60=rules.mapq_to_60.output.normal_bam_60
    output:
        normal_bam_60_bai=temp(rules.mapq_to_60.output.normal_bam_60 + ".bai")
    shell:
        """
        samtools index {input.normal_bam_60}
        """

rule call_germline_mutations:
    input:
        rules.make_haplotype_ref.output,
        hap_fa = rules.make_haplotype_ref.output.hap_fa,
        normal_bam=rules.mapq_to_60.output.normal_bam_60,
        normal_bam_bai=rules.index_60.output.normal_bam_60_bai
    output:
        out_vcf=temp(f'results/{patient}/calls/{sample}_{{gene}}_germline.vcf.gz')
    shell:
        """
        gatk HaplotypeCaller \
            --java-options "-Xmx8G" \
            -R {input.hap_fa} \
            -I {input.normal_bam} \
            -O {output.out_vcf}
        """

rule filter_germline_mutations:
    input:
        rules.make_haplotype_ref.output,
        hap_fa = rules.make_haplotype_ref.output.hap_fa,
        in_vcf = rules.call_germline_mutations.output.out_vcf
    output:
        out_vcf = f"results/{patient}/calls/{sample}_{{gene}}_germline_filtered.vcf"
    shell:
        """
        gatk VariantFiltration \
            -R {input.hap_fa} \
            -V {input.in_vcf} \
            -O {output.out_vcf} \
            -filter-name "my_filter" \
            -filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        """