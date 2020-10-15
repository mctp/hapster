import pysam
import re
import os
from pathlib import Path

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

# inputs
patient = config['patient']
normal = config['normal']
normal_bam = config['normal_bam']
tumor = config['tumor']
tumor_bam = config['tumor_bam']
germline_ref = config['germline_ref']
gff = config['gff']
original_normal_bam = config['normal_for_kmers']

rule all:
    input:
        mutect_filtered = f"results/{patient}/calls/{tumor}_{normal}_filtered.vcf",
        kmer_filtered = f"results/{patient}/calls/kmer_{tumor}_{normal}_filtered.vcf",
        csq_annotated = f"results/{patient}/calls/{tumor}_{normal}_annotated.txt"

rule germline_mapq_to_60:
    input:
        tumor_bam = tumor_bam,
        normal_bam = normal_bam
    output:
        tumor_bam_60=temp(f"temp/{tumor}/germline_{tumor}_realigned_deduped_60.bam"),
        normal_bam_60=temp(f"temp/{normal}/germline_{normal}_realigned_deduped_60.bam")
    run:
        tumor_bam = pysam.AlignmentFile(input.tumor_bam, 'rb')
        tumor_bam_60 = pysam.AlignmentFile(output.tumor_bam_60, 'wb', template = tumor_bam)
        for read in tumor_bam:
            read.mapq = 60
            tumor_bam_60.write(read)
        normal_bam = pysam.AlignmentFile(input.normal_bam, 'rb')
        normal_bam_60 = pysam.AlignmentFile(output.normal_bam_60, 'wb', template = normal_bam)
        for read in normal_bam:
            read.mapq = 60
            normal_bam_60.write(read)

rule germline_index_60:
    input:
        tumor_bam_60=rules.germline_mapq_to_60.output.tumor_bam_60,
        normal_bam_60=rules.germline_mapq_to_60.output.normal_bam_60
    output:
        tumor_bam_60_bai=temp(rules.germline_mapq_to_60.output.tumor_bam_60 + ".bai"),
        normal_bam_60_bai=temp(rules.germline_mapq_to_60.output.normal_bam_60 + ".bai")
    shell:
        """
        samtools index {input.tumor_bam_60}
        samtools index {input.normal_bam_60}
        """

rule run_mutect2:
    input:
        germ_fa = germline_ref,
        tumor_bam = rules.germline_mapq_to_60.output.tumor_bam_60,
        normal_bam = rules.germline_mapq_to_60.output.normal_bam_60,
        tumor_bam_bai = rules.germline_index_60.output.tumor_bam_60_bai,
        normal_bam_bai = rules.germline_index_60.output.normal_bam_60_bai
    output:
        out_vcf=temp(f'results/{patient}/calls/{tumor}_{normal}_raw.vcf.gz'),
        out_bam=temp(f'results/{patient}/calls/{tumor}_{normal}_raw.bam')
    shell:
        """
        /home/mumphrey/Projects/hla_pipeline/bin/gatk-4.1.2.0/gatk --java-options "-Xmx2g" Mutect2 \
            -R {input.germ_fa} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -tumor {tumor} \
            -normal {normal} \
            -O {output.out_vcf} \
            -bamout {output.out_bam} \
            --disable-read-filter MappingQualityAvailableReadFilter \
            --disable-read-filter MappingQualityNotZeroReadFilter \
            --disable-read-filter MappingQualityReadFilter
        """

rule filter_mutect2:
    input:
        germ_fa = germline_ref,
        in_vcf=rules.run_mutect2.output.out_vcf
    output:
        out_vcf = f"results/{patient}/calls/{tumor}_{normal}_filtered.vcf"
    params:
        out_vcf_gz = f"results/{patient}/calls/{tumor}_{normal}_filtered.vcf.gz"
    shell:
        """
        /home/mumphrey/Projects/hla_pipeline/bin/gatk-4.1.2.0/gatk --java-options "-Xmx2g" FilterMutectCalls \
            -V {input.in_vcf} \
            -O {params.out_vcf_gz} \
            -R {input.germ_fa} \
            --max-events-in-region 15
        gunzip {params.out_vcf_gz}
        """

#Filters out any variants that have kmer support in the normal sample
rule filter_normal_kmers:
    input:
        in_vcf = rules.filter_mutect2.output.out_vcf,
        normal_bam = normal_bam_for_kmers,
        germ_fa = germline_ref
    output:
        out_vcf = f"results/{patient}/calls/kmer_{tumor}_{normal}_filtered.vcf"
    shell:
        """
        python /home/mumphrey/Projects/hla_pipeline/scripts/filter_normal_kmers.py \
            {input.in_vcf} \
            {input.normal_bam} \
            {input.germ_fa} \
            25 \
            {output.out_vcf}
        """

rule run_csq:
    input:
        germ_fa = germline_refs,
        input_vcf=rules.filter_normal_kmers.output.out_vcf,
        input_gff= gff
    output:
        output_txt = f"results/{patient}/calls/{tumor}_{normal}_annotated.txt"
    shell:
    	"""
    	bcftools csq \
    	-f {input.germ_fa} \
    	-g {input.input_gff} \
    	{input.input_vcf} \
    	-Ot \
    	-o {output.output_txt} \
    	-l \
        -i 'FILTER="PASS"'
    	"""
