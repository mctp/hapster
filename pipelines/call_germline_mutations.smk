import pysam
import re
import os
from pathlib import Path

# hapster runtime
PD = Path(config['HAPSTER_DIR'])
NCORES = config['NCORES']

# inputs
patient = config['patient']
sample = config['sample']
haplotype = config['haplotype']

# hapster global references
genes = config['genes']
fq1s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_1.fq") for gene in genes]
fq2s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_2.fq") for gene in genes]
fasta_files = [str(PD / config['gene_prefix'] / "alts" / f"{gene}.fa") for gene in genes]

#Derived variables
with open(haplotype, 'r') as f:
    alleles = next(f).strip().split(',')[1:]
    haplotype_string = ' '.join(alleles)

rule all:
    input:
        hap_fa = f"results/{patient}/refs/{sample}_original.fa",
        deduped_bam = f"results/{patient}/alignments/{sample}/{sample}_haplotype_realigned.bam",
        bam_bai = f"results/{patient}/alignments/{sample}/{sample}_haplotype_realigned.bam.bai",
        germline_vcf = f"results/{patient}/calls/{sample}_germline_filtered.vcf"

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
        fa = [x for x in fasta_files],
        haplotype = haplotype
    output:
        temp_fa = temp(f"temp/{sample}/fa_full.fa"),
        hap_fa = temp(f"temp/{sample}/{sample}.fa"),
        hap_fa_amb = temp(f"temp/{sample}/{sample}.fa.amb"),
        hap_fa_ann = temp(f"temp/{sample}/{sample}.fa.ann"),
        hap_fa_bwt = temp(f"temp/{sample}/{sample}.fa.bwt"),
        hap_fa_pac = temp(f"temp/{sample}/{sample}.fa.pac"),
        hap_fa_fai = temp(f"temp/{sample}/{sample}.fa.fai"),
        hap_fa_sa = temp(f"temp/{sample}/{sample}.fa.sa"),
        hap_dict = temp(f"temp/{sample}/{sample}.dict"),
        hap_fa_saved = f"results/{patient}/refs/{sample}_original.fa"
    params:
        fasta_string = " ".join([x for x in fasta_files]),
        haplotype_string = haplotype_string
    shell:
        """
        for FASTA in "$(echo {params.fasta_string})"; do cat $FASTA >> {output.temp_fa}; done;
        for ALLELE in "$(echo {params.haplotype_string})"; do samtools faidx {output.temp_fa} $ALLELE >> {output.hap_fa}; done;
        bwa index {output.hap_fa}
        samtools faidx {output.hap_fa}
        picard CreateSequenceDictionary R={output.hap_fa}
        cp {output.hap_fa} {output.hap_fa_saved}
        """

rule realign_to_haplotype_ref:
    input:
        rules.make_haplotype_ref.output,
        hap_fa = rules.make_haplotype_ref.output.hap_fa,
        fq1 = [x for x in fq1s],
        fq2 = [x for x in fq2s]
    output:
        temp_fq1 = temp(f"temp/{sample}/{sample}_temp_1.fq"),
        temp_fq2 = temp(f"temp/{sample}/{sample}_temp_2.fq"),
        bam = temp(f"temp/{sample}/{sample}_realigned.bam"),
        deduped_bam = f"results/{patient}/alignments/{sample}/{sample}_haplotype_realigned.bam",
        bam_bai = f"results/{patient}/alignments/{sample}/{sample}_haplotype_realigned.bam.bai",
        metrics = temp(f"temp/{sample}/{sample}_realigned_deduped_metrics.txt")
    params:
        fq1 = " ".join([x for x in fq1s]),
        fq2 = " ".join([x for x in fq2s]),
        read_group = lambda w: f"@RG\\tSM:{sample}\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}"
    threads: NCORES
    shell:
        """
        for FQ in $(echo {params.fq1}); do if (( $(stat -c%s "$FQ") > 25 )); then cat $FQ >> {output.temp_fq1}; fi; done;
        for FQ in $(echo {params.fq2}); do if (( $(stat -c%s "$FQ") > 25 )); then cat $FQ >> {output.temp_fq2}; fi; done;
        bwa mem -t 4 {input.hap_fa} {output.temp_fq1} {output.temp_fq2} -R "{params.read_group}" | \
            samtools sort -@ 4 | \
            samtools view -@ 4 -hb > {output.bam}
        picard MarkDuplicates I={output.bam} O={output.deduped_bam} M={output.metrics}
        samtools index {output.deduped_bam}
        """

rule mapq_to_60:
    input:
        normal_bam=rules.realign_to_haplotype_ref.output.deduped_bam
    output:
        normal_bam_60=temp(f"temp/{sample}/{sample}_realigned_deduped_60.bam")
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
        out_vcf=temp(f'results/{patient}/calls/{sample}_germline.vcf.gz')
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
        out_vcf = f"results/{patient}/calls/{sample}_germline_filtered.vcf"
    shell:
        """
        gatk VariantFiltration \
            -R {input.hap_fa} \
            -V {input.in_vcf} \
            -O {output.out_vcf} \
            -filter-name "my_filter" \
            -filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        """
