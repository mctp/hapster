from pathlib import Path

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

#Inputs
germline_vcf = config['vcf']
ref = config['ref']
fq1 = config['fq1']
fq2 = config['fq2']

#Sample information
patient = config['patient']
sample = config['sample']

#Algorithm parameters
threads = config['threads']
genes = config['genes']

rule all:
    input:
        germline_ref = expand(f"results/{patient}/refs/{patient}_" + "{gene}_germline_imputed.fa", gene = genes),
        germline_realigned_bam = expand(f"results/{patient}/alignments/{sample}_" + "{gene}_germline_imputed.bam", gene = genes)

rule make_germline_ref:
    input:
        vcf = germline_vcf,
        hap_fa = ref
    output:
        germ_fa_raw = temp(f"results/{patient}/refs/{patient}_raw_germline_imputed.fa"),
        germ_fa = f"results/{patient}/refs/{patient}_germline_imputed.fa",
        germ_fa_amb = f"results/{patient}/refs/{patient}_germline_imputed.fa.amb",
        germ_fa_ann = f"results/{patient}/refs/{patient}_germline_imputed.fa.ann",
        germ_fa_bwt = f"results/{patient}/refs/{patient}_germline_imputed.fa.bwt",
        germ_fa_pac = f"results/{patient}/refs/{patient}_germline_imputed.fa.pac",
        germ_fa_fai = f"results/{patient}/refs/{patient}_germline_imputed.fa.fai",
        germ_fa_sa = f"results/{patient}/refs/{patient}_germline_imputed.fa.sa",
        germ_dict = f"results/{patient}/refs/{patient}_germline_imputed.dict"
    shell:
        """
        python /home/mumphrey/Projects/hla_pipeline/scripts/create_germline_ref.py {input.hap_fa} {input.vcf} {output.germ_fa}
        """

rule realign_to_germline_ref:
    input:
        rules.make_germline_ref.output,
        germ_fa = rules.make_germline_ref.output.germ_fa,
        fq1 = fq1,
        fq2 = fq2
    output:
        bam = temp(f"temp/{sample}/germline_{sample}_realigned.bam"),
        deduped_bam = f"results/{patient}/alignments/{sample}_germline_imputed.bam",
        bam_bai = f"results/{patient}/alignments/{sample}_germline_imputed.bam.bai",
        metrics = temp(f"temp/{sample}/germline_{sample}_realigned_deduped_metrics.txt")
    params:
        read_group = f"@RG\\tSM:{sample}\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}"
    threads: 4
    shell:
        """
        bwa mem -t 4 {input.germ_fa} {input.fq1} {input.fq2} -R "{params.read_group}" | \
            samtools sort -@ 4 | \
            samtools view -@ 4 -hb > {output.bam}
        picard MarkDuplicates I={output.bam} O={output.deduped_bam} M={output.metrics}
        samtools index {output.deduped_bam}
        """