from pathlib import Path

# hapster runtime
PD = Path(config['HAPSTER_DIR'])
NCORES = config['NCORES']

# inputs
germline_vcf = config['germline_vcf']
ref = config['ref']
patient = config['patient']
sample = config['sample']

# hapster global references
genes = config['genes']
gff = config['gff']
fq1s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_1.fq") for gene in genes]
fq2s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_2.fq") for gene in genes]

rule all:
    input:
        germline_ref = f"results/{patient}/refs/{patient}_germline_imputed.fa",
        germline_realigned_bam = f"results/{patient}/alignments/{sample}_germline_imputed.bam",
        germ_gff = f"results/{patient}/refs/{patient}_germline_imputed.gff",
        germ_gtf = f"results/{patient}/refs/{patient}_germline_imputed.gtf"

rule make_germline_ref:
    input:
        vcf = germline_vcf,
        hap_fa = ref,
        gff = gff
    output:
        germ_fa_raw = temp(f"results/{patient}/refs/{patient}_raw_germline_imputed.fa"),
        germ_fa = f"results/{patient}/refs/{patient}_germline_imputed.fa",
        germ_fa_amb = f"results/{patient}/refs/{patient}_germline_imputed.fa.amb",
        germ_fa_ann = f"results/{patient}/refs/{patient}_germline_imputed.fa.ann",
        germ_fa_bwt = f"results/{patient}/refs/{patient}_germline_imputed.fa.bwt",
        germ_fa_pac = f"results/{patient}/refs/{patient}_germline_imputed.fa.pac",
        germ_fa_fai = f"results/{patient}/refs/{patient}_germline_imputed.fa.fai",
        germ_fa_sa = f"results/{patient}/refs/{patient}_germline_imputed.fa.sa",
        germ_dict = f"results/{patient}/refs/{patient}_germline_imputed.dict",
        germ_gff = f"results/{patient}/refs/{patient}_germline_imputed.gff",
        germ_gtf = f"results/{patient}/refs/{patient}_germline_imputed.gtf"
    shell:
        """
        python {PD}/scripts/create_germline_ref_gff.py {input.hap_fa} {input.vcf} {output.germ_fa} {input.gff} {output.germ_gff}
        """

rule realign_to_germline_ref:
    input:
        rules.make_germline_ref.output,
        germ_fa = rules.make_germline_ref.output.germ_fa,
        fq1 = [x for x in fq1s],
        fq2 = [x for x in fq2s]
    output:
        temp_fq1 = temp(f"temp/{sample}/{sample}_temp_1.fq"),
        temp_fq2 = temp(f"temp/{sample}/{sample}_temp_2.fq"),
        bam = temp(f"temp/{sample}/germline_{sample}_realigned.bam"),
        deduped_bam = f"results/{patient}/alignments/{sample}_germline_imputed.bam",
        bam_bai = f"results/{patient}/alignments/{sample}_germline_imputed.bam.bai",
        metrics = temp(f"temp/{sample}/germline_{sample}_realigned_deduped_metrics.txt")
    params:
        fq1 = " ".join([x for x in fq1s]),
        fq2 = " ".join([x for x in fq2s]),
        read_group = f"@RG\\tSM:{sample}\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}"
    threads: 4
    shell:
        """
        for FQ in $(echo {params.fq1}); do if (( $(stat -c%s "$FQ") > 25 )); then cat $FQ >> {output.temp_fq1}; fi; done; 
        for FQ in $(echo {params.fq2}); do if (( $(stat -c%s "$FQ") > 25 )); then cat $FQ >> {output.temp_fq2}; fi; done;
        bwa mem -t 4 {input.germ_fa} {output.temp_fq1} {output.temp_fq2} -R "{params.read_group}" | \
            samtools sort -@ 4 | \
            samtools view -@ 4 -hb > {output.bam}
        picard MarkDuplicates I={output.bam} O={output.deduped_bam} M={output.metrics}
        samtools index {output.deduped_bam}
        """
