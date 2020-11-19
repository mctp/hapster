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
fq1s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_1.fq") for gene in genes]
fq2s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_2.fq") for gene in genes]
fasta_files = [str(PD / config['gene_prefix'] / "rna" / f"{gene}.fa") for gene in genes]

#Derived variables
with open(haplotype, 'r') as f:
    alleles = next(f).strip().split(',')[1:]
    haplotype_string = ' '.join(alleles)

rule all:
    input:
        hap_fa_saved = f"results/{patient}/refs/{sample}_transcripts.fa",
        germline_realigned_bam = f"results/{patient}/alignments/{sample}_germline_imputed.bam"

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
        hap_fa_saved = f"results/{patient}/refs/{sample}_transcripts.fa"
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

rule realign_to_germline_ref:
    input:
        rules.make_haplotype_ref.output,
        hap_fa = rules.make_haplotype_ref.output.hap_fa_saved,
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
