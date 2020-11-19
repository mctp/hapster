from pathlib import Path

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

# inputs
ref = config['ref']
patient = config['patient']
sample = config['sample']

# polytect global references
genes = config['genes']
gtf = str(PD / config['gene_prefix'] / "gff" / "alt_genes.gtf")
fq1s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_1.fq") for gene in genes]
fq2s = [str(PD / "results" / patient / "seqs" / sample / f"{sample}_{gene}_2.fq") for gene in genes]

rule all:
    input:
        star_index = f"results/{patient}/refs/{sample}_rna/SAindex",
        bam = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
        bam_bai = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam.bai"

rule make_rna_ref:
    input:
        hap_fa = ref,
        gtf = gtf
    output:
        star_index = f"results/{patient}/refs/{sample}_rna/SAindex"
    params:
        genome_dir = f"results/{patient}/refs/{sample}_rna"
    threads: NCORES
    shell:
        """
        mkdir -p {params.genome_dir}
        STAR \
            --runMode genomeGenerate \
            --genomeDir {params.genome_dir} \
            --genomeFastaFiles {input.hap_fa} \
            --runThreadN {threads} \
            --sjdbOverhang 125 \
            --sjdbScore 2 \
            --sjdbGTFfile {input.gtf}
        """

rule realign_to_germline_ref:
    input:
        rules.make_rna_ref.output,
        fq1 = [x for x in fq1s],
        fq2 = [x for x in fq2s],
        gtf = gtf
    output:
        temp_fq1 = temp(f"temp/{sample}/{sample}_temp_1.fq"),
        temp_fq2 = temp(f"temp/{sample}/{sample}_temp_2.fq"),
        bam = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
        bam_bai = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params:
        fq1 = " ".join([x for x in fq1s]),
        fq2 = " ".join([x for x in fq2s]),
        genome_dir = f"results/{patient}/refs/{sample}_rna",
        prefix = f"results/{patient}/alignments/{sample}."
    threads: NCORES
    shell:
        """
        for FQ in $(echo {params.fq1}); do if (( $(stat -c%s "$FQ") > 25 )); then cat $FQ >> {output.temp_fq1}; fi; done; 
        for FQ in $(echo {params.fq2}); do if (( $(stat -c%s "$FQ") > 25 )); then cat $FQ >> {output.temp_fq2}; fi; done;
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {output.temp_fq1} {output.temp_fq2} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --sjdbGTFfile {input.gtf}
        samtools index {output.bam}
        """
