from pathlib import Path

# hapster runtime
PD = Path(config['HAPSTER_DIR'])
NCORES = config['NCORES']

# inputs
patient = config['patient']
sample = config['sample']
ref = config['ref']
overhang = int(config['readlength']) - 1
gtf = config['gtf']

rule all:
    input:
        bam = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
        bam_bai = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam.bai"

rule make_rna_ref:
    input:
        hap_fa = ref,
        gtf = gtf
    output:
        star_index = temp(f"results/{patient}/refs/{patient}_STAR/SAindex")
    params:
        genome_dir = f"results/{patient}/refs/{patient}_STAR"
    threads: NCORES
    shell:
        """
        mkdir -p {params.genome_dir}
        STAR --runThreadN {NCORES} \
            --runMode genomeGenerate \
            --genomeDir {params.genome_dir} \
            --genomeFastaFiles {input.hap_fa} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {overhang}
        """

rule realign_to_germline_ref:
    input:
        rules.make_rna_ref.output,
        fq1 = f"results/{patient}/seqs/{sample}_1.fq",
        fq2 = f"results/{patient}/seqs/{sample}_2.fq",
        gtf = gtf
    output:
        bam = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam",
        bam_bai = f"results/{patient}/alignments/{sample}.Aligned.sortedByCoord.out.bam.bai"
    params:
        genome_dir = f"results/{patient}/refs/{patient}_STAR",
        prefix = f"results/{patient}/alignments/{sample}."
    threads: NCORES
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params.genome_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts
        samtools index {output.bam}
        """
