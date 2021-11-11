import os
import pysam
from pathlib import Path

# hapster runtime
PD = Path(config['HAPSTER_DIR'])
NCORES = config['NCORES']

# inputs
patient = config['patient']
sample = config['sample']
aligned_file = config['aligned_file']
cram_reference = config['cram_reference']
extraction_regions = config['extraction_regions']
read_length = int(config['readlength'])
overhang = read_length - 1

# hapster global references
kmers = PD / config['gene_prefix'] / "sim" / f"rna_kmers.txt"
regions = PD / config['gene_prefix'] / "sim" / "regions.txt"
extraction_reference = PD / config['gene_prefix'] / "fa" / f"extraction.fa"
extraction_reference_dir = PD / config['gene_prefix'] / "fa" / f"STAR_{read_length}"
gtf = PD / config['gene_prefix'] / "gff" / "alt_genes.gtf"

# derived variables
full_regions = ""
allele_regions = {}
with open(regions, 'r') as f:
    for line in f:
        line = line.strip().split()
        gene = line[0].replace(':', '')
        region = ' '.join(line[1:])
        full_regions += region + ' '
        allele_regions[gene] = region


# File for biobambam
# Changes depending on input format
bambam_format = Path(aligned_file).suffix[1:]
print(f"File format: {bambam_format}")
if bambam_format == "cram":
    bambam_file = f"temp/{sample}/{sample}_original.bam"
    bambam_format = "bam"
elif bambam_format == "bam" or bambam_format == "sam":
    bambam_file = f"temp/{sample}/{sample}_ra.bam"

rule all:
    input:
        no_blacklist_bam = f"results/{patient}/seqs/{sample}/{sample}_no_blacklist.bam",
        fq1 = f"results/{patient}/seqs/{sample}_1.fq",
        fq2 = f"results/{patient}/seqs/{sample}_2.fq"

# step 1 (optional): convert cram to bam
# this is required because bamtofastq breaks on some cram files
# TODO: support paired FASTQ file as input
rule cram_to_bam:
    input:
        cram = aligned_file,
        cram_reference = cram_reference,
        extraction_regions = extraction_regions
    output:
        bam = temp(f"temp/{sample}/{sample}_original.bam")
    threads: NCORES
    shell:
        "samtools view -hb -@ {threads} --reference {input.cram_reference} {input.cram} $(<{input.extraction_regions}) '*' > {output.bam}"

rule bam_regions:
    input:
        bam = aligned_file,
        extraction_regions = extraction_regions
    output:
        bam = temp(f"temp/{sample}/{sample}_ra.bam")
    threads: NCORES
    shell:
        "samtools view -hb -@ {threads} {input.bam} $(<{input.extraction_regions}) '*' > {output.bam}"

# step 2: convert bam to fastq
# Checkpoints used for dynamic output - we don't know how many pieces it will split into ahead of time
checkpoint unalign_reads_from_original:
    input:
        infile = bambam_file
    output:
        fastqs = directory(f"temp/bambamtemps/{sample}")
    params:
        out_dir = f"temp/bambamtemps/{sample}",
        fq1 = temp(f"temp/bambamtemps/{sample}/{sample}_original_1.fq"),
        fq2 = temp(f"temp/bambamtemps/{sample}/{sample}_original_2.fq"),
        singles = temp(f"temp/bambamtemps/{sample}/{sample}_original_singles.fq"),
        orphan1 = temp(f"temp/bambamtemps/{sample}/{sample}_original_orphans_1.fq"),
        orphan2 = temp(f"temp/bambamtemps/{sample}/{sample}_original_orphans_2.fq")
    shell:
        """
        mkdir -p {params.out_dir}
        bamtofastq \
            collate=1 \
            exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
            filename={input.infile} \
            gz=0 \
            inputformat={bambam_format} \
            F={params.fq1} \
            F2={params.fq2} \
            S={params.singles} \
            O={params.orphan1} \
            split=10K \
            O2={params.orphan2}
        """

# step 3: extract fq files that match HLA gene kmers
# Uses kmer extraction to try to grab all gene reads
# Nearly 100% sensitivity, but low specificity
# First pass extraction to save on compute time when doing bwa-postalt
# this runs in parallel for 10M unaligned chunks 
rule kmer_extraction:
    input:
        fq1 = f"temp/bambamtemps/{sample}/{sample}_original_1.fq_{{n}}",
        fq2 = f"temp/bambamtemps/{sample}/{sample}_original_2.fq_{{n}}"
    output:
        fq1_extracted = temp(f"temp/bambamtemps/{sample}/{sample}_extracted_1.fq_{{n}}"),
        fq2_extracted = temp(f"temp/bambamtemps/{sample}/{sample}_extracted_2.fq_{{n}}")
    threads: 1
    shell:
        """
        hpscan_cw \
            -s {kmers} \
            -1 {output.fq1_extracted} \
            -2 {output.fq2_extracted} \
            {input.fq1} \
            {input.fq2}
        """

def aggregate_fq(wildcards):
    checkpoint_output = checkpoints.unalign_reads_from_original.get(**wildcards).output[0]
    return sorted(expand("temp/bambamtemps/{sample}/{sample}_extracted_{i}.fq_{n}",
        sample = sample,
        i = wildcards.i,
        n = glob_wildcards(os.path.join(checkpoint_output, f"{sample}_original_{wildcards.i}.fq" + "_{n}")).n))

rule aggregate_kmer_extracted:
    input:
        infiles = aggregate_fq
    output:
        full_fq = f"temp/{patient}/seqs/{sample}_{{i}}.fq"
    shell:
        """
        cat {input.infiles} > {output.full_fq}
        """

rule remove_blacklisted:
    input:
        fq1 = f"temp/{patient}/seqs/{sample}_1.fq",
        fq2 = f"temp/{patient}/seqs/{sample}_2.fq",
        extraction_reference = extraction_reference
    output:
        extracted_sorted_bam = f"temp/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        extracted_sorted_bai = f"temp/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
        no_blacklist_bam = f"results/{patient}/seqs/{sample}/{sample}_no_blacklist.bam"
    params:
        prefix = f"temp/{sample}/{sample}."
    threads: NCORES
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {extraction_reference_dir} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            --outFilterMultimapNmax 2000 \
            --outReadsUnmapped Fastx \
            --sjdbGTFfile {gtf} \
            --sjdbOverhang {overhang} \
            --outFilterScoreMinOverLread 0 \
            --outFilterIntronStrands None \
            --outFilterMatchNminOverLread 0 \
            --outFilterMatchNmin 0 \
            --seedPerWindowNmax 50 \
            --winAnchorMultimapNmax 200 \
            --limitBAMsortRAM 50000000000
        samtools index {output.extracted_sorted_bam}
        samtools view -hb -@ {threads} {output.extracted_sorted_bam} {full_regions} |
            samtools sort -@ {threads} > {output.no_blacklist_bam}
        samtools index {output.no_blacklist_bam}
        """

checkpoint unalign_reads_from_extracted:
    input:
        infile = rules.remove_blacklisted.output.no_blacklist_bam
    output:
        fastqs = directory(f"temp/bambamtemps/f_{sample}"),
        fq1 = f"results/{patient}/seqs/{sample}_1.fq",
        fq2 = f"results/{patient}/seqs/{sample}_2.fq"
    params:
        out_dir = f"temp/bambamtemps/f_{sample}",
        fq1 = f"results/{patient}/seqs/{sample}_1.fq",
        fq2 = f"results/{patient}/seqs/{sample}_2.fq",
        singles = temp(f"temp/bambamtemps/f_{sample}/{sample}_original_singles.fq"),
        orphan1 = temp(f"temp/bambamtemps/f_{sample}/{sample}_original_orphans_1.fq"),
        orphan2 = temp(f"temp/bambamtemps/f_{sample}/{sample}_original_orphans_2.fq")
    shell:
        """
        mkdir -p {params.out_dir}
        bamtofastq \
            collate=1 \
            exclude=QCFAIL,SECONDARY,SUPPLEMENTARY \
            filename={input.infile} \
            gz=0 \
            inputformat=bam \
            F={params.fq1} \
            F2={params.fq2} \
            S={params.singles} \
            O={params.orphan1} \
            O2={params.orphan2}
        """