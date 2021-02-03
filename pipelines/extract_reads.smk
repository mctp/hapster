import os
import pysam
from pathlib import Path

# polytect runtime
PD = Path(config['POLYTECT_DIR'])
NCORES = config['NCORES']

# inputs
patient = config['patient']
sample = config['sample']
aligned_file = config['aligned_file']
cram_reference = config['cram_reference']
extraction_regions = config['extraction_regions']

# polytect global references
genes = config['genes']
kmers = PD / config['gene_prefix'] / "sim" / "kmers.txt"
regions = PD / config['gene_prefix'] / "sim" / "regions.txt"
complete_reference = PD / config['gene_prefix'] / "fa" / "complete.fa"
complete_reference_alt = PD / config['gene_prefix'] / "fa" / "complete.fa.alt"
extraction_reference = PD / config['gene_prefix'] / "fa" / "extraction.fa"
extraction_reference_alt = PD / config['gene_prefix'] / "fa" / "extraction.fa.alt"

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
        bam_for_kmer_filter = f"results/{patient}/alignments/{sample}/{sample}_extracted_with_blacklist.bam",
        fq1 = expand(f"results/{patient}/seqs/{sample}/{sample}_" + "{gene}_1.fq", gene = genes),
        fq2 = expand(f"results/{patient}/seqs/{sample}/{sample}_" + "{gene}_2.fq", gene = genes),
        bam = expand(f"results/{patient}/alignments/{sample}/{sample}_" + "{gene}_complete.bam", gene = genes)

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
        "samtools view -hb -@ {threads} --reference {input.cram_reference} {input.cram} $(<{input.extraction_regions})> {output.bam}"

rule bam_regions:
    input:
        bam = aligned_file,
        extraction_regions = extraction_regions
    output:
        bam = temp(f"temp/{sample}/{sample}_ra.bam")
    threads: NCORES
    shell:
        "samtools view -hb -@ {threads} {input.bam} $(<{input.extraction_regions})> {output.bam}"

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
        full_fq = f"temp/{sample}/{sample}_extracted_{{i}}.fq"
    shell:
        """
        cat {input.infiles} > {output.full_fq}
        """

# step 4: realign reads to the complete reference that contains
# all haplotypes, and perform postalt expansion
rule realign_reads_to_alt_ref:
    input:
        fq1 = f"temp/{sample}/{sample}_extracted_1.fq",
        fq2 = f"temp/{sample}/{sample}_extracted_2.fq"
    output:
        bam = temp(f"temp/{sample}/{sample}_alt_aligned.bam")
    params:
        read_group = f"@RG\\tSM:{sample}\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}"
    threads: NCORES
    shell:
        """
        bwa mem -t {threads} {complete_reference} {input.fq1} {input.fq2} -R "{params.read_group}" | \
            samtools sort -@ {threads} |
            samtools view -@ {threads} -hb > {output.bam}
        """

rule bwa_postalt:
    input:
        in_bam=rules.realign_reads_to_alt_ref.output.bam
    output:
        out_bam=temp(f'temp/{sample}/{sample}_postalt.bam'),
        out_bai=temp(f'temp/{sample}/{sample}_postalt.bam.bai')
    threads: NCORES        
    shell:
        """
        samtools view -h {input.in_bam} |
        k8 {PD}/bin/bwa-postalt.js {complete_reference_alt} |
        samtools view -hb |
        samtools sort -@ {threads} -o {output.out_bam}
        samtools index -@ {threads} {output.out_bam}
        """

# step 5: restrict reads to predefine regions
# Extracting from postalt gives high sensitivity and high specificity, but still
# leaves some reads that map to pseudogenes and highly similar blacklist genes
# these reads are further cleanuped
rule extract_all_gene_regions:
    input:
        in_bam = rules.bwa_postalt.output.out_bam,
        in_bai = rules.bwa_postalt.output.out_bai
    output:
        out_bam = temp(f"temp/{sample}/{sample}_gene_extracted.bam")
    shell:
        """
        samtools view -hb {input.in_bam} {full_regions} |\
        samtools sort -n -o {output.out_bam}
        """

# step 6: manual bamtofastq
# bwa-postalt.js creates many duplicate entries, and disrupts SAM
# flags so we have to extract pairs manually
rule unalign_reads_from_alt:
    input:
        in_bam = rules.extract_all_gene_regions.output.out_bam
    output:
        out_bam = f"results/{patient}/alignments/{sample}/{sample}_extracted_with_blacklist.bam",
        fq1 = temp(f"temp/{sample}/{sample}_gene_extracted_1.fq"),
        fq2 = temp(f"temp/{sample}/{sample}_gene_extracted_2.fq")
    run:
        bam = pysam.AlignmentFile(input.in_bam, 'rb')
        paired_bam = pysam.AlignmentFile(output.out_bam, 'wb', template=bam)
        pair = [False, False]
        cur_read = ''
        for read in bam:
            if read.query_name != cur_read:
                if type(pair[0]) != bool and type(pair[1]) != bool:
                    paired_bam.write(pair[0])
                    paired_bam.write(pair[1])
                pair = [False, False]
                cur_read = read.query_name
            if type(pair[0]) == bool and read.is_read1:
                pair[0] = read
            elif type(pair[1]) == bool and read.is_read2:
                pair[1] = read
        bam.close()
        paired_bam.close()
        subprocess.call(['bedtools', 'bamtofastq', '-i', output.out_bam, '-fq', output.fq1, '-fq2', output.fq2])

# step 7: remove reads which have a primary alignment to a blacklist gene
# Blacklist contains pseudogenes and highly similar genes that we aren't interested in
rule remove_blacklisted:
    input:
        fq1 = rules.unalign_reads_from_alt.output.fq1,
        fq2 = rules.unalign_reads_from_alt.output.fq2,
        extraction_reference = extraction_reference
    output:
        extracted_sorted_bam = temp(f"temp/{sample}/{sample}_extraction.bam"),
        extracted_sorted_bai = temp(f"temp/{sample}/{sample}_extraction.bam.bai"),
        no_blacklist_bam = temp(f"temp/{sample}/{sample}_extracted_no_blacklist.bam")
    threads: NCORES
    shell:
        """
        bwa mem -t {threads} {input.extraction_reference} {input.fq1} {input.fq2} | \
            samtools sort -@ {threads} > {output.extracted_sorted_bam}
        samtools index {output.extracted_sorted_bam}
        samtools view -hb -@ {threads} {output.extracted_sorted_bam} {full_regions} |
            samtools sort -@ {threads} -n > {output.no_blacklist_bam}
        """
# bwa-postalt.js creates many duplicate entries, and disrupts SAM
# flags so we have to extract pairs manually
rule unalign_reads_from_extraction:
    input:
        in_bam = rules.remove_blacklisted.output.no_blacklist_bam
    output:
        out_bam = temp(f"temp/{sample}/f_{sample}_no_blacklist_pairs.bam"),
        fq1 = temp(f"temp/{sample}/f_{sample}_extracted_no_blacklist_1.fq"),
        fq2 = temp(f"temp/{sample}/f_{sample}_extracted_no_blacklist_2.fq")
    run:
        # TODO: turn this into a function and re-use with rule unalign_reads_from_alt
        bam = pysam.AlignmentFile(input.in_bam, 'rb')
        paired_bam = pysam.AlignmentFile(output.out_bam, 'wb', template=bam)
        pair = [False, False]
        cur_read = ''
        for read in bam:
            if read.query_name != cur_read:
                if type(pair[0]) != bool and type(pair[1]) != bool:
                    paired_bam.write(pair[0])
                    paired_bam.write(pair[1])
                pair = [False, False]
                cur_read = read.query_name
            if type(pair[0]) == bool and read.is_read1:
                pair[0] = read
            elif type(pair[1]) == bool and read.is_read2:
                pair[1] = read
        bam.close()
        paired_bam.close()
        subprocess.call(['bedtools', 'bamtofastq', '-i', output.out_bam, '-fq', output.fq1, '-fq2', output.fq2])

rule align_to_alt_ref_final:
    input:
        fq1 = rules.unalign_reads_from_extraction.output.fq1,
        fq2 = rules.unalign_reads_from_extraction.output.fq2
    output:
        bam = temp(f"temp/{sample}/f_{sample}_alt_aligned_final.bam")
    params:
        read_group = f"@RG\\tSM:{sample}\\tID:{sample}\\tPL:ILLUMINA\\tLB:{sample}"
    threads: NCORES
    shell:
        """
        bwa mem -t {threads} {complete_reference} {input.fq1} {input.fq2} -R "{params.read_group}" | \
            samtools sort -@ {threads} |
            samtools view -@ {threads} -hb > {output.bam}
        """

rule bwa_postalt_final:
    input:
        in_bam=rules.align_to_alt_ref_final.output.bam
    output:
        out_bam=temp(f'temp/{sample}/f_{sample}_postalt_final.bam'),
        out_bai=temp(f'temp/{sample}/f_{sample}_postalt_final.bam.bai')
    shell:
        """
        samtools view -h {input.in_bam} |
        k8 {PD}/bin/bwa-postalt.js {complete_reference_alt} |
        samtools view -hb |
        samtools sort -o {output.out_bam}
        samtools index {output.out_bam}
        """

rule extract_single_gene_region:
    input:
        in_bam = rules.bwa_postalt_final.output.out_bam,
        in_bai = rules.bwa_postalt_final.output.out_bai
    output:
        out_bam = f"results/{patient}/alignments/{sample}/{sample}_{{gene}}_complete.bam"
    params:
        regions = lambda w: allele_regions[w.gene]
    shell:
        """
        samtools view -hb {input.in_bam} {params.regions} |\
        samtools sort -n -o {output.out_bam}
        """

rule unalign_reads:
    input:
        in_bam = rules.extract_single_gene_region.output.out_bam
    output:
        out_bam = temp(f"temp/{sample}/f_{sample}_{{gene}}_pairs.bam"),
        fq1 = f"results/{patient}/seqs/{sample}/{sample}_{{gene}}_1.fq",
        fq2 = f"results/{patient}/seqs/{sample}/{sample}_{{gene}}_2.fq"
    run:
        # TODO: turn this into a function and re-use with rule unalign_reads_from_alt
        bam = pysam.AlignmentFile(input.in_bam, 'rb')
        paired_bam = pysam.AlignmentFile(output.out_bam, 'wb', template=bam)
        pair = [False, False]
        cur_read = ''
        for read in bam:
            if read.query_name != cur_read:
                if type(pair[0]) != bool and type(pair[1]) != bool:
                    paired_bam.write(pair[0])
                    paired_bam.write(pair[1])
                pair = [False, False]
                cur_read = read.query_name
            if type(pair[0]) == bool and read.is_read1:
                pair[0] = read
            elif type(pair[1]) == bool and read.is_read2:
                pair[1] = read
        bam.close()
        paired_bam.close()
        subprocess.call(['bedtools', 'bamtofastq', '-i', output.out_bam, '-fq', output.fq1, '-fq2', output.fq2])
