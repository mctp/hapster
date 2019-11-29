"""
Takes a set of gene names and sequences and creates reference files for the mutation calling pipeline
"""
import re
import subprocess
from Bio import SeqIO
from pathlib import Path

primary_ref = config['primary_ref']
analysis_ref = config['analysis_ref']
gene_prefix = config['gene_prefix']
genes = config['genes']

# create gff/fa dictionaries
gff_files = {}
fasta_files = {}
for gene in genes:
    gff_files[gene] = "%s/gff/%s.gff" % (gene_prefix, gene)
    fasta_files[gene] = "%s/alts/%s.fa" % (gene_prefix, gene)
blacklist_fasta = Path(gene_prefix) / "alts" / "blacklist.fa"

rule all:
    input:
        expand("{gene_prefix}/sim/kmers.txt", gene_prefix=gene_prefix),
        expand("{gene_prefix}/sim/regions.txt", gene_prefix=gene_prefix),
        expand("{gene_prefix}/gff/alt_genes.gff", gene_prefix=gene_prefix),
        # expand("{gene_prefix}/fa/extraction.fa", gene_prefix=gene_prefix),
        # expand("build/refs/finished_single_refs_{gene}.txt", gene = genes),
        expand("{gene_prefix}/fa/complete.fa", gene_prefix=gene_prefix),
        expand("{gene_prefix}/fa/complete.fa.alt", gene_prefix=gene_prefix)

# Makes each record its own fasta to provide easy access when generating
# germline refs later. Probably better to eventually just use a fasta index
# for the original fastas rather than writing so many files.
rule make_single_refs:
    input:
        fasta = lambda w: fasta_files[w.gene]
    output:
        finished = temp("build/refs/finished_single_refs_{gene}.txt")
    run:
        output_dir = Path(gene_prefix) / "fa/single_refs"
        output_dir.mkdir(parents=True, exist_ok=True)
        for record in SeqIO.parse(input.fasta, "fasta"):
            seq_filename = output_dir / (re.sub('[^0-9a-zA-Z]+', '_', record.id) + ".fa")
            with open(seq_filename, 'w') as of:
                of.write(f">{record.id}\n")
                of.write(f"{str(record.seq)}\n")
        with open(output.finished, 'w') as of:
            pass

# Makes a list of regions for each gene that captures all nucleotides of all alleles
# Used later with samtools to extract reads across all alleles for a given gene
rule make_allele_regions:
    input:
        fastas = [fasta for fasta in fasta_files.values()]
    output:
        regions = "{gene_prefix}/sim/regions.txt"
    run:
        with open(output.regions, 'w') as of:
            for gene, fasta in fasta_files.items():
                of.write(f"{gene}: ")
                for record in SeqIO.parse(fasta, "fasta"):
                    of.write(f"{record.id}:1-{len(record.seq)} ")
                of.write("\n")

# Makes a single gff for later use with the somatic mutation consequence caller
rule consolidate_gffs:
    input:
        gffs = [x for x in gff_files.values()]
    output:
        full_gff = "{gene_prefix}/gff/alt_genes.gff"
    shell:
        """
        cat {input.gffs} >> {output.full_gff}
        """

# Makes a single gff for later use with the somatic mutation consequence caller
rule consolidate_fastas:
    input:
        fasta_files = [x for x in fasta_files.values()]
    output:
        consolidated_fasta = temp("build/refs/alt_genes.fa")
    shell:
        """
        cat {input} >> {output.consolidated_fasta}
        """

# Makes an alt index for our new fasta for use with bwa-postalt.js
rule make_alt_index:
    input:
        fasta = rules.consolidate_fastas.output.consolidated_fasta
    output:
        alt_index_raw = temp("build/refs/alt_genes_raw.fa.alt"),
        alt_index = temp("build/refs/alt_genes.fa.alt")
    threads: 12
    run:
        alt_index_raw = open(output.alt_index_raw, 'w')
        subprocess.call(['minimap2', '-r7k', '-a', primary_ref, input.fasta], stdout = alt_index_raw)
        alt_index_raw.close()
        alt_index_raw = open(output.alt_index_raw, 'r')
        alt_index = open(output.alt_index, 'w')
        for line in alt_index_raw:
            if not line.startswith('@'):
                line = line.strip().split()
                line[9] = '*'
                new_line = '\t'.join(line[:12] + line[13:])
                alt_index.write(f"{new_line}\n")
        alt_index.close()

# Creates a complete indexed reference for alt allele aware alignment
rule make_alt_ref:
    input:
        analysis = analysis_ref,
        analysis_alt = analysis_ref + ".alt",
        alt_fa = rules.consolidate_fastas.output.consolidated_fasta,
        alt_index = rules.make_alt_index.output.alt_index
    output:
        complete_fasta = "{gene_prefix}/fa/complete.fa",
        complete_index = "{gene_prefix}/fa/complete.fa.alt"
    threads: 4
    shell:
        """
        cat {input.analysis} > {output.complete_fasta}
        cat {input.alt_fa} >> {output.complete_fasta}
        cat {input.analysis_alt} > {output.complete_index}
        cat {input.alt_index} >> {output.complete_index}
        bwa index {output.complete_fasta}
        samtools faidx {output.complete_fasta}
        picard CreateSequenceDictionary R={output.complete_fasta}
        """

# Creates a reference containing only our whitelist and blacklist alleles, no Grch38
# Used to separate blacklist reads from whitelist reads
rule make_extraction_ref:
    input:
        main_fasta = rules.consolidate_fastas.output.consolidated_fasta,
        blacklist_fasta = blacklist_fasta
    output:
        fasta_raw = temp("build/refs/{gene_prefix}/extraction_raw.fa"),
        fasta = "{gene_prefix}/fa/extraction.fa"
    shell:
        """
        cat {input} >> {output.fasta_raw}
        picard NormalizeFasta I={output.fasta_raw} O={output.fasta} LINE_LENGTH=70
        bwa index {output.fasta}
        samtools faidx {output.fasta}
        picard CreateSequenceDictionary R={output.fasta}
        """

# Creates all kmers from all alt alleles for use in kmer read extraction
rule make_kmers:
    input:
        fastas = rules.consolidate_fastas.output.consolidated_fasta  
    output:
        kmers = "{gene_prefix}/sim/kmers.txt"
    shell:
        """
        python scripts/gen_kmers.py {input.fastas} {output.kmers} 30 1
        """