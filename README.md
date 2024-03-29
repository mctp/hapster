# Hapster

## Installation

Hapster is developed and tested on Ubuntu 18.04, but should work on modern Linux distributions.

### Pre-requisites

- Python 3.7
- Miniconda3 (https://docs.conda.io/en/latest/miniconda.html)
- gcc tool chain to compile programs
- unzip
- R >= 3.2
- R packages:
  - rtracklayer
  - Biostrings
  - GenomicFeatures
  - BSgenome
  - foreach 
  - readr
  - dplyr
  - tibble
  - stringr

### Setup Hapster base directory

Clone repository from github:
```bash
git clone https://github.com/mctp/hapster
```

Modify path to get access to `hapster` command:
```bash
cd hapster
PATH=$PWD/bin:$PATH
```

### Setup External Dependencies

Hapster depends on a number of dependencies managed by conda and select ones which are not: GATK4, hpseq, biobambam2 and bwa.kit.

- bwa.kit (==0.7.15) No other version supported.
  (https://sourceforge.net/projects/bio-bwa/files/bwakit/)  
- gatk4 (==4.1.2.0) No other version supported.
  (https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip)  
- hpseq
  (https://github.com/mcieslik-mctp/hpseq)  
- biobambam2
  (https://gitlab.com/german.tischler/biobambam2/-/releases)  

These dependencies are included in the provided tools resources bundle.

```bash
wget https://storage.googleapis.com/mctp-open-share/hapster/tools-1.0.1.tar.gz --directory-prefix=resources
hapster setup_tools
```

### Running

Before you can use any of the hapster commands, you need to activate the Conda environment and setup path.

```bash
conda activate hapster
```

### Setup references

```bash
conda activate hapster
wget https://storage.googleapis.com/mctp-open-share/hapster/hs-hg38-1.0.2.tar.gz --directory-prefix=resources
tar --strip 1 -xf resources/hs-hg38-1.0.2.tar.gz -C refs
hapster make_refs hs-hg38-hla
```
### Simulate likelihood matrices

If using the haplotyping module, you must first simulate matrices that represent the probability of reads generated from one allele aligning to another allele. This can be done using the make_matrices pipeline and requires the following inputs:
 - gene: gene name
 - protocol: wgs or wes, depending on the experimental setup  
 - min_insert_length: min insert size to simulate, recommended as mean insert length - 2 * SD from experiment being simulated, or if this value is less than your read_length then just set equal to read_length
 - max_insert_length: max insert size to simulate, recommended as mean insert length + 2 * SD from experiment being simulated  
 - read_length: length of reads in the experiment being simulated  
 - n_reads: number of reads to simulate, recommend at least 2000  
 - nm: max nm score to consider a "good" alignment, recommend 1  
 - capture_targets: required if simulating a WES experiment, list of capture probes  
 - similarity: percent sequence similarity between probe and insert to consider it captured, recommend .82  
  
This module expects a set of reference files that can be created with the make_refs module.  
  
NOTE: This is just an example command. There are no defaults for min_insert_length and max_insert_length, or read_length. These must be determined from the sequencing experiment being simulated.
  ```
  # sample command for whole genome sequencing
  # make_matrices [gene] [protocol] [min_insert_length] [max_insert_length] [read_length] [n_reads] [nm] [capture_targets] [similarity]
  hapster make_matrices hs-hg38-hla wgs 151 674 151 2000 1
  ```

### Extract reads specific to the region of interest
To improve run time for later parts of the algorithm, we extract reads specific to the region of interest. Since we expect the base alignment to be poor, we don't know for sure that all of our reads are at the annotated coordinates for our genes, so we can't just extract based on location. Requires the following inputs:
 - gene: gene name
 - patient: patient ID
 - sample: sample ID
 - aligned_file: filepath to alignment we are extracting reads from
 - extraction_regions: BED file containing the regions to be extracted from the original aligned file
 - cram_reference: filepath to reference used to generate a cram file, only required if extracting reads from a .cram file

This module expects a set of reference files that can be created with the make_refs module.
```
# hapster extract_reads [gene] [patient] [sample] [aligned_file] [extraction_regions] <cram_reference>
hapster extract_reads hs-hg38-hla patient1 sample1 /path/to/sample1.bam /path/to/extraction/regions.bed
```
### Infer haplotype
If haplotyping has already been done, this step can be skipped. Otherwise, the command can be run with the following inputs:
 - gene: gene name
 - patient: patient ID
 - sample: sample ID
 - nm: max nm score to consider a "good" alignment, must match what was used in make_matrices
 - procotol: wgs or wes depending on experimental setup

This module expects a set of reference files that can be created with the make_refs module, and a set of extracted BAMs that can be produced with the extract_reads module.
```
# hapster infer_haplotype [gene] [patient] [sample] [nm] [protocol]
hapster infer_haplotype hs-hg38-hla patient1 sample1 1 wgs
```

### Call germline mutations
Calls germline mutations relative to a haplotype reference. This can be either produced with the infer_haplotype module, or provided based on external haplotyping. The command can be run with the following inputs:
 - gene: gene name
 - patient: patient ID
 - sample: sample ID
 - haplotype: filepath to CSV containing haplotype information

This module expects a set of reference files that can be created with the make_refs module, and a set of extracted BAMs that can be produced with the extract_reads module.
```
# hapster germline_mutations [gene] [patient] [sample] [haplotype]
hapster germline_mutations hs-hg38-hla patient1 sample1 /path/to/haplotype.csv
```

### Realign to imputed germline reference
Realigns DNA sequencing reads to the new haplotype reference with imputed novel germline variants. This needs to be done for both the tumor and the normal. The command can be run with the following inputs:
 - gene: gene name
 - patient: patient ID
 - sample: sample ID
 - germline_vcf: The VCF created by the germline_mutations function
 - haplotype_ref: The haplotype reference fasta created by the germline_mutations function

This module expects a set of reference files that can be created with the make_refs module, and a set of extracted BAMs that can be produced with the extract_reads module.
```
# hapster dna_realign [gene] [patient] [sample] [germline_vcf] [haplotype_ref]
hapster dna_realign hs-hg38-hla patient1 sample1 /path/to/germline.vcf /path/to/hapotype.fa
```
### Call somatic mutations
Calls somatic mutations using mutect2 and the haplotype reference realigned sequencing data.The command can be run with the following inputs:
 - gene: gene name
 - patient: patient ID
 - normal_sample: Normal sample ID
 - normal_alignment: path to the normal alignment generated by dna_realign
 - tumor_sample: Tumor sample ID
 - tumor_alignment: path to the tumor alignment generated by dna_realign
 - germline_ref: path to the germline imputed reference created by dna_realign
 - gff: path to the GFF file for all alternate genes created by make_refs
 - normal_for_kmers: path to original normal sequencing data used to look for variant supporting kmers in unaligned reads

```
# hapster somatic_mutations [gene] [patient] [normal_sample] [normal_alignment] [tumor_sample] [tumor_alignment] [germline_ref] [gff] [normal_for_kmers]
hapster somatic_mutations hs-hg38-hla patient1 sample1 /path/to/sample1.bam sample2 /path/to/sample2.bam /path/to/germline_imputed_ref.fa /path/to/alt_genes.gff /path/to/normal.bam
```
