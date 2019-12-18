# Polytect

## Installation

Polytect is developed and tested on Ubuntu 18.04, but should work on modern Linux distributions.

### Pre-requisites

- Python 3.7
- Miniconda3 (https://docs.conda.io/en/latest/miniconda.html)
- gcc tool chain to compile programs
- unzip

### Setup Polytect base directory

Clone repository from github:
```bash
git clone https://github.com/mctp/polytect
```

Modify path to get access to `polytect` command:
```bash
cd polytect
PATH=$PWD/bin:$PATH
```

### Running

Before you can use any of the polytect commands, you need to activate the Conda environment and setup path.

```bash
conda activate polytect
```

### Setup External Dependencies

Polytect depends on a number of dependencies managed by conda and select ones which are not: GATK4, hpseq, biobambam2 and bwa.kit.

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
wget https://storage.googleapis.com/mctp-open-share/polytect/tools-1.0.0.tar.gz --directory-prefix=resources
polytect setup_tools
```

### Setup references

```bash
conda activate polytect
wget https://storage.googleapis.com/mctp-open-share/polytect/hs-hg38-1.0.1.tar.gz --directory-prefix=resources
tar --strip 1 -xf resources/hs-hg38-1.0.1.tar.gz -C refs
polytect make_refs hs-hg38-hla
```
### Simulate likelihood matrices

If using the haplotyping module, you must first simulate matrices that represent the probability of reads generated from one allele aligning to another allele. This can be done using the make_matrices pipeline and requires the following inputs:
  protocol: wgs or wes, depending on the experimental setup
  min_insert_length: min insert size to simulate, recommended as mean insert length - 2 * SD from experiment being simulated
  max_insert_length: max insert size to simulate, recommended as mean insert length + 2 * SD from experiment being simulated
  read_length: length of reads in the experiment being simulated
  n_reads: number of reads to simulate, recommend at least 2000
  nm: max nm score to consider a "good" alignment, recommend 1
  capture_targets: required if simulating a WES experiment, list of capture probes
  similarity: percent sequence similarity between probe and insert to consider it captured, recommend .82
  
  ```
  # sample command for whole genome sequencing
  # make_matrices [gene] [protocol] [min_insert_length] [max_insert_length] [read_length] [n_reads] [nm]
  make matrices hs-hg38-hla wgs 125 325 151 2000 1
  ```
