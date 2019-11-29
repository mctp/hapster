# Polytect

## Installation

Polytect is developed and tested on Ubuntu 18.04, but should work on modern Linux distributions.

### Pre-requisites

- Python 3.7
- Miniconda3 (https://docs.conda.io/en/latest/miniconda.html)
- gcc tool chain to compile programs

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

### Setup External Dependencies

Polytect depends on a number of dependencies managed by conda and select one which are not GATK4, hpseq, biobambam2 and bwa.kit.

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
conda activate polytect
```

### Setup references

```bash
wget https://storage.googleapis.com/mctp-open-share/polytect/hs-hg38-1.0.0.tar.gz --directory-prefix=resources
tar --strip 1 -xf resources/hs-hg38-1.0.0.tar.gz -C refs
polytect setup_refs
```

## Running

### Acticate the Conda environment and setup path.

```bash
conda activate polytect
```

### Buidling Re

```bash
polytect make_refs
```
