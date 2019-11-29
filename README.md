# Polytect

## Installation

Clone repository
```
git clone https://github.com/mctp/polytect
```

### Pre-requisites

- Python 3.7
- Miniconda3

### Setup External Dependencies

Polytect depends on a number of dependencies managed by conda and select one which are not GATK4, hpseq, biobambam2 and bwa.kit.

- bwa.kit (==0.7.15)
  No other version supported.
  (https://sourceforge.net/projects/bio-bwa/files/bwakit/)  
- gatk4 (==4.1.2.0)
  No other version supported.
  (https://github.com/broadinstitute/gatk/releases/download/4.1.2.0/gatk-4.1.2.0.zip)  
- hpseq
  (https://github.com/mcieslik-mctp/hpseq)  
- biobambam2
  (https://gitlab.com/german.tischler/biobambam2/-/releases)  

These dependencies are included in the provide tools resource bundle.

```
wget <TODO>
tar xf resources/tools-1.0.0.tar.gz -C build
bash build/tools-1.0.0/setup-tools.sh
```

### Setup references

```bash
wget <TODO>
tar --strip 1 -xf resources/hs-hg38-1.0.0.tar.gz -C refs
```

## Running

### Acticate the Conda environment and setup path.

```bash
conda activate polytect
PATH=$POLYTECT_DIR/bin:$PATH
```

### Buidling Re

```bash
polytect make_refs
```
