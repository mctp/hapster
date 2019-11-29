# Polytect

## Installation

Clone repository
```
git clone <TODO>
```

### Setup External Dependencies

Polytect depends on GATK4, hpseq, biobambam2 and bwa.kit

- bwa.kit (>= 0.7.15):
  (https://sourceforge.net/projects/bio-bwa/files/bwakit/)  
- gatk4 (==4.1.2.0)
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

```
wget <TODO>
tar --strip 1 -xf resources/hs-hg38-1.0.0.tar.gz -C refs
```

##
```
```


## Running
```
conda activate polytect
PATH=$POLYTECT_DIR/bin:$PATH
```
