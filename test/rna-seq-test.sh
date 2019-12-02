bbduk.sh in=mctp_SI_25272_CCW1MANXX_8_1.fq.gz in2=mctp_SI_25272_CCW1MANXX_8_2.fq.gz ref=refs/hs-hg38/hla/transcripts.fa \
         outm=match_1.fq outm2=match_2.fq k=21 \
         removeifeitherbad=f \
         threads=4 -Xmx12000m overwrite=true

STAR \
    --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeSAindexNbases 10 \
    --genomeDir star_MO_2503 \
    --genomeFastaFiles MO_2503.fa 

STAR \
    --runThreadN 4 \
    --twopassMode Basic \
    --genomeDir star_MO_2503 \
    --readFilesIn match_1.fq match_2.fq \
    --scoreGenomicLengthLog2scale 0 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000

samtools view -F 256 -b -o MO_2503.bam Aligned.out.sam
samtools sort -o MO_2503_sort.bam MO_2503.bam
samtools index MO_2503_sort.bam
