import pysam

configfile: "/home/mumphrey/Projects/hla_pipeline/config/call_somatic_mutations_config.yaml"

#Inputs
germline_refs = config['germline_refs']
tumor_bams = config['tumor_bams']
normal_bams = config['normal_bams']
gffs = config['gffs']
normal_bam_for_kmers = config['normal_bam_for_kmers']

#Sample information
patient = config['patient']
normal = config['normal']
tumor = config['tumor']

#Algorithm parameters
threads = config['threads']
genes = config['genes']

rule all:
    input:
        mutect_filtered = expand(f"results/{patient}/calls/{tumor}_{normal}_" + "{gene}_filtered.vcf", gene = genes),
        kmer_filtered = expand(f"results/{patient}/calls/kmer_{tumor}_{normal}_" + "{gene}_filtered.vcf", gene = genes),
        csq_annotated = expand(f"results/{patient}/calls/{tumor}_{normal}_" + "{gene}_annotated.txt", gene = genes)

rule germline_mapq_to_60:
    input:
        tumor_bam = lambda w: tumor_bams[w.gene],
        normal_bam = lambda w: normal_bams[w.gene]
    output:
        tumor_bam_60=temp(f"temp/{tumor}/germline_{tumor}_{{gene}}_realigned_deduped_60.bam"),
        normal_bam_60=temp(f"temp/{normal}/germline_{normal}_{{gene}}_realigned_deduped_60.bam")
    run:
        tumor_bam = pysam.AlignmentFile(input.tumor_bam, 'rb')
        tumor_bam_60 = pysam.AlignmentFile(output.tumor_bam_60, 'wb', template = tumor_bam)
        for read in tumor_bam:
            read.mapq = 60
            tumor_bam_60.write(read)
        normal_bam = pysam.AlignmentFile(input.normal_bam, 'rb')
        normal_bam_60 = pysam.AlignmentFile(output.normal_bam_60, 'wb', template = normal_bam)
        for read in normal_bam:
            read.mapq = 60
            normal_bam_60.write(read)

rule germline_index_60:
    input:
        tumor_bam_60=rules.germline_mapq_to_60.output.tumor_bam_60,
        normal_bam_60=rules.germline_mapq_to_60.output.normal_bam_60
    output:
        tumor_bam_60_bai=temp(rules.germline_mapq_to_60.output.tumor_bam_60 + ".bai"),
        normal_bam_60_bai=temp(rules.germline_mapq_to_60.output.normal_bam_60 + ".bai")
    shell:
        """
        samtools index {input.tumor_bam_60}
        samtools index {input.normal_bam_60}
        """

rule run_mutect2:
    input:
        germ_fa = lambda w: germline_refs[w.gene],
        tumor_bam = rules.germline_mapq_to_60.output.tumor_bam_60,
        normal_bam = rules.germline_mapq_to_60.output.normal_bam_60,
        tumor_bam_bai = rules.germline_index_60.output.tumor_bam_60_bai,
        normal_bam_bai = rules.germline_index_60.output.normal_bam_60_bai
    output:
        out_vcf=temp(f'results/{patient}/calls/{tumor}_{normal}_{{gene}}_raw.vcf.gz'),
        out_bam=temp(f'results/{patient}/calls/{tumor}_{normal}_{{gene}}_raw.bam')
    shell:
        """
        /home/mumphrey/Projects/hla_pipeline/bin/gatk-4.1.2.0/gatk --java-options "-Xmx2g" Mutect2 \
            -R {input.germ_fa} \
            -I {input.tumor_bam} \
            -I {input.normal_bam} \
            -tumor {tumor} \
            -normal {normal} \
            -O {output.out_vcf} \
            -bamout {output.out_bam} \
            --disable-read-filter MappingQualityAvailableReadFilter \
            --disable-read-filter MappingQualityNotZeroReadFilter \
            --disable-read-filter MappingQualityReadFilter
        """

rule filter_mutect2:
    input:
        germ_fa = lambda w: germline_refs[w.gene],
        in_vcf=rules.run_mutect2.output.out_vcf
    output:
        out_vcf = f"results/{patient}/calls/{tumor}_{normal}_{{gene}}_filtered.vcf"
    params:
        out_vcf_gz = f"results/{patient}/calls/{tumor}_{normal}_{{gene}}_filtered.vcf.gz"
    shell:
        """
        /home/mumphrey/Projects/hla_pipeline/bin/gatk-4.1.2.0/gatk --java-options "-Xmx2g" FilterMutectCalls \
            -V {input.in_vcf} \
            -O {params.out_vcf_gz} \
            -R {input.germ_fa} \
            --max-events-in-region 15
        gunzip {params.out_vcf_gz}
        """

#Filters out any variants that have kmer support in the normal sample
rule filter_normal_kmers:
    input:
        in_vcf = rules.filter_mutect2.output.out_vcf,
        normal_bam = normal_bam_for_kmers,
        germ_fa = lambda w: germline_refs[w.gene],
    output:
        out_vcf = f"results/{patient}/calls/kmer_{tumor}_{normal}_{{gene}}_filtered.vcf"
    shell:
        """
        python /home/mumphrey/Projects/hla_pipeline/scripts/filter_normal_kmers.py \
            {input.in_vcf} \
            {input.normal_bam} \
            {input.germ_fa} \
            25 \
            {output.out_vcf}
        """

rule run_csq:
    input:
        germ_fa = lambda w: germline_refs[w.gene],
        input_vcf=rules.filter_normal_kmers.output.out_vcf,
        input_gff=lambda w: gffs[w.gene]
    output:
        output_txt = f"results/{patient}/calls/{tumor}_{normal}_{{gene}}_annotated.txt"
    shell:
    	"""
    	bcftools csq \
    	-f {input.germ_fa} \
    	-g {input.input_gff} \
    	{input.input_vcf} \
    	-Ot \
    	-o {output.output_txt} \
    	-l \
        -i 'FILTER="PASS"'
    	"""
