#!/usr/bin/env python3

# containers
porechop_container = 'shub://TomHarrop/ont-containers:porechop_0.2.4'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'
minimap_container = 'docker://staphb/minimap2:2.24'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'
kraken_container = 'docker://staphb/kraken2:2.1.2-no-db'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# RULES #
#########

rule target:
    input:
        # nanopore mapping
        expand('output/nanopore-{mode}/minimap2-{mode}/samtools_coverage.out', mode=["hac", 'sup']),
        expand('output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1_trimmed_fastqc.zip', mode=["hac", 'sup']),
        expand('output/nanopore-{mode}/minimap2-{mode}/mapping_stats.out', mode=["hac", 'sup']),
        # non-hyp reads
            # illumina
        'output/illumina/extract-non-hyp/non-hyp-reads-repaired_r1.fq',
        'output/illumina/extract-non-hyp/non-hyp-reads_fastqc.zip',
        'output/illumina/kraken/kraken_all_report.txt',
            # nanopore
        'output/nanopore-sup/extract-non-hyp/non-hyp-reads.fq.gz',
        expand('output/nanopore-{mode}/extract-non-hyp/non-hyp-fasta-faidx.out', mode=["hac", "sup"]),
        #'output/nanopore/minimap2_viral_contigs/mapping_stats.out'

###########################################
## map non-hyp nanopore to viral contigs ##
###########################################

rule samtools_flagstat_viral_contigs:
    input:
        'output/minimap2_initial_viral_contigs/minimap2.sam'
    output:
        'output/minimap2_initial_viral_contigs/mapping_stats.out'
    log:
        'output/logs/samtools_flagstat_viral_contigs.log'
    shell:
        'samtools flagstat '
        '{input} > {output} 2> {log}'

rule minimap2_viral_contigs:
    input:
        genome = 'data/initial_viral_contigs.fasta',
        reads = 'output/illumina/extract-non-hyp/non-hyp-reads.fq.gz'
    output:
        sam = 'output/minimap2_initial_viral_contigs/minimap2.sam'
    threads:
        20
    singularity:
        minimap_container
    log:
        'output/logs/minimap2_initial_viral.log'
    shell:
        'minimap2 '
        '-ax map-ont '
        '{input.genome} '
        '{input.reads} '
        '-t {threads} '
        '> {output} '
        '2> {log}'

############################
##  extract non-hyp reads ##
############################

rule fastqc_illumina_nonhyp:
    input:
        'output/illumina/extract-non-hyp/non-hyp-reads.fq.gz'
    output:
        'output/illumina/extract-non-hyp/non-hyp-reads_fastqc.zip'
    params:
        outdir = directory('output/illumina/extract-non-hyp')
    singularity:
        fastqc_container
    shell:
        'mkdir -p {params.outdir} ; '
        'fastqc --outdir {params.outdir} {input}'

rule gzip_reads:
    input:
        non_hyp_reads = 'output/{reads}/extract-non-hyp/non-hyp-reads.fq'
    output:
        non_hyp_reads = 'output/{reads}/extract-non-hyp/non-hyp-reads.fq.gz'
    shell:
        'gzip {input.non_hyp_reads}'

rule nanopore_non_hyp_faidx:
    input:
        'output/nanopore-{guppy}/extract-non-hyp/non-hyp.fasta'
    output:
        'output/nanopore-{guppy}/extract-non-hyp/non-hyp-fasta-faidx.out'
    log:
        'output/logs/nanopore_non_hyp_faidx-{guppy}.log'
    shell:
        'samtools faidx {input} -o {output} 2> {log}'

rule nanopore_non_hyp_fasta:
    input:
        'output/nanopore-{guppy}/extract-non-hyp/non-hyp-reads.fq'
    output:
        'output/nanopore-{guppy}/extract-non-hyp/non-hyp.fasta'
    shell:
        "sed -n '1~4s/^@/>/p;2~4p' {input} > {output}"

rule repair_illumina:
    input:
        r1 = 'output/illumina/extract-non-hyp/non-hyp-reads_r1.fq',
        r2 = 'output/illumina/extract-non-hyp/non-hyp-reads_r2.fq'
    output:
        r1 = 'output/illumina/extract-non-hyp/non-hyp-reads-repaired_r1.fq',
        r2 = 'output/illumina/extract-non-hyp/non-hyp-reads-repaired_r2.fq',
        singletons = 'output/illumina/extract-non-hyp/non-hyp-reads-repaired_singletons.fq'
    log:
        'output/logs/repair_illumina.log'
    singularity:
        bbduk_container
    shell:
        'repair.sh '
        'in1={input.r1} in2={input.r2} '
        'out1={output.r1} out2={output.r2} outs={output.singletons} '
        'repair '
        '2> {log}'

rule extract_non_hyp_reads_illumina_r1_r2: # should contain unmapped also I think?
    input:
        fil_bam = 'output/illumina/extract-non-hyp/non-Mh.bam'
    output:
        r1 = 'output/illumina/extract-non-hyp/non-hyp-reads_r1.fq',
        r2 = 'output/illumina/extract-non-hyp/non-hyp-reads_r2.fq'
    threads:
        20
    log:
        'output/logs/extract_non_hyp_illumina_r1_r2_reads.log'
    shell:
        'samtools bam2fq -1 {output.r1} -2 {output.r2} --threads {threads} {input.fil_bam} 2> {log}'

rule extract_non_hyp_reads: # should contain unmapped also I think?
    input:
        fil_bam = 'output/{reads}/extract-non-hyp/non-Mh.bam'
    output:
        non_hyp_reads = 'output/{reads}/extract-non-hyp/non-hyp-reads.fq'
    threads:
        20
    log:
        'output/logs/extract_non_hyp_{reads}_reads.log'
    shell:
        'samtools bam2fq --threads {threads} {input.fil_bam} > {output.non_hyp_reads} 2> {log}'

rule filter_bam:
    input:
        fil_bed = 'output/{reads}/extract-non-hyp/non-Mh.bed',
        bam = 'data/reads/Mh_{reads}.bam'
    output:
        fil_bam = 'output/{reads}/extract-non-hyp/non-Mh.bam'
    log:
        'output/logs/filter_bam_{reads}.log'
    threads:
        20
    shell:
        'samtools view --threads {threads} -h '
        '-L {input.fil_bed} '
        '{input.bam} > {output.fil_bam} '
        '2> {log}'

rule bed_extract_non_hyp:
    input:
        bed = 'output/{reads}/extract-non-hyp/Mh.bed'
    output:
        fil_bed = 'output/{reads}/extract-non-hyp/non-Mh.bed'
    params:
        pattern = "'scaffold_1|scaffold_2|scaffold_3|scaffold_4|scaffold_5|scaffold_6|scaffold_7|scaffold_8|scaffold_9|scaffold_10|scaffold_11|scaffold_12'"
    log:
        'output/logs/bed_extract_non_hyp_{reads}.log'
    shell:
        'grep -E --word-regexp --invert-match {params.pattern} {input.bed} > {output.fil_bed} 2> {log}' # look for lines that don't match ids in that file

rule make_bed_file:
    input:
        'data/reads/Mh_{reads}.bam'
    output:
        'output/{reads}/extract-non-hyp/Mh.bed'
    log:
        'output/logs/make_bed_file_{reads}.log'
    shell:
        'bedtools bamtobed '
        '-i {input} > {output} '
        '2> {log}'

##################################
# nanopore basecalling & mapping #
##################################

rule samtools_coverage: # looks good!
    input:
        sorted_bam = 'output/nanopore-{mode}/minimap2-{mode}/sorted.bam'
    output:
        depth_out = 'output/nanopore-{mode}/minimap2-{mode}/samtools_coverage.out'
    log:
        'output/logs/samtools_coverage_{mode}.log'
    threads:
        20
    shell:
        'samtools coverage '
        '{input.sorted_bam} '
        '-o {output.depth_out} '
        '2> {log}'

rule samtools_flagstat:
    input:
        'output/nanopore-{mode}/minimap2-{mode}/sorted.bam'
    output:
        'output/nanopore-{mode}/minimap2-{mode}/mapping_stats.out'
    log:
        'output/logs/samtools_flagstat-{mode}.log'
    shell:
        'samtools flagstat '
        '{input} > {output} 2> {log}'

rule samtools_sort:
    input:
        sam = 'output/nanopore-{mode}/minimap2-{mode}/mapped_nanopore.sam'
    output:
        sorted_bam = 'output/nanopore-{mode}/minimap2-{mode}/sorted.bam',
        bam_link = 'data/reads/Mh_nanopore-{mode}.bam'
    params:
        bam_full_path = '/Volumes/archive/deardenlab/sarahinwood/mh_projects/mh-virus-nanopore/output/nanopore-{mode}/minimap2-{mode}/sorted.bam'
    threads:
        20
    log:
        'output/logs/samtools_sort-{mode}.log'
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '--threads {threads} '
        '2> {log} || exit 1 ; '
        'ln -s {params.bam_full_path} {output.bam_link}'

rule minimap2:
    input:
        genome = 'data/Mh_genome.fa',
        reads = 'output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1_trimmed.fq.gz'
    output:
        sam = 'output/nanopore-{mode}/minimap2-{mode}/mapped_nanopore.sam'
    threads:
        20
    singularity:
        minimap_container
    log:
        'output/logs/minimap2_{mode}.log'
    shell:
        'minimap2 '
        '-ax map-ont '
        '{input.genome} '
        '{input.reads} '
        '-t {threads} '
        '> {output} '
        '2> {log}'

rule fastqc_nanopore:
    input:
        'output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1_trimmed.fq.gz'
    output:
        'output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1_trimmed_fastqc.zip'
    params:
        outdir = directory('output/nanopore/joined')
    threads:
        10
    singularity:
        fastqc_container
    shell:
        'mkdir -p {params.outdir} ; '
        'fastqc --outdir {params.outdir} {input} -t {threads}'

rule porechop_adapters:
    input:
        'output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1.fq.gz'
    output:
        'output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1_trimmed.fq.gz'
    params:
        out = 'output/nanopore-{mode}/joined/hyp-nanopore-{mode}_1_trimmed.fq'
    log:
        'output/logs/porechop_adapters_{mode}.log'
    threads:
        40
    singularity:
        porechop_container
    shell:
        'porechop '
        '-i {input} '
        '-o {params.out} '
        '--threads {threads} '
        '2> {log} || exit 1 ; '
        'gzip {params.out}'

rule guppy_basecalling_sup:
    input:
        'data/reads/hyperodae-virus-nanopore/1/'
    output:
        joined_fq = 'output/nanopore-sup/joined/hyp-nanopore-sup_1.fq.gz'
    params:
        wd = 'output/nanopore-sup/guppy_sup',
        config_file = '/Volumes/archive/scratch/deardenlab/Jgilligan/ont-guppy/data/dna_r9.4.1_450bps_sup.cfg'
    log:
        'output/logs/guppy_basecalling_sup.log'
    threads:
        40
    shell:
        "bin/ont-guppy/bin/guppy_basecaller "
        "--device 'cuda:0' " # on bcc3
        "--input_path {input} "
        "--save_path {params.wd} "
        "-c {params.config_file}  " # config file sets basecalling mode
        "--recursive "
        "--compress_fastq "
        "2> {log} || exit 1 ; "
        "cat output/nanopore/guppy_sup/pass/*.fastq.gz > {output.joined_fq}"

rule guppy_basecalling:
    input:
        'data/reads/hyperodae-virus-nanopore/1/'
    output:
        joined_fq = 'output/nanopore-hac/joined/hyp-nanopore-hac_1.fq.gz'
    params:
        wd = 'output/nanopore/guppy_hac'
    log:
        'output/logs/guppy_basecalling_hac.log'
    threads:
        40
    shell:
        "bin/ont-guppy/bin/guppy_basecaller " # version 6.0.0
        "--device 'cuda:0' " # on bcc3
        "--input_path {input} "
        "--save_path {params.wd} "
        "--flowcell 'FLO-MIN106' "
        "--kit 'SQK-LSK110' "
        "--recursive "
        "--compress_fastq "
        "2> {log} || exit 1 ; "
        "cat output/nanopore/guppy_hac/pass/*.fastq.gz > {output.joined_fq}"

# config file - Data trimming:
#trim_strategy = dna
#trim_threshold = 2.5
#trim_min_events = 3
#min_qscore = 9.0
