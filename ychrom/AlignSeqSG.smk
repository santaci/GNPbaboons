# This one is originated from the AlignSeqV2.smk in the folder  /hpcnfs/scratch/EO/baboonproj/cindyproj/alignment

shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
# shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

import os

# For fancy configuration file also for the cluster look at https://github.com/crazyhottommy/pyflow-cellranger/blob/master/Snakefile
# Name of the config file
# 2 change following line / and script
configfile: "config_SG.yaml"

# HOME FOLDER
home_folder=config['home_folder']
out_folder=config['OUT_DIR']
#FASTA REF
REF=config['REF_FASTA']

# Load the link for the ouput
FASTQ_DIR=config['FASTQ_DIR']
BAM_DIR=out_folder+config['BAM_DIR']
BAM_DIR_vw=out_folder+config['BAM_DIR_view']
FASTA_DIR=out_folder+config['FASTA_DIR']

# Load the CLUSTER and SAMPLES .json
# CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES_SAMPLE = json.load(open(config['SAMPLES_JSON']))

# Extract the sammples for which you have the FASTQ
SAMPLES = sorted(FILES_SAMPLE.keys())

# Assign the folder for download log Snakemake creates the folder if they do not exist

FASTA=[]
for sample in SAMPLES:
    FASTA.append(FASTA_DIR+sample+"_"+config['angsd_chr2extract']+".fa")

# Empty list
TARGETS = []

# Append all the file that will be the output
TARGETS.append(FASTA)

# Define the containter
container:config['singularity_img']

# Define a rule that is local to all
localrules: all
rule all:
    input: TARGETS

# Make the rule for the mapping
rule bwa_mem2_mem:
    input:
        reads=[FASTQ_DIR+"{sample}_R1.fastq.gz", FASTQ_DIR+"{sample}_R2.fastq.gz"]
    output:
        BAM_DIR+"{sample}.sorted.bam"
    log:
       home_folder+"logs/bwa_mem2/{sample}.log"
    params:
        index=REF,
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'coordinate' (default) or 'queryname'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 4
    wrapper:
        "v0.69.0/bio/bwa-mem2/mem"

# Define the rule for the MarkDupicate with picard
rule mark_duplicates:
    input:
        BAM_DIR+"{sample}.sorted.bam"
    output:
        bam=BAM_DIR+"{sample}.sorted.markdup.bam",
        metrics=BAM_DIR+"dedup/{sample}.metrics.txt"
    log:
        home_folder+"logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=false"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=40000
    wrapper:
        "v0.69.0/bio/picard/markduplicates"

#Add a flagstat

rule samtools_view:
    input:
        BAM_DIR+"{sample}.sorted.markdup.bam"
    output:
        BAM_DIR_vw+"{sample}.filtered.bam"
    params:
        "-b -q 10 -F 1292 -f 2" # optional params string
    wrapper:
        "v0.69.0/bio/samtools/view"

# Add a flagstat

rule samtools_index:
    input:
        BAM_DIR_vw+"{sample}.filtered.bam"
    output:
        BAM_DIR_vw+"{sample}.filtered.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "v0.69.0/bio/samtools/index"


rule prepare_bam4fasta:
    input:
        bai=BAM_DIR_vw+"{sample}.filtered.bam.bai",
        bam=BAM_DIR_vw+"{sample}.filtered.bam"
    conda:
        config['angsd_env_file']+".yml"
    log:
        log1=home_folder+"logs/rsync/{sample}bam.log",
        log2=home_folder+"logs/rsync/{sample}bai.log"
    output:
        baio=BAM_DIR_vw+"{sample}.filtered_index.bam.bai",
        bamo=BAM_DIR_vw+"{sample}.filtered_index.bam"
    shell:
        """
rsync -v {input.bam} {output.bamo} > {log.log1}
rsync -v {input.bai} {output.baio} > {log.log2}
        """

rule generate_fasta:
    input:
        bam=BAM_DIR_vw+"{sample}.filtered_index.bam"
    output:
        FASTA_DIR+"{sample}_"+config['angsd_chr2extract']+".fa"
    conda:
        config['angsd_env_file']+".yml"
    params:
        minMapQ=config['angsd_minmapq'],
        minq=config['angsd_minq'],
        setMinDepth=config['angsd_mindepth'],
        fasta=REF,
        fastaout_tmp=FASTA_DIR+"{sample}_tmp",
        chr2extract=config['angsd_chr2extract']
    # log:
    shell:
        """
angsd \
-i {input.bam} \
-minMapQ {params.minMapQ} \
-minQ {params.minq} \
-doCounts 1 \
-setMinDepth {params.setMinDepth} \
-doFasta 3 \
-ref {params.fasta} \
-out {params.fastaout_tmp} \
-r {params.chr2extract}:

gzip -d {params.fastaout_tmp}.fa.gz

sed -e "s/V00654.1/{wildcards.sample}/g" {params.fastaout_tmp}.fa > {output}
        """
