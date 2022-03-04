shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")

import os

# Name of the config file
configfile: "config_alignv3.yaml"

# HOME FOLDER
home_folder=config['home_folder']
out_folder=config['OUT_DIR']
#FASTA REF
REF=config['REF_FASTA']

QC_DIR=out_folder+config['QC_DIR']
FASTQ_DIR_seprun=config['FASTQ_DIR_seprun']
SAM_DIR=out_folder+config['SAM_DIR']
BAM_DIR=out_folder+config['BAM_DIR']
BAM_DIR_vw=out_folder+config['BAM_DIR_view']
MPILEUP_DIR=out_folder+config['MPILEUP_DIR']
FASTA_DIR=out_folder+config['FASTA_DIR']

# Load the SAMPLES .json
FILES_SAMPLE = json.load(open(config['SAMPLES_JSON']))

# Extract the sammples for which you have the FASTQ
SAMPLES = sorted(FILES_SAMPLE.keys())


outrule=[]
for id in SAMPLES:
    for el in FILES_SAMPLE[id]:
            outrule.append(BAM_DIR+el+".sorted.bam")

run=[]
for id in SAMPLES:
    for el in FILES_SAMPLE[id]:
            run.append(el)


BAM=[]
for sample in SAMPLES:
    BAM.append(BAM_DIR+sample+".merge.sorted.markdup.bam")


FASTA=[]
for sample in SAMPLES:
    FASTA.append(FASTA_DIR+sample+"_"+config['angsd_chr2extract']+".fa")

TARGETS = []

TARGETS.append(FASTA)

# Define the containter
container:config['singularity_img']

# Define a rule that is local to all
localrules: all
rule all:
    input: TARGETS

# Map with bwa-mem 2 and sort with samtools

rule bwa_mem2_mem:
    input:
        reads=[FASTQ_DIR_seprun+"{el}R1.fastq.gz",FASTQ_DIR_seprun+"{el}R2.fastq.gz"]
    output:
        BAM_DIR+"{el}.sorted.bam"
    log:
       home_folder+"logs/bwa_mem2/{el}.log"
    params:
        index=REF,
        extra=r"-R '@RG\tID:{el}\tSM:{el}'",
        sort="samtools",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate", # Can be 'coordinate' (default) or 'queryname'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: 4
    wrapper:
        "v0.69.0/bio/bwa-mem2/mem"


# Function to get bams from the same sample which have to be merged

def get_bam2merge(wildcards):
    bam_2_merge=[]
    for i in FILES_SAMPLE[wildcards.sample].keys():
        bam_2_merge.append(BAM_DIR+FILES_SAMPLE[wildcards.sample][i]+".sorted.bam")
    return bam_2_merge

# Merge bams from the same sample

rule samtools_merge:
    input:
        get_bam2merge
    output:
        BAM_DIR+"{sample}.merge.sort.bam"
    params:
        "-cp" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        4    # This value - 1 will be sent to -@
    wrapper:
        "0.70.0/bio/samtools/merge"


# # This work as well..
# rule samtools_merge:
#     input:
#         expand([BAM_DIR+"{obj}.sorted.bam"], obj=run)
#     output:
#         BAM_DIR+"{sample}.merge.sort.bam"
#     params:
#         "" # optional additional parameters as string
#     threads:  # Samtools takes additional threads through its option -@
#         8     # This value - 1 will be sent to -@
#     wrapper:
#         "0.70.0/bio/samtools/merge"


# Mark duplicates with Picard

rule mark_duplicates:
    input:
        BAM_DIR+"{sample}.merge.sort.bam"
    output:
        bam=BAM_DIR+"{sample}.merge.sorted.markdup.bam",
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

# Filter out duplicated reads and min quality=10

rule samtools_view:
    input:
        BAM_DIR+"{sample}.merge.sorted.markdup.bam"
    output:
        BAM_DIR_vw+"{sample}.filtered.bam"
    params:
        "-b -q 10 -F 1292 -f 2" # optional params string
    wrapper:
        "v0.69.0/bio/samtools/view"

# index the bam file

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
        home_folder+config['angsd_env_file']+".yml"
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

# Get the consensus fasta from the bam file

rule generate_fasta:
    input:
        bam=BAM_DIR_vw+"{sample}.filtered_index.bam"
    output:
        FASTA_DIR+"{sample}_"+config['angsd_chr2extract']+".fa"
    conda:
        home_folder+config['angsd_env_file']+".yml"
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

sed -e "s/{params.chr2extract}/{wildcards.sample}/g" {params.fastaout_tmp}.fa > {output}
        """
