import os
import yaml

# https://github.com/lh3/psmc
VCFUTILS="/home/krishang/software/bcftools/bcftools-1.9/bin/vcfutils.pl"
SAMTOOLS="/home/krishang/software/samtools/samtools-1.9/bin/samtools"
BCFTOOLS="/home/krishang/software/bcftools/bcftools-1.9/bin/bcftools"

PSMC_DIR="/home/krishang/software/psmc/psmc"
FQ2PSMCFA=os.path.join(PSMC_DIR,"utils/fq2psmcfa")
PSMC=os.path.join(PSMC_DIR,"psmc")
PSMC_PLOT=os.path.join(PSMC_DIR,"utils/psmc_plot.pl")

FA="/home/users/cindy/local/baboons/work/bams/GCA_000264685.2_Panu_3.0_genomic.fna"

MIN_MQ=30
MIN_BQ=30

## OUTMAIN = "results"
OUTMAIN="/home/users/cindy/local/baboons/work/bams/results/"

BEDS = {
    # "baserepeat": "/home/users/cindy/local/baboons/work/psmc/include.bed",
    "baserepeat": "/home/users/cindy/local/baboons/work/psmc/genome.bed"
}

allele_support = ["1", "2", "3"]

#configfile: "config_roger.yaml"
configfile: "config.yaml"

wildcard_constraints:
    t = "|".join(allele_support),
    sample = "|".join(config["samples"].keys()),
    bed = "|".join(BEDS.keys()),

rule all:
    input:
        expand(os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}.psmc"),
               sample=config["samples"].keys(),
               bed=BEDS.keys(),
               t=allele_support
        ),
        expand(os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.pdf"),
               bed=BEDS.keys(),
               t=allele_support
        ),
    shell:
        "ln -s -f {OUTMAIN} ."



rule gen_bcftools_genome_wide:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp(os.path.join(OUTMAIN, "psmc_input", "{sample}.bcf.gz"))
    params:
        names = "{sample}"
    threads: 1
    shell:
        """
{BCFTOOLS} mpileup  -B -Q {MIN_BQ}  -q {MIN_MQ} --threads {threads} -O u --fasta-ref {FA} -R regions.list --per-sample-mF -a FORMAT/AD,FORMAT/DP {input}  | {BCFTOOLS} call -Ob -o {output} --threads {threads} -c

"""

#    samtools mpileup -C50 -uf ref.fa aln.bam | bcftools view -c - \
#       | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz

# Here option -d sets and minimum read depth and -D sets the maximum. It is
# recommended to set -d to a third of the average depth and -D to twice.

rule gen_fq_mindepthx:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{sample}.bcf.gz"),
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.fq.gz"),
    params:
        mindepth=lambda wildcards: config["depths"][wildcards.sample][0], # 30/3
        maxdepth=lambda wildcards: config["depths"][wildcards.sample][1],  # 30*2
        B = lambda wildcards: BEDS[wildcards.bed]
    shell:
        """
        {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        """
        #{BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f /home/leopard/users/krishang/psmc/scripts/rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        
	#{BCFTOOLS} view -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i 'GT=="het" || GT=="hom"' |  {VCFUTILS} vcf2fq | gzip > {output}
        #{BCFTOOLS} view -i 'sum(INFO/DP)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i 'GT=="het" || GT=="hom"' | awk -f /home/leopard/users/krishang/psmc/scripts/rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        ## '(GT="hom" & FMT/DP<{wildcards.dp}) | (GT="het" & (FMT/AD[:0]<{wildcards.ad} | FMT/AD[:1]<{wildcards.ad} | FMT/DP<{wildcards.dp}))'

        # """
        # {BCFTOOLS} view -i 'sum(INFO/DP4)>={params.mindepth}' -T {params.B} -V mnps,indels -Ou {input} |  {BCFTOOLS} view -i '(GT=="het" && (INFO/DP4[0]+INFO/DP4[1])>={wildcards.t} && (INFO/DP4[2]+INFO/DP4[3])>={wildcards.t} ) || GT=="hom"' | awk -f scripts/rm_indels.awk |  {VCFUTILS} vcf2fq -d {params.mindepth} -D {params.maxdepth} | gzip > {output}
        # """


rule gen_psmcfa:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.fq.gz"),
    output:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    shell:
        """
        {FQ2PSMCFA} -q 20 {input} > {output}
        """

rule run_psmc:
    input:
        os.path.join(OUTMAIN, "psmc_input", "{bed}", "{t}", "{sample}.psmcfa"),
    output:
        os.path.join(OUTMAIN, "psmc_output", "{bed}", "{t}", "{sample}.psmc")
    shell:
        """{PSMC} -N25 -t15 -r5 -p "4+25*2+4+6" -o {output} {input}"""


rule plot:
    input:
        expand(os.path.join(OUTMAIN, "psmc_output", "{{bed}}", "{{t}}", "{sample}.psmc"),
               sample=config["samples"].keys())
    output:
        pdf=os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.pdf"),
        png=os.path.join(OUTMAIN, "psmc_output_figures", "{bed}_{t}.png")
    params:
        base = lambda wildcards, output: output.pdf[:-4],
        names = ",".join(config["samples"].keys())
    shell: """
    {PSMC_PLOT} -w 4 -T "{wildcards.bed} {wildcards.t}" -M {params.names} -p -P "right bottom" -g 11 -x 10000 -u 0.9e-8 {params.base} {input}
    pdftoppm  -png -r 300 -singlefile {output.pdf} {params.base}
    """
