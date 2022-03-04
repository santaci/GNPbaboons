#!/bin/bash

# This file is baboon

#SBATCH --partition=htc
#SBATCH --time=120:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=baboons
#SBATCH --mail-type=ALL
#SBATCH --mem=120000

#module load bwa
module load bcftools
module load java-jdk/13.0.1
#module load samtools/1.3.1

GATK=$DATA/gatk-4.1.8.1/gatk
REF=$DATA/baboon/X18057/papAnu_3/soft_masked/GCA_000264685.2_Panu_3.0_genomic.fna
NAME=bf186_paired


gatk_precall () {

#INDEL REALIGNMENT
#java -jar $GATK -T RealignerTargetCreator -R ${REF} -I ${NAME}.bam -o realigner_${NAME}.intervals
#java -jar $GATK -T IndelRealigner -R ${REF} -I ${NAME}.bam -targetIntervals realigner_${NAME}.intervals -o indel_realigned_${NAME}.bam
#java -jar $PICARD AddOrReplaceReadGroups I=indel_realigned_${NAME}.bam O=indel_${NAME}.bam SORT_ORDER=coordinate RGLB=${NAME} RGPL=Illumina RGSM=${NAME} RGPU=${NAME} RGCN=BCM
#mv indel_${NAME}.bam indel_realigned_${NAME}.bam
#samtools index indel_realigned_${NAME}.bam

#java -jar $GATK -T BaseRecalibrator -R ${REF} -I indel_realigned_${NAME}.bam -knownSites papAnu_3/actual_sorted.snps.filtered.removed.vcf -o recal_${NAME}.table
#java -jar $GATK -T PrintReads -R ${REF} -I indel_realigned_${NAME}.bam -BQSR recal_${NAME}.table -o recal_${NAME}.bam
#samtools index recal_${NAME}.bam
#java -jar $GATK -T BaseRecalibrator -R ${REF} -I indel_realigned_${NAME}.bam -knownSites papAnu_3/actual_sorted.snps.filtered.removed.vcf -BQSR recal_${NAME}.table -o after_recal_${NAME}.table
#java -jar $GATK -T AnalyzeCovariates -R ${REF} -before recal_${NAME}.table -after after_recal_${NAME}.table -plots recal_${NAME}_plots.pdf -csv plots_${NAME}.csv
#
}


#ACTUAL CALLING

#for name in $(cat bams.id); do

name=$1

#$GATK HaplotypeCaller -R $REF -I $DATA/baboon/X18057/rogers/${name}.bam -ERC GVCF -O $SCRATCH/${name}_raw.g.vcf 


#JOINT CALLING AND VARIANT FILTERING
#$GATK --java-options -Xmx200G GenomicsDBImport -R $REF \
#-L chr1 \
#-L chr2 \
#-L chr3 \
#-L chr4 \
#-L chr5 \
#-L chr6 \
#-L chr7 \
#-L chr8 \
#-L chr9 \
#-L chr10 \
#-L chr11 \
#-L chr12 \
#-L chr13 \
#-L chr14 \
#-L chr15 \
#-L chr16 \
#-L chr17 \
#-L chr18 \
#-L chr19 \
#-L chr20 \
#-V $SCRATCH/softmask_4_raw/realignedBam.LIV5_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.30977_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.30877_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.L142_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/34449.dups_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/34474.dups_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/34472.dups_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.16098_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.16066_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.97124_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.97074_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.28755_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.28697_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.28547_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/realignedBam.30388_raw.g.vcf \
#-V $SCRATCH/softmask_4_raw/bf186.panu3sm.g.vcf \
#-V $SCRATCH/softmask_4_raw/gelada.dups_raw.g.vcf \
#--genomicsdb-workspace-path $SCRATCH/autosome_db \
#--tmp-dir $SCRATCH/tmp

#$GATK --java-options -Xmx200G GenotypeGVCFs -R $REF \
#-V gendb://$SCRATCH/autosome_db \
#-L chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20 \
#-O $SCRATCH/moz_rogers_4.vcf

# VARIANT FILTERING
$GATK --java-options -Xmx200G SelectVariants -R $REF -V $SCRATCH/moz_rogers_4.vcf -select-type SNP -O $SCRATCH/moz_rogers_4.snp.vcf
$GATK --java-options -Xmx200G VariantFiltration -R $REF -V $SCRATCH/moz_rogers_4.snp.vcf -filter "QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum< -12.5 || ReadPosRankSum < -8.0" --filter-name "hard_snp_filter" -O $SCRATCH/moz_rogers_4.filtered.vcf
