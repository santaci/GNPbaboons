#!/bin/bash

# This file is dstat

#SBATCH --partition=htc
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=2
#SBATCH --job-name=dstat
#SBATCH --mail-type=ALL

# use submission environment
#module load python/anaconda2/4.2.0

#for ((i=1,j=8;i<=14;i+=8,j=i+7)); do
#for ((i=1,j=21;i<=210;i+=21,j=i+20)); do

for ((i=1;i<=6;i++)); do
while read pop; do

# Individual Baboons
#python /home/scro2664/bin/popstats.py --file babs --pops $pop --informative >> dstat_indivs_C.out 

# Individual Baboons
#python2 ~/bin/popstats.py -p babs.tped -f babs_pops.tfam --pops $pop --informative >> dstat_all_pops.out

# Variable B in D calculation = D(Gelada,B,Moz,Ursinus)
#python2 ~/bin/popstats.py -p babs.tped -f babs_pops.tfam --pops $pop --informative >> dstat_var_B_mozur.out

#Other variable B in D calculations
#python2 ~/bin/popstats.py -p babs.tped -f babs.tfam --pops $pop --informative >> dstat_var_B_others.out

#H3 is variable and H1 is Ursinus indivs
#python2 ~/bin/popstats.py -p babs.tped -f babs.tfam --pops $pop --informative >> dstat_h3_var_ursinus.out

#H3 is variable and H1 is Ursinus indivs - USING BCFTOOLS NO REPEATMASKER
#python2 ~/bin/popstats.py -p /steveData/cindy/autosome/vcfs/bcfmerge_baboons_d.tped -f /steveData/cindy/autosome/vcfs/bcfmerge_baboons_dclust.fam --pops $pop --informative >> /steveData/cindy/autosome/vcfs/bcf_dstat_h3_var_ursinus.out

#H3 is variable H1 is Ursinus indivs - USING GATK - X CHROMOSOME
#python2 ~/bin/popstats.py --not23 -p /steveData/cindy/autosome/vcfs/gatk/moz_rogers_4_chrX.tped -f /steveData/cindy/autosome/vcfs/gatk/moz_rogers_snps4.filtered_d.tfam --pops $pop --informative >> /steveData/cindy/autosome/vcfs/gatk/x_gatk_dstat_h3_var_ursinus.out

# H3 is variable, H1 is union Ursinus - USING GATK - X CHROMOSOME & AUTOSOME
#python2 ~/bin/popstats.py --not23 -p /steveData/cindy/autosome/vcfs/gatk/moz_rogers_4_chrX.tped -f /steveData/cindy/autosome/vcfs/gatk/moz_rogers_snps4.filtered.fam --pops $pop --informative >> /steveData/cindy/autosome/vcfs/gatk/union_x_gatk_dstat_h3_var_ursinus.out
#python2 ~/bin/popstats.py --not23 -p /steveData/cindy/autosome/vcfs/gatk/moz_rogers_snps4.filtered_d.tped -f /steveData/cindy/autosome/vcfs/gatk/moz_rogers_snps4.filtered.fam --pops $pop --informative >> /steveData/cindy/autosome/vcfs/gatk/union_gatk_dstat_h3_var_ursinus.out


#H3 is variable H1 is Anubis indivs (not 30877)
python2 ~/bin/popstats.py --not23 -p /steveData/cindy/autosome/vcfs/gatk/moz_rogers_snps4.filtered_d.tped -f /steveData/cindy/autosome/vcfs/gatk/moz_rogers_anubis.tfam --pops $pop --informative >> /steveData/cindy/autosome/vcfs/gatk/gatk_dstat_h3_L142_anubis.out
#python2 ~/bin/popstats.py --not23 -p /steveData/cindy/autosome/vcfs/gatk/moz_rogers_4_chrX.tped -f /steveData/cindy/autosome/vcfs/gatk/moz_rogers_anubis.tfam --pops $pop --informative >> /steveData/cindy/autosome/vcfs/gatk/x_gatk_dstat_h3_var_anubis.out


#done < <(cat h3_var_ursinus_union.list | sed -n ${i}p) &

#done < <(cat anu_var_ursinus.list | sed -n ${i}p) &
done < <(cat anu_cont_ursinus.list | sed -n ${i}p) &

#done < <(cat h3_var_ursinus.list | sed -n ${i}p) &

#done < <(cat other_var_B.list | sed -n ${i}p) &


#done < <(cat vary_B_mozur.list | sed -n ${i}p) &


#done < <(cat baboons_comp_pop.list | sed -n ${i},${j}p) &

done

wait
