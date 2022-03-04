angsd=/home/users/cindy/software/angsd/angsd
whitelist=include.bed
bams=bams.list
out=dstat_lowCov

$angsd -GL 2 -out $out -nThreads 20 -doAbbababa 1 -doCounts 1 -bam $bams -minmapQ 30 -minq 30 -sites $whitelist -blockSize 5000000
Rscript ~/software/angsd/R/jackKnife.R file=dstat_lowCov.abbababa indNames=baboon_samples.list outfile=lowCov_dstats


