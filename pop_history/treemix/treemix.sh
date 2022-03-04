for m in {1..5}; do
for i in {1..100}
do

treemix -i gatk_baboons.treemix.frq.gz \
-m $m \
-o ../treemix/gatk_baboons.${m}.${i} \
-root Gelada -bootstrap \
-global \
-k 500
#grep 'Exiting' ../treemix/gatk_baboons.${i}.llik >> ../treemix/mig.${i}.llik
#grep 'SEED' ../treemix/gatk_${i}_log >> ../treemix/mig.${i}.llik

done &
done
