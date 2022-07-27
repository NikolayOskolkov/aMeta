ssh nikolay@rackham.uppmax.uu.se
cd /proj/uppstore2018095/private/NBIS_Demo/HOPS
conda activate HOPS
ml bioinfo-tools cutadapt java
for i in {1..10}
do
echo Working with sample${i}
cutadapt -a AGATCGGAAGAG --minimum-length 30 -o /proj/uppstore2018095/private/NBIS_Demo/HOPS/data/sample${i}.trimmed.fastq.gz /proj/uppstore2018095/private/NBIS_Demo/HOPS/data_untrimmed/sample${i}.fastq.gz
done
hops -Xmx1000G -input /proj/uppstore2018095/private/NBIS_Demo/HOPS/data/*trimmed.fastq.gz -output /proj/uppstore2018095/private/NBIS_Demo/HOPS/results -m full -c /proj/uppstore2018095/private/NBIS_Demo/HOPS/HOPS_Config_Default_HOPS.txt
/proj/uppstore2018095/private/NBIS_Demo/rma-tabuliser/rma-tabuliser/rma-tabuliser -d /proj/uppstore2018095/private/NBIS_Demo/HOPS/results/malt -r 'S' -n


#score HOPS
cd /home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome
/usr/bin/Rscript HOPS_Score.R


