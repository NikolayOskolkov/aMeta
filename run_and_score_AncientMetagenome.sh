ml bioinfo-tools perl samtools java gnuparallel
cd /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow
conda activate ancient_microbiome_workflow
#snakemake -j 20 results/MULTIQC_BEFORE_TRIMMING/multiqc_report.html results/MULTIQC_AFTER_TRIMMING/multiqc_report.html results/KRAKENUNIQ/{simulation_s1,simulation_s2}/krakenuniq.output.filtered results/BOWTIE2/{simulation_s1,simulation_s2}/AlignedToPathogenome.bam results/MAPDAMAGE/{simulation_s1,simulation_s2} results/KRAKENUNIQ/{simulation_s1,simulation_s2}/taxonomy.krona.html results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix.txt results/MALT_DB/maltDB.dat results/MALT/{simulation_s1,simulation_s2}.trimmed.rma6 results/MALT_ABUNDANCE_MATRIX/malt_abundance_matrix.txt results/AUTHENTICATION/{simulation_s1,simulation_s2}
snakemake -j 20 results/MULTIQC_BEFORE_TRIMMING/multiqc_report.html results/MULTIQC_AFTER_TRIMMING/multiqc_report.html results/KRAKENUNIQ/{sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10}/krakenuniq.output.filtered results/BOWTIE2/{sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10}/AlignedToPathogenome.bam results/MAPDAMAGE/{sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10} results/KRAKENUNIQ/{sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10}/taxonomy.krona.html results/KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix.txt results/MALT_DB/maltDB.dat results/MALT/{sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10}.trimmed.rma6 results/MALT_ABUNDANCE_MATRIX/malt_abundance_matrix.txt results/AUTHENTICATION/{sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8,sample9,sample10}
conda deactivate
conda activate HOPS
/proj/uppstore2018095/private/NBIS_Demo/rma-tabuliser/rma-tabuliser/rma-tabuliser -d /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/results/MALT -r 'S' -n


#score AncientMetagenome
cd /home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome
for i in $(ls -d */); do /usr/bin/Rscript /home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/AncientMetagenome_Score.R simulation_s1.trimmed.rma6 ~/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/$i; done

cd /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/results/AUTHENTICATION/sample9
cp /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/results/MALT/sample9.trimmed.rma6 .
for i in $(ls -d */); do Rscript /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/scripts/AncientMetagenome_Score.R sample9.trimmed.rma6 /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/results/AUTHENTICATION/sample9/$i; done >> sample9_AncientMetagenome_scores.txt


