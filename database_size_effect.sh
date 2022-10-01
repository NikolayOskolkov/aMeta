################################################## DATABASE SIZE EFFECT: SIMULATIONS ###############################################
cd /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims
cp /proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/results/CUTADAPT_ADAPTER_TRIMMING/* /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/data
cd /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/data
find $PWD -name '*.fastq.gz' > /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/sample.list

############################### RefSeq Microbial DB ##############################################
threads=10
DB=/crex/proj/uppstore2018095/private/db/Kraken2-DB/Complete_genome_DB/db
for i in 2 3 4 5 6 7 8 9 10
do
/proj/uppstore2018095/private/NBIS_Demo/KrakenUniq_v0.7.3/KrakenUniq/krakenuniq --report-file sample${i}.trimmed.fastq.gz_krakenuniq.output --fastq-input /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/data/sample${i}.trimmed.fastq.gz --db ${DB} --threads ${threads} --output sample${i}.trimmed.fastq.gz_classified_sequences.krakenuniq --only-classified-output --preload
Rscript filter_krakenuniq.R sample${i}.trimmed.fastq.gz_krakenuniq.output
done

############################### Kraken1 Standard DB ##############################################
threads=10
DB=/proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/temp/DBDIR_Kraken1_Standard_DB_indexed_for_KrakenUniq

for i in 1 2 3 4 5 6 7 8 9 10
do
/proj/uppstore2018095/private/NBIS_Demo/KrakenUniq_v0.7.3/KrakenUniq/krakenuniq --report-file sample${i}.trimmed.fastq.gz_krakenuniq.output_Kraken1_Standard --fastq-input /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/data/sample${i}.trimmed.fastq.gz --db ${DB} --threads ${threads} --output sample${i}.trimmed.fastq.gz_classified_sequences.krakenuniq_Kraken1_Standard --only-classified-output --preload
Rscript filter_krakenuniq.R sample${i}.trimmed.fastq.gz_krakenuniq.output_Kraken1_Standard
done

############################### Microbial NT DB ##############################################
threads=20
DB=/proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/temp/DBDIR_KrakenUniq_MicrobialNT

for i in 1 2 3 4 5 6 7 8 9 10
do
/proj/uppstore2018095/private/NBIS_Demo/KrakenUniq_v0.7.3/KrakenUniq/krakenuniq --report-file sample${i}.trimmed.fastq.gz_krakenuniq.output_Microbial_NT --fastq-input /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/data/sample${i}.trimmed.fastq.gz --db ${DB} --threads ${threads} --output sample${i}.trimmed.fastq.gz_classified_sequences.krakenuniq_Microbial_NT --only-classified-output --preload
Rscript filter_krakenuniq.R sample${i}.trimmed.fastq.gz_krakenuniq.output_Microbial_NT
done

############################### Full NT DB ##############################################
threads=10
DB=/proj/uppstore2018095/private/NBIS_Demo/DBDIR_KrakenUniq_Full_NT

for i in 1 2 3 4 5 6 7 8 9 10
do
/proj/uppstore2018095/private/NBIS_Demo/KrakenUniq_v0.7.3/KrakenUniq/krakenuniq --report-file sample${i}.trimmed.fastq.gz_krakenuniq.output_Full_NT --fastq-input /proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/data/sample${i}.trimmed.fastq.gz --db ${DB} --threads ${threads} --output sample${i}.trimmed.fastq.gz_classified_sequences.krakenuniq_Full_NT --only-classified-output --preload
Rscript filter_krakenuniq.R sample${i}.trimmed.fastq.gz_krakenuniq.output_Full_NT
done


