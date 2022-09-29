######################################### MAPPING MODERN METAGENOMIC READS TO Y.PESTIS ##############################################
# first, download the Y.pestis ref genome from here https://www.ncbi.nlm.nih.gov/genome/153?genome_assembly_id=299265
# second, download the modern infant stool metagenomic sample G69146_pe_1.fastq.gz from here https://diabimmune.broadinstitute.org/diabimmune/three-country-cohort/resources/metagenomic-sequence-data
module load bioinfo-tools cutadapt bowtie2 perl fastqc samtools
cd /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes
cutadapt -a CTGTCTCTTATA --minimum-length 30 -o G69146_pe_1.trimmed.fastq.gz G69146_pe_1.fastq.gz -j 4
# third, download RefSeq bacterial reference genomes
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
grep 'Complete Genome' assembly_summary.txt > assembly_summary_complete_latest_reference_genomes.txt
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > assembly_summary_complete_latest_reference_genomes_paths.txt

mkdir BacterialGenomes
for i in $(cat assembly_summary_complete_latest_reference_genomes_paths.txt)
do 
wget -P BacterialGenomes ${i}/$(basename ${i})_genomic.fna.gz
done

#remove Y.pestis RefSeq reference genomes from database
grep 'Yersinia pestis' assembly_summary_complete_latest_reference_genomes.txt | cut -f20 | cut -f10 -d '/' > Y.pestis_RefSeq_genomes.txt
for i in $(cat Y.pestis_RefSeq_genomes.txt)
do
echo $i
mv /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/BacterialGenomes/${i} /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/Y.pestis_RefSeq_genomes
done

#select 10, 100, 1000 and 10 000 random bacteria
cd BacterialGenomes
ls | shuf -n 10 > list_10_random_bacteria.txt
cat Y.pestis_plus_hg38.fasta.gz $(cat list_10_random_bacteria.txt | tr '\n' ' ') > Y.pestis_plus_hg38_plus_10_random_bacteria.fna.gz

ls | shuf -n 100 > list_100_random_bacteria.txt
cat Y.pestis_plus_hg38.fasta.gz $(cat list_100_random_bacteria.txt | tr '\n' ' ') > Y.pestis_plus_hg38_plus_100_random_bacteria.fna.gz

ls | shuf -n 1000 > list_1000_random_bacteria.txt
cat Y.pestis_plus_hg38.fasta.gz $(cat list_1000_random_bacteria.txt | tr '\n' ' ') > Y.pestis_plus_hg38_plus_1000_random_bacteria.fna.gz

ls | shuf -n 10000 > list_10000_random_bacteria.txt
cat Y.pestis_plus_hg38.fasta.gz $(cat list_10000_random_bacteria.txt | tr '\n' ' ') > Y.pestis_plus_hg38_plus_10000_random_bacteria.fna.gz

ls *fna.gz > /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/list_28896_RefSeq_bacteria.txt
cat Y.pestis_plus_hg38.fasta.gz $(cat /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/list_28896_RefSeq_bacteria.txt | tr '\n' ' ') > /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/Y.pestis_plus_hg38_plus_28896_RefSeq_bacteria.fna.gz


# Alignment to Y.pestis alone
samtools faidx Y.pestis_C092_ref_genome.fasta
bowtie2-build --large-index Y.pestis_C092_ref_genome.fasta Y.pestis_C092_ref_genome.fasta --threads 4
time bowtie2 --large-index -x /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/Y.pestis_C092_ref_genome.fasta --end-to-end --threads 4 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 4 - | samtools sort > G69146_pe_1.trimmed.fastq.gz.aligned_to_Y.pestis.bam
samtools index G69146_pe_1.trimmed.fastq.gz.aligned_to_Y.pestis.bam
samtools view G69146_pe_1.trimmed.fastq.gz.aligned_to_Y.pestis.bam NC_003143.1 | wc -l
# 22490

# Alignment to Y.pestis + hg38
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
cat Y.pestis_alone/Y.pestis_C092_ref_genome.fasta hg38.fasta > Y.pestis_plus_hg38.fasta
bowtie2-build --large-index Y.pestis_plus_hg38.fasta Y.pestis_plus_hg38.fasta --threads 4
time bowtie2 --large-index -x /proj/uppstore2018095/private/NBIS_Demo/misc/Y.pestis/genomes/Y.pestis_plus_hg38.fasta --end-to-end --threads 4 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 4 - | samtools sort > G69146_pe_1.trimmed.fastq.gz.aligned_to_Y.pestis_plus_hg38.bam
samtools index G69146_pe_1.trimmed.fastq.gz.aligned_to_Y.pestis_plus_hg38.bam
samtools view G69146_pe_1.trimmed.fastq.gz.aligned_to_Y.pestis_plus_hg38.bam NC_003143.1 | wc -l
# 22490

# Alignment to Y.pestis + hg38 + 10 random bacteria
bowtie2-build --large-index Y.pestis_plus_hg38_plus_10_random_bacteria.fna.gz Y.pestis_plus_hg38_plus_10_random_bacteria.fna.gz --threads 4
bowtie2 --large-index -x genomes/Y.pestis_plus_hg38_plus_10_random_bacteria/Y.pestis_plus_hg38_plus_10_random_bacteria.fna.gz --end-to-end --threads 4 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 4 - | samtools sort > G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_10_random_bacteria_bowtie2.bam
samtools index G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_10_random_bacteria_bowtie2.bam
samtools view G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_10_random_bacteria_bowtie2.bam NC_003143.1 | wc -l
# 8507

# Alignment to Y.pestis + hg38 + 100 random bacteria
bowtie2-build --large-index Y.pestis_plus_hg38_plus_100_random_bacteria.fna.gz Y.pestis_plus_hg38_plus_100_random_bacteria.fna.gz --threads 20
bowtie2 --large-index -x genomes/Y.pestis_plus_hg38_plus_100_random_bacteria/Y.pestis_plus_hg38_plus_100_random_bacteria.fna.gz --end-to-end --threads 20 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 4 - | samtools sort > G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_100_random_bacteria_bowtie2.bam
samtools index G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_100_random_bacteria_bowtie2.bam
samtools view G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_100_random_bacteria_bowtie2.bam NC_003143.1 | wc -l
# 1017

# Alignment to Y.pestis + hg38 + 1000 random bacteria
bowtie2-build --large-index Y.pestis_plus_hg38_plus_1000_random_bacteria.fna.gz Y.pestis_plus_hg38_plus_1000_random_bacteria.fna.gz --threads 20
bowtie2 --large-index -x genomes/Y.pestis_plus_hg38_plus_1000_random_bacteria/Y.pestis_plus_hg38_plus_1000_random_bacteria.fna.gz --end-to-end --threads 20 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 4 - | samtools sort > G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_1000_random_bacteria_bowtie2.bam
samtools index G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_1000_random_bacteria_bowtie2.bam
samtools view G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_1000_random_bacteria_bowtie2.bam NC_003143.1 | wc -l
# 197

# Alignment to Y.pestis + hg38 + 10000 random bacteria
bowtie2-build --large-index Y.pestis_plus_hg38_plus_10000_random_bacteria.fna.gz Y.pestis_plus_hg38_plus_10000_random_bacteria.fna.gz --threads 20
bowtie2 --large-index -x genomes/Y.pestis_plus_hg38_plus_10000_random_bacteria/Y.pestis_plus_hg38_plus_10000_random_bacteria.fna.gz --end-to-end --threads 20 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 20 - | samtools sort > G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_10000_random_bacteria_bowtie2.bam
samtools index G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_10000_random_bacteria_bowtie2.bam
samtools view G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_10000_random_bacteria_bowtie2.bam NC_003143.1 | wc -l
# 23


# Alignment to Y.pestis + hg38 + 28896 bacteria
bowtie2-build --large-index Y.pestis_plus_hg38_plus_28896_RefSeq_bacteria.fna.gz Y.pestis_plus_hg38_plus_28896_RefSeq_bacteria.fna.gz --threads 20
bowtie2 --large-index -x genomes/Y.pestis_plus_hg38_plus_28896_RefSeq_bacteria/Y.pestis_plus_hg38_plus_28896_RefSeq_bacteria.fna.gz --end-to-end --threads 20 --very-sensitive -U G69146_pe_1.trimmed.fastq.gz | samtools view -bS -q 1 -h -@ 20 - | samtools sort > G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_28896_RefSeq_bacteria_bowtie2.bam
samtools index G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_28896_RefSeq_bacteria_bowtie2.bam
samtools view G69146_pe_1.trimmed.fastq.gz_Aligned_to_Ypestis_plus_hg38_plus_28896_RefSeq_bacteria_bowtie2.bam NC_003143.1 | wc -l

