setwd("/proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims")

###################################################### RefSeq Complete Genomes ########################################################################
IoU_vector<-vector()
for(i in seq(1,10,1))
{
#print(paste0("Working with sample ",i))
detected<-readLines(paste0("sample",i,".trimmed.fastq.gz_krakenuniq.output.species_detected"))
ground_truth<-readLines(paste0("/proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/sim_ground_truth_lists/sample",i,".trimmed.fastq.gz_ground_truth"))
IoU<-length(intersect(detected,ground_truth))/length(union(detected,ground_truth))
print(IoU)
IoU_vector<-append(IoU_vector,IoU)
#print(intersect(detected,ground_truth))
#print(union(detected,ground_truth))
}
IoU_vector
mean(IoU_vector)
#0.1674132
sd(IoU_vector)
#0.0483214


###################################################### Kraken1 Standard ########################################################################
IoU_vector<-vector()
for(i in seq(1,10,1))
{
#print(paste0("Working with sample ",i))
detected<-readLines(paste0("sample",i,".trimmed.fastq.gz_krakenuniq.output_Kraken1_Standard.species_detected"))
ground_truth<-readLines(paste0("/proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/sim_ground_truth_lists/sample",i,".trimmed.fastq.gz_ground_truth"))
IoU<-length(intersect(detected,ground_truth))/length(union(detected,ground_truth))
print(IoU)
IoU_vector<-append(IoU_vector,IoU)
#print(intersect(detected,ground_truth))
#print(union(detected,ground_truth))
}
IoU_vector
mean(IoU_vector)
#0.5156419
sd(IoU_vector)
#0.1287115


###################################################### Microbial NT ########################################################################
IoU_vector<-vector()
for(i in seq(1,10,1))
{
#print(paste0("Working with sample ",i))
detected<-readLines(paste0("sample",i,".trimmed.fastq.gz_krakenuniq.output_Microbial_NT.species_detected"))
ground_truth<-readLines(paste0("/proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/sim_ground_truth_lists/sample",i,".trimmed.fastq.gz_ground_truth"))
IoU<-length(intersect(detected,ground_truth))/length(union(detected,ground_truth))
print(IoU)
IoU_vector<-append(IoU_vector,IoU)
#print(intersect(detected,ground_truth))
#print(union(detected,ground_truth))
}
IoU_vector
mean(IoU_vector)
#0.7295603
sd(IoU_vector)
#0.2007677


###################################################### Full NT ########################################################################
IoU_vector<-vector()
for(i in seq(1,10,1))
{
#print(paste0("Working with sample ",i))
detected<-readLines(paste0("sample",i,".trimmed.fastq.gz_krakenuniq.output_Full_NT.species_detected"))
ground_truth<-readLines(paste0("/proj/uppstore2018095/private/NBIS_Demo/misc/DatabaseSizeEffect_Sims/sim_ground_truth_lists/sample",i,".trimmed.fastq.gz_ground_truth"))
IoU<-length(intersect(detected,ground_truth))/length(union(detected,ground_truth))
print(IoU)
IoU_vector<-append(IoU_vector,IoU)
#print(intersect(detected,ground_truth))
#print(union(detected,ground_truth))
}
IoU_vector
mean(IoU_vector)
#0.735958
sd(IoU_vector)
#0.1742338

################################################## PLOT DATABSE SIZE EFFECT #######################################################################
library("ggplot2")

df_plot<-data.frame(JACCARD = c(0.1674132, 0.5156419, 0.7295603, 0.735958), 
                    DATABASE_CHAR = c(65e+9, 78e+9, 110e+9, 233e+9), 
                    JACCARD_SD = c(0.0483214, 0.1287115, 0.2007677, 0.1742338), 
                    LABELS = c("Microbial NCBI RefSeq complete genomes","Kraken standard", 
                               "Microbial NCBI NT","Full NCBI NT"))
df_plot

ggplot(df_plot,aes(x=DATABASE_CHAR,y=JACCARD)) + geom_line(size=1.5, color = "blue")  + geom_point(size=5) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 25), 
        legend.title = element_blank(), axis.text = element_text(size = 15), axis.title = element_text(size = 25)) +
  xlab("Database size (characters)") + 
  ylab("Jaccard similarity") + 
  scale_color_brewer(palette="Dark2") + 
  geom_errorbar(aes(ymin=JACCARD-JACCARD_SD, ymax=JACCARD+JACCARD_SD)) + 
  ylim(c(0,1)) + xlim(c(0.5e+11,2.5e+11)) + 
  annotate("text", x = 1.15e+11, y = 0.18, label = "Microbial NCBI RefSeq complete genomes", size = 7) + 
  annotate("text", x = 0.98e+11, y = 0.48, label = "Kraken standard", size = 7) + 
  annotate("text", x = 1.32e+11, y = 0.77, label = "Microbial NCBI NT", size = 7) + 
  annotate("text", x = 2.15e+11, y = 0.77, label = "Full NCBI NT", size = 7)


