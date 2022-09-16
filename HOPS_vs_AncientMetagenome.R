setwd("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome")

#true_list<-readLines("/proj/uppstore2018095/private/NBIS_Demo/HOPS/HOPS_github_clone/Resources/custom_list_short_names.txt")
#true_list<-readLines("/proj/uppstore2018095/private/NBIS_Demo/HOPS/HOPS_github_clone/Resources/custom_list_short_names_extended.txt")

#df_hops<-read.delim("count_table.tsv",header=TRUE,row.names=1)
#df_am<-read.delim("/proj/sllstore2017093/b2015028/b2015028_nobackup/nikolay/ancient_microbiome_workflow/results/MALT/count_table.tsv",
#                  header=TRUE,row.names=1)

###################################################### HOPS CUTOFF ##############################################################
true_list<-readLines("custom_list_short_names_extended.txt")
df_hops<-read.delim("HOPS_count_table_name.tsv",header=TRUE,row.names=1)

df_am<-read.delim("AncientMetagenome_count_table_name.tsv",header=TRUE,row.names=1)
am_list<-rownames(df_am)
IoU_am<-length(intersect(true_list,am_list))/length(union(true_list,am_list))
TP_am<-sum(am_list%in%true_list)
FP_am<-sum(!am_list%in%true_list)
FN_am<-sum(!true_list%in%am_list)
F1_am<-(2*TP_am)/(2*TP_am+FP_am+FN_am)

F1_hops_vector<-vector()
IoU_hops_vector<-vector()
cutoffs_hops<-seq(from=0,to=800,by=50)
for(i in cutoffs_hops)
{
df_hops<-df_hops[rowMeans(df_hops)>i,]
hops_list<-rownames(df_hops)

IoU_hops<-length(intersect(true_list,hops_list))/length(union(true_list,hops_list))
IoU_hops_vector<-append(IoU_hops_vector,IoU_hops)

TP_hops<-sum(hops_list%in%true_list)
FP_hops<-sum(!hops_list%in%true_list)
FN_hops<-sum(!true_list%in%hops_list)
F1_hops<-(2*TP_hops)/(2*TP_hops+FP_hops+FN_hops)
F1_hops_vector<-append(F1_hops_vector,F1_hops)
}


##################################################### KRAKENUNIQ CUTOFF #########################################################
F1_krakenuniq_vector<-vector()
IoU_krakenuniq_vector<-vector()
cutoffs_krakenuniq<-seq(from=0,to=800,by=50)
for(i in cutoffs_krakenuniq)
{
  my_path<-"/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/krakenuniq_abundance_loop/"
  df_krakenuniq<-read.delim(paste0(my_path,"krakenuniq_abundance_matrix_",i,".txt"),header=TRUE,row.names=1)
  
  krakenuniq_list<-rownames(df_krakenuniq)
  
  IoU_krakenuniq<-length(intersect(true_list,krakenuniq_list))/length(union(true_list,krakenuniq_list))
  IoU_krakenuniq_vector<-append(IoU_krakenuniq_vector,IoU_krakenuniq)
  
  TP_krakenuniq<-sum(krakenuniq_list%in%true_list)
  FP_krakenuniq<-sum(!krakenuniq_list%in%true_list)
  FN_krakenuniq<-sum(!true_list%in%krakenuniq_list)
  F1_krakenuniq<-(2*TP_krakenuniq)/(2*TP_krakenuniq+FP_krakenuniq+FN_krakenuniq)
  F1_krakenuniq_vector<-append(F1_krakenuniq_vector,F1_krakenuniq)
}



#################################### PLOT HOPS VS KRAKENUNIQ SENSITIVITY VS SPECIFICITY ##########################################
par(mfrow=c(1,2))
plot(cutoffs_hops,rep(IoU_am,length(cutoffs_hops)),lty=2,ylim=c(0,1),type='l',
     xlab="READ NUMBER CUTOFF",ylab="INTERSECTION OVER UNION",col='darkgreen',lwd=2,main="Overlap with ground truth")
lines(cutoffs_hops,IoU_hops_vector,type='o',col='red')
lines(cutoffs_krakenuniq,IoU_krakenuniq_vector,type='o',col='blue')
legend("bottomright", legend=c("HOPS","KrakenUniq","aMeta"),col=c("red","blue","darkgreen"),
       lty=c(1,1,2),lwd=c(1,1,2),cex=1.5,inset=0.02)

plot(cutoffs_hops,rep(F1_am,length(cutoffs_hops)),lty=2,ylim=c(0,1),type='l',
     xlab="READ NUMBER CUTOFF",ylab="F1 SCORE",col='darkgreen',lwd=2,main="Microbial detection sensitivity vs. specificity")
lines(cutoffs_hops,F1_hops_vector,type='o',col='red')
lines(cutoffs_krakenuniq,F1_krakenuniq_vector,type='o',col='blue')
legend("bottomright", legend=c("HOPS","KrakenUniq","aMeta"),col=c("red","blue","darkgreen"),
       lty=c(1,1,2),lwd=c(1,1,2),cex=1.5,inset=0.02)

library("ggplot2")
library("ggpubr")
iou<-data.frame(IOU=c(IoU_krakenuniq_vector,IoU_hops_vector,IoU_am,IoU_am),CUTOFF=c(cutoffs_krakenuniq,cutoffs_hops,0,800),
                METHOD=c(rep("KrakenUniq",length(IoU_krakenuniq_vector)),rep("HOPS",length(IoU_hops_vector)),
                         rep("aMeta",2)))
f1<-data.frame(F1=c(F1_krakenuniq_vector,F1_hops_vector,F1_am,F1_am),CUTOFF=c(cutoffs_krakenuniq,cutoffs_hops,0,800),
               METHOD=c(rep("KrakenUniq",length(IoU_krakenuniq_vector)),rep("HOPS",length(IoU_hops_vector)),
                        rep("aMeta",2)))
iou_plot<-ggplot(iou,aes(x=CUTOFF,y=IOU,color=METHOD)) + 
  geom_line(aes(color=METHOD,linetype=METHOD),size=1.5) + geom_point(aes(color=METHOD,shape=METHOD),size=2.5) + 
  scale_linetype_manual(values=c("dashed","solid", "solid")) +
  scale_shape_manual(values=c(NA,19,19)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 25), 
        legend.title = element_blank(), axis.text = element_text(size = 15), axis.title = element_text(size = 25)) +
  xlab("Read number threshold") + 
  ylab("Jaccard similarity") + scale_color_brewer(palette="Dark2") +
  scale_y_continuous(breaks = round(seq(0, 0.9, by = 0.1),1))
f1_plot<-ggplot(f1,aes(x=CUTOFF,y=F1,color=METHOD)) + 
  geom_line(aes(color=METHOD,linetype=METHOD),size=1.5) + geom_point(aes(color=METHOD,shape=METHOD),size=2.5) +
  scale_linetype_manual(values=c("dashed","solid", "solid")) +
  scale_shape_manual(values=c(NA,19,19)) +
  theme(legend.position = "bottom", legend.text = element_text(size = 25), 
        legend.title = element_blank(), axis.text = element_text(size = 15), axis.title = element_text(size = 25)) +
  xlab("Read number threshold") + 
  ylab("F1 score") + scale_color_brewer(palette="Dark2") +
  scale_y_continuous(breaks = round(seq(0, 0.9, by = 0.1),1))
ggarrange(iou_plot,f1_plot,labels = c("A", "B"), ncol = 2, nrow = 1, 
          font.label = list(size = 30),common.legend = TRUE, legend = "bottom")


