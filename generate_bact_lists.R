setwd("/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel")

ancient<-read.delim("list_ancient",header=FALSE,sep="\t")
modern<-read.delim("list_modern",header=FALSE,sep="\t")

rand_vect_cont <- function(N, M, sd = 1) {
  vec <- abs(rnorm(N, M/N, sd))
  vec / sum(vec) * M
}

set.seed(123)
ground_truth_list<-list()
for(i in 1:10)
{
  ancient_present<-ancient[sample(1:dim(ancient)[1],10),]
  ancient_present$V2<-rand_vect_cont(10,1)
  ancient_absent<-ancient[!as.character(ancient$V1)%in%as.character(ancient_present$V1),]
  ancient_absent$V2<-0
  ancient_merged<-rbind(ancient_present,ancient_absent)
  
  modern_present<-modern[sample(1:dim(modern)[1],10),]
  modern_present$V2<-rand_vect_cont(10,1)
  modern_absent<-modern[!as.character(modern$V1)%in%as.character(modern_present$V1),]
  modern_absent$V2<-0
  modern_merged<-rbind(modern_present,modern_absent)
  
  ground_truth_list[[i]]<-rbind(ancient_merged,modern_merged)
  names(ground_truth_list[[i]])[2]<-paste0("Sample",i)
  print(paste0("Generated ",i," samples"))
}
ground_truth_matrix<-Reduce(function(dtf1,dtf2) merge(dtf1,dtf2,by="V1",all=TRUE),ground_truth_list)
rownames(ground_truth_matrix)<-ground_truth_matrix$V1
ground_truth_matrix$V1<-NULL
ground_truth_matrix
write.table(ground_truth_matrix,file="ground_truth_all_microbes.txt",col.names=FALSE,row.names=TRUE,quote=FALSE,sep="\t")

ground_truth_microbes<-rownames(ground_truth_matrix)
ground_truth_microbes<-gsub("\\.fa","",ground_truth_microbes)
ground_truth_microbes
write.table(ground_truth_microbes,file="ground_truth_microbes_list.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

ground_truth_matrix_ancient<-ground_truth_matrix[match(as.character(ancient$V1),rownames(ground_truth_matrix)),]
ground_truth_matrix_ancient
write.table(ground_truth_matrix_ancient,file="ground_truth_ancient_microbes.txt",col.names=FALSE,row.names=TRUE,quote=FALSE,sep="\t")
ground_truth_matrix_modern<-ground_truth_matrix[match(as.character(modern$V1),rownames(ground_truth_matrix)),]
ground_truth_matrix_modern
write.table(ground_truth_matrix_modern,file="ground_truth_modern_microbes.txt",col.names=FALSE,row.names=TRUE,quote=FALSE,sep="\t")


####################################### PLOT GROUND TRUTH AND COMPARE WITH KRAKENUNIQ ABUNDANCE MATRIX ###################################
library("pheatmap")
setwd("/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel")
ground_truth_matrix<-read.delim("ground_truth_all_microbes.txt",row.names=1,header=FALSE,sep="\t")
ground_truth_plot<-ground_truth_matrix
rownames(ground_truth_plot)<-gsub("\\.fa","",rownames(ground_truth_plot))
colnames(ground_truth_plot)<-paste0("Sample",seq(1,10,1))
pheatmap(ground_truth_plot,display_numbers=TRUE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,main="Ground truth: microbial species")

ground_truth_plot_binary<-ground_truth_plot
ground_truth_plot_binary[ground_truth_plot_binary>0.01]<-1
ground_truth_plot_binary[ground_truth_plot_binary<=0.01]<-0
pheatmap(ground_truth_plot_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Binary ground truth: microbial species")

krakenuniq<-read.delim("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/krakenuniq_abundance_matrix_10samples.txt",
                       header=TRUE,row.names=1,check.names=FALSE,sep="\t")
krakenuniq<-subset(krakenuniq,select=paste0("sample",seq(1,10,1)))
colnames(krakenuniq)<-paste0("Sample",seq(1,10,1))
krakenuniq_binary<-krakenuniq
krakenuniq_binary[krakenuniq_binary>0]<-1
#krakenuniq_binary[krakenuniq_binary<=10]<-0
#krakenuniq_binary[krakenuniq_binary>10]<-1
rownames(krakenuniq_binary)<-gsub(" ","_",rownames(krakenuniq_binary))
rownames(krakenuniq_binary)[rownames(krakenuniq_binary)=="Yersinia_pestis"]<-"Yersinia_pestis_CO92"
rownames(krakenuniq_binary)[rownames(krakenuniq_binary)=="Pseudomonas_sp._C27(2019)"]<-"Pseudomonas_caeni"
krakenuniq_binary<-krakenuniq_binary[match(rownames(ground_truth_plot_binary),rownames(krakenuniq_binary)),]
rownames(krakenuniq_binary)<-rownames(ground_truth_plot_binary)
krakenuniq_binary[is.na(krakenuniq_binary)]<-0
pheatmap(krakenuniq_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Binary KrakenUniq abundance: microbial species")


confusion_matrix_krakenuniq<-table(as.vector(t(ground_truth_plot_binary)),as.vector(t(krakenuniq_binary)))
confusion_matrix_krakenuniq
acc_krakenuniq<-sum(diag(confusion_matrix_krakenuniq))/sum(confusion_matrix_krakenuniq)
acc_krakenuniq
fisher.test(table(as.vector(t(ground_truth_plot_binary)),as.vector(t(krakenuniq_binary))))


hops<-read.delim("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/hops_abundance_matrix_10samples.txt",
                 header=TRUE,row.names=1,check.names=FALSE,sep="\t")
hops<-subset(hops,select=paste0("sample",seq(1,10,1),".trimmed"))
colnames(hops)<-paste0("Sample",seq(1,10,1))
hops_binary<-hops
hops_binary[hops_binary<=300]<-0
hops_binary[hops_binary>300]<-1
rownames(hops_binary)<-gsub(" ","_",rownames(hops_binary))
rownames(hops_binary)[rownames(hops_binary)=="Yersinia_pestis"]<-"Yersinia_pestis_CO92"
rownames(hops_binary)[rownames(hops_binary)=="Pseudomonas_sp._C27"]<-"Pseudomonas_caeni"
hops_binary<-hops_binary[match(rownames(ground_truth_plot_binary),rownames(hops_binary)),]
rownames(hops_binary)<-rownames(ground_truth_plot_binary)
hops_binary[is.na(hops_binary)]<-0
pheatmap(hops_binary,display_numbers=FALSE,fontsize=12,cluster_cols=FALSE,cluster_rows=FALSE,
         main="Binary HOPS abundance: microbial species")

confusion_matrix_hops<-table(as.vector(t(ground_truth_plot_binary)),as.vector(t(hops_binary)))
confusion_matrix_hops
acc_hops<-sum(diag(confusion_matrix_hops))/sum(confusion_matrix_hops)
acc_hops
fisher.test(table(as.vector(t(ground_truth_plot_binary)),as.vector(t(hops_binary))))


for(i in 1:length(rownames(ground_truth_plot_binary)))
{
  if(length(colnames(ground_truth_plot_binary[i,])[ground_truth_plot_binary[i,]==0 & krakenuniq_binary[i,]==0 & hops_binary[i,]==1])>0)
  {
    print(rownames(ground_truth_plot_binary)[i])
    print(colnames(ground_truth_plot_binary[i,])[ground_truth_plot_binary[i,]==0 & krakenuniq_binary[i,]==0 & hops_binary[i,]==1])
  }
}


for(i in 1:length(rownames(ground_truth_plot_binary)))
{
  if(length(colnames(ground_truth_plot_binary[i,])[ground_truth_plot_binary[i,]==0 & krakenuniq_binary[i,]==1 & hops_binary[i,]==0])>0)
  {
    print(rownames(ground_truth_plot_binary)[i])
    print(colnames(ground_truth_plot_binary[i,])[ground_truth_plot_binary[i,]==0 & krakenuniq_binary[i,]==1 & hops_binary[i,]==0])
  }
}


################################################### KRAKENUNIQ VS HOPS DETECTION ACCURACY ########################################################
acc_hops_vector<-vector()
cutoffs_krakenuniq<-seq(from=0,to=800,by=50)
for(i in cutoffs_krakenuniq)
{
  hops<-read.delim("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/hops_abundance_matrix_10samples.txt",
                   header=TRUE,row.names=1,check.names=FALSE,sep="\t")
  hops<-subset(hops,select=paste0("sample",seq(1,10,1),".trimmed"))
  colnames(hops)<-paste0("Sample",seq(1,10,1))
  hops_binary<-hops
  hops_binary[hops_binary<=i]<-0
  hops_binary[hops_binary>i]<-1
  rownames(hops_binary)<-gsub(" ","_",rownames(hops_binary))
  rownames(hops_binary)[rownames(hops_binary)=="Yersinia_pestis"]<-"Yersinia_pestis_CO92"
  rownames(hops_binary)[rownames(hops_binary)=="Pseudomonas_sp._C27"]<-"Pseudomonas_caeni"
  hops_binary<-hops_binary[match(rownames(ground_truth_plot_binary),rownames(hops_binary)),]
  rownames(hops_binary)<-rownames(ground_truth_plot_binary)
  hops_binary[is.na(hops_binary)]<-0
  
  confusion_matrix_hops<-table(as.vector(t(ground_truth_plot_binary)),as.vector(t(hops_binary)))
  acc_hops_vector<-append(acc_hops_vector,sum(diag(confusion_matrix_hops))/sum(confusion_matrix_hops))
}
plot(cutoffs_krakenuniq,rep(acc_krakenuniq,length(cutoffs_krakenuniq)),lty=2,ylim=c(0,1),type='l',
     xlab="READ NUMBER CUTOFF",ylab="ACCURACY",col='darkgreen',lwd=2,main="Overlap with ground truth")
lines(cutoffs_krakenuniq,acc_hops_vector,type='o',col='red')
legend("bottomright", legend=c("HOPS","KrakenUniq"),col=c("red","darkgreen"),
       lty=c(1,2),lwd=c(1,2),cex=1.5,inset=0.02)
