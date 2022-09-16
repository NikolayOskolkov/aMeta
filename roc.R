setwd("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/")
library("ROCit")

###################################################### SCORES 2 SAMPLES ####################################################################
df_2<-read.delim("AncientMetagenome_scores_2samples.txt",header=FALSE,sep="\t")
df_2<-na.omit(df_2)
df_2

roc_obj_AncientMetagenome_2<-rocit(as.numeric(df_2$V2),as.numeric(df_2$V4))
roc_obj_HOPS_2<-rocit(as.numeric(df_2$V3),as.numeric(df_2$V4))

par(cex.axis=1.2,cex.lab=1.5)
par(mar=c(4,5,1,3),mgp=c(2.2,1,0))

plot(roc_obj_AncientMetagenome_2$TPR~roc_obj_AncientMetagenome_2$FPR, col="blue",lwd=2,xlab="1-Specificity (FPR)",ylab="Sensitivity (TPR)",type='o')
lines(roc_obj_HOPS_2$TPR~roc_obj_HOPS_2$FPR, col="red",lwd=2,type='o')
lines(c(0,1),c(0,1),lwd=2,lty=2,col="darkgreen")
legend("bottomright",legend=c("HOPS","AncientMetagenome", "Random"),col=c("red","blue","darkgreen"), lty=c(1,1,2),lwd=2,inset=0.02,cex=1.5)


###################################################### SCORES 10 SAMPLES ####################################################################
df_10<-read.delim("AncientMetagenome_scores_10samples.txt",header=FALSE,sep="\t")
df_10<-na.omit(df_10)
table(df_10$V2)
df_10

df_10_hops<-read.delim("HOPS_scores_10samples.txt",header=FALSE,sep="\t")
df_10_hops<-na.omit(df_10_hops)
table(df_10_hops$V2)
df_10_hops

roc_obj_AncientMetagenome_10<-rocit(as.numeric(df_10$V2),as.numeric(df_10$V3),method="binormal")
roc_obj_HOPS_10<-rocit(as.numeric(df_10_hops$V2),as.numeric(df_10_hops$V3),method="binormal")
#roc_obj_HOPS_10$TPR[roc_obj_HOPS_10$TPR<roc_obj_HOPS_10$FPR]<-roc_obj_HOPS_10$FPR[roc_obj_HOPS_10$TPR<roc_obj_HOPS_10$FPR]

par(mfrow=c(1,1))
par(cex.axis=1.2,cex.lab=1.5)
par(mar=c(4,5,1,3),mgp=c(2.2,1,0))

plot(roc_obj_AncientMetagenome_10$TPR~roc_obj_AncientMetagenome_10$FPR, col="blue",lwd=2,xlab="1-Specificity (FPR)",ylab="Sensitivity (TPR)",type='l')
lines(roc_obj_HOPS_10$TPR~roc_obj_HOPS_10$FPR, col="red",lwd=2,type='l')
lines(c(0,1),c(0,1),lwd=3,lty=2,col="darkgreen")
legend("bottomright",legend=c(paste0("HOPS: ROC AUC = ",round(roc_obj_HOPS_10$AUC,2)),
                              paste0("AncientMetagenome: ROC AUC = ",round(roc_obj_AncientMetagenome_10$AUC,2)), "Random"),
       col=c("red","blue","darkgreen"), lty=c(1,1,2),lwd=2,inset=0.02,cex=1.5)

library("ggplot2")
df2plot<-data.frame(TPR=c(roc_obj_AncientMetagenome_10$TPR,roc_obj_HOPS_10$TPR,0,1),
                    FPR=c(roc_obj_AncientMetagenome_10$FPR,roc_obj_HOPS_10$FPR,0,1),
                    METHOD=c(rep("aMeta",length(roc_obj_AncientMetagenome_10$TPR)),
                             rep("HOPS",length(roc_obj_HOPS_10$TPR)),
                             rep("Random",2)))
my_color <- RColorBrewer::brewer.pal(8, "Dark2")
ggplot(data=df2plot,aes(y=TPR,x=FPR)) + geom_line(aes(color=METHOD,linetype=METHOD), size=1.5) + scale_color_brewer(palette="Dark2") + 
  scale_linetype_manual(values=c("solid", "solid","dashed")) + 
  geom_text(x=0.75, y=0.10, label="aMeta ROC AUC = 0.94", size = 8, color = my_color[1]) + 
  geom_text(x=0.75, y=0.05, label="HOPS ROC AUC = 0.81", size = 8, color = my_color[2]) + 
  ylab("Sensitivity (TPR)") + 
  xlab("1-Specificity (FPR)") + 
  theme(legend.position = "bottom", legend.text = element_text(size = 25), legend.title = element_blank(), 
                                      axis.text = element_text(size = 15), axis.title = element_text(size = 25))


######################################################### GROUND TRUTH #####################################################################
gt_ancient<-read.delim("/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/ground_truth_ancient_microbes.txt",
                       header=FALSE,row.names=1,check.names=FALSE,sep="\t")
gt_ancient_list<-gsub(" CO92","",gsub("_"," ",gsub("\\.fa","",rownames(gt_ancient))))
gt_ancient_df<-data.frame(MICROBE=gt_ancient_list,GROUND_TRUTH=1)
gt_modern<-read.delim("/home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel/ground_truth_modern_microbes.txt",
                       header=FALSE,row.names=1,check.names=FALSE,sep="\t")
gt_modern_list<-gsub(" CO92","",gsub("_"," ",gsub("\\.fa","",rownames(gt_modern))))
gt_modern_df<-data.frame(MICROBE=gt_modern_list,GROUND_TRUTH=0)

gt_df<-rbind(gt_ancient_df,gt_modern_df)
gt_df

hops_scores_10samples<-read.delim("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/HOPS_scores_10samples1.txt",
                                  header=FALSE,sep="\t")
hops_scores_10samples$V3<-gt_df$GROUND_TRUTH[match(as.character(hops_scores_10samples$V1),as.character(gt_df$MICROBE))]
head(hops_scores_10samples)
write.table(hops_scores_10samples,file="HOPS_scores_10samples.txt",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
