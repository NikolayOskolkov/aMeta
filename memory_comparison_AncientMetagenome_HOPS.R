##################################################### RESOURCES: 2 SAMPLES #############################################################
mem_HOPS<-c(0,69,129,243,564,711,711,711,711,213,100,0)
t_HOPS<-c(0,10,25,40,50,75,100,125,150,160,165,179)
plot(mem_HOPS~t_HOPS,type='o',col="red",xlab="TIME (MIN)",ylab="MEMORY USAGE (GB)",xaxt="n",
     main="Memory usage comparison between AncientMetagenome and HOPS",xlim=c(0,180),pch=19)
xtick<-seq(0, 180, by=10)
axis(side=1,at=xtick,labels=TRUE)


mem_AncientMetagenome<-c(0,110,110,150,110,110,5,5,50,10,10,50,50,0)
t_AncientMetagenome<-c(0,25,35,50,65,75,100,105,110,115,120,125,130,150)
lines(mem_AncientMetagenome~t_AncientMetagenome,type='o',col="blue",pch=19)



legend("topleft", legend=c("HOPS","AncientMetagenome"),col=c("red","blue"),lty=1,lwd=1,cex=1.5,inset=0.02)


##################################################### RESOURCES: 10 SAMPLES ############################################################
setwd("/home/nikolay/WABI/A_Gotherstrom/Manuscript/Method_Paper/HOPS_vs_AncientMetagenome/")

#CPU USAGE
df_cpu_AncientMetagenome<-as.numeric(readLines("cpu_usage_AncientMetagenome_10samples.txt"))
df_cpu_AncientMetagenome<-df_cpu_AncientMetagenome[1:562]
plot(df_cpu_AncientMetagenome,type="l",xlab="TIME (MIN)",ylab="CPU USAGE (%)",col="blue",
     main="CPU USAGE COMPARISON BETWEEN HOPS AND ANCIENT_METAGENOME",xlim=c(0,700))

df_cpu_HOPS<-as.numeric(readLines("cpu_usage_HOPS_10samples_max1000G_20threads.txt"))
df_cpu_HOPS<-df_cpu_HOPS-min(df_cpu_HOPS)
lines(df_cpu_HOPS,type="l",col="red",lwd=2)

df_cpu_HOPS_1thread<-as.numeric(readLines("cpu_usage_HOPS_10samples_max1000G_1thread.txt"))
df_cpu_HOPS_1thread<-df_cpu_HOPS_1thread-min(df_cpu_HOPS_1thread)
lines(df_cpu_HOPS_1thread,type="l",col="darkorange",lwd=2)

legend("topright", legend=c("HOPS 20 threads","HOPS 1 thread","AncientMetagenome 20 threads"),col=c("red","darkorange","blue"),
       lty=c(1,1,1),lwd=c(2,2,1),cex=1.5,inset=0.02)


#MEMORY USAGE
df_mem_AncientMetagenome<-as.numeric(readLines("memory_usage_AncientMetagenome_10samples.txt"))
df_mem_AncientMetagenome<-df_mem_AncientMetagenome[1:562]
df_mem_AncientMetagenome<-(df_mem_AncientMetagenome/100)*512
plot(df_mem_AncientMetagenome,type="l",xlab="TIME (MIN)",ylab="MEMORY (RAM) USAGE (GB)",col="blue",ylim=c(0,850),xlim=c(0,700),
     main="MEMORY USAGE COMPARISON BETWEEN HOPS AND ANCIENT_METAGENOME")

df_mem_HOPS<-as.numeric(readLines("memory_usage_HOPS_max1000G_10samples_20threads.txt"))
df_mem_HOPS<-(df_mem_HOPS/100)*1024
lines(df_mem_HOPS,type="l",col="red",lwd=2)

df_mem_HOPS_1thread<-as.numeric(readLines("memory_usage_HOPS_10samples_max1000G_1thread.txt"))
df_mem_HOPS_1thread<-(df_mem_HOPS_1thread/100)*1024
lines(df_mem_HOPS_1thread,type="l",col="darkorange",lwd=2)

legend("topright", legend=c("HOPS 20 threads","HOPS 1 thread","AncientMetagenome 20 threads"),col=c("red","darkorange","blue"),
       lty=c(1,1,1),lwd=c(2,2,1),cex=1.5,inset=0.02)



library("ggplot2")
df_plot<-data.frame(MEMORY=c(df_mem_AncientMetagenome,df_mem_HOPS,df_mem_HOPS_1thread),
                    TIME=c(seq(length(df_mem_AncientMetagenome)),seq(length(df_mem_HOPS)),seq(length(df_mem_HOPS_1thread))),
                    METHOD=c(rep("aMeta 20 threads",length(df_mem_AncientMetagenome)),
                             rep("HOPS 20 threads",length(df_mem_HOPS)),
                             rep("HOPS 1 thread",length(df_mem_HOPS_1thread))))
head(df_plot)

ggplot(df_plot,aes(x=TIME,y=MEMORY,color=METHOD)) + 
        geom_line(aes(color=METHOD,linetype=METHOD),size=1.5)  +
        scale_linetype_manual(values=c("solid","solid", "solid")) +
        theme(legend.position = "bottom", legend.text = element_text(size = 25), 
              legend.title = element_blank(), axis.text = element_text(size = 15), axis.title = element_text(size = 25)) +
        xlab("Time (min)") + 
        ylab("Memory (RAM) usage (GB)") + scale_color_brewer(palette="Dark2")
