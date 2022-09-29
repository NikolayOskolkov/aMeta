library("ggplot2")
df_plot<-data.frame(MAPPED_READS=c(22490,22490,8507,1017,197,23), GENOME_NUMBER=c(1,2,12,102,1002,10002))
head(df_plot)

base_breaks <- function(n = 10){function(x) {axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)}}


ggplot(df_plot,aes(x=GENOME_NUMBER,y=MAPPED_READS)) + geom_line(size=1.5, color = "blue")  + geom_point(size=5) + 
  theme(legend.position = "bottom", legend.text = element_text(size = 25), 
        legend.title = element_blank(), axis.text = element_text(size = 15), axis.title = element_text(size = 25)) +
  xlab("Number of reference genomes in database") + 
  ylab("Number of reads mapped uniquely to Y.pestis") + 
  scale_color_brewer(palette="Dark2") + scale_x_log10()  + scale_y_log10(breaks = base_breaks()) + 
  annotate("text", x = 1, y = 18000, label = "Y.pestis", size = 7) + 
  annotate("text", x = 5, y = 25000, label = "Y.pestis + human", size = 7) + 
  annotate("text", x = 90, y = 10000, label = "Y.pestis + human + 10 random bacteria", size = 7) + 
  annotate("text", x = 800, y = 1200, label = "Y.pestis + human + 100 random bacteria", size = 7) + 
  annotate("text", x = 100, y = 200, label = "Y.pestis + human + 1000 random bacteria", size = 7) + 
  annotate("text", x = 1000, y = 20, label = "Y.pestis + human + 10000 random bacteria", size = 7)
  
