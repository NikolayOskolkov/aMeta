#SCRIPT FOR COMPUTING HOPS SCORE PER MICROBE

sample_num<-10

MaltExtract_output_path<-"/proj/uppstore2018095/private/NBIS_Demo/HOPS/results_max1000G_16threads/maltExtract"
RMA6<-paste0("sample",sample_num,".trimmed.rma6")

count_table<-read.delim("count_table.tsv",header=TRUE,row.names=1,check.names=FALSE,sep="\t")
microbe_list<-rownames(count_table)[count_table[,paste0("sample",sample_num,".trimmed")]>200]

microbe_list_gt<-readLines("/proj/uppstore2018095/private/NBIS_Demo/HOPS/HOPS_github_clone/Resources/custom_list_short_names_extended.txt")
microbe_list<-microbe_list[microbe_list%in%microbe_list_gt]

microbe_list<-gsub(" ","_",microbe_list)

for(i in microbe_list)
{

total_score<-0

#DAMAGE PATTERN
dam<-suppressWarnings(read.table(paste0(MaltExtract_output_path,"/default/damageMismatch/",RMA6,"_damageMismatch.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char=''))
if(round(dam[i,'C>T_1'],4)>0){total_score<-total_score+1}
if(round(dam[i,'G>A_20'],4)>0){total_score<-total_score+1}


#EVENNES OF COVERAGE
rd<-suppressWarnings(read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_alignmentDist.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char=''))
if(rd[i,'TotalAlignmentsOnReference']>200){total_score<-total_score+3}


#EDIT DISTANCE FOR ALL READS
df<-suppressWarnings(read.delim(paste0(MaltExtract_output_path,"/default/editDistance/",RMA6,"_editDistance.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t"))
df$higher<-NULL
if(as.numeric(df[i,1])>as.numeric(df[i,2]) & as.numeric(df[i,2])>as.numeric(df[i,3]) & as.numeric(df[i,3])>as.numeric(df[i,4]) & as.numeric(df[i,4])>as.numeric(df[i,5])){total_score<-total_score+1}

#EDIT DISTANCE FOR ANCIENT READS
df<-suppressWarnings(read.delim(paste0(MaltExtract_output_path,"/ancient/editDistance/",RMA6,"_editDistance.txt"),header=TRUE,check.names=FALSE,row.names=1,sep="\t"))
df$higher<-NULL
if(as.numeric(df[i,2])>as.numeric(df[i,3]) & as.numeric(df[i,3])>as.numeric(df[i,4]) & as.numeric(df[i,4])>as.numeric(df[i,5])){total_score<-total_score+1}


#PMD SCORES
if(round(dam[i,'C>T_1'],4)>0.05 & round(dam[i,'G>A_20'],4)>0.05){total_score<-total_score+1}


#READ LENGTH DISTRIBUTION
rd<-suppressWarnings(read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_readLengthDist.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char=''))
if(sum(as.numeric(rd[i,])[1:16])/(sum(as.numeric(rd[i,]))+1)>0.9){total_score<-total_score+1}


#DEPTH OF COVERAGE
rd<-suppressWarnings(read.table(paste0(MaltExtract_output_path,"/default/readDist/",RMA6,"_alignmentDist.txt"),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char=''))
if(rd[i,'TotalAlignmentsOnReference']>200){total_score<-total_score+1}



print(paste0(gsub("_"," ",i)," score: ",total_score))
}
