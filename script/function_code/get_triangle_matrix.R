#####test the influence of blood samples
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggplot2)))

#### parameters of Monocle
option_list = list(
  make_option(c("-i", "--input_dir"), type="character",
              help="Input directory of expression matrix by 10xGenomics", metavar="character"),
  make_option(c("-r", "--runname"), type="character", 
              help="Runname for cmCluster heatmap", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
); 

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
input_info1 <- opt$input_dir       ## input information
f_prefix <- opt$runname            ## runname for cmCluster heatmap
outputpath <- opt$output_file     ## path and prefix of output files

da_path <- input_info1
setwd(da_path)

tsne_filelist <- list.files(da_path)
tsne_filelist <- tsne_filelist[which(str_count(tsne_filelist,"evaluation.tab") != 0)]

temp <- read.table(tsne_filelist[1],header=T,sep="\t")

all_para <- temp[,c("macro_Precision","macro_Recall","macro_F1.score","macro_TP","macro_TN")]

for(i in seq(2,length(tsne_filelist))){
	cat("start",i,"time of merge file...\n")
	temp <- read.table(tsne_filelist[i],header=T,sep="\t")
	temp_para <- temp[,c("macro_Precision","macro_Recall","macro_F1.score","macro_TP","macro_TN")]
	cat("get information from",tsne_filelist[i],"success!\n")
	name_para <- union(rownames(all_para),rownames(temp))
	if(dim(all_para)[1] < dim(temp_para)[1]){
		na_sub_para <- setdiff(name_para,rownames(all_para))
		cat("choose 1:add old parameter list from",dim(all_para)[1],"to match new of",length(name_para),"\n")
		na_sub <- matrix(nrow=length(na_sub_para),ncol=dim(all_para)[2])
		rownames(na_sub) <- na_sub_para
		colnames(na_sub) <- colnames(all_para)
		all_para <- rbind(all_para,na_sub)
	}else{
		na_sub_para <- setdiff(name_para,rownames(temp_para))
		cat("choose 2:add new parameter list from",dim(temp_para)[1],"to match old of",length(name_para),"\n")
		na_sub <- matrix(nrow=length(na_sub_para),ncol=dim(temp_para)[2])
		rownames(na_sub) <- na_sub_para
		colnames(na_sub) <- colnames(temp_para)
		temp_para <- rbind(temp_para,na_sub)
	}
	all_para <- cbind(all_para,temp_para[rownames(all_para),])
	cat("finish merge file!\n")
	cat("\n")
}

paralist <- read.table(paste0(outputpath,"/filelist.txt"),sep="\t",header=T)
rownames(paralist) <- rownames(read.table(paste0("1",f_prefix,".evaluation.tab"),sep="\t",header=T))
rownames(all_para) <- paralist[rownames(all_para),1]

namelist <- colnames(all_para)
all_pre <- all_para[,which(str_count(namelist,"Precision") != 0)]
all_rec <- all_para[,which(str_count(namelist,"Recall") != 0)]
all_f1 <- all_para[,which(str_count(namelist,"F1") != 0)]
all_tp <- all_para[,which(str_count(namelist,"TP") != 0)]
all_fp <- all_para[,which(str_count(namelist,"TN") != 0)]
colnames(all_pre) <- tsne_filelist
colnames(all_rec) <- tsne_filelist
colnames(all_f1) <- tsne_filelist
colnames(all_tp) <- tsne_filelist
colnames(all_fp) <- tsne_filelist
### sort by columns
no_nas <- vector()
for(j in 1:dim(all_pre)[2]){
	#print(j)
	no_nas <- c(no_nas,length(which(all_pre[,j]>=0)))
}
no_nas <- rank(no_nas)
names(no_nas) <- colnames(all_pre)
all_pre <- all_pre[,names(sort(no_nas,decreasing = T))]
all_rec <- all_rec[,names(sort(no_nas,decreasing = T))]
all_f1 <- all_f1[,names(sort(no_nas,decreasing = T))]
all_tp <- all_tp[,names(sort(no_nas,decreasing = T))]
all_fp <- all_fp[,names(sort(no_nas,decreasing = T))]
### sort by rows
no_nas <- vector()
for(j in 1:dim(all_pre)[1]){
	#print(j)
	no_nas <- c(no_nas,length(which(all_pre[j,]>=0)))
}
no_nas <- rank(no_nas)
names(no_nas) <- rownames(all_pre)
all_pre <- all_pre[names(sort(no_nas,decreasing = T)),]
all_rec <- all_rec[names(sort(no_nas,decreasing = T)),]
all_f1 <- all_f1[names(sort(no_nas,decreasing = T)),]
all_tp <- all_tp[names(sort(no_nas,decreasing = T)),]
all_fp <- all_fp[names(sort(no_nas,decreasing = T)),]
print("finish parameter get!")
write.table(all_pre,paste0(outputpath,"/pre.txt"),sep="\t")
write.table(all_rec,paste0(outputpath,"/rec.txt"),sep="\t")
write.table(all_f1,paste0(outputpath,"/f1.txt"),sep="\t")
write.table(all_tp,paste0(outputpath,"/tp.txt"),sep="\t")
write.table(all_fp,paste0(outputpath,"/fp.txt"),sep="\t")
##

print("deal with the f1-information to heatmap...")
f1_ori <- all_f1
for(i in 1:dim(f1_ori)[1]){ f1_ori[i,which(is.na(f1_ori[i,])==T)] <- 0}
para_com <- vector()
for(e in rownames(f1_ori)){
	temp <- unlist(strsplit(e,split="\\."))
	the_para <- unlist(strsplit(temp[2],split="_"))[2]
	the_para <- gsub("pc","-",the_para)
	the_para <- gsub("k","-",the_para)
	the_para <- gsub("r","-",the_para)
	if(temp[3] == "tab"){
		para_com <- rbind(para_com,c(e,unlist(strsplit(the_para,split="-"))[2:4],"0"))
	}else{
		para_com <- rbind(para_com,c(e,unlist(strsplit(the_para,split="-"))[2:4],temp[3]))
	}
}
colnames(para_com) <- c("file_name","pc","k","r1","r2")
rownames(para_com) <- para_com[,1]
para_com <- data.frame(para_com)
para_com$r <- as.numeric(paste(para_com$r1,para_com$r2,sep="."))
para_com <- para_com[rownames(f1_ori),]

library(pheatmap)
colnames(f1_ori) <- seq(1,dim(f1_ori)[2])
rownames(f1_ori) <- paste0("pc",para_com[,"pc"],"_k",para_com[,"k"],"_r",para_com[,"r"])

annotation_row <- para_com[,c("pc","k","r")]
rownames(annotation_row) <- rownames(f1_ori)

pdf(file=paste0(outputpath,"/heatmap_parameter.pdf"))
pheatmap(f1_ori,color = colorRampPalette(c("white", "firebrick3"))(50),
	cluster_row = FALSE,cluster_col = FALSE,
	show_colnames = FALSE,show_rownames = FALSE,
	annotation_row =annotation_row)
dev.off()
###

plot_info <- cbind(para_com[,c("pc","k","r")],f1_ori)
write.table(plot_info,paste0(outputpath,"/heatmap_info_parameter.csv"),sep=",",quote=F)