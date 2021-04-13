#####test the influence of blood samples
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggplot2)))

#### parameters of Monocle
option_list = list(
  make_option(c("-i", "--input_dir"), type="character",
              help="Input directory of expression matrix by 10xGenomics", metavar="character"),
  make_option(c("-r", "--runname"), type="character", default = "None",
              help="runname for this cell detection", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
);

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
da_path1 <- opt$input_dir       ## input information
sample_name <- opt$runname            ## runname
re_path <- opt$output_file     ## path and prefix of output files

#read in all of the tsne cluster information as a matrix contain 500 thousand cell and 10 parameters
print("start seprate cell detect...")
print("try to merge tsne information...")
tsne_filelist <- list.files(da_path1)
tsne_filelist <- tsne_filelist[intersect(which(str_count(tsne_filelist,".tab") != 0),which(str_count(tsne_filelist,"tsne_") != 0))]
write.table(tsne_filelist,paste0(re_path,"/filelist.txt"),sep="\t",quote=F)
cell_barcode <- rownames(read.table(paste0(da_path1,"/",tsne_filelist[1]),header=T,row.names=1,sep="\t"))
tsne_ori <- as.vector(read.table(paste0(da_path1,"/",tsne_filelist[1]),header=T,row.names=1,sep="\t")[,"ident"])
my_num <- length(table(tsne_ori))
my_para <- 1
for(i in 2:length(tsne_filelist)){
	temp <- as.vector(read.table(paste0(da_path1,"/",tsne_filelist[i]),header=T,row.names=1,sep="\t")[cell_barcode,"ident"])
	tsne_ori <- cbind(tsne_ori,temp)
	temp_num <- length(table(temp))
	if(temp_num < my_num){
		my_num <- temp_num
		my_para <- i
	}
	print(paste0("finish read and paste file:",tsne_filelist[i]))
}
print(paste0("and the standard file is:",tsne_filelist[my_para])) #"and the standard file is:pc.tsne_PC11.csv"
stand_cell <- tsne_ori[,my_para]
names(stand_cell) <- cell_barcode
tsne_ori <- data.frame(tsne_ori)
rownames(tsne_ori) <- cell_barcode
write.table(tsne_ori,paste0(re_path,"/",sample_name,"-tsne-ori.tab"),quote=F,sep="\t")

print("try to correspond cluster to common standard...")
tsne <- tsne_ori
for(j in 1:dim(tsne_ori)[2]){
	clusters <- as.vector(levels(tsne_ori[,j]))
	ori_ident <- tsne_ori[,j]
	for(e in 1:length(clusters)){
		sub_loc <- which(ori_ident==clusters[e])
		sub_cell <- rownames(tsne_ori)[sub_loc]
		cor_stand_clu <- stand_cell[sub_cell]
		stand_tab <- table(cor_stand_clu)
		cor_cluname <- names(which.max(stand_tab))
		tsne[sub_loc,j] <- cor_cluname
		print(paste0("transform from ",clusters[e]," to ",cor_cluname," in col ",j))
	}
}

clu_times <- rep(0,dim(tsne)[1])
clu_max <- rep(0,dim(tsne)[1])
for(i in 1:dim(tsne)[1]){
	temp <- table(as.vector(as.matrix(tsne[i,])))
	clu_times[i] <- length(temp)
	clu_max[i] <- max(temp)
}
print("finish cell cluster calculate!")
tsne$cluster_times <- clu_times
tsne$cluster_max <- clu_max
write.table(tsne,paste0(re_path,"/",sample_name,"-tsne-correspond.tab"),quote=F,sep="\t")

unstable_loc <- intersect(which(clu_times>1),which(clu_max<7)) #pc:44855/544905=0.0823
write.table(tsne[unstable_loc,],paste0(re_path,"/",sample_name,"-tsne-unstable.tab"),quote=F,sep="\t")






