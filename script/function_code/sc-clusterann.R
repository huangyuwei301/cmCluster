#####test the influence of blood samples
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(rhdf5)))

#### parameters of Monocle
option_list = list(
  make_option(c("-i", "--input_dir"), type="character",
              help="Input directory of expression matrix by 10xGenomics", metavar="character"),
  make_option(c("-r", "--runname"), type="character", default = "None",
              help="runname for task of agree cell preparation", metavar="character"),
  make_option(c("-n", "--name"), type="character", 
              help="The name of input seurat RData", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
  
); 

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
input_info1 <- opt$input_dir       ## input information
f_prefix <- opt$runname
da_name <- opt$name
re_path <- opt$output_file     ## path and prefix of output files

cor_filename <- paste0(re_path,"/",f_prefix,"-tsne-correspond.tab")
tsne_all <- read.table(cor_filename,header=T,sep="\t")

agree_loc <- which(tsne_all[,"cluster_times"]==1)
agree_cell <- rownames(tsne_all)[agree_loc] #tsne_all[agree_loc,1]
write.table(tsne_all[agree_loc,],paste0(re_path,"/",f_prefix,".tsne.tab"),quote=F,sep="\t")

load(paste0(input_info1,"/",da_name))
clean_data <- FindVariableGenes(object = clean_data, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = FALSE)
clean_data1 <- SubsetData(clean_data, cells=agree_cell) #cells.use = agree_cell)
colnames(clean_data1@data) <- agree_cell
print(clean_data1)
rm(clean_data)
gc()
tryCatch({
	cat("load data",paste0(input_info1,"/",da_name),"and start convert data structure...","\n")
	anndata_fileobject <- Convert(clean_data1,to='anndata')
	cat("transform to anndata_fileobject!")
	},error=function(e){
		cat("ERROR :",conditionMessage(e),"\n")
		stop("Seurat Convert failed!\n",call.=FALSE)
		}
)
tryCatch({
	cat("start writing anndata data...","\n")
	re_filename = paste0(re_path,"/",f_prefix,".anndata_advance.h5")
	h5createFile(re_filename)
	},error=function(e){
		cat("ERROR :",conditionMessage(e),"\n")
		stop("Output H5 file create failed!\n",call.=FALSE)
		}
)
tryCatch({
	###build the h5 data storage structure
	#X : the expression matrix for genes and cells
	#it's a sparse matrix, so use special discription to decrease storage place
	h5createGroup(re_filename,"X")
	#obs : the annotation of observations 
	h5createGroup(re_filename,"obs")
	#var : the annotation of variable or features
	h5createGroup(re_filename,"var")
	#obsm : 
	h5createGroup(re_filename,"obsm")

	###get the corresponding data for the storage structure
	#data:value;indices:line numbers;indptr:column offset
	#shape:the dimension of matrix
	#feature:identify genes
	#barcode:identify cells
	print(attributes(attributes(clean_data1)))
	print(attributes(anndata_fileobject))
	anndata_fileobject
	indptr = anndata_fileobject$X["indptr"]
	indices = anndata_fileobject$X["indices"]
	data = anndata_fileobject$X["data"]
	shape = c(anndata_fileobject$X["shape"][[2]],anndata_fileobject$X["shape"][[1]])
	features = rownames(clean_data1@data) #rownames(anndata_fileobject$var)
	barcodes = colnames(clean_data1@data) #rownames(anndata_fileobject$obs)
	print("finish raw matrix load...")

	n_genes = anndata_fileobject$obs["n_genes"]
	n_counts = anndata_fileobject$obs["n_counts"]
	print("finish gene and count load...")
	orig_ident = anndata_fileobject$obs["orig_ident"]
	print("finish ident load...")
	percent_mito = anndata_fileobject$obs["percent_mito"]
	print("finish mito load...")
	#res_0_6 = anndata_fileobject$obs["res_0_6"]
	print("finish feature load...")

	gene_mean = anndata_fileobject$var["gene.mean"]
	print("finish var load...")
	#gene_dispersion = anndata_fileobject$var["gene.dispersion"]
	#gene_dispersion_scaled = anndata_fileobject$var["gene.dispersion.scaled"]

	X_pca = anndata_fileobject$obsm["X_pca"]
	print("finish pca load...")
	},error=function(e){
		cat("ERROR : can not find corresponding values in input data matrix :\n",conditionMessage(e),"\n","please check input data!")
		stop("H5 file structure build failed!\n",call.=FALSE)
		}
)
tryCatch({
	###write in the data according to the data structure
	h5write(indptr, re_filename,"X/indptr")
	h5write(indices,re_filename,"X/indices")
	h5createDataset(re_filename,"X/data",length(data),storage.mode="double",chunk=floor(length(data)/5),level=6)
	h5write(data,re_filename,"X/data")
	h5write(shape,re_filename,"X/shape")
	h5write(features,re_filename,"X/features")
	h5write(barcodes,re_filename,"X/barcodes")
	###added this two column to fit scanpy read in file structure
	h5write(features,re_filename,"X/gene_names")
	h5write(features,re_filename,"X/genes")

	h5write(n_genes,re_filename,"obs/n_genes")
	h5write(n_counts,re_filename,"obs/n_counts")
	h5write(orig_ident,re_filename,"obs/orig_ident")
	h5write(percent_mito,re_filename,"obs/percent_mito")
	#h5write(res_0_6,re_filename,"obs/res_0_6")

	h5write(gene_mean,re_filename,"var/gene_mean")
	#h5write(gene_dispersion,re_filename,"var/gene_dispersion")
	#h5write(gene_dispersion_scaled,re_filename,"var/gene_dispersion_scaled")

	h5write(X_pca,re_filename,"obsm/X_pca")

	cat("task done!","\n")
	},error=function(e){
		cat("ERROR :",conditionMessage(e),"\n")
		stop("write H5 file failed!\n",call.=FALSE)
		}
)
tryCatch({
	cat("show the structure of .h5 result file:","\n")
	#show the result file structure of output files
	h5ls(re_filename)
	},error=function(e){
		cat("ERROR :",conditionMessage(e),"\n")
		stop("",call.=FALSE)
		}
)
cat("create copy data to original RData...")
save(clean_data1, file=paste0(re_path,"/",f_prefix,".seurat.RData"))

tryCatch({
	cat("and create batch file...\n")
		#get batch information from cell barcode
		batch <- do.call(rbind,strsplit(clean_data1@cell.names,"-"))
		if(dim(batch)[2] > 1){
			cat("contain different batch and sample!\n")
			batch_name <- batch[,2]
			batch_file <- cbind(clean_data1@cell.names,batch[,2],batch_name)
		}else{
			cat("contain one batch and sample!\n")
			temp <- rep("1",dim(batch)[1])
			batch_name <- temp
			batch_file <- cbind(clean_data1@cell.names,temp,batch_name)
		}
		colnames(batch_file) <- c("rn","batch","batch_name")
		write.table(batch_file,paste0(re_path,"/",f_prefix,".group.txt"),row.names=F,sep="\t")
	},error=function(e){
		cat("ERROR:",conditionMessage(e),"\n")
		stop("write batch file failed!\n",call.=FALSE)
	}
)
print("finish data transform!")
