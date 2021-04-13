#################################
#### scRNA-seq Analysis #########
#### Seurat2            #########
#### Date: 20191113     #########
#################################
## install cellranger R package
##source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
#### library
suppressWarnings(suppressMessages(library(Seurat)))
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(monocle)))

#### parameters of Monocle
option_list = list(
  make_option(c("-i", "--input_dir"), type="character",
              help="Input directory of expression matrix by 10xGenomics", metavar="character"),
  make_option(c("-r", "--runname"), type="character", default = "None",
              help="runname for disagree preparation and size factor calculate", metavar="character"),
  make_option(c("-n", "--name"), type="character", 
              help="The name of input seurat RData", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
  
); 

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
input_info1 <- opt$input_dir       ## input information
f_prefix <- opt$runname            ## Library id
da_name <- opt$name
re_path <- opt$output_file     ## path and prefix of output files

cor_filename <- paste0(re_path,"/",f_prefix,"-tsne-correspond.tab")
tsne_all <- read.table(cor_filename,header=T,sep="\t")
clu_times <- as.vector(tsne_all$cluster_times)

disagree_loc <- which(tsne_all[,"cluster_times"]!=1)
disagree_cell <- rownames(tsne_all)[disagree_loc] #tsne_all[disagree_loc,1]
print("finish read in and select disagree cells!")

load(paste0(input_info1,"/",da_name))
clean_data <- FindVariableGenes(object = clean_data, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = FALSE)
clean_data1 <- SubsetData(clean_data, cells=disagree_cell) #cells.use = disagree_cell)
clean_data1@raw.data <- clean_data1@raw.data[,rownames(clean_data1@meta.data)]  #[,disagree_cell]
print(clean_data1)
save(clean_data1, file=paste0(re_path,"/",f_prefix,".disagree.RData"))
print("finish save disagree RData!")
rm(clean_data)
gc()

#import the mydata1 seurat object
print("try to convert seurat data format into monocle...")
my_CDS1 <- importCDS(clean_data1,import_all=TRUE)
print("finish read in clean_data1!")

myCDS1 <- newCellDataSet(as(exprs(my_CDS1),"sparseMatrix"),phenoData = new("AnnotatedDataFrame", data = pData(my_CDS1)),featureData = new("AnnotatedDataFrame", data = fData(my_CDS1)),lowerDetectionLimit = 0.5,expressionFamily=negbinomial.size())
print("create cell data!")

#estimate size factors and dispersions of the cell expression data set
print("try to estimate size factor...")
myCDS1 <- estimateSizeFactors(myCDS1)

print("try to save sizefactor result...")
raw_sizefactor <- sizeFactors(myCDS1)
save(raw_sizefactor, file=paste0(re_path,"/",f_prefix,".sizefactor.RData"))

