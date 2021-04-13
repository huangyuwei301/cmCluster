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

#### parameters of seurat
option_list = list(
  make_option(c("-i", "--input_dir"), type="character",
              help="Input directory of expression matrix by 10xGenomics", metavar="character"),
  make_option(c("-s", "--species"), type="character", default = "human",
              help="Library id for multiple GEMs when 'cellranger aggr' is used.", metavar="character"),
  make_option(c("-p", "--parameter_list"), type="character", default = "None",
              help="parameter during seurat cell filter and clustering", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
); 

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
input_info <- opt$input_dir       ## input information
outputpath <- opt$output_file     ## path and prefix of output files

sam_species <- opt$species
my_para <- opt$parameter_list

#####################################################################################
###### Load the dataset
#####################################################################################
input.data <- Read10X(data.dir = input_info)
clean_data <- CreateSeuratObject(raw.data = input.data, min.cells = 3, min.genes = 200, project = "10X")
print("finish read data!")

### QC and selecting cells for further analysis
if(sam_species == "human"){
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = clean_data@data), ignore.case=TRUE, value = TRUE)
}else{
	mito.genes <- grep(pattern = "^mt-", x = rownames(x = clean_data@data), ignore.case=TRUE, value = TRUE)
}
percent.mito <- Matrix::colSums(clean_data@raw.data[mito.genes, ])/Matrix::colSums(clean_data@raw.data)
clean_data <- AddMetaData(object = clean_data, metadata = percent.mito, col.name = "percent.mito")
print(paste0("para_list is", my_para))
if(my_para == "None"){
	nGene_low <- unname(quantile(clean_data@meta.data$nGene,0.001))
	nGene_high <- unname(quantile(clean_data@meta.data$nGene,0.999))
	mt_high <- unname(quantile(clean_data@meta.data$percent.mito,0.98))
	mt_low <- -Inf
	my_reso <- 0.6
}else{
	my_para <- str_trim(unlist(str_split(my_para,",")))
	if(my_para[1] == "default" && my_para[2] == "default"){
		nGene_low <- unname(quantile(clean_data@meta.data$nGene,0.001))
		nGene_high <- unname(quantile(clean_data@meta.data$nGene,0.999))
	}else{
		nGene_low <- as.numeric(my_para[1])
		nGene_high <- as.numeric(my_para[2])
	}
	
	if(my_para[3] == "default" && my_para[4] == "default"){
		mt_low <- -Inf
		mt_high <- unname(quantile(clean_data@meta.data$percent.mito,0.98))
	}else{
		mt_low <- as.numeric(my_para[3])
		mt_high <- as.numeric(my_para[4])
	}
	
	if(my_para[5] == "default"){
		my_reso <- 0.6
	}else{
		my_reso <- as.numeric(my_para[5])
	}
}
clean_data <- FilterCells(object = clean_data, subset.names = c("nGene", "percent.mito"), low.thresholds = c(nGene_low, -Inf), high.thresholds = c(nGene_high, mt_high))
print("finish filter standard list get!")

### Normalizing the data
clean_data <- NormalizeData(object = clean_data, normalization.method = "LogNormalize")
print("finish normalize!")

### Detection of variable genes across the single cells
clean_data <- FindVariableGenes(object = clean_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = FALSE)
print("finish find variable genes!")

### Scaling the data and removing unwanted sources of variation
clean_data <- ScaleData(object = clean_data, vars.to.regress = c("nUMI", "percent.mito"))
print("finish scale data!")

### Perform linear dimensional reduction (PCA)
clean_data <- RunPCA(object = clean_data, pc.genes = clean_data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

print("try to plot elbow plot...")
pdf(file=paste0(outputpath,".standard_deviations_of_the_principle_components.pdf"))
### visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line)
tmp <- PCElbowPlot(object = clean_data)
tmp
fwrite(tmp$data, file=paste0(outputpath,".PC_selection.txt"),sep="\t")
dev.off()
print("finish plot!")

print("try to save dataset...")
save(clean_data, file=paste0(outputpath,".pca.clean_data.RData"))






