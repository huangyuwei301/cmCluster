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
  make_option(c("-w", "--windowSize"), type="double", default = 7,
              help="gene list for tSNE plot", metavar="double"),
  make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
  
); 

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
input_info <- opt$input_dir       ## input information
outputpath <- opt$output_file     ## path and prefix of output files

sam_species <- opt$species
my_win <- opt$windowSize
tsne_genes <- opt$genes

get_default_pc <- function(e_result){
        ###input the result matrix for elbow plot
        print("calculate properly pc...")
        the_pc <- 1
        the_deta <- max(e_result[,2])-min(e_result[,2])
        the_cutoff <- the_deta*0.05
        print(paste("the cutoff in this dataset will be",the_cutoff))
        for(i in 2:dim(e_result)[1]){
                print(paste0("for pc ",the_pc," with values ",e_result[i,2],", and the deta will be",the_deta))
                if(the_deta >= the_cutoff){
                        the_pc <- i
                        the_deta <- e_result[(i-1),2] - e_result[i,2]
                }else{
                        break
                }
        }
        return(the_pc)
}

#####################################################################################
###### Load the dataset
#####################################################################################
print("try to load merge clean data...")
load(paste0(input_info,".pca.clean_data.RData"))

## get default pc from elbow plot result
tmp <- read.table(paste0(input_info,".PC_selection.txt"),header=T,sep="\t")
pc_center <- get_default_pc(tmp)

## Cluster the cells
reso_list <- seq(from=0.8-floor(my_win/2)*0.1,to=0.8+floor(my_win/2)*0.1,by=0.1)
knei_list <- seq(from=30-floor(my_win/2),to=30+floor(my_win/2),by=1)
pc_list <- seq(from=pc_center - floor(my_win/2),to=pc_center + floor(my_win/2),by=1)
for(i in pc_list){
    for(j in knei_list){
        for(k in reso_list){
            print(paste("for ",i," PCs,", j," neighbors,",k," resolutions, try to cluster by tsne..."))
            clean_data <- FindClusters(object = clean_data, reduction.type = "pca", dims.use = 1:i, resolution = k, k.param = j, print.output = 0, save.SNN = TRUE)
            PCA_output <- as.data.table(clean_data@ident,keep.rownames=TRUE)
            colnames(PCA_output) <- c("Barcode","Group")
            #fwrite(PCA_output, file=paste0(outputpath,".PCA_group.tab"),sep="\t")

            #####################################################################################
            ###### Run Non-linear dimensional reduction (tSNE)
            #####################################################################################
            if(length(clean_data@cell.names)>80){
                clean_data <- RunTSNE(object = clean_data, dims.use = 1:i, do.fast = TRUE)
            }else{
                clean_data <- RunTSNE(object = clean_data, dims.use = 1:i,  perplexity=10)
            }
            if(length(levels(clean_data@ident))!=1){
                print("try to save tsne result...")
                current.cluster.ids <- levels(clean_data@ident)
                new.cluster.ids <- paste0("cluster_",as.numeric(current.cluster.ids)+1)
                clean_data@ident <- plyr::mapvalues(x = clean_data@ident, from = current.cluster.ids, to = new.cluster.ids)
                pdf(file=paste0(outputpath,".tsne_pc",i,"k",j,"r",k,".pdf"))
                tsne_plot <- TSNEPlot(object = clean_data, do.label = TRUE, pt.size = 0.5)
                dev.off()
                fwrite(as.data.table(tsne_plot$data,keep.rownames=TRUE), file=paste0(str_trim(outputpath),".tsne_pc",i,"k",j,"r",k,".tab"),sep="\t")
            }else{
                stop("No clustering succeeded")
            }
        }
    }
}
