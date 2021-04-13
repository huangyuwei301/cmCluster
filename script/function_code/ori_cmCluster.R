#####test the influence of blood samples
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(ggplot2)))
##### library
suppressWarnings(suppressMessages(library(optparse)))
suppressWarnings(suppressMessages(library(SingleCellExperiment)))
suppressWarnings(suppressMessages(library(cellassign)))
suppressWarnings(suppressMessages(library(pheatmap)))

#### parameters of Monocle
option_list = list(
make_option(c("-i", "--input_dir"), type="character",
              help="Input directory of expression matrix by 10xGenomics", metavar="character"),
make_option(c("-r", "--runname"), type="character", default = "None",
              help="runname of cmcluster", metavar="character"),
make_option(c("-o", "--output_file"), type="character", 
              help="Absolute path of output directory and prefix of output, exp: {output_path}/{output_prefix}", metavar="character")
); 

opt <- parse_args(OptionParser(option_list=option_list))

Args <- commandArgs()
## parameters
da_path1 <- opt$input_dir       ## input information
sample_name <- opt$runname
re_path <- opt$output_file     ## path and prefix of output files

cor_filename <- paste0(sample_name,"-tsne-correspond.tab")
cluann_filename <- paste0(sample_name,".annotation.tab")

cell_marker_filename <- "user_marker.txt" 
exp_name <- paste0(sample_name,".disagree.RData") 
si_filename <- paste0(sample_name,".sizefactor.RData") 
tsne_filename <- paste0(sample_name,".tsne.tab")
#output data
ann_filename <- paste0(sample_name,".cellassign.annotation.tab")
ann_tsne_filename <- paste0(sample_name,".cellassign.ann.tsne.tab")
ce_po_filename <- paste0(sample_name,".cellassign.cell_cell_score.csv")
cl_po_filename <- paste0(sample_name,".cellassign.cluster_cell_score.csv")


#calculate statistic parameter
getParameter <- function(true_label,pred_label){
    #calculate the average accuracy,precision and recall
    my_clu <- names(table(true_label))
    para <- vector()
    for(k in 1:length(my_clu)){
    #true_DEG,true_non_DEG,pre_DEG,pre_non_DEG){
        temp_true <- true_label
        temp_true[which(true_label == my_clu[k])] <- 1
        temp_true[which(true_label != my_clu[k])] <- 0
        temp_pred <- pred_label
        temp_pred[which(pred_label == my_clu[k])] <- 1
        temp_pred[which(pred_label != my_clu[k])] <- 0
        TP = length(intersect(which(temp_true==1),which(temp_pred==1))) #length(intersect(true_DEG,pre_DEG))        #True Positive
        FP = length(intersect(which(temp_true==0),which(temp_pred==1)))#length(intersect(true_non_DEG,pre_DEG))    #False Positive
        TN = length(intersect(which(temp_true==0),which(temp_pred==0)))#length(intersect(true_non_DEG,pre_non_DEG))#True Negative
        FN = length(intersect(which(temp_true==1),which(temp_pred==0)))#length(intersect(true_DEG,pre_non_DEG))    #False Negative

        #print(TP)
        #print(FP)
        #print(TN)
        #print(FN)

        if(TP == 0 && FP == 0){
            Precision = 0                      
        }else{
            Precision = TP/(TP+FP)
        }
        if(TP == 0 && FN == 0){
            Recall = 0                        
        }else{
            Recall = TP/(TP+FN)
        }
        #Accuracy = (TP+TN)/(TP+FP+TN+FN)          
        if((Precision+Recall) != 0){
            F1 = 2*Precision*Recall/(Precision+Recall) 
        }else{
            F1 = 0
        }
        Support = length(which(temp_true == 1))    

        para = rbind(para,c(TP,FP,TN,FN,Precision,Recall,F1,Support))
    }
    colnames(para) <- c("TP","FP","TN","FN","Precision","Recall","F1-score","Support")
    rownames(para) <- my_clu
    
    if(length(my_clu) > 1){
    total_value <- apply(para[,c(1:4,8)],2,sum)
    #print(total_value)
    
    macro_avg <- c(total_value[1:4],apply(para[,5:7],2,mean),total_value[5])
    names(macro_avg) <- c("TP","FP","TN","FN","Precision","Recall","F1-score","Support")
    #print(macro_avg)
    
    if(total_value[1]==0 && total_value[2]==0){
        mPrecision=0
    }else{
        mPrecision=total_value[1]/(total_value[1]+total_value[2])
    }
    if(total_value[1]==0 && total_value[4]){
        mRecall=0
    }else{
        mRecall=total_value[1]/(total_value[1]+total_value[4])
    }
    if((mPrecision+mRecall) != 0){
        mF1 = 2*mPrecision*mRecall/(mPrecision+mRecall)
    }else{
        mF1 = 0
    }
    
    micro_avg <- c(total_value[1:4],mPrecision,mRecall,mF1,total_value[5])
    names(micro_avg) <- c("TP","FP","TN","FN","Precision","Recall","F1-score","Support")
    #print(micro_avg)
    }else{
    macro_avg <- para
    names(macro_avg) <- c("TP","FP","TN","FN","Precision","Recall","F1-score","Support")
    micro_avg <- para
    names(micro_avg) <- c("TP","FP","TN","FN","Precision","Recall","F1-score","Support")
    }
    
    para_table <- rbind(rbind(para,macro_avg),micro_avg)
    rownames(para_table) <- c(rownames(para),"macro_avg","micro_avg")
    
    print("finish statistic evaluate!")
    return(para_table)
}

print("read in corresponding matrix...")
cor_matrix_ori <- read.table(paste0(da_path1,"/",cor_filename),header=T,row.names=1,sep="\t")
cat("finish read in correspond matrix!\n")
cor_matrix <- cor_matrix_ori

print("read in cluster annotation result for agree cells...")
clu_ann <- read.table(paste0(da_path1,"/",cluann_filename),header=T,row.names=1,sep="\t")

print("read in ori data for cellassign...")
cat("try to load Seurat data...\n")
load(paste0(da_path1,"/",exp_name))
cat("try to read in gene marker file and get cellassign input gene list...\n")
cell_markers_file <- read.table(paste0(da_path1,"/",cell_marker_filename),sep="\t",header=T)
cell_marker <- list()
for(i in 1:dim(cell_markers_file)[1]){
    cell_marker[[as.character(cell_markers_file[i,1])]] <- strsplit(as.character(cell_markers_file[i,2]),split=",")[[1]]
}
cell_marker_ori <- marker_list_to_mat(cell_marker,include_other=FALSE)
cat("try to load sizefactor file...\n")
load(paste0(da_path1,'/',si_filename))


##### annotate all of the disagree cell may be use in the parameter
print("annotate the disagree cell maybe use in the parameter calculation!") 
print("prepare true annotation vector...")
###read in expression matrix
cluster_times <- table(cor_matrix$cluster_times)
dis_loc <- vector()
i <- length(cluster_times) 
while(length(dis_loc)<1000 && i!=1){
    cat("get disagree cell from ",cluster_times[i],"\n")
    dis_loc <- c(dis_loc,which(cor_matrix$cluster_times==names(cluster_times)[i]) )
    i <- i-1

}
dis_cell <- rownames(cor_matrix)[dis_loc]
da_exp_ori <- as.matrix(clean_data1@data[,dis_cell])

###read in gene markers for each cell type and convert it to cellassign marker matrix
#the default marker list is from cellassign article supplement table2:HGSC and FL combined
cat("try to read gene markers and convert to cellassign format...\n")
cell_marker_cellassign <- cell_marker_ori[intersect(rownames(da_exp_ori),rownames(cell_marker_ori)),]
cat("there maybe",length(rownames(cell_marker_cellassign)),"marker genes for",length(colnames(da_exp_ori)),"cells!")

###covert count matrix to object SingleCellExperiment
cat("get count data and create cellassign required count data format...\n")
da_exp_count <- as.matrix(clean_data1@raw.data)[rownames(da_exp_ori),colnames(da_exp_ori)]
da_exp_sce <- SingleCellExperiment(assays = list(counts = da_exp_count,logcounts = log2(da_exp_count+1)),
        colData=data.frame(cell_names = colnames(da_exp_ori)),
        rowData=data.frame(gene_names = rownames(da_exp_ori)))

###read in sizefactor result from monocle
cat("read in size factors...\n")
s <- raw_sizefactor[colnames(da_exp_ori)]
rm(da_exp_count)

###do cellassign
cat("start cellassign...\n")
fit <- cellassign(exprs_obj = da_exp_sce[rownames(cell_marker_cellassign),],
        marker_gene_info = cell_marker_cellassign,
        s = s,
        learning_rate = 1e-2,
        shrinkage = TRUE,
        verbose = FALSE)
res <- cellprobs(fit)
rownames(res) <- colnames(da_exp_ori)
cat("cellassign seccess!\n")

###label cell by cellassign possibility matrix and record the possibility
cat("deal with the result and write into files...\n")
ce_ce_label <- c()
for(i in 1:dim(res)[1]){
    temp <- colnames(res)[which(res[i,]==max(res[i,]))]
    ce_ce_label <- rbind(ce_ce_label,c(temp,max(res[i,])))
}
rownames(ce_ce_label) <- colnames(da_exp_ori)
colnames(ce_ce_label) <- c("Cell_type","probability")
rm(da_exp_ori)
cat('finish calculate cell type probability!\n')
print("finish true annotaion!")


run_times <- 1
myround <- dim(cor_matrix)[2] - 3

while(run_times < myround){
    print("calculate the cluster_times for corresponding matrix...")
    if(dim(cor_matrix)[2] == dim(cor_matrix_ori)[2]){
        cat("load the cluster times for all parameter...\n")
        cluster_times <- table(cor_matrix$cluster_times)
        print(cluster_times)
    }else{
        cat("recalculate cluster times after delete parameter...\n")
        for(i in 1:dim(cor_matrix)[1]){
            temp <- table(as.vector(as.matrix(cor_matrix[i,1:(dim(cor_matrix)[2]-2)])))
            cor_matrix[i,(dim(cor_matrix)[2]-1)] <- length(temp) #clu_times
            cor_matrix[i,(dim(cor_matrix)[2])] <- max(temp) #clu_max
        }
        cluster_times <- table(cor_matrix$cluster_times)
        print(cluster_times)
        cat("finish recalculate cluster times for new parameter!\n")
    }

    max_times_cell <- rownames(cor_matrix)[which(cor_matrix$cluster_times==max(names(cluster_times)))]
    print(paste0("the max cluster times is:",max(names(cluster_times)),", and the number of cells is:",length(max_times_cell)))


    mytimes <- 2
    times_sort <- sort(names(cluster_times),decreasing=T)
    if(length(times_sort)==1){
        print("finish calculate for all the rest parameter has agree clusters!")
        break;
    }
    while(length(max_times_cell) < 11 && mytimes < length(times_sort) ){
        more_cell <- rownames(cor_matrix)[which(cor_matrix$cluster_times==times_sort[mytimes])]
        max_times_cell <- c(max_times_cell,more_cell)
        print(paste0("add cluster times is:",times_sort[mytimes],", and the number of cells is:",length(max_times_cell)))
        mytimes <- mytimes+1
    }

    print("prepare predict annotation matrix...")
    cat("read in cluster annotation result!\n")
    clu_ann_matrix <- as.matrix(cor_matrix[max_times_cell,1:(dim(cor_matrix)[2]-2)])
    for(i in 1:dim(clu_ann_matrix)[2]){
        for(j in 1:dim(clu_ann)[1]){
            clu_ann_matrix[which(clu_ann_matrix[,i] == as.character(clu_ann[j,1])),i] <- as.character(clu_ann[j,2])
        }
    }
    cat("finish transform cluster to cell type!\n")

    print("start to calculate the accuracy, precision and recall...")
    eva_value <- vector()
    for(i in 1:dim(clu_ann_matrix)[2]){
        pre_type <- as.vector(data.frame(strsplit(clu_ann_matrix[,i],"_"))[1,])
        tru_type <- as.vector(ce_ce_label[rownames(clu_ann_matrix),"Cell_type"])
        print(paste0("parameter ",i," has the precision matrix:"))
        temp <- getParameter(tru_type,pre_type)
        print(temp)
        eva_value <- rbind(eva_value,c(as.vector(temp["macro_avg",]),as.vector(temp["micro_avg",])))
    }
    rownames(eva_value) <- colnames(clu_ann_matrix)
    colnames(eva_value) <- c("macro_TP","macro_FP","macro_TN","macro_FN","macro_Precision","macro_Recall","macro_F1-score","macro_Support","micro_TP","micro_FP","micro_TN","micro_FN","micro_Precision","micro_Recall","micro_F1-score","micro_Support")
    write.table(eva_value,paste0(re_path,"/",run_times,sample_name,".evaluation.tab"),sep="\t")

    #the result of most cluster is more important, so we take the macro-f-score as the standard to select worst parameter
    #
    #here the cluster annotation only has part of all the cluster, so the predict cell type label may not be complete
    #two way to solve it: 1.delete the NA label first, 2.take the NA label as the wrong predict
    #
    #and the cellassign may annotate some cell as other cell type, so the true label is much more strict
    #two way to solve it: 1.give the most samilar cell type as true label, 2.take any cell type as the wrong predict
    #
    #all of these problem will effect the selection of worst parameter

    #which.max(eva_value$macro_F1-score) #only return first min value location, when there are more one mins this will not be properly
    #so check F1 score first, then precision, then TP and finally TN
    worst_para_list <- which(eva_value[,7]==min(eva_value[,7])) #f1 score
    if(length(worst_para_list)>1){
        sub_eva <- eva_value[worst_para_list,]
        worst_para_list <- intersect(worst_para_list,which(eva_value[,5]==min(sub_eva[,5]))) #precision
        if(length(worst_para_list)>1){
            sub_eva <- eva_value[worst_para_list,]
            worst_para_list <- intersect(worst_para_list,which(eva_value[,1]==min(sub_eva[,1]))) #TP
            if(length(worst_para_list)>1){
                sub_eva <- eva_value[worst_para_list,]
                worst_para_list <- intersect(worst_para_list,which(eva_value[,3]==min(sub_eva[,3]))) #TN
            }
        }
    }
    worst_para <- rownames(eva_value)[worst_para_list]
    print("the result of worst parameter is:")
    print(worst_para)
    cor_matrix <- cor_matrix[,setdiff(colnames(cor_matrix),worst_para[1])]
    print("current parameter list is:")
    print(colnames(cor_matrix))
    run_times <- run_times + 1
}





