# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:54:20 2019

@author: user
"""

#import numpy as np
import pandas as pd
import tables
from scipy import sparse
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
import os
import argparse

'''
    use gene marker to detect cell type of each cluster
    score gene marker for the expression count of cells, and get cluster score by average cell score
    plot heatmap to visualize the score and decide cell type for cluster
'''

def _read_other_information(filename, info_loc, datatype):
    """
    Read hdf5 file from Cell Ranger v2 or earlier versions.
    info_loc is a string for corresponding structure path
    """
    with tables.open_file(str(filename), 'r') as f:
        try:
            #f.list_nodes(f.root)
            children = [x._v_name for x in f.list_nodes(info_loc)]
            print('the variable names are:',children)
            dsets = {}
            for node in f.walk_nodes(info_loc, datatype):
                print('read in data from:',node)
                dsets[node.name] = node.read()
            
            return dsets
        except KeyError:
            raise Exception('File is missing one or more required datasets.')

def get_data(file_pathway,filename):
    f = open(file_pathway+'/'+filename)
    filecontain = f.readlines()
    print('finish read in data:',file_pathway+'/'+filename)
    f.close()
    return filecontain

def split_file_get_infor(string_list,the_spliter):
    #split the file with a spliter and get all information from the original files
    all_information_list = [string_list[i].split(the_spliter) for i in range(len(string_list))]
    print('finish split string list',len(string_list),'by spliter',the_spliter,'and get all information you need')
    return all_information_list

###define parameters for cluster annotation
#necessary: input_dir,name
#optional: data_expr,input_marker,output_dir
parser = argparse.ArgumentParser(description='parameters of GSVA')
parser.add_argument('-i','--input_dir',required=True,help='Absolute path of input file')
parser.add_argument('-n','--name',required=True,help='project name of cluster analysis for cells')
parser.add_argument('-im','--input_marker',default='user_markers.txt',help='file name of user input gene markers for different cells')
parser.add_argument('-o','--output_dir',help='Absolute path of output file')

args = parser.parse_args()

da_path = args.input_dir 
print('read data from ',da_path)
pro_name = args.name
tsne_file = pro_name+'.tsne.tab'
print('and the sample ID is',pro_name,', the tsne cluster file name is',tsne_file)
exp_file = pro_name+'.anndata_advance.h5'
print('the input expression matrix is',exp_file)
marker_file = args.input_marker
print('the user input marker is',marker_file)
if args.output_dir == 'None':
    re_path = da_path
    print('and all of result file will be put in',re_path)
else:
    re_path = args.output_dir
    print('and all of result file will be put in',re_path)


###start analysis...
###read in data
os.chdir(da_path)
print('change work path to',os.getcwd())

###read in corresponding expression matrix information of seurat result files
print('try to read expression matrix...')
X_ori = _read_other_information(da_path+'/'+exp_file,'/X','Array')
#trans the sparse matrix to DataFrame
X = pd.DataFrame(sparse.csc_matrix((X_ori['data'], X_ori['indices'], X_ori['indptr']), shape=X_ori['shape']).todense())
features = [str(X_ori['features'][i],'utf-8') for i in range(0,len(X_ori['features']))]
barcodes = [str(X_ori['barcodes'][i],'utf-8') for i in range(0,len(X_ori['barcodes']))]
#give the rownames by features(genes) and colnames by barcodes(cells)
X.index = features
X.columns = barcodes

#read in tsne cluster result
print('try to read cluster information...')
tsne_cluster = split_file_get_infor(get_data(da_path,tsne_file),'\t')
temp_clu = [tsne_cluster[i][3] for i in range(1,len(tsne_cluster))]
temp_bar = [tsne_cluster[i][0] for i in range(1,len(tsne_cluster))]
tsne = pd.DataFrame(temp_clu,index=temp_bar,columns=['tsne'])

###read in user marker
print('try to read user marker gene list...')
user_marker = split_file_get_infor(get_data(da_path,marker_file),'\t')

###get all marker genes
print('get marker gene list of sequencing matrix...')
all_marker_genes = []
for i in range(1,len(user_marker)):
    #temp = user_marker[i][1].split('\n')[0].split(',')
    temp = user_marker[i][1].strip().split(',')
    for j in range(0,len(temp)):
        if temp[j] not in all_marker_genes and temp[j] in X.index:
            all_marker_genes.append(temp[j])

print('extract the sub matrix of marker genes...')
#get expression matrix of marker genes        
X_sub = X.loc[all_marker_genes]#[list(set(all_marker_genes).intersection(set(X.index)))]

###score marker genes for each cell and correct by detection rate
print('score marker genes and correct...')
cell_type = X_sub
for i in cell_type.columns:
    temp = cell_type[i]
    marker_index = temp[temp > 0].index
    if len(marker_index) > 0:
        marker_weight = 1/len(marker_index)
        #cell_type[i][marker_index] = marker_weight
        cell_type[i][marker_index] = cell_type[i][marker_index]*marker_weight

print('calculate average cluster score for each marker...')
#define the average score matrix 
cluster_marker = pd.DataFrame(index=cell_type.index)
#get cluster type for each cell
cluster_type = Counter(tsne['tsne'])
#calculate average cell score for each cluster
for e in cluster_type.keys():
    temp = cell_type[tsne[tsne.tsne==e].index].mean(axis=1)
    cluster_marker[e] = temp

print('calculate average cluster score for each cell type...')
#define the average score cluster matrix
cluster_celltype = pd.DataFrame(columns=['cell_ann']) #pd.DataFrame(columns=cluster_marker.columns)
#get marker genes for each cell types
#marker_genes = []
for i in range(1,len(user_marker)):
    temp = user_marker[i][1].strip().split(',')
    for j in range(0,len(temp)):
        #cluster_celltype.append([user_marker[i][0],temp[j]])
        cluster_celltype.loc[temp[j]] = {'cell_ann':user_marker[i][0]}

#calculate average gene score for each cluster to difine their cell type
cell_ann = Counter(cluster_celltype['cell_ann'])
cluster_cell = pd.DataFrame(columns=cluster_marker.columns)
for e in cell_ann.keys():
    index_cell_marker = cluster_celltype[cluster_celltype.cell_ann==e].index
    index_cell_marker_in = [i for i in index_cell_marker if i in cluster_marker.index]
    if len(index_cell_marker_in) > 0:
        temp = cluster_marker.loc[index_cell_marker_in].mean(axis=0)
        #temp = cluster_marker.loc[cluster_celltype[cluster_celltype.cell_ann==e].index].mean(axis=0)
        cluster_cell.loc[e] = temp

#annotate cluster by average cluster score matrix: cluster_cell
print('annotate cluster by cluster score...')
cluster_cell_ann = pd.DataFrame(columns=["id","Cluster","Cell_type"],index=[list(range(0,len(cluster_cell.columns)))])
i=0
for e in cluster_cell.columns:
    cluster_cell_ann.loc[i]['id'] = str(int(e.split('_')[1]) - 1)
    cluster_cell_ann.loc[i]['Cluster'] = e
    cluster_cell_ann.loc[i]['Cell_type'] = cluster_cell[e].idxmax() #define cluster max score as the cell type
    i=i+1

#sort cluster annotation dataframe by cell type
for i in cluster_cell_ann.index:
    cluster_cell_ann.loc[i]['id'] = int(cluster_cell_ann.loc[i]['id'])
#cluster_cell_ann = cluster_cell_ann.sort_values(by=['id'],axis=0,ascending=True)
#marker the duplicated cell type, add -n behind the cell type
Cell_type = Counter(cluster_cell_ann['Cell_type'])
for e in Cell_type.keys():
    one_cell_type_index = cluster_cell_ann[cluster_cell_ann.Cell_type==e].index
    if len(one_cell_type_index)>1:
        print('mark',e,'duplicate!')
        one_cell_type = cluster_cell_ann.loc[one_cell_type_index]
        one_cell_type = one_cell_type.sort_values(by=['id'],axis=0,ascending=True)
        i=1
        for j in one_cell_type.index:
            cluster_cell_ann.loc[j,'Cell_type'] = e+'_'+str(i)
            i=i+1

###
print('redefine tsne cluster file by replacing cluster name of cell type...')
for i in range(1,len(tsne_cluster)):
    new_clu = list(cluster_cell_ann[cluster_cell_ann.Cluster==tsne_cluster[i][3]]['Cell_type'])[0]
    tsne_cluster[i][3] = new_clu

print('plot and save result file...')
### plot and write result files
#sort cluster by their cell type
cluster_cell_ann = cluster_cell_ann.sort_values(by = ['Cell_type'],axis = 0,ascending = True)
#marker average score
fig = sns.heatmap(cluster_marker.loc[all_marker_genes][cluster_cell_ann['Cluster']],cmap='RdBu_r')#,center=0,vmax= 1, vmin= -1)
plt.savefig(re_path+'/'+pro_name+'.cluster_marker_score.pdf',bbox_inches='tight',pad_inches=0.0)
plt.close()
cluster_marker.loc[all_marker_genes][cluster_cell_ann['Cluster']].to_csv(path_or_buf=re_path+'/'+pro_name+'.cluster_marker_score.csv')
#marker average score
fig = sns.heatmap(cluster_cell[cluster_cell_ann['Cluster']],cmap='RdBu_r')#,center=0,vmax= 1, vmin= -1)
plt.savefig(re_path+'/'+pro_name+'.cluster_cell_score.pdf',bbox_inches='tight',pad_inches=0.0)
plt.close()
cluster_cell[cluster_cell_ann['Cluster']].to_csv(path_or_buf=re_path+'/'+pro_name+'.cluster_cell_score.csv')
#annotation corresponding file
cluster_cell_ann.to_csv(path_or_buf=re_path+'/'+pro_name+'.annotation.tab', sep='\t', index=False)
#annotation tsne file
f = open(re_path+'/'+pro_name+'.ann.tsne.tab','w')
for line in tsne_cluster:
    f.write('\t'.join(line))
print('finish read in data:',re_path+'/'+pro_name+'.ann.tsne.tab')
f.close()

print('finish cluster annotation!')
