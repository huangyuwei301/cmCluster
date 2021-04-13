import sys,getopt
import os, re, sys

def run(runName,sourcePath,codePath,singularityPath,inputInfo,outputPath,windowSize,qsubmode):
    if not os.path.exists("{od}/{rn}".format(od=outputPath,rn=runName)):
        print("Dir {od}/{rn} does not exist. Trying to make it...".format(od=outputPath,rn=runName))
        if os.system("mkdir -p {od}/{rn}".format(od=outputPath,rn=runName)) != 0:
            sys.exit("Failed!")
    qsub_sh="{od}/{rn}/run_cmCluster.{rn}.sh".format(od=outputPath,rn=runName)
    qsubfile=open(qsub_sh,"w")
    qsubfile.write('''
#!/bin/sh
set -e
#run cmCluster for cell ranger output file

###copy inputdata and turn into simple file format for read10X in Seurat
echo -n `date` && echo "try to copy data..."
    mkdir -p {outputPath}/{rn}
    mkdir -p {outputPath}/{rn}/raw_data
    mkdir -p {outputPath}/{rn}/shell
    cp {inputInfo}/*features.tsv* {outputPath}/{rn}/raw_data/
    cp {inputInfo}/*barcodes.tsv* {outputPath}/{rn}/raw_data/
    cp {inputInfo}/*matrix.mtx* {outputPath}/{rn}/raw_data/
    cp {sourcePath}/user_marker.txt {outputPath}/{rn}
    cd {outputPath}/{rn}/raw_data
    gunzip -fq {outputPath}/{rn}/raw_data/*
    if [[ ! -s {outputPath}/{rn}/raw_data/genes.tsv ]];then
        mv {outputPath}/{rn}/raw_data/*features.tsv {outputPath}/{rn}/raw_data/genes.tsv
    fi
    if [[ ! -s {outputPath}/{rn}/raw_data/barcodes.tsv ]];then
        mv {outputPath}/{rn}/raw_data/*barcodes.tsv {outputPath}/{rn}/raw_data/barcodes.tsv
    fi
    if [[ ! -s {outputPath}/{rn}/raw_data/matrix.mtx ]];then
        mv {outputPath}/{rn}/raw_data/*matrix.mtx {outputPath}/{rn}/raw_data/matrix.mtx
    fi
    
    '''.format(inputInfo=inputInfo,outputPath=outputPath,sourcePath=sourcePath,rn=runName))
    
    ######1.data preprocess by seurat(filter&normalization&pca)
    qsubfile.write('''
echo -n `date` && echo "seurat started preprocess for cellranger result..."
    cd {sourcePath}
    python3 -c 'from function_work.sc_seurat import seurat; seurat(runName="{rn}",pipelinelog="{od}/{rn}/shell/seurat_preprocess.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.seurat.preprocess.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}/raw_data", parameter_list="None" ,qsubmode="{qsubmode}")'
    echo -n `date` && echo "seurat preprocess done!"
    
    '''.format(rn=runName,singularityPath=singularityPath,od=outputPath,qsubmode=qsubmode,sourcePath=sourcePath,codePath=codePath))

    ######2.cluster with different parameters
    qsubfile.write('''
echo -n `date` && echo "seurat started clustering for normalized data..."
    cd {sourcePath}
    mkdir -p {od}/{rn}/cluster
    python3 -c 'from function_work.sc_seurat import seurat_cluster; seurat_cluster(runName="{rn}",pipelinelog="{od}/{rn}/shell/seurat_cluster.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.seurat.cluster.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}", windowSize={windowSize} ,qsubmode="{qsubmode}")'
    echo -n `date` && echo "seurat clustering done!"
    
    '''.format(rn=runName,singularityPath=singularityPath,od=outputPath,windowSize=windowSize,qsubmode=qsubmode,sourcePath=sourcePath,codePath=codePath))

    ######3.get agree&noise cell and prepare data for annotation
    qsubfile.write('''
echo -n `date` && echo "cmCluster started for cell detection and annotation..."
    cd {sourcePath}
    echo `date` "cmCluster get agree&noise cell..."
    python3 -c 'from function_work.sc_splitcell import cellDetection; cellDetection(runName="{rn}",pipelinelog="{od}/{rn}/shell/cell_detection.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.cell.detection.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}/cluster", qsubmode="{qsubmode}")'
    echo `date` "cmCluster start to prapare input data for agree cell to annotate..."
    python3 -c 'from function_work.sc_annotation import agreePrepare; agreePrepare(runName="{rn}",pipelinelog="{od}/{rn}/shell/agree_cell.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.agree.cell.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}", qsubmode="{qsubmode}")'
    echo `date` "cmCluster start to prapare input data for noise cell to annotate..."
    python3 -c 'from function_work.sc_annotation import noisePrepare; noisePrepare(runName="{rn}",pipelinelog="{od}/{rn}/shell/noise_cell.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.noise.cell.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}", qsubmode="{qsubmode}")'
    echo `date` "cmCluster start to annotate cluster ..."
    python3 -c 'from function_work.sc_annotation import cluAnnotation; cluAnnotation(runName="{rn}",pipelinelog="{od}/{rn}/shell/cluster_annotation.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.cluster.annotation.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}", qsubmode="{qsubmode}")'
    echo `date` "cmCluster finish cell detection and annotation!"
    
    '''.format(rn=runName,singularityPath=singularityPath,od=outputPath,windowSize=windowSize,qsubmode=qsubmode,sourcePath=sourcePath,codePath=codePath))

    ######4.evaluate cluster result
    qsubfile.write('''
echo `date` "cmCluster starts evaluation..."
    cd {sourcePath}
    echo `date` "cmCluster starts calculate accuracy and sort parameter ..."
    mkdir -p {od}/{rn}/evaluate
    python3 -c 'from function_work.sc_cmCluster import cmCluster; cmCluster(runName="{rn}",pipelinelog="{od}/{rn}/shell/cmCluster.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.cmCluster.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}", qsubmode="{qsubmode}")'
    echo `date` "heatmap starts ..."
    mkdir -p {od}/{rn}/statistic
    mv {od}/{rn}/filelist.txt {od}/{rn}/statistic
    python3 -c 'from function_work.sc_cmCluster import cmHeatmap; cmHeatmap(runName="{rn}",pipelinelog="{od}/{rn}/shell/heatmap.{rn}.log", tasklogprefix="{od}/{rn}/shell/qsub.heatmap.{rn}", singularityPath="{singularityPath}", codePath="{codePath}", outputPath="{od}/{rn}",inputPath="{od}/{rn}", qsubmode="{qsubmode}")'
    echo `date` "cmCluster finishes evaluation!"
    '''.format(rn=runName,singularityPath=singularityPath,od=outputPath,windowSize=windowSize,qsubmode=qsubmode,sourcePath=sourcePath,codePath=codePath))

    print("shell file to run cmCluster will be see as {od}/{rn}/run_cmCluster.{rn}.sh".format(od=outputPath,rn=runName))
    qsubfile.close()

def run_wrap(runName,inputInfo,outputPath,windowSize):
    sourcePath=sys.path[0]+"/script"            ###"/cmCluster/script"
    singularityPath=sys.path[0]+"/singularity"  ###"/cmCluster/singularity"
    codePath=sys.path[0]+"/script/function_code"
    thread=4
    qsubmode="local"
    print("input data will be collect form: "+inputInfo)
    print("and output result will be put in: "+outputPath)
    run(runName,sourcePath,codePath,singularityPath,inputInfo,outputPath,windowSize,qsubmode)

def main(argv):
    try:
        opts,args = getopt.getopt(argv,"hi:o:d:w:")
    except getopt.GetoptError:
        print('usage of scrnaseq_sh_create.py \n -i <absolute directory of input file> \n -o <absolute output directory> \n -d <output runname&runName> \n -w <slide window size of cmCluster>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('usage of scrnaseq_sh_create.py \n -i <absolute directory of input file> \n -o <absolute output directory> \n -d <output runname&runName> \n -w <slide window size of cmCluster>')
            sys.exit()
        elif opt in("-i"):
            inputInfo=arg
        elif opt in("-o"):
            outputPath=arg
        elif opt in("-d"):
            runName=arg
        elif opt in("-w"):
            windowSize=arg
    run_wrap(runName,inputInfo,outputPath,windowSize)

if __name__=='__main__':
    main(sys.argv[1:])
    sourcePath=sys.path[0]+"/script" ###"/cmCluster/script"
    thread=4
    qsubmode="local"
