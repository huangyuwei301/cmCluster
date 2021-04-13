import luigi
import os, re, sys
from utils import *

def seurat(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            parameter_list="None",
            qsubmode="local"):

    print("[Task  Added] seurat on {}...".format(runName))
    rn = runName
    file_path = codePath

    if not os.path.exists(inputPath):
        sys.exit("{} does not exist!".format(inputPath))
    inputPath = os.popen("readlink -f {}".format(inputPath)).read().strip()
    input_prefix = "/" + inputPath.split("/")[1]

    if not os.path.exists(outputPath):
        print("Dir {} does not exist. Trying to make it...".format(outputPath))
        if os.system("mkdir -p {}".format(outputPath)) != 0:
            sys.exit("Failed!")

    od = os.popen("readlink -f {}".format(outputPath)).read().strip()
    output_prefix = "/" + outputPath.split("/")[1]

    outputPathPrefix = od + "/" + rn
    mountstring = mountgen([input_prefix, file_path, output_prefix])

    commandline = "Rscript {file_path}/ori_seurat.R -i {inputPath} -o {outputPathPrefix} -s human -p {para}".format(
        file_path=file_path, inputPath=inputPath, outputPathPrefix=outputPathPrefix, para=parameter_list)
    
    qsub_sh = "{od}/shell/qsub.seurat.{taskid}.sh".format(od=od, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts seurat ..."
    singularity exec -B {mountstring} {singularityPath}/advance.sif {commandline}
    echo `date` "{rn} finishes seurat"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] seurat on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} seurat failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True

def seurat_cluster(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            windowSize=7,
            qsubmode="local"):

    print("[Task  Added] seurat cluster on {}...".format(runName))
    rn = runName
    file_path = codePath
    outputOld = outputPath
    outputPath = "{}/cluster".format(outputPath)

    if not os.path.exists(inputPath):
        sys.exit("{} does not exist!".format(inputPath))
    inputPath = os.popen("readlink -f {}".format(inputPath)).read().strip()
    input_prefix = "/" + inputPath.split("/")[1]
    inputPathPrefix = inputPath + "/" + rn

    if not os.path.exists(outputPath):
        print("Dir {} does not exist. Trying to make it...".format(outputPath))
        if os.system("mkdir -p {}".format(outputPath)) != 0:
            sys.exit("Failed!")

    od = os.popen("readlink -f {}".format(outputPath)).read().strip()
    output_prefix = "/" + outputPath.split("/")[1]

    outputPathPrefix = od + "/" + rn
    mountstring = mountgen([input_prefix, file_path, output_prefix])


    commandline = "Rscript {file_path}/ori_seurat_cluster.R -i {inputPathPrefix} -o {outputPathPrefix} -s human -w {windowSize}".format(
        file_path=file_path, inputPathPrefix=inputPathPrefix, outputPathPrefix=outputPathPrefix, windowSize=windowSize)
    
    qsub_sh = "{outputOld}/shell/qsub.seurat.{taskid}.sh".format(outputOld=outputOld, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts seurat cluster ..."
    singularity exec -B {mountstring} {singularityPath}/advance.sif {commandline}
    echo `date` "{rn} finishes seurat cluster"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] seurat cluster on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} seurat_cluster failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True

