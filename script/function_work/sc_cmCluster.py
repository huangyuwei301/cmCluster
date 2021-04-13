import luigi
import os, re, sys
from utils import *

def cmCluster(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            qsubmode="local"):

    print("[Task  Added] cmCluster on {}...".format(runName))
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
    outputPathPrefix = od + "/evaluate"

    mountstring = mountgen([input_prefix, file_path, output_prefix])

    commandline = "Rscript {file_path}/ori_cmCluster.R -i {inputPath} -o {outputPathPrefix} -r {rn}".format(
        file_path=file_path, inputPath=inputPath, outputPathPrefix=outputPathPrefix, rn=rn)
    
    qsub_sh = "{od}/shell/qsub.cmcluster.{taskid}.sh".format(od=od, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts cmcluster ..."
    singularity exec -B {mountstring} {singularityPath}/cellassign.sif {commandline}
    echo `date` "{rn} finishes cmcluster"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] cmCluster on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} cmCluster failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True



def cmHeatmap(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            qsubmode="local"):

    print("[Task  Added] heatmap on {}...".format(runName))
    rn = runName
    file_path = codePath

    if not os.path.exists(inputPath):
        sys.exit("{} does not exist!".format(inputPath))
    inputPath = os.popen("readlink -f {}".format(inputPath)).read().strip()
    input_prefix = "/" + inputPath.split("/")[1]
    inputPathPrefix = inputPath + "/evaluate"

    if not os.path.exists(outputPath):
        print("Dir {} does not exist. Trying to make it...".format(outputPath))
        if os.system("mkdir -p {}".format(outputPath)) != 0:
            sys.exit("Failed!")

    od = os.popen("readlink -f {}".format(outputPath)).read().strip()
    output_prefix = "/" + outputPath.split("/")[1]
    outputPathPrefix =  "{}/statistic".format(od)

    mountstring = mountgen([input_prefix, file_path, output_prefix])

    commandline = "Rscript {file_path}/get_triangle_matrix.R -i {inputPathPrefix} -o {outputPathPrefix} -r {rn}".format(
        file_path=file_path, inputPathPrefix=inputPathPrefix, outputPathPrefix=outputPathPrefix, rn=rn)
    
    qsub_sh = "{od}/shell/qsub.heatmap.{taskid}.sh".format(od=od, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts heatmap for cmcluster ..."
    singularity exec -B {mountstring} {singularityPath}/advance.sif {commandline}
    echo `date` "{rn} finishes heatmap for cmcluster!"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] heatmap on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} heatmap failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True


