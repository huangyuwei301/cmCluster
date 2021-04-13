import luigi
import os, re, sys
from utils import *


def agreePrepare(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            qsubmode="local"):

    print("[Task  Added] prepare agree cells on {}...".format(runName))
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

    mountstring = mountgen([input_prefix, file_path, output_prefix])

    commandline = "Rscript {file_path}/sc-clusterann.R -i {inputPath} -o {od} -r {rn} -n {rn}.pca.clean_data.RData".format(
        file_path=file_path, inputPath=inputPath, od=od, rn=rn)
    
    qsub_sh = "{od}/shell/qsub.agreePrepare.{taskid}.sh".format(od=od, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts agree preparation ..."
    singularity exec -B {mountstring} {singularityPath}/advance.sif {commandline}
    echo `date` "{rn} finishes agree preparation"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] prepare agree cells on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} agreePrepare failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True


def noisePrepare(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            qsubmode="local"):

    print("[Task  Added] prepare noise cells on {}...".format(runName))
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

    mountstring = mountgen([input_prefix, file_path, output_prefix])

    commandline = "Rscript {file_path}/sc-clustersf.R -i {inputPath} -o {od} -r {rn} -n {rn}.pca.clean_data.RData".format(
        file_path=file_path, inputPath=inputPath, od=od, rn=rn)
    
    qsub_sh = "{od}/shell/qsub.noisePrepare.{taskid}.sh".format(od=od, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts noise preparation and size factor calculation ..."
    singularity exec -B {mountstring} {singularityPath}/scrnaseq.sif {commandline}
    echo `date` "{rn} finishes noise preparation"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] prepare noise cells on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} noisePrepare failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True


def cluAnnotation(runName,
            pipelinelog,
            tasklogprefix,
            outputPath,
            inputPath,
            singularityPath,
            codePath,
            qsubmode="local"):

    print("[Task  Added] cluster annotation on {}...".format(runName))
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

    mountstring = mountgen([input_prefix, file_path, output_prefix])


    commandline = "python3 {file_path}/ori_cluster_annotation.py -i {inputPath} -o {od} -n {rn} -im user_marker.txt".format(
        file_path=file_path, inputPath=inputPath, od=od, rn=rn)
    
    qsub_sh = "{od}/shell/qsub.cluAnnotation.{taskid}.sh".format(od=od, taskid=runName)
    qsubfile = open(qsub_sh, "w")
    qsubfile.write('''#!/bin/sh
    set -e
    ls {inputPath} {od} > /dev/null

    echo `date` "{rn} starts noise preparation and size factor calculation ..."
    singularity exec -B {mountstring} {singularityPath}/advance.sif {commandline}
    echo `date` "{rn} finishes noise preparation"

    sleep 10 && echo "RUN COMPLETED!" && sync
    '''.format(inputPath=inputPath, od=od, rn=rn, mountstring=mountstring, singularityPath=singularityPath, commandline=commandline))

    qsubfile.close()
    os.system("chmod +x {}".format(qsub_sh))
    success = qsub(qsub_sh, 1, tasklogprefix, mode=qsubmode)

    if success:
        # self.is_complete = 1
        print("[Task Done] cluster annotation on {}.".format(rn))
    else:
        with open(pipelinelog, 'a') as pipelinelog:
            pipelinelog.write("{} cluAnnotation failed!\n".format(rn))
        print("[Task Failed] analysis on {} ...".format(rn))
        return False

    return True


