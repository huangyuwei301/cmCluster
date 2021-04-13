#-*- encoding: UTF-8 -*-
import os
import sys
import re
import luigi
#from importlib import reload
#reload(sys)
#sys.setdefaultencoding('utf8')

class advanced(luigi.Config):
	shared = luigi.BoolParameter(default=True)

STYLE = {
	'fore':
		{
			'black': 30,
			'red'	: 31,
			'green'	: 32,
			'yellow'   : 33,
			'blue'	 : 34,
			'purple'   : 35,
			'cyan'	 : 36,
			'white'	: 37,
		},

	'back' :
		{
			'black'	 : 40,
			'red'	   : 41,
			'green'	 : 42,
			'yellow'	: 43,
			'blue'	  : 44,
			'purple'	: 45,
			'cyan'	  : 46,
			'white'	 : 47,
		},

	'mode' :
		{
			'mormal'	: 0,
			'bold'	  : 1,
			'underline' : 4,
			'blink'	 : 5,
			'invert'	: 7,
			'hide'	  : 8,
		},

	'default' :
		{
			'end' : 0,
		},
}



def qsub(subfile, thread, jobname, mode="local"):
	
	#if mode == "qsub":
		
	#	if queue().queuename != '':
	#		qstring = " -q {}".format(queue().queuename)
	#	else:
	#		qstring = " "

	#	if queue().nodename != '':
	#		lstring = " -l nodes={nodename}:ppn={t}".format(nodename=queue().nodename, t=thread)
	#	else:
	#		lstring = " -l nodes=1:ppn={t}".format(t=thread)
		
	#	print("    Submitting: sbatch -W {qstring} {lstring} {subfile} > {jobname}.log 2> {jobname}.err".format(qstring=qstring, lstring=lstring, subfile=subfile, jobname=jobname))
	#	sub = os.system("sbatch -W {qstring} {lstring} {subfile} > {jobname}.log 2> {jobname}.err".format(qstring=qstring, lstring=lstring, subfile=subfile, jobname=jobname))
	if mode == "qsub":
		if queue().queuename != '':
			qstring = " -p {}".format(queue().queuename)
		else:
			qstring = " "

		if queue().nodename != '':
			lstring = "--nodelist={nodename} --cpus-per-task={t}".format(nodename=queue().nodename, t=thread)
		else:
			lstring = "-N 1 --cpus-per-task={t}".format(t=thread)

		print("    Submitting: sbatch -W {qstring} {lstring} -o {jobname}.log -e {jobname}.err {subfile}".format(
			qstring=qstring, lstring=lstring, subfile=subfile, jobname=jobname))
		sub = os.system(
			"sbatch -W {qstring} {lstring} -o {jobname}.log -e {jobname}.err {subfile}".format(qstring=qstring,lstring=lstring,subfile=subfile,jobname=jobname))

		'''print("    Submitting: salloc {qstring} {lstring} /bin/sh {subfile} > {jobname}.log 2> {jobname}.err ".format(
			qstring=qstring, lstring=lstring, subfile=subfile, jobname=jobname))
		sub = os.system(
			"salloc {qstring} {lstring} /bin/sh {subfile} > {jobname}.log 2> {jobname}.err".format(qstring=qstring,lstring=lstring,subfile=subfile,jobname=jobname))'''
	elif mode == "local":
		
		print("    Submitting: sh {subfile} > {jobname}.log 2> {jobname}.err".format(subfile=subfile, jobname=jobname))
		sub = os.system("sh {subfile} > {jobname}.log 2> {jobname}.err".format(subfile=subfile, jobname=jobname))
	# print("sub value:",sub)	
	if sub != 0:
		# sys.exit("[Failed] qsub: {}". format(subfile))
		# input()
		return False
	success = False
	with open("{}.log".format(jobname), 'r') as statusfile:
		if mode == "local":
			statusstring = statusfile.readlines()[-1].strip()
		else:
			statusstring = statusfile.readlines()[-1].strip()
		if statusstring == 'RUN COMPLETED!':
			success = True

	return success


def mountgen(list, shared=True):
	mountdict = {}
	mountstring = ""
	shared = advanced().shared
	for i in list:
		if i == "":
			continue
		realpath = os.path.abspath(i)
		prefix = "/" + realpath.split("/")[1]
		mountdict[prefix] = "yes"
	
	for i in mountdict.keys():
		if shared:
			mountstring += " {i}:{i} ".format(i=i) #mountstring += " -v {i}:{i}:shared ".format(i=i)
		else:
			mountstring += " {i}:{i} ".format(i=i) #mountstring += " -v {i}:{i} ".format(i=i)
	return mountstring

