#!/usr/bin/python
##################################################################
# Script to copy MT2trees from T3 SE to the user's home          #
# -> MT2trees can be skimmed according to specified cuts.        #
#                                                                #  
# Pascal Nef                             September 18th, 2011    #  
# ################################################################
import sys, os, commands, shlex
from optparse import OptionParser
from sys import argv,exit

def MT2PostProcessing():
	# parse arguments
	usage = """
	usage: %prog [options] file"
	
	This script will copy MT2trees from the T3 SE to the user's home directory 
	using the T3 batch system.
	-> use the --skim option to skim the MT2trees according to your 
	   preferred selection cuts. 

	Argument:
	file=lfn of the sample on the T3 SE
	-> in case the specified path contains subdirectories, the script will process
	   the all subdirectories contraining MT2trees

	Options:
	 --skim=shlib: set path to shlib of MT2trees that was used to generate the MT2trees. 
	               make sure the shlib corresponds to the correct MT2tag!
	               in case this option is set, the MT2trees will be skimmed according 
		       to the cuts specified in 'MT2treeSkimming.C'
	 --verbose:    set this option to get more prints
	 --dryRun:     if this option is set, no jobs to the T3 batch system
	               will be submitted. 

	Example:
	./MT2PostProcessing.py --skim=~/MT2Analysis/Code/MT2AnalysisCode_V01-00-00/MT2Code/shlib/libDiLeptonAnalysis.so [file]
	with [file]=/store/user//pnef/SUSY/MassTrees/MT2_V01-00-00/20110915_MSSM_MC_nocuts/

	"""
	parser = OptionParser(usage)

	parser.add_option("--skim"   ,dest="shlib",      help="shlib: path to shlib used to generate the MT2trees")
	parser.add_option("--verbose",dest="verbose",    help="set verbose", action="store_true")
	parser.add_option("--dryRun", dest="dryRun",     help="dry Run: don't do anything", action="store_true")

	global options, args, MT2tag
	(options, args) = parser.parse_args()

	if(options.verbose): print options
	
	if(len(args)!=1):
		print "exactly one argument must be given!"
		print "try --help option"
		exit(-1)

	FILE   = args[0]

	MT2tag = FILE[FILE.find("/MT2_V")+1:FILE.find("/",FILE.find("/MT2_V")+1)]
	if(len(MT2tag)!=13):
		print "error parsing MT2_tag from argument"
		exit(-1)
	else:
		print "setting MT2tag= " +MT2tag
	
	# get list of files
	t3_se_dir       ="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/"
	command   ="srmls "+ t3_se_dir + FILE+ " | grep -v .root | awk '{print $2}' | awk -F trivcat '{print $2}'"
	status, output  = commands.getstatusoutput(command)
	if options.verbose: 
		print output  
	subdirs =output.splitlines()
	if(len(subdirs) >1): subdirs.pop(0)  # remove base dir from list in case there are sub-directories 

	for subdir in subdirs:
		command   ="srmls "+ t3_se_dir + subdir +" | grep .root | awk '{print $2}' | awk -F trivcat '{print $2}'"
		status, files  = commands.getstatusoutput(command)
		if(options.verbose):
			print files
			print ""
		filelist =""
		for file in files.splitlines():
			filelist += file + " " 
		os.system("rm -rf "+logfile(subdir))
		cmd = "qsub -q short.q -N "+sample(subdir)+" -o "+logfile(subdir)+" -j y  ./MT2PostProcessing.csh "
		cmd = cmd + " "+ MT2tag + " "  + production(subdir) +  " "  + sample(subdir) 
		if(options.shlib!= None):
			cmd = cmd+ " " + options.shlib
		else:
			cmd = cmd+ " " + "none"
		cmd = cmd + " "+ filelist
		print cmd
		if not options.dryRun:
			os.system(cmd)
		print "---------------------------------------------------------------------"


def logfile(subdir):
	logname = subdir.split(MT2tag)[1]
	logname = logname.replace("/","_")
	logname = logname.strip("_")
	logname = "logs/"+logname+".log"
	return logname

def production(subdir):
	production = subdir.split(MT2tag)[1]
	production = production.strip("/")
	production = production.split("/")[0]
	return production

def sample(subdir):
	sample = subdir.split(production(subdir))[1]
	sample = sample.strip("/")
	return sample

if __name__ == "__main__":
	MT2PostProcessing()
