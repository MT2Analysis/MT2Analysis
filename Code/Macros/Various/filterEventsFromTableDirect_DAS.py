#!/bin/env python

import os
from sys import argv, exit
from optparse import OptionParser

START_FILE=0


#### Filter on allowed datasets
ALLOWED_DATASETS=""

usage = argv[0]+"""[options] tableFile.txt 

This tool reads a generic run/lumi table and run a filter on files.
The table must be in the format:

**************************************************************************************************
*   Run  * Lumi *      Event   *   HT   *  MET  *  MT2  *  NJets * NJets40 * minDphi * hemiDphi *
**************************************************************************************************
* 163332 *   778 *    519897214 * 965.83 * 787.0 * 618.8 *      5 *       5 *   1.195 *    0.870 *

or similar. In case more than one file is found, it's run on all the files, and then duplicates are
rejected.

The datasets in which the files are looked into is hardcoded: 
"""+str(ALLOWED_DATASETS)+"""

In case you want to specify a different set of datasets, use the --dataset option

BEWARE: with DAS no AND|OR is supported, so you can only search a single dataset string (* wildcard allowed)

This tool can be run only WHERE THE DATASETS ARE DIRECTLY AVAILABLE. Running it from
lxplus is suggested. You have to check by yourself that the datasets are available directly
(dcap, xrootd, rfio...)
"""

parser = OptionParser(usage = usage)
parser.add_option("--dataset",
                  action="store", dest="datasets", default="",
                  help="Name of dataset to search in (* wildcard accepted, only one name)")
parser.add_option("--reco",
                  action="store_true", dest="reco", default=False,
                  help="Store events in RECO format, if possible (default is AOD)")
parser.add_option("--run",
                  action="store_true", dest="run", default=False,
                  help="Directly run the created CFG")



(options, args) = parser.parse_args()
if len(args)<1:
    parser.print_help()
    exit(1)

if options.datasets!="":
    ALLOWED_DATASETS = options.datasets

if ALLOWED_DATASETS=="":
    print "No dataset is given as input, please use --dataset DATASETNAME"
    exit(1)

PWD=os.getenv("PWD")
if not os.path.isfile(PWD+"/cli"):
    print "DAS Command Line utility not here, downloading"
    os.system("wget --no-check-certificate https://cmsweb.cern.ch/das/cli && chmod +x cli")
    

print "\n----------------- "
print "I will search in the following dataset:"
print ALLOWED_DATASETS
print "----------------- \n"

TABLE=args[0]
OUTPUT=TABLE.replace(".txt",".root")
CFG=TABLE.replace(".txt",".py")

RUN_LUMI=[]
FILES = []

####### parsing the table file
tfile = open(args[0])
run_i = 0
lumi_i=0
event_i=0
for line in tfile:
    line = line.strip()
    if len(line.strip(" ") )== 0: continue
    if line[0]=="#": continue
    if line.find("**")!=-1: continue
    sline = line.split("*")
    for i in range(0, len(sline)): sline[i] = sline[i].strip(" ")
    if "Run" in sline or "Row" in sline:
        for i in range(0, len(sline)):
            if sline[i].find("Run")!=-1: run_i=i
            elif sline[i].find("Lumi")!=-1: lumi_i=i
            elif sline[i].find("Event")!=-1: event_i=i
    elif len(sline[run_i]) !=0:
        RUN_LUMI.append( [sline[run_i], sline[lumi_i], sline[event_i] ])

if len(RUN_LUMI)==0:
    print "No Run, Lumi information found, please check your input file"
    exit(1)


### preparing the dbs query
DST = ALLOWED_DATASETS

### run/lumi loop
filelist=""
sizeConv = {"kB":1024, "MB":1024*1024, "GB":1024*1024*1024, "TB":1024*1024*1024*1024}
N_FILE=0
for rl in RUN_LUMI:
    if N_FILE< START_FILE: continue
    N_FILE +=1
    cmd="./cli --limit=0 --query='file dataset="+DST+" run="+str(rl[0])+" lumi="+str(rl[1])+"'"
    out = os.popen(cmd)
    found=False
    
    for l in out:
        l.strip().strip("\n")
        if l=="" or l.find("#")!=-1 or len(l)<2: continue
        file=l.split()[0]

        if file.find("/store")!=-1:
            filelist+="""'"""+file+"""',\n"""
            found=True
            print N_FILE,len(FILES),"Found file run=",str(rl[0])," lumi:",str(rl[1])," FILE=",file

    if not found: print "NO FILE FOUND FOR LUMI="+str(rl[1])+" AND RUN="+str(rl[0])

    
runevt=""
for r in RUN_LUMI:
        runevt+="'"+r[0]+":"+r[2]+"',\n"

format="process.AODSIMEventContent.outputCommands"
if options.reco: format="process.RECOSIMEventContent.outputCommands"


txt="""import FWCore.ParameterSet.Config as cms

process = cms.Process('FILTER')
process.load('FastSimulation.Configuration.EventContent_cff')
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(-1)
)

process.source = cms.Source(
'PoolSource',
fileNames = cms.untracked.vstring(
"""+filelist+"""
),

eventsToProcess = cms.untracked.VEventRange(
"""+runevt+"""
),
secondaryFileNames = cms.untracked.vstring(),
noEventSort = cms.untracked.bool(True),
#duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)

process.display = cms.OutputModule('PoolOutputModule',
outputCommands = """+format+""",
#outputCommands = cms.untracked.vstring('keep *'),
#process.RECOSIMEventContent,
fileName = cms.untracked.string('"""+OUTPUT+"""')
)


process.outpath = cms.EndPath(process.display)
"""

cfg = open(CFG,"w")
cfg.write(txt)
cfg.close()

        
if os.getenv("CMSSW_BASE")=="":
    print "CMSSW env not set up, set it up and then run: cmsRun "+CFG
    exit(0)
elif options.run:
    print "\nNow running: cmsRun "+CFG
    os.system("cmsRun "+CFG)
else:
    print "\nNow run: cmsRun "+CFG
    
