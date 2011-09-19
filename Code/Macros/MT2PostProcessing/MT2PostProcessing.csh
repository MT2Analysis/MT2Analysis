#!/bin/tcsh -f
##################################################################
# Script to skim MT2trees on T3 batch system                     #
# Pascal Nef                             September 18th, 2011    #  
# ################################################################

# Avoid wild cards
set noglob
############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

##### MONITORING/DEBUG INFORMATION ###############################
set DATE_START=`date +%s`
echo "Job started at " `date`

#### Check n agruments ###########################################
if ($#argv < 5) then
	echo "usage "$0" MT2tag production sample shlib file(s)"
	exit(-1)
endif

set NARGS      = $#argv
set MT2TAG     = $1
set PRODUCTION = $2
set SAMPLE     = $3
set SHLIB      = $4
set HOME=`pwd`
set OUTDIR=/shome/`whoami`/MT2Analysis/MT2trees/$MT2TAG/$PRODUCTION/
set SKIMSCRIPT=$HOME"/MT2treeSkimming.C"
set TOPWORKDIR=/scratch/`whoami`
set WORKDIR=/scratch/`whoami`/$PRODUCTION/$SAMPLE
set MERGED=/scratch/`whoami`/$PRODUCTION/$SAMPLE/$SAMPLE.root
mkdir -pv $WORKDIR

########## Set the CMSSW environment (for ROOT, mainly...)
source $VO_CMS_SW_DIR/cmsset_default.csh
setenv SCRAM_ARCH slc5_amd64_gcc434
cd /shome/pnef/SUSY/CMSSW_4_1_3/src/
cmsenv
# Fix problem with DCache and ROOT
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
cd -

#### copy scaples from SE to scratch ###########################################
set lfiles=()
foreach i ( $argv[5-$NARGS] )
	echo "copying file " $i
	set file = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/"$i
	dccp "$file" "$WORKDIR"
   	set lfiles = ($lfiles $file)
end

# Now merge
rm -vf $MERGED
hadd -f $MERGED $lfiles
if ( $? != 0 ) then
  echo Merging failed
  rm -vf $MERGED $lfiles
  exit -1
endif

#### start skimming ##########################################################
if ( $SHLIB != "none" ) then
	echo " -------------->do skimming <---------------------- "
	cd $WORKDIR
	cp  $SKIMSCRIPT .
	mkdir -pv skimmed
	set file=$SAMPLE.root
	echo root -l -b -q  'MT2treeSkimming.C("'$file'", "'$SHLIB'")'	
	root -l -b -q  'MT2treeSkimming.C("'$file'", "'$SHLIB'")'	
	if (-e skimmed/$file) then
		echo "skimming terminated successfully..."
	else
		echo "skimming failed!!! "
		rm -r $WORKDIR 
		exit(1)
	endif
	echo "copying skimmed MT2tree"
	mkdir -pv $OUTDIR/skimmed/skimlogs
	cp -v skimmed/$SAMPLE.root $OUTDIR/skimmed
	if (-e $SAMPLE.root.skim.log) then
		cp -v $SAMPLE.root.skim.log $OUTDIR/skimmed/skimlogs
	else
		echo "$SAMPLE.root.skim.log not found..."
	endif 
else
	echo "--------------> no skimming performed <-------------"
	mkdir -pv $OUTDIR	
	cp -v $MERGED $OUTDIR
endif
	
# clean up ################################################################
rm -r $WORKDIR 
echo -------------------- ALL DONE --------------------------------------
###########################################################################
set DATE_END=`date +%s`
@ RUNTIME = $DATE_END - $DATE_START
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0

