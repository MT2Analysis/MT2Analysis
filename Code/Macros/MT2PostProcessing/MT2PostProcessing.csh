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
set PREFIX     = $5
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
cd /shome/pnef/CMSSW/CMSSW_4_2_8/src
cmsenv
# Fix problem with DCache and ROOT
setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:$LD_LIBRARY_PATH
cd -

#### copy scaples from SE to scratch ###########################################
set lfiles=()
foreach i ( $argv[6-$NARGS] )
	echo "found file " $i
	set file = "dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/"$i
   	set lfiles = ($lfiles $file)
end

#### start skimming ##########################################################
if ( $SHLIB != "none" ) then
	echo " -------------->do skimming <---------------------- "
	echo "with root version " `which root`
	cd $WORKDIR
	cp  $SKIMSCRIPT .
	mkdir -pv $PREFIX
	set lnewfiles=()
	foreach i ($lfiles)
		echo "processing $i ++++++++++++++++++++++++++++++++++++++++++++++++++++"
		echo "copying file $i to $WORKDIR"
		dccp "$i" "$WORKDIR"
		set file = `echo $i | awk -F '/' '{print $(NF)}'`
		root -l -b -q  'MT2treeSkimming.C("'$file'", "'$SHLIB'", "'$PREFIX'")'	
		set lnewfiles = ($lnewfiles $PREFIX/$file)
		if (-e $PREFIX/$file) then
			echo "skimming terminated successfully..."
			rm -fv $file
		else
			echo "skimming failed!!! "
			rm -r $WORKDIR 
			exit(1)
		endif
	end
	echo "merging skimmed MT2trees ++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	rm -vf $SAMPLE.root
	echo "files: $lnewfiles"
	hadd -f $SAMPLE.root $lnewfiles
	if ( $? != 0 ) then
	  echo Merging failed
	  rm -vf $SAMPLE.root $lnewfiles
	  exit -1
	endif
	echo "copying merged skimmed MT2tree"
	mkdir -pv $OUTDIR/$PREFIX/skimlogs
	cp -v $SAMPLE.root $OUTDIR/$PREFIX
	if (-e $SAMPLE.root.skim.log) then
		cp -v $SAMPLE.root.skim.log $OUTDIR/$PREFIX/skimlogs
	else
		echo "$SAMPLE.root.skim.log not found..."
	endif 
else
	# Now merge from SE
	rm -vf $MERGED
	hadd -f $MERGED $lfiles
	if ( $? != 0 ) then
	  echo Merging failed
	  rm -vf $MERGED $lfiles
	  exit -1
	endif
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

