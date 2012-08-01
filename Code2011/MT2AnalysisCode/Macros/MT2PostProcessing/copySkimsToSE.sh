#!/bin/bash

if [ $# != 1 ]; then
    echo "USAGE: $0 SHORT_DIR (like: MT2_V01-05-01/20120116_MC_HT600_data_nocuts_PU2_GR_R_42_V23_AK5PF_UP)"
    exit
fi


DIR=$1
BASEDIR="/shome/leo/MT2Analysis/MT2trees/"
DIR=`echo ${DIR} | sed "s:${BASEDIR}::g"` 
DEST="/store/user/leo/SUSY/Skims/"


echo "Looking in ${BASEDIR}${DIR}"
for i in `ls ${BASEDIR}${DIR}/*.root | sed "s:${BASEDIR}::g"`; do
    echo "From: " ${BASEDIR}/$i;
    echo "To: " $SE_PSI${DEST}$i;
    srmcp -2 file:///${BASEDIR}/$i $SE_PSI${DEST}$i;
done

for i in `ls ${BASEDIR}${DIR}/*/* | sed "s:${BASEDIR}::g"`; do
    echo $i;
    echo "To: " $SE_PSI${DEST}$i;
    srmcp -2 file:///${BASEDIR}/$i $SE_PSI${DEST}$i;
done

