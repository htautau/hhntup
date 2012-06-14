#!/bin/bash

TAG=1011
MERGE=1015
OUT=data_p${TAG}_p${MERGE}.txt
GRL=../../../../grl/2012/data12_8TeV.periodAllYear_DetStatus-v45-pro13_CoolRunQuery-00-04-08_Higgs_tautau_lh.xml

rm -f $OUT

allds=`ami list data --grl $GRL --project data12_8TeV --type NTUP_TAUMEDIUM --parent-type AOD p${TAG} | sed "s/$/_p${MERGE}\//g"`

for ds in $allds
do
    mds=`dq2-ls $ds`
    if [ -n "$mds" ]
    then
        echo $mds
        echo $mds >> $OUT
    fi
done
