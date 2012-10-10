#! /bin/bash

source /cluster/data10/software/root-5.32-patches-64/bin/thisroot.sh
cp $ROOTSYS/etc/HistFactorySchema.dtd ./config

rm -f hist2workspace.log

for mass in $(seq 115 5 150)
do
    (hist2workspace ./config/hh_combination_${mass}.xml 2>&1) | tee --append hist2workspace.log
    for category in ggf boosted vbf
    do
        (hist2workspace ./config/hh_combination_${category}_${mass}.xml 2>&1) | tee --append hist2workspace.log
    done
done
