#! /bin/bash

source /atlas/software/bleedingedge/root-5.32-patches-64/bin/thisroot.sh

for mass in $(seq 100 5 150)
do
    hist2workspace config/hh_combination_${mass}.xml
    for category in ggf boosted vbf
    do
        hist2workspace config/hh_combination_${category}_${mass}.xml
    done
done
