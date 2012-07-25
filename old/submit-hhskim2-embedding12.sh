#!/bin/bash

if [ $1 -eq 'test' ]
then

grid-batch --dataset embedding12_ztautau_hh_p1130_v1_skim \
--metadata datasets.cfg \
--student HHSkim2 /global/endw/mc12_8TeV/skimtest

else

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 1 \
embedding12_ztautau_hh_p1130_v1_skim \
--official --voms=atlas:/atlas/phys-higgs/Role=production

fi
