#!/bin/bash

if [ $1 -eq 'test' ]
then

grid-batch --dataset embedding11_ztautau_p851_v1_skim --metadata datasets.cfg \
--student HHSkim2 /global/endw/embedding11_7TeV/test

else

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 1 \
embedding11_ztautau_p851_v1_skim \
--official --voms=atlas:/atlas/phys-higgs/Role=production

fi
