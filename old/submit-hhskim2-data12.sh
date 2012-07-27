#!/bin/bash

if [ $1 -eq 'test' ]
then

grid-batch --dataset data12_p1130_hadhad_v2_skim --metadata datasets.cfg \
--student HHSkim2 \
/global/oneil/hadhad/test/group.phys-higgs.HHSkim.data12_8TeV.00203680.physics_JetTauEtmiss.merge.NTUP_TAU.f446_m1148_p1130.v2.120720151724

else

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 3 \
data12_p1130_hadhad_v3_skim \
--official --voms=atlas:/atlas/phys-higgs/Role=production

fi
