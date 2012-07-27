#!/bin/bash

if [ $1 -eq 'test' ]
then

grid-batch --dataset data11_hadhad --metadata datasets.cfg --student HHSkim2 \
/global/endw/data11_7TeV/higgs_tautau_hh_skim_p851/user.NoelDawe.HTauSkim.data11_7TeV.00191933.physics_JetTauEtmiss.merge.NTUP_TAUMEDIUM.f415_m1025_p851.v4.120204172411

else

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 5 \
data11_p851_hadhad_v4_skim \
--antiMatch "*cutflow.p" \
--official --voms=atlas:/atlas/phys-higgs/Role=production

fi
