#!/bin/bash

grid-batch --dataset Ztautau_hadhad --metadata datasets.cfg --anti-match "*TPileupReweighting.prw.root*" --student HHSkim2 /global/endw/mc11_7TeV/higgs_tautau_hh_reskim_p851/user.NoelDawe.HHSkim.mc11_7TeV.107675.AlpgenJimmyZtautauNp5_pt20.merge.NTUP_TAUMEDIUM.e835_s1299_s1300_r3043_r2993_p851.v7.120501172721
grid-batch --dataset data11_hadhad --metadata datasets.cfg --student HHSkim2 /global/endw/data11_7TeV/higgs_tautau_hh_skim_p851/user.NoelDawe.HTauSkim.data11_7TeV.00191933.physics_JetTauEtmiss.merge.NTUP_TAUMEDIUM.f415_m1025_p851.v4.120204172411
