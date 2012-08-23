#! /bin/bash
dq2-ls data12_8TeV.*.physics_Egamma.merge.NTUP_TAU.*_p1130_*/ | sort > data12_Egamma_p1130.txt 
dq2-ls data12_8TeV.*.physics_Muons.merge.NTUP_TAU.*_p1130_*/ | sort > data12_Muons_p1130.txt
dq2-ls data12_8TeV.*.physics_JetTauEtmiss.merge.NTUP_TAU.*_p1130_*/ | sort > data12_JetTauEtmiss_p1130.txt