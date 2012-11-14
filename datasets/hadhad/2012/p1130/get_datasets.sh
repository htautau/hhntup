#!/bin/bash
ami list data --project data12_8TeV --type NTUP_TAU p1130_p | sed 's/$/\//g' > data_ami.txt
dq2-ls data12*physics_JetTauEtmiss*NTUP_TAU.*p1130_p*/ | sort > data_dq2.txt
