#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 3 \
mc12_p1130_hadhad_v3_skim \
--antiMatch "*TPileupReweighting.prw.root*" \
--official --voms=atlas:/atlas/phys-higgs/Role=production
