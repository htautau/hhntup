#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 1 \
mc12_p1011_hadhad_v1_skim \
--antiMatch "*TPileupReweighting.prw.root*" \
--official --voms=atlas:/atlas/phys-higgs/Role=production
