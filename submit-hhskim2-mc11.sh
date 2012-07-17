#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 2 \
mc11_p851_hadhad_v1_skim \
--antiMatch "*TPileupReweighting.prw.root*" \
--official --voms=atlas:/atlas/phys-higgs/Role=production
