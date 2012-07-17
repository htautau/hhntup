#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 5 \
data11_p851_hadhad_v4_skim \
--antiMatch "*cutflow.p" \
--official --voms=atlas:/atlas/phys-higgs/Role=production
