#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 3 \
data12_p1130_hadhad_v3_skim \
--official --voms=atlas:/atlas/phys-higgs/Role=production
