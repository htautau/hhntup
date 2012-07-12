#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 5 \
data12_p1011_hadhad_v1_skim \
--official --voms=atlas:/atlas/phys-higgs/Role=production
