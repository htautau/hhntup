#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim.py -v 3 \
mc12_p1130_hadhad --outputs TPileupReweighting.prw.root  \
--official --voms=atlas:/atlas/phys-higgs/Role=production
