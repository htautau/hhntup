#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim.py -v 3 \
data12_p1130_hadhad \
--official --voms=atlas:/atlas/phys-higgs/Role=production
