#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim.py -v 1 \
data12_p1011_hadhad \
--official --voms=atlas:/atlas/phys-higgs/Role=production
