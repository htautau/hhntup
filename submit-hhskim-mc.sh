#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim.py -v 1 \
--outputs TPileupReweighting.prw.root mc11_p851_hadhad \
--official --voms=atlas:/atlas/phys-higgs/Role=production
