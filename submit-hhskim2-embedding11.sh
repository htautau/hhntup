#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim2.py -v 1 \
embedding11_ztautau_hh_p851_v1_skim embedding11_ztautau_incl_p851_v1_skim \
--official --voms=atlas:/atlas/phys-higgs/Role=production
