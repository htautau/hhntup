#!/bin/bash

grid-submit -u group.phys-higgs -m datasets.cfg -s HHSkim.py -v 1 \
embedding11_ztautau_hh_p851 embedding11_ztautau_incl_p851 \
--official --voms=atlas:/atlas/phys-higgs/Role=production \
--site ANALY_SFU_bugaboo
