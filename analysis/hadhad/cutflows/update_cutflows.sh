#!/bin/bash

./cutflow.py --rst --db datasets_hh --short > cutflow.txt
./cutflow_mc_bkg.py --rst --db datasets_hh > cutflow_mc_bkg.txt
./cutflow_mc_sig.py --rst --db datasets_hh > cutflow_mc_sig.txt
