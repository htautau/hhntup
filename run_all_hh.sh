#!/bin/bash

./run_data_hh.py --year 11 --nsplit 40
./run_data_hh.py --year 12 --nsplit 80
./run_mc_hh.py --nominal-only --year 11
./run_mc_hh.py --nominal-only --year 12
