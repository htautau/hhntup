#!/bin/bash

#./run_data_hh.py --year 12 --nsplit 200
#./run_data_hh.py --year 11 --nsplit 40
./run_embed_hh.py --year 12
#./run_embed_hh.py --nominal-only --year 11
./run_mc_hh.py --year 12
#./run_mc_hh.py --nominal-only --year 11
