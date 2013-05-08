#!/bin/bash

./run_mc_hh.py --year 12 --redo-mmc
#./run_mc_hh.py --nominal-only --year 11
./run_data_hh.py --year 12 --nsplit 200 --redo-mmc
#./run_data_hh.py --year 11 --nsplit 40
./run_embed_hh.py --year 12 --redo-mmc
#./run_embed_hh.py --nominal-only --year 11
