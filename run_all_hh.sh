#!/bin/bash

./run_embed_hh.py --nominal-only --year 12 --redo-mmc
./run_mc_hh.py --nominal-only --year 12 --redo-mmc
./run_data_hh.py --year 12 --nsplit 200 --redo-mmc
./run_embed_hh.py --systematics-only --year 12 --redo-mmc
./run_mc_hh.py --systematics-only --year 12 --redo-mmc

#./run_embed_hh.py --nominal-only --year 11 --redo-mmc
#./run_mc_hh.py --nominal-only --year 11 --redo-mmc
#./run_data_hh.py --year 11 --nsplit 40 --redo-mmc
#./run_embed_hh.py --systematics-only --year 11 --redo-mmc
#./run_mc_hh.py --systematics-only --year 11 --redo-mmc
