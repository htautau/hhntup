#!/bin/bash

#./run_data_hh.py --year 12 --nsplit 50 --local --redo-selection

#./run_embed_hh.py --nominal-only --year 12 --nsplit 30 --local --redo-selection
#./run_mc_hh.py --nominal-only --year 12 --local --redo-selection

#./run_embed_hh.py --systematics-only --year 12 --nsplit 30 --local --redo-selection
#./run_mc_hh.py --systematics-only --year 12 --local --redo-selection

#./run_embed_hh.py --nominal-only --year 11 --redo-mmc
#./run_mc_hh.py --nominal-only --year 11 --redo-mmc
#./run_data_hh.py --year 11 --nsplit 40 --redo-mmc
#./run_embed_hh.py --systematics-only --year 11 --redo-mmc
#./run_mc_hh.py --systematics-only --year 11 --redo-mmc

./run_mc_hh.py --nominal-only --year 12 --local --redo-selection \
McAtNloJimmy_AUET2CT10_ggH125_tautau.mc12a \
PowhegJimmy_AUET2CT10_ggH125_tautauInclusive.mc12a \
