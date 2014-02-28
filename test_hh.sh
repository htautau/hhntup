#!/bin/bash

#source grid.setup && python run --output-path . -s hhskim.py -n 1 \
#--db datasets_hh --nice 10 \
#AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp3.mc12a $@

source grid.setup && python run --profile --output-path . -s hhskim.py -n 2 --db datasets_hh --nice 10 --suffix TES_TRUE_FINAL_UP AlpJimZtautauNp3_pt20.mc11c --syst-terms MET_SCALESOFTTERMS_UP --local
