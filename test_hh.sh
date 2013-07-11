#!/bin/bash

source grid.setup && python run --output-path . -s hhskim.py -n 1 \
--db datasets_hh --nice 10 \
AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp3.mc12a $@

#source grid.setup && python run --output-path . -s hhskim.py -n 1 \
#--db datasets_hh --nice 10 \
#data12-Muons-207528 $@
