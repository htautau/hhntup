#!/bin/bash

source grid.setup && python run --output-path . -s HHProcessor.py -n 1 \
--db datasets_hh --nice 10 \
AlpgenJimmy_AUET2CTEQ6L1_WtaunuNp3.mc12a --redo-mmc
