#!/bin/bash

#grid-batch --dataset data --metadata datasets.yml --student HTauProcessor ../testdata/data
#grid-batch --dataset VBFH130hh --metadata SM_datasets.yml --student HTauProcessor ../../testdata/VBFH130_tautauhh
grid-batch --dataset VBFH100hh --metadata datasets.yml --student HTauProcessor --events 1000 ../../testdata/Ztautau
