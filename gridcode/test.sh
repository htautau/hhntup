#!/bin/bash

#grid-batch --dataset data --metadata rel17_datasets.yml --student HTauProcessor --events 1000 ../../testdata/rel17data

#grid-batch --dataset data --metadata datasets.cfg --student HTauProcessor --events 1000 ../../testdata/data
#grid-batch --dataset VBFH130hh --metadata SM_datasets.yml --student HTauProcessor ../../testdata/VBFH130_tautauhh
grid-batch --dataset VBFH100hh --metadata datasets.cfg --student HTauProcessor --events 1000 ../../testdata/Ztautau
