#!/bin/bash

#grid-batch --dataset data --metadata htau_datasets.yml --student HTauProcessor ../testdata/data

grid-batch --dataset VBFH100hh --metadata datasets.yml --student HTauProcessor ../testdata/Ztautau

#grid-batch --dataset Ztautau25 --metadata tauid_datasets.yml --student TauIDProcessor ../testdata/Ztautau
