#!/bin/bash

grid-batch --dataset data --metadata rel17_datasets.yml --student HTauSkim ../../testdata/rel17data
grid-batch --dataset Ztautau --metadata rel17_datasets.yml --student HTauSkim ../../testdata/r17ztautau
