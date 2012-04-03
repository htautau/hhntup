#!/bin/bash

batch -n 10 --dataset VBFH130hh_SM --metadata SM_datasets.yml --student HTauProcessor VBFH130hh_SM
batch -n 10 --dataset VBFH130hh --metadata datasets.yml --student HTauProcessor VBFH130hh
