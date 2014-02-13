#!/bin/bash

nohup ./skim --local-test --yall embed11_hadhad > embed11_hadhad.log &
nohup ./skim --local-test --yall mc11_hadhad > mc11_hadhad.log & 
nohup ./skim --local-test --yall data11_hadhad > data11_hadhad.log &

nohup ./skim --local-test --yall embed12_hadhad > embed12_hadhad.log &
nohup ./skim --local-test --yall mc12_hadhad > mc12_hadhad.log &
nohup ./skim --local-test --yall data12_hadhad > data12_hadhad.log &
