#!/bin/bash

TAG=p1444

ami list data --project data12_8TeV --type NTUP_TAU %physics_HadDelayed%$TAG | sed 's/$/\//g' > data_ami.txt
dq2-ls data12*physics_HadDelayed*NTUP_TAU.*$TAG/ | sort > data_dq2.txt
