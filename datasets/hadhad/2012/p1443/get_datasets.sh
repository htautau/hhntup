#!/bin/bash

TAG=p1443

ami list data --project data12_8TeV --type NTUP_TAU %$TAG | sed 's/$/\//g' > data_ami.txt
dq2-ls data12*physics_JetTauEtmiss*NTUP_TAU.*$TAG/ | sort > data_dq2.txt
