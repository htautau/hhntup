#!/bin/bash

for cont in `cat data11.txt`
do
    echo $cont
    for ds in `dq2-list-datasets-container $cont`
    do
        echo '\t'$ds
        dq2-register-datasets-container 'user.NoelDawe.data11_7TeV.physics_JetTauEtmiss.NTUP_TAUMEDIUM.p741/' $ds 
    done
done
