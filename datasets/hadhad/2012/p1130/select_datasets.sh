#!/bin/bash
# select runs in the current GRL
grl runs ../../../../grl/2012/current.xml > /tmp/grlruns.txt
grep -f /tmp/grlruns.txt data_dq2.txt | grep -v v100 > data.txt
