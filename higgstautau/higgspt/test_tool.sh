#!/bin/sh
exec=bin/test_ggF_tool.exe
mkdir -p bin
rm -f $exec
c++ $($ROOTSYS/bin/root-config --cflags --glibs) -lTreePlayer -I. test/test_ggF_tool.cxx -o $exec && ./$exec

#
####
