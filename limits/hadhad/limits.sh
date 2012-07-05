#! /bin/bash
source /atlas/software/bleedingedge/root-5.32-patches-64/bin/thisroot.sh


#//main
#void runAsymptoticsCLs(const char* infile,
#		       const char* workspaceName,
#		       const char* modelConfigName,
#		       const char* dataName,
#		       const char* asimovDataName,
#		       string folder,
#		       string mass,
#		       double CL);

rm -f bandPlotData_combined.txt
rm -f bandPlotData_ggf.txt
rm -f bandPlotData_boosted.txt
rm -f bandPlotData_vbf.txt

for mass in $(seq 100 5 150)
do
    echo "*************************************"
    echo "Calculating limits for mass point ${mass}"
    
    for category in ggf boosted vbf
    do
        echo "--------------------------------------"
        echo "Category: ${category}"
        root -l -b -q runAsymptoticsCLs.C+"(\
            \"./results/hh_${category}_${mass}_combined_STATONLY_model.root\",\
            \"combined\",\
            \"ModelConfig\",\
            \"asimovData\",\
            \"asimovData_0\",\
            \"hadhad\",\
            \"${mass}\",\
            0.95)" > log${mass}_${category}.txt
        echo ${mass} >> bandPlotData_${category}.txt
        grep -A 6 -h "Correct bands" log${mass}_${category}.txt >> bandPlotData_${category}.txt
    done
    
    echo "--------------------------------------"
    echo "Combination"
    root -l -b -q runAsymptoticsCLs.C+"(\
        \"./results/hh_${mass}_combined_STATONLY_model.root\",\
        \"combined\",\
        \"ModelConfig\",\
        \"asimovData\",\
        \"asimovData_0\",\
        \"hadhad\",\
        \"${mass}\",\
        0.95)" > log${mass}_combined.txt
    echo ${mass} >> bandPlotData_combined.txt
    grep -A 6 -h "Correct bands" log${mass}_combined.txt >> bandPlotData_combined.txt
done
