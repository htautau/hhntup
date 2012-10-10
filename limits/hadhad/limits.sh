#! /bin/bash
source /cluster/data10/software/root-5.32-patches-64/bin/thisroot.sh


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

TAG=00-01-02
WORKSPACE=FULL_SYSTEMATICS

echo "Using tag runAsymptoticsCLs-"${TAG}

for mass in $(seq 115 5 150)
do
    echo "*************************************"
    echo "Calculating limits for mass point ${mass}"
    
    for category in ggf boosted vbf
    do
        echo "--------------------------------------"
        echo "Category: ${category}"
        echo "--------------------------------------"

        if [ "$TAG" == "old" ]
        then
            echo "Running the old limit script"
            (root -l -b -q ../common/runAsymptoticsCLs-old/runAsymptoticsCLs.C+"(\
                \"./results/hh_${category}_${mass}_combined_${WORKSPACE}_model.root\",\
                \"combined\",\
                \"ModelConfig\",\
                \"asimovData\",\
                \"asimovData_0\",\
                \"conditionalGlobs_0\",\
                \"nominalGlobs\",\
                \"${mass}\",\
                ${mass},\
                0.95)" 2>&1) | tee log${mass}_${category}.txt
        else
            (root -l -b -q ../common/runAsymptoticsCLs-${TAG}/runAsymptoticsCLs.C+"(\
                \"./results/hh_${category}_${mass}_combined_${WORKSPACE}_model.root\",\
                \"combined\",\
                \"ModelConfig\",\
                \"asimovData\",\
                \"asimovData_0\",\
                \"hadhad\",\
                \"${mass}\",\
                0.95)" 2>&1) | tee log${mass}_${category}.txt
            echo ${mass} >> bandPlotData_${category}.txt
            grep -A 6 -h "Correct bands" log${mass}_${category}.txt >> bandPlotData_${category}.txt
        fi
    done
    
    echo "--------------------------------------"
    echo "Combination"
    echo "--------------------------------------"
    
    if [ "$TAG" == "old" ]
    then
        echo "Running the old limit script"
        (root -l -b -q ../common/runAsymptoticsCLs-old/runAsymptoticsCLs.C+"(\
            \"./results/hh_${mass}_combined_${WORKSPACE}_model.root\",\
            \"combined\",\
            \"ModelConfig\",\
            \"asimovData\",\
            \"asimovData_0\",\
            \"conditionalGlobs_0\",\
            \"nominalGlobs\",\
            \"${mass}\",\
            ${mass},\
            0.95)" 2>&1) | tee log${mass}_combined.txt
    else
        (root -l -b -q ../common/runAsymptoticsCLs-${TAG}/runAsymptoticsCLs.C+"(\
            \"./results/hh_${mass}_combined_${WORKSPACE}_model.root\",\
            \"combined\",\
            \"ModelConfig\",\
            \"asimovData\",\
            \"asimovData_0\",\
            \"hadhad\",\
            \"${mass}\",\
            0.95)" 2>&1) | tee log${mass}_combined.txt
        echo ${mass} >> bandPlotData_combined.txt
        grep -A 6 -h "Correct bands" log${mass}_combined.txt >> bandPlotData_combined.txt
    fi
done
