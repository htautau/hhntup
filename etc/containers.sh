#
# source this script to extract a list of datasets from their container
# and add that list to a new container
#
newcontainer="group.phys-higgs.HHSkim2.data12_8TeV.periodsCD/"
#
# which containers do we want to open?
#
export datasets=`dq2-ls group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*00207*NTUP_TAU*p1130.v3.v3/`
export datasets=$datasets" "`dq2-ls group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*00208*NTUP_TAU*p1130.v3.v3/`
export datasets=$datasets" "`dq2-ls group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*002069*NTUP_TAU*p1130.v3.v3/`
export datasets=$datasets" "`dq2-ls group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*002067*NTUP_TAU*p1130.v3.v3/`
export datasets=$datasets" "`dq2-ls group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*002068*NTUP_TAU*p1130.v3.v3/`
#
# loop over dataset containers, extract datasets from each
#
export listofsets=""
for dataset in $datasets
do
 listofsets=$listofsets" "`dq2-list-datasets-container $dataset`
done
echo "--------------------------------------------------------"
echo $listofsets
echo "adding to container " $newcontainer
dq2-register-datasets-container $newcontainer $listofsets