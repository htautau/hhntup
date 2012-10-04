#! /bin/bash

# 2012 Official Tau D3PDs
##################################################
grid-submit -s LHSkim.py -u group.phys-higgs -l small --official --destSE=TOKYO-LCG2_PHYS-HIGGS,NIKHEF-ELPROD_PHYS-HIGGS -v 6 -m datasets.cfg mc12_p1130_lephad

grid-submit -s LHSkim.py -u group.phys-higgs -l small --official --destSE=TOKYO-LCG2_PHYS-HIGGS,NIKHEF-ELPROD_PHYS-HIGGS -v 6 -m datasets.cfg embedded_p1130_lephad

grid-submit -s LHSkim.py -u group.phys-higgs -l small --official --destSE=TOKYO-LCG2_PHYS-HIGGS,NIKHEF-ELPROD_PHYS-HIGGS -v 6 -m datasets.cfg data12_p1130_muhad
grid-submit -s LHSkim.py -u group.phys-higgs -l small --official --destSE=TOKYO-LCG2_PHYS-HIGGS,NIKHEF-ELPROD_PHYS-HIGGS -v 6 -m datasets.cfg data12_p1130_ehad
grid-submit -s LHSkim.py -u group.phys-higgs -l small --official --destSE=TOKYO-LCG2_PHYS-HIGGS,NIKHEF-ELPROD_PHYS-HIGGS -v 6 -m datasets.cfg data12_p1130_leptau

# 2012 Testing Tau D3PDs
##################################################
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=TOKYO-LCG2_PHYS-HIGGS -v 2 -m datasets.cfg mc11c_p851_lephad_signal
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=TOKYO-LCG2_PHYS-HIGGS -v 2 -m datasets.cfg mc12_p1011_lephad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=UNI-FREIBURG_PHYS-HIGGS -v 2 -m datasets.cfg mc12_p1011_SVTF_lephad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=UNI-FREIBURG_PHYS-HIGGS -v 2 -m datasets.cfg embedded_p851_lephad
#grid-submit -s LHSkim.py -u group.phys-higgs -l small --official --destSE=TOKYO-LCG2_PHYS-HIGGS,UNI-FREIBURG_PHYS-HIGGS,NIKHEF-ELPROD_PHYS-HIGGS -v 2 -m datasets.cfg embedded_p851_lephad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=NIKHEF-ELPROD_PHYS-HIGGS -v 2 -m datasets.cfg data12_p1011_muhad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=NIKHEF-ELPROD_PHYS-HIGGS -v 2 -m datasets.cfg data12_p1011_ehad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=NIKHEF-ELPROD_PHYS-HIGGS -v 2 -m datasets.cfg data12_p1011_leptau
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=TOKYO-LCG2_PHYS-HIGGS -v 2 -m datasets.cfg data12_SVTF_muhad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=NIKHEF-ELPROD_PHYS-HIGGS -v 2 -m datasets.cfg data12_SVTF_ehad
# grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=UNI-FREIBURG_PHYS-HIGGS -v 2 -m datasets.cfg data12_SVTF_leptau

# 2011 Official Tau D3PDs
##################################################
#grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=NIKHEF-ELPROD_PHYS-HIGGS -v 1 -m datasets.cfg data11_p851_leptau
#grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=NIKHEF-ELPROD_PHYS-HIGGS -v 1 -m datasets.cfg data11_p851_muhad
#grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=UNI-FREIBURG_PHYS-HIGGS -v 1 -m datasets.cfg data11_p851_ehad
#grid-submit -s LHSkim.py -u mtm -p phys-higgs -l small --destSE=TOKYO-LCG2_PHYS-HIGGS -v 1 -m datasets.cfg mc11c_p851_lephad