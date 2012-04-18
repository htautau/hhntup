#!/bin/bash

/cluster/data05/michel/NoelUtilities/atlastools/scripts/batch --nproc 16 --metadata datasets.cfg --student LHSkim data_muhad
/cluster/data05/michel/NoelUtilities/atlastools/scripts/batch --nproc 16 --metadata datasets.cfg --student LHSkim data_ehad
/cluster/data05/michel/NoelUtilities/atlastools/scripts/batch --nproc 16 --metadata datasets.cfg --student LHSkim Ztautau_lephad
#grid-batch --dataset Ztautau_lephad --metadata datasets.cfg --student muLHSkim /global/mtm/data/D3PD/MC/testSkim
#grid-batch --dataset data_ehad --metadata datasets.cfg --student muLHSkim /global/mtm/data/D3PD/DATA/
#grid-batch --dataset Ztautau_lephad --metadata datasets.cfg --student eHSkim /global/mtm/data/D3PD/MC
