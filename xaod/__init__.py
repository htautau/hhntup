# Set up ROOT and RootCore:
import os
import sys

BASE_DIR = os.getenv('DIR_HIGGSTAUTAU_SETUP')
if not BASE_DIR:
    sys.exit('You did not source setup.sh!')

CACHE_DIR = os.path.join(BASE_DIR, 'cache')

import ROOT
ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
# Initialize the xAOD infrastructure: 
RC = ROOT.xAOD.Init()
if not RC.isSuccess():
    raise RuntimeError('Error in the xAOD initialization')

ROOT.xAOD.AuxContainerBase()

# Set up the input files:
ftemp = ROOT.TFile(os.path.join(CACHE_DIR, 'xaod_struct.root'))
ttemp = ROOT.xAOD.MakeTransientTree(ftemp)
ftemp.Close()

TOOLS = ROOT.asg.ToolStore()
