import os
import ROOT

try:
    from externaltools import PileupReweighting
except ImportError:
    ROOT.gSystem.Load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  'libPileupReweighting.so'))
from ROOT import Root

TPileupReweighting = Root.TPileupReweighting
