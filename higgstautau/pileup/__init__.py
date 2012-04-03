import os
import ROOT

ROOT.gSystem.Load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  'libPileupReweighting.so'))
from ROOT import Root
# what's with the stupid T in the name guys? Come on... No need to propagate
# dumb ROOT naming conventions...
PileupReweighting = Root.TPileupReweighting
