import os
import ROOT

ROOT.gSystem.Load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  'libTauDiscriminant.so'))
from ROOT import TauID

MethodBDT = TauID.MethodBDT
Error = TauID.Error
Types = TauID.Types
