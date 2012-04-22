import os
import ROOT

#ROOT.gSystem.Load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
#                  'libTauTriggerCorrections.so'))
from externaltools import TauTriggerCorrections

from ROOT import TauTriggerCorrections

__all__ = []

tool = TauTriggerCorrections()

tool.loadInputFile("../root/triggerSF_wmcpara_EF_tau20_medium1.root", "1P3P","BDTm")

"""
The scale factor for the ditau triggered event is the product of the scale
factor of each tau.

In an MC event if a reco tau matches a true hadronic tau within dR<.2
then use the scale factor otherwise use the fake rate.
"""

def get_scale_factor(*args):

    return tool.getSF(*args)
