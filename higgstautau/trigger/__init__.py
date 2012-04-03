import os
import ROOT

ROOT.gSystem.Load(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  'libTrigRootAnalysis.so'))
from ROOT import D3PD

__all__ = [
    'get_trigger_config',
    'update_trigger_config'
]


def get_trigger_config():

    return D3PD.TrigConfigSvcD3PD()


def update_trigger_config(tool, name, file):

    configTreeName = "%sMeta/TrigConfTree"% name
    configTree = file.Get(configTreeName)
    if not configTree:
        raise RuntimeError("Could not find %s in %s"% (configTreeName, file.GetName()))
    tool.SetConfigTree(configTree)

