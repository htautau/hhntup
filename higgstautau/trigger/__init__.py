import os
import ROOT

from externaltools import TrigRootAnalysis
from ROOT import D3PD


__all__ = [
    'get_trigger_config',
    'update_trigger_config'
]


def get_trigger_config():

    return D3PD.TrigConfigSvcD3PD()


def update_trigger_config(tool, name, file, tree):

    configTreeName = "%sMeta/TrigConfTree"% name
    configTree = file.Get(configTreeName)
    if not configTree:
        raise RuntimeError("Could not find %s in %s"% (configTreeName, file.GetName()))
    tool.SetConfigTree(configTree)

