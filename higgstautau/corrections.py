import re
from .units import GeV
from externaltools import ggFReweighting
ggFResources = ggFReweighting.RESOURCE_PATH
from ROOT import ggFReweighting

#################################################
#ggF reweighting
#################################################


ggF_tool = dict([
    (mass, ggFReweighting("PowHeg", mass, "Mean", ggFResources, "mc11"))
    for mass in range(100, 155, 5)])

DS_PATTERN = re.compile('ggH(?P<mass>\d+)')


def reweight_ggf(event, dataname):
    """
    Reweight the ggF samples
    """

    match = re.search(DS_PATTERN, dataname)
    if not match:
        return 1.

    # Get corresponding ggF tool setting
    try:
        ggFTool = ggF_tool[int(match.group('mass'))]
    except KeyError:
        raise ValueError("Higgs mass %d is not handled by ggF tool..."
                % mass)

    # Find the Higgs particle in mc
    pt = 0
    higgs = None
    for mc in event.mc:
        if mc.pdgId == 25:
            # get the Higgs just before the decay
            higgs = mc.last_self
            break
    if higgs is None:
        raise RuntimeError('could not find the Higgs!')
    pt = higgs.pt / GeV
    if pt > 0:
        weight = ggFTool.getWeight(pt)
    else:
        weight = 1.
    return weight
