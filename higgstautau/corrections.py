import re
from atlastools.units import GeV
from externaltools import ggFReweighting

ggFResources = ggFReweighting.RESOURCE_PATH

#################################################
#ggF reweighting
#################################################

from ROOT import ggFReweighting

ggF_tool = {}
ggF_tool[100] = ggFReweighting("PowHeg", 100, "Mean", ggFResources, "mc11")
ggF_tool[105] = ggFReweighting("PowHeg", 105, "Mean", ggFResources, "mc11")
ggF_tool[110] = ggFReweighting("PowHeg", 110, "Mean", ggFResources, "mc11")
ggF_tool[115] = ggFReweighting("PowHeg", 115, "Mean", ggFResources, "mc11")
ggF_tool[120] = ggFReweighting("PowHeg", 120, "Mean", ggFResources, "mc11")
ggF_tool[125] = ggFReweighting("PowHeg", 125, "Mean", ggFResources, "mc11")
ggF_tool[130] = ggFReweighting("PowHeg", 130, "Mean", ggFResources, "mc11")
ggF_tool[135] = ggFReweighting("PowHeg", 135, "Mean", ggFResources, "mc11")
ggF_tool[140] = ggFReweighting("PowHeg", 140, "Mean", ggFResources, "mc11")
ggF_tool[145] = ggFReweighting("PowHeg", 145, "Mean", ggFResources, "mc11")
ggF_tool[150] = ggFReweighting("PowHeg", 150, "Mean", ggFResources, "mc11")

DS_PATTERN = re.compile('PowHegPythia_ggH(?P<mass>\d+)_tautau')


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
    for mc in event.mc:
        if mc.pdgId == 25 and mc.status != 3:
            pt = mc.pt / GeV

    if pt > 0:
        return ggFTool.getWeight(pt)
    else:
        return 1.
