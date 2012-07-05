from ..Configuration import *
def ABCDweight(sample, nJets, tree, channel, control = False):
    """
    Returns the normalization for QCD background to control region D
    """

    if not control: return 1.0

    VBFCategory = None
    if channel == 'mu': VBFCategory = muVBFCategory
    if channel == 'e': VBFCategory = eVBFCategory

    if sample.name.find('data') == -1: return 1.0
    if nJets == 0: return 0.0345532453359
    if nJets == 1: return 0.0147802451474
    if nJets >= 2 and VBFCategory(tree): return 0.00579902076536
    if nJets >= 2 and not VBFCategory(tree): return 0.00878085894136
    return 1.0