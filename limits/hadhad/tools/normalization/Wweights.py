from ..Configuration import *
def Wweight(OS, sample, nJets, tree, channel):
    """
    Return the W background estimation normalization weight
    """

    if sample.name.find('JimmyW') == -1: return 1.0

    VBFCategory = None
    if channel == 'mu': VBFCategory = muVBFCategory
    if channel == 'e': VBFCategory = eVBFCategory

    if OS:
        if nJets == 0: return 0.544206624064
        if nJets == 1: return 0.589419229752
        if nJets >= 2 and VBFCategory(tree):
            if tree.tau_numTrack == 1: return 0.611231803365
            if tree.tau_numTrack == 3: return 1.04092163384
        if nJets >= 2 and not VBFCategory(tree): return 0.66623248012
    else:
        if nJets == 0: return 0.688403294784
        if nJets == 1: return 0.720450849261
        if nJets >= 2 and VBFCategory(tree):
            if tree.tau_numTrack == 1: return 0.654933519175
            if tree.tau_numTrack == 3: return 0.920139540838
        if nJets >= 2 and not VBFCategory(tree): return 0.716120267197
    return 1.0