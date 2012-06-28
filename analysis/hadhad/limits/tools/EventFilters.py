from CutFlow import CutFlow

def BDTresponse(tree, variables, BDT):
    """
    Returns BDT response for a specified BDT with a specified input variable list
    """

    for key in variables.variables.iterkeys():
        variables.SetValue(key, getattr(tree, key))

    return BDT.EvaluateMVA("BDT")


cutFlowVBF = CutFlow()
cutFlowVBF.add('Total')
#cutFlowVBF.add('Muon Isolation')
#cutFlowVBF.add('Opposite Signs')
cutFlowVBF.add('MET > 20 GeV')
cutFlowVBF.add('mT < 30 GeV')
cutFlowVBF.add('Njets >= 2')
cutFlowVBF.add('Jet eta product < 0')
cutFlowVBF.add('leptons between jets')
cutFlowVBF.add('Pseudorapidity gap > 3.0')
cutFlowVBF.add('j1j2 mass > 300 GeV')


def muLH_VBF_Cuts(tree, cutFlow = None):
    """
    Returns the VBF cuts decsion
    """

    if cutFlow: cutFlow.increment('Total')

    #Muon Isolation
    # if not tree.muon_isolated: return False
    # if cutFlow: cutFlow.increment('Muon Isolation')
        
    #Opposite charge requirement
    # if not (tree.muon_charge*tree.tau_charge == -1): return False
    # if cutFlow: cutFlow.increment('Opposite Signs')

    #MET cut
    if not tree.MET > 20000: return False
    if cutFlow: cutFlow.increment('MET > 20 GeV')

    # Transverse mass
    if not tree.mass_transverse_met_muon < 30000: return False
    if cutFlow: cutFlow.increment('mT < 30 GeV')

    #Jet multiplicity, VBF cut
    if not tree.numJets >= 2: return False
    if cutFlow: cutFlow.increment('Njets >= 2')

    #Eta product
    if not tree.eta_product_j1_j2 < 0: return False
    if cutFlow: cutFlow.increment('Jet eta product < 0')

    #Leptons between the jets
    if not (tree.tau_centrality_j1_j2 > 1/math.e and tree.muon_centrality_j1_j2 > 1/math.e): return False
    if cutFlow: cutFlow.increment('leptons between jets')

    #Delta Eta of the jets
    if not tree.eta_delta_j1_j2 > 3.0: return False
    if cutFlow: cutFlow.increment('Pseudorapidity gap > 3.0')

    #Jet mass
    if not tree.mass_j1_j2 > 300000: return False
    if cutFlow: cutFlow.increment('j1j2 mass > 300 GeV')

    return True


cutFlowggF = CutFlow()
cutFlowggF.add('Total')
# cutFlowggF.add('Muon Isolation')
# cutFlowggF.add('Opposite Signs')
cutFlowggF.add('MET > 20 GeV')
cutFlowggF.add('mT < 30 GeV')


def muLH_ggF_Cuts(tree, cutFlow = None):
    """
    Returns the ggF cuts decision
    """

    if cutFlow: cutFlow.increment('Total')

    #Muon Isolation
    # if not tree.muon_isolated: return False
    # if cutFlow: cutFlow.increment('Muon Isolation')
        
    #Opposite charge requirement
    # if not (tree.muon_charge*tree.tau_charge == -1): return False
    # if cutFlow: cutFlow.increment('Opposite Signs')

    #MET cut
    if not tree.MET > 20000: return False
    if cutFlow: cutFlow.increment('MET > 20 GeV')

    # Transverse mass
    if not tree.mass_transverse_met_muon < 30000: return False
    if cutFlow: cutFlow.increment('mT < 30 GeV')

    return True
