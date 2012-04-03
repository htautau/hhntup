from atlastools import pdg

def get_tau_initial_final_states(event):
    """
    Get all taus and their decay products
    """
    states = []
    for mc in event.mc:
        if mc.pdgId in (pdg.tau_plus, pdg.tau_minus):
            init_state = mc
            final_state = mc.final_state()
            # some decays are not fully stored in the D3PDs
            # ignore them...
            if len(final_state) > 1:
                states.append((init_state, final_state))
    return states

def get_VBF_partons(event):
    
    return []
