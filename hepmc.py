from atlastools import pdg

def get_tau_final_states(event):

    final_states = []
    for mc in event.mc:
        if mc.pdgId in (pdg.tau_plus, pdg.tau_minus):
            final_states.append(mc.final_state())
    return final_states

def get_VBF_partons(event):
    
    return []
