

def get_tau_trigger_obj_idx(config, event, trigger):
    """
    Return list of trigger tau indicies associated with
    RoIs for this trigger
    """
    chainID = config.GetChainId(trigger)

    trigger_idx = []
    for i, id in enumerate(event.trig_Nav_chain_ChainId):
        if id == chainID + 10000:
            for roi_idx in event.trig_Nav_chain_RoIIndex[i]:
                if len(event.trig_RoI_EF_tau_Analysis__TauJetContainer) == 0:
                    continue
                for k, status in enumerate(event.trig_RoI_EF_tau_Analysis__TauJetContainerStatus[roi_idx]):
                    if status != 1:
                        continue
                    idx = event.trig_RoI_EF_tau_Analysis__TauJetContainer[roi_idx][k]
                    trigger_idx.append(idx)
    # remove duplicates
    return list(set(trigger_idx))
