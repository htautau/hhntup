import copy

# common systematics
SYSTEMATICS_COMMON = [
#     ('JES_UP',),
    ('JES_Statistical_UP',),
    ('JES_Modelling_UP',),
    ('JES_Detector_UP',),
    ('JES_Mixed_UP',),
    ('JES_EtaModelling_UP',),
    ('JES_EtaMethod_UP',),
    ('JES_PURho_UP',),
    ('JES_PUPt_UP',),
    ('JES_PUNPV_UP',),
    ('JES_PUMu_UP',),
    ('JES_FlavComp_UP',),
    ('JES_FlavResp_UP',),
    ('JES_BJet_UP',),
    ('JES_NonClosure_UP',),
    ('JVF_UP',),
    ('TES_UP',),
#     ('JES_DOWN',),
    ('JES_Statistical_DOWN',),
    ('JES_Modelling_DOWN',),
    ('JES_Detector_DOWN',),
    ('JES_Mixed_DOWN',),
    ('JES_EtaModelling_DOWN',),
    ('JES_EtaMethod_DOWN',),
    ('JES_PURho_DOWN',),
    ('JES_PUPt_DOWN',),
    ('JES_PUNPV_DOWN',),
    ('JES_PUMu_DOWN',),
    ('JES_FlavComp_DOWN',),
    ('JES_FlavResp_DOWN',),
    ('JES_BJet_DOWN',),
    ('JES_NonClosure_DOWN',),
    ('JVF_DOWN',),
    ('TES_DOWN',),
    ('JER_UP',),
] # ADD MORE HERE

# hadhad-only systematics
SYSTEMATICS_HADHAD = [
    ('TAUBDT_UP',),
    ('TAUBDT_DOWN',),

    ## Extra TES systematics
    ('TES_UP',),
    ('TES_DOWN',),
    ('TES_EOP_UP',),
    ('TES_EOP_DOWN',),
    ('TES_CTB_UP',),
    ('TES_CTB_DOWN',),
    ('TES_Bias_UP',),
    ('TES_Bias_DOWN',),
    ('TES_EM_UP',),
    ('TES_EM_DOWN',),
    ('TES_LCW_UP',),
    ('TES_LCW_DOWN',),
    ('TES_PU_UP',),
    ('TES_PU_DOWN',),
    ('TES_OTHERS_UP',),
    ('TES_OTHERS_DOWN',),
] # ADD MORE HERE

# lepton-hadron-only systematics
SYSTEMATICS_LEPHAD = [
    ## JES systematics
    #('ATLAS_JES_Detector_DOWN',),
    #('ATLAS_JES_Detector_UP',),
    
    ('ATLAS_JES_2012_Modelling1_DOWN',),
    ('ATLAS_JES_2012_Modelling1_UP',),
    
    ('ATLAS_JES_FlavComp_DOWN',),
    ('ATLAS_JES_FlavComp_UP',),
    
    ('ATLAS_JES_FlavResp_DOWN',),
    ('ATLAS_JES_FlavResp_UP',),

    ('ATLAS_JES_EtaModelling_DOWN',),
    ('ATLAS_JES_EtaModelling_UP',),

    #('ATLAS_JER_DOWN',),
    #('ATLAS_JER_UP',),

    ('ATLAS_JVF_DOWN',),
    ('ATLAS_JVF_UP',),

    ## TES systematics
    ('ATLAS_TES_2012_DOWN',),
    ('ATLAS_TES_2012_UP',),

    ## MET systematics
    #('ATLAS_MET_RESOSOFT_DOWN',),
    #('ATLAS_MET_RESOSOFT_UP',),
    #('ATLAS_MET_SCALESOFT_DOWN',),
    #('ATLAS_MET_SCALESOFT_UP',),

    ## Electron systematics
    #('ATLAS_EL_ES_DOWN',),
    #('ATLAS_EL_ES_UP',),
    #('ATLAS_EL_RES_DOWN',),
    #('ATLAS_EL_RES_UP',),

    ## Muon systematics
    #('ATLAS_MU_ES_DOWN',),
    #('ATLAS_MU_ES_UP',),

    # ('ATLAS_ANA_LH12_01jet_Wlnu_DETA_UP',),
    # ('ATLAS_ANA_LH12_01jet_Wlnu_DETA_DOWN',),
    # ('ATLAS_ANA_LH12_01jet_Wlnu_PTRAT_UP',),
    # ('ATLAS_ANA_LH12_01jet_Wlnu_PTRAT_DOWN',),

    # ('ATLAS_ANA_EMB_ISOL_UP',),
    # ('ATLAS_ANA_EMB_ISOL_DOWN',),

    # ('ATLAS_ANA_LH12_SR_FF_UP',),
    # ('ATLAS_ANA_LH12_SR_FF_DOWN',),

    # ('ATLAS_ANA_LH12_SR_FFRW_UP',),
    # ('ATLAS_ANA_LH12_SR_FFRW_DOWN',),
    
] # ADD MORE HERE

SYSTEMATICS = {
    'HADHAD': SYSTEMATICS_COMMON + SYSTEMATICS_HADHAD,
    'LEPHAD': SYSTEMATICS_LEPHAD,
}

def iter_systematics(channel, include_nominal=False):

    if include_nominal:
        yield 'NOMINAL'
    for sys_variations in SYSTEMATICS[channel.upper()]:
        yield set(sys_variations)
