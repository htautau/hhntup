import copy


# common systematics
SYSTEMATICS_COMMON = [
    ('JES_UP',),
    ('TES_UP',),
    ('JES_DOWN',),
    ('TES_DOWN',),
    ('JER_UP',),
] # ADD MORE HERE

# hadhad-only systematics
SYSTEMATICS_HADHAD = [
    ('TAUBDT_UP',),
    ('TAUBDT_DOWN',),
] # ADD MORE HERE

# lepton-hadron-only systematics
SYSTEMATICS_LEPHAD = [
    ## JES systematics
    ('ATLAS_JES_Detector_DOWN',),
    ('ATLAS_JES_Detector_UP',),
    ('ATLAS_JES_EtaModelling_DOWN',),
    ('ATLAS_JES_EtaModelling_UP',),
    ('ATLAS_JES_Modelling_DOWN',),
    ('ATLAS_JES_Modelling_UP',),

    ('ATLAS_JER_DOWN',),
    ('ATLAS_JER_UP',),

    ('ATLAS_JVF_DOWN',),
    ('ATLAS_JVF_UP',),

    ## TES systematics
    ('ATLAS_TAU_ES_DOWN',),
    ('ATLAS_TAU_ES_UP',),

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
