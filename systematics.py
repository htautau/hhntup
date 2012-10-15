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
] # ADD MORE HERE

SYSTEMATICS = {
    'HADHAD': SYSTEMATICS_COMMON + SYSTEMATICS_HADHAD,
    'LEPHAD': SYSTEMATICS_COMMON + SYSTEMATICS_LEPHAD,
}

def iter_systematics(channel, include_nominal=False):

    if include_nominal:
        yield 'NOMINAL'
    for sys_variations in SYSTEMATICS[channel.upper()]:
        yield set(sys_variations)
