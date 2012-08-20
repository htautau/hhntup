import copy


# common systematics
SYSTEMATICS_COMMON = [
    ('JES_UP', 'TES_UP'),
    ('JES_DOWN', 'TES_DOWN'),
    ('JER_UP',),
] # ADD MORE HERE

# hadhad-only systematics
SYSTEMATICS_HADHAD = [
    ('TAUBDT_UP',),
    ('TAUBDT_DOWN',),
] # ADD MORE HERE

# electron-hadron-only systematics
SYSTEMATICS_EHAD = [
] # ADD MORE HERE

# muon-hadron-only systematics
SYSTEMATICS_MUHAD = [
] # ADD MORE HERE

SYSTEMATICS = {
    'HADHAD': SYSTEMATICS_COMMON + SYSTEMATICS_HADHAD,
    'EHAD': SYSTEMATICS_COMMON + SYSTEMATICS_EHAD,
    'MUHAD': SYSTEMATICS_COMMON + SYSTEMATICS_MUHAD,
}

def iter_systematics(channel, include_nominal=False):

    if include_nominal:
        yield 'NOMINAL'
    for sys_variations in SYSTEMATICS[channel.upper()]:
        yield set(sys_variations)
