import copy


# common systematics
SYSTEMATICS_COMMON = {
    'Jets': {
        'JES': ('UP', 'DOWN'),
        'JER': ('UP',),
    },
    'Taus': {
        'TES': ('UP', 'DOWN'),
    }
} # ADD MORE HERE

# hadhad-only systematics
SYSTEMATICS_HADHAD = {
} # ADD MORE HERE

# electron-hadron-only systematics
SYSTEMATICS_EHAD = {
} # ADD MORE HERE

# muon-hadron-only systematics
SYSTEMATICS_MUHAD = {
} # ADD MORE HERE

SYSTEMATICS_HADHAD.update(SYSTEMATICS_COMMON)
SYSTEMATICS_EHAD.update(SYSTEMATICS_COMMON)
SYSTEMATICS_MUHAD.update(SYSTEMATICS_COMMON)

SYSTEMATICS = {
    'HADHAD': SYSTEMATICS_HADHAD,
    'EHAD': SYSTEMATICS_EHAD,
    'MUHAD': SYSTEMATICS_MUHAD,
}
