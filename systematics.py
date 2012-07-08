# common systematics
SYSTEMATICS_COMMON = [
    ('Jets', 'JES_UP'),
    ('Jets', 'JES_DOWN'),
    ('Jets', 'JER_UP'),
    # ('Jets', 'JER_DOWN'), NOT USED!
    ('Taus', 'TES_UP'),
    ('Taus', 'TES_DOWN'),
] # ADD MORE HERE

# hadhad-only systematics
SYSTEMATICS_HADHAD = [
] # ADD MORE HERE

# electron-hadron-only systematics
SYSTEMATICS_EHAD = [
] # ADD MORE HERE

# muon-hadron-only systematics
SYSTEMATICS_MUHAD = [
] # ADD MORE HERE

SYSTEMATICS = {
        'HADHAD': SYSTEMATICS_HADHAD + SYSTEMATICS_COMMON,
        'EHAD': SYSTEMATICS_EHAD + SYSTEMATICS_COMMON,
        'MUHAD': SYSTEMATICS_MUHAD + SYSTEMATICS_COMMON,
}
