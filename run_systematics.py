import cluster

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


def run(channel, *args, **kwargs):

    for sys_type, sys_term in SYSTEMATICS[channel.upper()]:

        print
        print '============== Running %s systematics ==============' % sys_term
        print

        suffix = '--suffix %s' % sys_term
        syst = '--syst-type Systematics.%s --syst-term Systematics.%s.%s' % (
                sys_type, sys_type, sys_term)

        cluster.run(*args,
                args=suffix.split(),
                student_args=syst.split(),
                **kwargs)
