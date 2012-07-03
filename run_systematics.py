import cluster

# common systematics
SYSTEMATICS_COMMON = [
    ('Jets', 'JES_UP'),
    ('Jets', 'JES_DOWN'),
]

# hadhad-only systematics
SYSTEMATICS_HADHAD = [
]

# electron-hadron-only systematics
SYSTEMATICS_EHAD = [
]

# muon-hadron-only systematics
SYSTEMATICS_MUHAD = [
]

SYSTEMATICS = {
        'HADHAD': SYSTEMATICS_HADHAD + SYSTEMATICS_COMMON,
        'EHAD': SYSTEMATICS_EHAD + SYSTEMATICS_COMMON,
        'MUHAD': SYSTEMATICS_MUHAD + SYSTEMATICS_COMMON,
}


def run(channel, *args, **kwargs):

    for sys_type, sys_term in SYSTEMATICS[channel.upper()]:

        suffix = '--suffix %s' % sys_term
        syst = '--syst-type Systematics.%s --syst-term Systematics.%s.%s' % (
                sys_type, sys_type, sys_term)

        cluster.run(*args,
                args=suffix.split(),
                student_args=syst.split(),
                **kwargs)
