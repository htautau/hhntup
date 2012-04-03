[__many__]
    container = string
    grl = string(default='')
    tree = string(default=tauPerf)
    type = option(DATA, MC, default=MC)
    class = option(SIGNAL, BACKGROUND, default=BACKGROUND)
    label = option(TAU, ELEC, MUON, JET, default=TAU)
    weight = float(default=1.)
