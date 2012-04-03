[__many__]
    container = string(default='')
    local = string(default='')
    grl = string(default='')
    tree = string(default=tau)
    type = option(DATA, MC, default=MC)
    class = option(SIGNAL, BACKGROUND, default=BACKGROUND)
    label = option(TAU, ELEC, MUON, JET, default=TAU)
    weight = float(default=1.)
