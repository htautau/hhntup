
is_visible = lambda fourvect: (
    fourvect.Et() > 10 * GeV and abs(fourvect.Eta()) < 2.5)


def closest_reco_object(objects, thing, dR=0.2):

    closest_object = None
    closest_dR = 1111
    for other in objects:
        dr = utils.dR(other.eta, other.phi, thing.Eta(), thing.Phi())
        if dr < dR and dr < closest_dR:
            closest_object = other
            closest_dR = dr
    return closest_object, closest_dR
