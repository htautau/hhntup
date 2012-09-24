from atlastools import utils
from math import sin


def count_tracks(tau, event):
    """
    2011 track recounting
    """
    nOuterKtTrack = 0
    threshold = 4.

    for trk in event.tracks:

        dR = utils.dR(tau.eta, tau.phi, trk.eta, trk.phi)

        if (dR > 0.2 and dR < 0.6 and trk.pt / 1000.0 > 0.5
            and abs(trk.d0_wrtPV) < 1.0
            and abs(trk.z0_wrtPV * sin(trk.theta)) < 1.5
            and (trk.nPixHits + trk.nPixHoles) > 1
            and (trk.nPixHits + trk.nPixHoles + trk.nSCTHits + trk.nSCTHoles) > 6):
            iCheckKtTrack = 0.

            for j in xrange(tau.track_atTJVA_n):
                dR1 = utils.dR(tau.track_atTJVA_eta[j], tau.track_atTJVA_phi[j],
                               trk.eta, trk.phi)
                ptdR1 = tau.track_atTJVA_pt[j] * dR1 / trk.pt

                if ptdR1 > iCheckKtTrack:
                    iCheckKtTrack = ptdR1

            if iCheckKtTrack < threshold:
                nOuterKtTrack += 1

    ntrack_core = tau.track_atTJVA_n
    ntrack_full = ntrack_core + nOuterKtTrack
    return ntrack_core, ntrack_full
