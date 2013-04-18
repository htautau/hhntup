from atlastools import utils
from atlastools import datasets
from math import sin, tan, log


def count_tracks(tau, event, year, datatype):

    year = year % 1000
    if year == 11:
        return count_tracks_2011(tau, event)
    elif year == 12:
        if datatype in (datasets.DATA, datasets.EMBED):
            return count_tracks_2012_data(tau, event)
        return count_tracks_2012(tau, event)
    raise ValueError('No track recounting defined for year %d' % year)


def count_tracks_2011(tau, event):
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

    return tau.track_atTJVA_n + nOuterKtTrack


def count_tracks_2012(tau, event):
    """
    2012 track recounting
    """
    nOuterKtTrack1 = 0
    tau_index = tau.index

    for trk in event.tracks:
        sinth  = sin(trk.atTJVA_theta[tau_index])
        trkpt  = sinth / abs(trk.atTJVA_qoverp[tau_index])
        trketa = -log(tan(0.5 * trk.atTJVA_theta[tau_index]))
        trkphi = trk.atTJVA_phi[tau_index]
        dR = utils.dR(tau.eta, tau.phi, trketa, trkphi)
        if dR > 0.6:
            continue
        if (trkpt / 1000.0 > 0.5
            and abs(trk.atTJVA_d0[tau_index]) < 1.0
            and abs(trk.atTJVA_z0[tau_index] * sinth) < 1.5
            and trk.nBLHits > 0
            and trk.nPixHits + trk.nPixelDeadSensors > 1
            and trk.nPixHits + trk.nPixelDeadSensors + trk.nSCTHits + trk.nSCTDeadSensors > 6):
            iCheckKtTrack = 0
            for j in xrange(tau.track_atTJVA_n):
                dR1 = utils.dR(tau.track_atTJVA_eta[j], tau.track_atTJVA_phi[j],
                        trketa, trkphi)
                ptdR1 = tau.track_atTJVA_pt[j] * dR1 / trkpt
                if ptdR1 < 4.0 and dR > 0.2:
                    iCheckKtTrack += 1
                if ptdR1 < 4.0 and dR < 0.2 and trkpt / 1000.0 < 1.0:
                    iCheckKtTrack += 1
            if iCheckKtTrack>0:
                nOuterKtTrack1 += 1
    return nOuterKtTrack1 + tau.track_atTJVA_n


def count_tracks_2012_data(tau, event):

    return tau.out_track_n_extended
