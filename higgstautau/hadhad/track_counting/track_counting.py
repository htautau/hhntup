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

    return tau.track_atTJVA_n + nOuterKtTrack


def count_tracks_new(tau, event):
    """
    New 2011/2012 track recounting
    """

    nOuterKtTrack1 = 0
    """
    for trk in event.tracks:
        sinth  = TMath::Sin((*trk_atTJVA_theta)[i][indexTau]);
        trkpt  = sinth/TMath::Abs((*trk_atTJVA_qoverp)[i][indexTau]);
        trketa = -TMath::Log(TMath::Tan(0.5*(*trk_atTJVA_theta)[i][indexTau]));
        trkphi = (*trk_atTJVA_phi)[i][indexTau];
        dR = deltaR((*tau_phi)[indexTau],trkphi,(*tau_eta)[indexTau],trketa);
        if dR > 0.6:
            continue
        if (trkpt / 1000.0 > 0.5
            and abs((tau.trk_atTJVA_d0)[i][indexTau])<1.0
            and abs((*trk_atTJVA_z0)[i][indexTau]*sinth)<1.5
            and (*trk_nPixHits)[i]>1 && (*trk_nBLHits)[i]>0 && (*trk_nPixHits)[i]+(*trk_nSCTHits)[i]>6)
        {
            iCheckKtTrack = 0
            for (int j=0; j<(*tau_track_atTJVA_n)[indexTau]; j++) {
                dR1 = deltaR((*tau_track_atTJVA_phi)[indexTau][j],trkphi,(*tau_track_atTJVA_eta)[indexTau][j],trketa);
                ptdR1 = (*tau_track_atTJVA_pt)[indexTau][j]*dR1/trkpt;
                if ptdR1 < 4.0 and dR > 0.2:
                    iCheckKtTrack += 1
                if ptdR1 < 4.0 and dR < 0.2 and trkpt / 1000.0 < 1.0:
                    iCheckKtTrack += 1
            if iCheckKtTrack > 0:
                nOuterKtTrack1 += 1
    """
    return tau.track_atTJVA_n + nOuterKtTrack1
