Event Selection
===============

* GRL
* Trigger
* At least one type-1 vertex with at least 4 associated tracks
* ``larError <= 1``
* LAr hole: veto event if for any jet with ``jet_pt > 20`` GeV,
  ``jet_pt > 20 * (1-BCH_CORR_JET)/(1-BCH_CORR_CELL)|`` in data
  (Only applied to data/MC runs 180614 to 184169)
  ``-0.2 < jet_eta < 1.6 and -0.988 < jet_phi < -0.392``
* Jet cleaning: veto if there is a "loose-minus" bad jet with ``jet_pt > 20`` GeV and ``abs(jet_eta) < 4.5``
* Electron veto: ``pt > 15`` GeV, not in cracks, mediumPP,
  ``OQ & 1446 == 0``, ``author == 1, 3``, ``abs(charge) == 1``
* Muon veto: loose, staco, ``pt > 10`` GeV, ``eta < 2.5``, with good track

Require that the event has at least 2 preselected taus after each selection

* ``tau_author == 1, 3``
* ``tau_numTrack > 0``
* ``tau_muonVeto == 0``
* ``tau_EleBDTLoose == 0``
* lead track not in crack: ``!(1.37 < abs(eta) < 1.52)``
* lead track not in LAr hole for runs 180614 - 184169 (data and MC):
  ``!((-0.1 < eta < 1.55) && (-0.9 < phi < -0.5))``
* ``tau_JetBDTSigLoose == 1``
* Select two taus matching di-tau trigger RoIs
* Lead ``tau_pt > 35`` GeV and sublead ``tau_pt > 25`` GeV


