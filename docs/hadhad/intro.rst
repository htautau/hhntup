Introduction
============

Plans
-----

* Use boosted decision trees (BDTs) to separate signal from background
* Perform analysis separately in ggF, boosted and VBF categories
* Estimate QCD from same-sign (SS) data
* Determine normalization of QCD and Ztautau from simultaneous fit with
  BDT distribution in low MMC mass region (80 - 110 GeV)

Outstanding Problems
--------------------

* What is the excess in data at pi/2 in alpha (3D angle) between the two taus?
  Drop cos alpha if we cannot resolve this.
  Using dR for now instead of cos alpha.
* Do fitting of QCD and Ztautau with 2D track distribution in low MMC mass bin
  [80, 110]. Separately in each category. Try fitting without selection on
  tracks or charge and use nonOS data - nonOS MC as QCD.
* Use shape of nonOS data - nonOS MC for QCD in OS region. Should smooth out
  QCD distributions with more stats. Possible bias in QCD background model since
  this will be the only sample with 2-track, 4-track etc taus.
  This is now implemented for the VBF category.
* Develop method of estimating systematics on the BDT score *WORK IN PROGRESS*
* Treat JES and TES as correlated (separately in signal and background) *DONE*
* Try embedding samples!!! Do they fix the angular distributions in VBF?

Ideas
-----

* Use Higgs pT as an input variable to the BDT.
* Use Higgs pT to define the boosted category instead of #jet
* Need theory uncert based on Higgs pT
* Use sum pT.
* Require "Higgs pT" > 100 GeV or X in the boosted category.
* Try defining separate categories on final BDT output and calculating limits in
  each subcategory separately using the MMC output.
  Swagato suggested this. gammagamma analysis does this.
  OR try calculating limits using 2D BDT vs MMC distributions.
* Look at sideband above 180GeV to check shape and norm
* Use SS for VBF norm but !OS for QCD in plots after norm to smooth it out
  Check that BDTJetScore is still well modeled.
* Next skims: use baseline trigger and make separate skim with new triggers
  E- and H7-
* Add new runs to skim


Notes
-----

* When shifting TES and BDT score, medium and tight must be recalculated.
  MET term must also be recalculated...

