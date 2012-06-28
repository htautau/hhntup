Introduction
============

Plans:

* Use boosted decision trees (BDTs) to separate signal from background
* Perform analysis separately in ggF, boosted and VBF categories
* Estimate QCD from same-sign (SS) data
* Determine normalization of QCD and Ztautau from simultaneous fit with
  BDT distribution in low MMC mass region (80 - 110 GeV)


Status:

* Pileup reweighting *DONE*
* Define event selection *DONE*
* Define categories *DONE*
* Trigger emulation *DONE*
* Fake rate scale factor correction (will have updated factors from Marcus) *DONE*
* Trigger efficiency scale factor correction *TODO*
* Background estimation: determine QCD and Ztautau normalizations *WORK IN PROGRESS*
* !! Apply New jet calibration (ApplyJetCalibration-00-01-06?)
* Final BDT training and show output in low mass region (and sideband?) *TODO*
* Systematics *TODO*
* Determine limits on Higgs cross sections *TODO*
* Document analysis here! *WORK IN PROGRESS*


Outstanding Problems:

* What is the excess in data at pi/2 in alpha (3D angle) between the two taus?
  Drop cos alpha if we cannot resolve this.
* Do fitting of QCD and Ztautau with 2D track distribution in low MMC mass bin
  [80, 110]. Separately in each category. Try fitting without selection on
  tracks or charge and use nonOS data - nonOS MC as QCD.
* Use shape of nonOS data - nonOS MC for QCD in OS region. Should smooth out
  QCD distributions with more stats. Possible bias in QCD background model since
  this will be the only sample with 2-track, 4-track etc taus.
* Develop method of estimating systematics on the BDT score

Ideas
=====

* Try dR<3.2 cut on the taus to easily remove some QCD
* Use dR>0.8 or 1.0 to fix low mass modeling
* Use Higgs pT (see Michel's analysis). Bottom of ranking.
* Use sum pT. Bottom of ranking
* Split 0/1-jet category into boosted (1 jet above 150GeV) and non-boosted to
  improve sensitivity. DONE. Large improvement.
* Reduce number of bins in VBF BDT fit.
* Try defining separate categories on final BDT output and calculating limits in
  each subcategory separately using the MMC output.
  Swagato suggested this. gammagamma analysis does this.
  OR try calculating limits using 2D BDT vs MMC distributions.
