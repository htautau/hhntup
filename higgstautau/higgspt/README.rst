Higgs pT and Njets reweighing of Powheg ggF
===========================================

Download the reweighing histogram here::

  https://dgillber.web.cern.ch/dgillber/Higgs_pT_reweigh/Higgs_pT_reweighing_April2014.root

Also, a few slides with motivation is available here::

  https://dgillber.web.cern.ch/dgillber/Higgs_pT_reweigh/ggF_jet_bin_and_Higgs_pT_April10.pdf

The procedure outlined below reweighs the default Powheg ggF
samples used at ATLAS to simultaneously match the shape of the Higgs pT
spectrum of HRes 2.1, with dynamic scale, and the truth jet multiplicity of
some of the best available predictions to date (0/1 jet bins at NNLO and ≥2-jet
bin at NLO). Two ingredients are needed as input: The Higgs pT, in GeV, from
the MC truth record of Powheg+Pythia8 ggF. Find it by requiring PDGID=25 and
status=62.  The number of anti kt R = 0.4 truth jets with pT>25 GeV, not
originating from the decay products of the Higgs boson.  Start from the
AntiKt4Truth collection. Reject any jet with pT<25 GeV.  Reject any jet withing
ΔR < 0.4 of any electron, tau, photon or parton (directly) produced in the
Higgs decay.  Count the number of remaining jets.  Cross check: Doing this
before any selection, you should see about 57%, 30% and 13% of the events
having 0, 1, and ≥2 jets, respectively.  Calculate the event weight like this::

  // In init
  TFile *f = TFile::Open("Higgs_pT_reweighing_April2014.root");
  m_ggF_weight = (TH1F*)f->Get("Reweigh_Powheg_To_HRes2Dynamic_2jetsWeight15");

  // in the event loop, after finding Ntrutjets and H_pT as described above
  double ggF_weight = m_ggF_weight->Interpolate(H_pT); // H_pT in GeV
  if (Ntruthjets>=2) ggF_weight *= 1.5;

After applying this weight (but before applying any selection), the Higgs pT
spectrum will match the "best available prediction to date" (HRes2.1, with loop
mass effects and dynamical scale), and the fraction of events in the Ntruthjets
bins should be 61%, 25%, 14% corresponding to cross sections of 11.8 pb, 4.8 pb
and 2.7 pb based on the Yellow Report ggF cross section of 19.27 pb.


Gluon-gluon fusion jet bin uncertainty tool
===========================================

Download the tool here::

   https://dgillber.web.cern.ch/dgillber/ggF_uncertainty_tool/ggF_cross_section_uncertainty.tar.gz

Unpack and test the tool by::

  tar xzf ggF_cross_section_uncertainty.tar.gz

  # make sure to have ROOT setup, then test it by:
  ./test_tool.sh

A short example on how to use the tool is availalbe on page 13 of this talk::

  https://indico.cern.ch/event/305830/contribution/1/material/slides/0.pdf

Note that the number of truth jets above 25 GeV shoudl be used as input to the
tool. This should be taken from AntiKt4Truth for ATLAS and jets reconstructed
from Higgs decay products shoudl be removed (overlap removal with ΔR = 0.4).
