The scale factor should be applied to "SF29(leading tau) x SF20(sub-leading tau)".
(we should not apply like "SF29(leading tau) x SF29(sub-leading tau)", because either of
tau also fire the tau20Ti_medium1 as well.

So far, I do not use these variable (sorry I don't know how to use it...)

My example is

indxTau[0] … index of leading tau after ID
indxTau[1] … index of sub-leading tau after ID

    int indexEF1  = -1; double mindREF1 = 100.0;
    int indexEF2  = -1; double mindREF2 = 100.0;
    for (int j=0; j<trig_EF_tau_n; j++) {
      double dR1 = deltaR((*tau_phi)[indxTau[0]],(*trig_EF_tau_phi)[j],(*tau_eta)[indxTau[0]],(*trig_EF_tau_eta)[j]);
      double dR2 = deltaR((*tau_phi)[indxTau[1]],(*trig_EF_tau_phi)[j],(*tau_eta)[indxTau[1]],(*trig_EF_tau_eta)[j]);
      if (dR1<mindREF1) { mindREF1 = dR1; indexEF1 = j; }
      if (dR2<mindREF2) { mindREF2 = dR2; indexEF2 = j; }
    }

    double weightTrigSF = -1.0;
    if (mindREF1<0.2 && mindREF2<0.2) {
      if (indexEF1==indexEF2) {
        cout << "Strange trigger matching index!! " << indexEF1 << " " << indexEF2 << endl;
        abort();
      }

      std::string period = "periodB";
      if (RunNumber<201557) { period = "periodA"; }

      std::string nprong1 = "3p";
      if ((*tau_numTrack)[indxTau[0]]==1) { nprong1 = "1p"; }

      std::string nprong2 = "3p";
      if ((*tau_numTrack)[indxTau[1]]==1) { nprong2 = "1p"; }

      std::string id1 = "BDTm";
      if ((*tau_JetBDTSigTight)[indxTau[0]]==1) { id1 = "BDTt"; }

      std::string id2 = "BDTm";
      if ((*tau_JetBDTSigTight)[indxTau[1]]==1) { id2 = "BDTt"; }

      if ((*trig_EF_tau_EF_tau29Ti_medium1)[indexEF1]==1 && (*trig_EF_tau_EF_tau20Ti_medium1)[indexEF1]) {
        double sf29Ti = ttc29Ti->getSF((*tau_pt)[indxTau[0]], 0, period, nprong1, id1, "EVl");
        double sf20Ti = ttc20Ti->getSF((*tau_pt)[indxTau[1]], 0, period, nprong2, id2, "EVl");
        weightTrigSF = sf29Ti*sf20Ti;
      }
      else if ((*trig_EF_tau_EF_tau29Ti_medium1)[indexEF2]==1 && (*trig_EF_tau_EF_tau20Ti_medium1)[indexEF2]) {
        double sf29Ti = ttc29Ti->getSF((*tau_pt)[indxTau[1]], 0, period, nprong2, id2, "EVl");
        double sf20Ti = ttc20Ti->getSF((*tau_pt)[indxTau[0]], 0, period, nprong1, id1, "EVl");
        cout << "Scale factor " << sf29Ti << " " << sf20Ti << " " << (*tau_pt)[indxTau[1]] << " " << (*tau_pt)[indxTau[0]] << endl;
        cout << "Very strange!!! " << indexEF1 << " "  << indexEF2 << " " 
             << (*trig_EF_tau_EF_tau20Ti_medium1)[indexEF1] << " " << (*trig_EF_tau_EF_tau29Ti_medium1)[indexEF1] << " " 
             << (*trig_EF_tau_EF_tau20Ti_medium1)[indexEF2] << " " << (*trig_EF_tau_EF_tau29Ti_medium1)[indexEF2] << endl;
        weightTrigSF = sf29Ti*sf20Ti;
      }
      else {
        cout << "Very strange!!! " << indexEF1 << " "  << indexEF2 << " " 
             << (*trig_EF_tau_EF_tau20Ti_medium1)[indexEF1] << " " << (*trig_EF_tau_EF_tau29Ti_medium1)[indexEF1] << " " 
             << (*trig_EF_tau_EF_tau20Ti_medium1)[indexEF2] << " " << (*trig_EF_tau_EF_tau29Ti_medium1)[indexEF2] << endl;
        abort();
      }
    }
    else {
      cout << "Strange trigger matching!! " << endl;
      abort();
    }

    if (weightTrigSF<0.0) {
      cout << "could not find a trigger SF." << endl;
      abort();
    }


However, this logic is sometime broken, in case that leading tau fire tau20Ti_medium1 and
sub-leading tau fire tau29Ti_medium1 but offline pT of the sub-leading tau is less than 35GeV,
because the trigger SF does not support below 35GeV for tau29Ti_medium1.

That fraction is less than 0.5%,  thus, in the such a case, I just ignore to set the scale factor.
(that is, SF=1)
