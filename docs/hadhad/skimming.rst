Skimming
========

Skimming is performed in two stages. Stage I performs a loose selection
and stage II performs event cleaning and further selection.

The latest versions of the skims are:

MC11:

* Stage I: **group.phys-higgs.HHSkim.mc11_7TeV*p851.v1/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.mc11_7TeV*p851.v1.v3/**

MC12:

* Stage I: **group.phys-higgs.HHSkim.mc12_8TeV*NTUP_TAU*p1130.v4/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.mc12_8TeV*NTUP_TAU*p1130.v4.v1/**

Data11:

* Stage I: **user.NoelDawe.HTauSkim.data11_7TeV*p851.v4/**
* Stage II: **group.phys-higgs.HHSkim2.HTauSkim.data11_7TeV*p851.v4.v6/**

Data12:

* Stage I: **group.phys-higgs.HHSkim.data12_8TeV*NTUP_TAU*p1130.v3/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*NTUP_TAU*p1130.v3.v3/**

Embedding11:

* Stage I: **group.phys-higgs.HHSkim.period*.DESD_SGLMU.pro10.embedding-02-41.Ztautau_*_rereco_p851_EXT0.v1/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.period*.DESD_SGLMU.pro10.embedding-02-41.Ztautau_*_rereco_p851_EXT0.v1.v2/**

Embedding12:

* Stage I: **group.phys-higgs.HHSkim.period*.DESD_SGLMU.pro13.embedding-02-42.Ztautau_hh_isol_mfsim_rereco_p1131_EXT0.v1/**
  Notes: Some files were on offline sites.
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.period*.DESD_SGLMU.pro13.embedding-02-42.Ztautau_hh_isol_mfsim_rereco_p1131_EXT0.v1.v1/**


These skims may still be running on the grid. Check
`Noel's <http://panda.cern.ch/server/pandamon/query?ui=user&name=Edmund%20Dawe%20ptu-382>`_ or 
`Dugan's <http://panda.cern.ch/server/pandamon/query?ui=user&name=Dugan%20ONeil%20xba-044>`_
panda page for jobs still running.

The current skims will always be backed up on SFU-LCG2_LOCALGROUPDISK

Stage I
-------

Stage I greatly reduces the dataset sizes
by selecting events that pass the di-tau trigger (data and MC) and that have at
least two loose tau candidates with a pT greater than 18 GeV (data only).

A histogram named "cutflow" is also saved.
The first bin is the total number of events and for MC the second
bin is the weighted number of events (the sum of mc_event_weight over all events
in the original sample).

Data
~~~~

1) 2011 triggers::

      if 177986 <= event.RunNumber <= 187815: # Periods B-K
         return event.EF_tau29_medium1_tau20_medium1
      elif 188902 <= event.RunNumber <= 191933: # Periods L-M
         return event.EF_tau29T_medium1_tau20T_medium1

   2012 triggers::

      event.EF_tau29Ti_medium1_tau20Ti_medium1 or event.EF_2tau38T_medium1

2) Two loose taus:

   * ``tau_author!=2 && tau_pT > 18 GeV && tau_numTrack > 0``
   * ``tau_JetBDTLoose==1 || tau_tauLlhLoose==1``

   Note: the pT>18GeV cut was chosen to be able to fluctuate the TES.

MC
~~

1) Triggers: The di-tau triggers are first emulated and then the same trigger
   requirements as above are imposed. Branches are created to hold the emulated
   trigger decisions:
	  
   * ``bool EF_tau29_medium1_tau20_medium1_EMULATED``
   * ``bool EF_tau29T_medium1_tau20T_medium1_EMULATED``
	
   During the emulation procedure offline taus are matched to the trigger
   objects. These matches are saved in additional branches:

   * ``vector<int> tau_trigger_match_index``: index of the matching EF tau, -1 if not matched.
   * ``vector<int> tau_trigger_match_thresh``: the threshold of the matched
     trigger in GeV, 0 if not matched

No other selection is performed on MC.

Stage I for MC also creates the templates needed for pileup reweighting. These
are saved in ``*.TPileupReweighting.prw.root`` in each dataset. To perform
pileup reweighting you can either use the default templates that are provided by
the pileup reweighting package or you can ``hadd`` together all of these
``*.TPileupReweighting.prw.root`` provided by stage I into one file and use
that::

   hadd TPileupReweighting.prw.root */*TPileupReweighting.prw.root

Stage II
--------

Stage II performs all event selection, cleaning and object selection. These
additional branches are also created for convenience:

* ``vector<bool> tau_selected``: true if this tau was and false otherwise.
* ``float pileup_weight``: The output of the pileup reweighting tool

Since the trigger matching is not performed on the data during stage I, this is
performed in stage II and the same branches created during stage I for MC are
now created for the data in stage II:

* ``vector<int> tau_trigger_match_index``: index of the matching EF tau, -1 if not matched.
* ``vector<int> tau_trigger_match_thresh``: the threshold of the
  matched trigger in GeV, 0 if not matched


Issues / Ideas
--------------

* Fix trigger threshold association in 2012 MC skims
* Run skim II directly on D3PDs to get absolute cut-flows
* Also keep track of weighted cut-flow (``mc_event_weight``, pileup weight)
* Remove 2tau38 trigger in 2012 skims
* Implement new trigger emulation in 2012 skims
* Set trigger scale factors, efficiency scale factors, fake rate scale factors
  in the skims. 
* Separate skims for systematics?
