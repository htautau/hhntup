Skimming
========

Skimming is performed in two stages. Stage I performs a loose selection
and stage II performs event cleaning and further selection.

The latest versions of the skims are:

MC11:

* Stage I: **group.phys-higgs.HHSkim.mc11_7TeV*p851.v1/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.mc11_7TeV*p851.v1.v1/**

MC12:

(some datasets still missing)
* Stage I: **group.phys-higgs.HHSkim.mc12_8TeV*p1011.v1/**
* Stage II: On the way... 

Data11:

* Stage I: **user.NoelDawe.HTauSkim.data11_7TeV*p851.v4/**
* Stage II: **group.phys-higgs.HHSkim2.HTauSkim.data11_7TeV*p851.v4.v3/**

Data12:

* Stage I: **group.phys-higgs.HHSkim.data12_8TeV*p1011_p1015.v1/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.data12_8TeV*p1011_p1015.v1.v5/**

Embedding11:

* Stage I: **group.phys-higgs.HHSkim.period*.DESD_SGLMU.pro10.embedding-02-41.Ztautau_*_rereco_p851_EXT0.v1/**
* Stage II: **group.phys-higgs.HHSkim2.HHSkim.period*.DESD_SGLMU.pro10.embedding-02-41.Ztautau_*_rereco_p851_EXT0.v1.v1/**

Embedding12:

On the way...

These skims may still be running on the grid. Check
`my pandamon page <http://panda.cern.ch/server/pandamon/query?ui=user&name=Edmund%20Dawe%20ptu-382>`_ or 'Dugans page <http://panda.cern.ch/server/pandamon/query?ui=user&name=Dugan%20ONeil%20xba-044>'
for jobs still running.

The current skims will always be backed up on SFU-LCG2_LOCALGROUPDISK

Stage I
-------

Stage I greatly reduces the dataset sizes
by selecting events that pass the di-tau trigger (data and MC) and that have at
least two loose tau candidates with a pT greater than 18 GeV (data only).

A histogram named "cutflow" is also saved.
The first bin is the total number of events and for MC the second
bin is the weighted number of events (from mc_event_weight).

The stage I script is here:
`HHSkim.py <https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/higgspy_svn/trunk/HHSkim.py>`_
and grid jobs are launched with the ``submit-hhskim2-*.sh`` scripts.

Data
~~~~

1) Triggers::

    if 177986 <= event.RunNumber <= 187815: # Periods B-K
  		return event.EF_tau29_medium1_tau20_medium1
    elif 188902 <= event.RunNumber <= 191933: # Periods L-M
  		return event.EF_tau29T_medium1_tau20T_medium1

2) Two loose taus:

   * ``tau_author!=2 && tau_pT > 18 GeV && tau_numTrack > 0``
   * ``tau_JetBDTLoose==1 || tau_tauLlhLoose==1``

   Note: the pT>18GeV cut was chosen to be able to fluctuate the TES.

MC
~~

1) Triggers: The di-tau triggers are first emulated and then the same trigger
   requirements as above are imposed. Branches are created to hold the emulated
   trigger decisions:
	  
   * ``EF_tau29_medium1_tau20_medium1_EMULATED``
   * ``EF_tau29T_medium1_tau20T_medium1_EMULATED``
	
   During the emulation procedure offline taus are matched to the trigger
   objects. These matches are saved in additional branches:

   * ``tau_trigger_match_index``: index of the matching EF tau, -1 if not matched.
   * ``tau_trigger_match_thresh``: 20 or 29 depending on the threshold of the
     matched trigger, 0 if not matched

No other selection is performed on MC.

Stage II
--------

Stage II performs all event selection, cleaning and object selection. These
additional branches are also created for convenience:

* ``tau_selected``: true if this is a selected tau and false otherwise.
  This is after pT cuts of 35/25 GeV (after the trigger matching) 
* ``pileup_weight``: The output of the pileup reweighting tool

Since the trigger matching is not performed on the data during stage I, this is
performed in stage II and the same branches created during stage I for MC are
now created for the data in stage II:

* ``tau_trigger_match_index``: index of the matching EF tau, -1 if not matched.
* ``tau_trigger_match_thresh``: 20 or 29 depending on the threshold of the
  matched trigger, 0 if not matched

The stage II script is here:
`HHSkim2.py <https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/higgspy_svn/trunk/HHSkim2.py>`_
and grid jobs are launched with the ``submit-hhskim-*.sh`` scripts.
