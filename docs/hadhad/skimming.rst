Skimming
========

The latest versions of the skims are:

Data
----

* group.phys-higgs.hhskim.data11_7TeV*physics_JetTauEtmiss*v9/
* group.phys-higgs.hhskim.data12_8TeV*physics_JetTauEtmiss*v9/

MC
--

* group.phys-higgs.hhskim.mc11_7TeV*p851.v9/
* group.phys-higgs.hhskim.mc12_8TeV*p1130.v9/

MC with no trigger requirement
------------------------------

(same selection as embedded samples)

* group.phys-higgs.hhskim.mc12_8TeV*p1130.no_trigger.v9/

Embedded
--------

* group.phys-higgs.hhskim.data12_8TeV.period*.DESD_ZMUMU.pro13.embedding-02-48*v9/
* group.phys-higgs.hhskim.period*.DESD_SGLMU.pro10.embedding-02-41*v9/

These skims may still be running on the grid. Check
`Noel's <http://panda.cern.ch/server/pandamon/query?ui=user&name=Edmund%20Dawe%20ptu-382>`_ or 
`Dugan's <http://panda.cern.ch/server/pandamon/query?ui=user&name=Dugan%20ONeil%20xba-044>`_
panda page for jobs still running.

Skims are DaTRI'd to TOKYO-LCG2_PHYS-HIGGS and NIKHEF-ELPROD_PHYS-HIGGS
automatically as soon as the jobs finish.

Notes on skim v9
----------------

* The 2012 v9 skims were produced with the wrong MMC tag. The next version will use tag
  7. MMC 9 is for 2011 but tag 7 should be used for 2012 until we have a new tag
  for 2012. Strange tagging scheme going on here...

* At least one file is corrupt (missing tauMeta/TrigConfTree)::

  group.phys-higgs.hhskim.mc11_7TeV.107671.AlpgenJimmyZtautauNp1_pt20.e835_s1299_s1300_r3043_r2993_p851.v9.121105031615/
  group.phys-higgs.156433_052542.107671._00037.hhskim.mc11_p851_hadhad.root

* A small portion of data11 must have been recorded as "complete" by panda when
  it was in fact not done. This can be picked up in the next skim.

* The total lumi in skim v9 of data11 is 4604.52 pb-1 and of data12 is 14130.8
  pb-1

* The ggF reweighting wasn't applied so is 1. for all events in mc11 ggF
  samples. You must apply this yourself. ggF weights will be set properly in the
  next skim.

Additional Branches Created by the Skim
---------------------------------------

These additional branches are also created for convenience:

* ``vector<bool> tau_selected``: true if this tau was and false otherwise.
* ``float pileup_weight``: The output of the pileup reweighting tool
* ``vector<int> tau_trigger_match_index``: index of the matching EF tau, -1 if
  not matched.
* ``vector<int> tau_trigger_match_thresh``: the threshold of the
  matched trigger in GeV, 0 if not matched
* ``int number_of_good_vertices``
* ``vector<int> tau_numTrack_recounted`` 
* ``float ggf_weight``
* ``vector<float> tau_efficiency_scale_factor``
* ``vector<float> tau_efficiency_scale_factor_high``
* ``vector<float> tau_efficiency_scale_factor_low``
* ``vector<float> tau_fakerate_scale_factor``
* ``vector<float> tau_fakerate_scale_factor_high``
* ``vector<float> tau_fakerate_scale_factor_low``
* ``vector<float> tau_trigger_scale_factor``
* ``vector<float> tau_trigger_scale_factor_high``
* ``vector<float> tau_trigger_scale_factor_low``
* ``float tau_MMC_mass``
* ``TLorentzVector tau_MMC_resonance``
* ``float tau_MMC_resonance_pt``
* ``float tau_MMC_MET``
* ``float tau_MMC_MET_x``
* ``float tau_MMC_MET_y``
* ``float tau_MMC_MET_phi``
* ``float tau_collinear_mass``
* ``vector<float> tau_collinear_momentum_fraction``

Skims of MC also create the templates needed for pileup reweighting. These
are saved in ``*.TPileupReweighting.prw.root`` in each dataset. To perform
pileup reweighting you can either use the default templates that are provided by
the pileup reweighting package (used to create the ``pileup_weight`` branch in
the skims) or you can ``hadd`` together all of these
``*.TPileupReweighting.prw.root`` into one file and use that::

   hadd TPileupReweighting.prw.root */*TPileupReweighting.prw.root

Issues / Ideas
--------------

* Separate skims for systematics?

2011 Package Versions
---------------------

.. literalinclude:: ../../../externaltools/bundles/2011.lst

2012 Package Versions
---------------------

.. literalinclude:: ../../../externaltools/bundles/2012.lst

Skimmer Source
--------------

.. literalinclude:: ../../hhskim.py
   :emphasize-lines: 127-259
