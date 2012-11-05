Skimming
========

The latest versions of the skims are:

Data
----

* group.phys-higgs.hhskim.data11_7TeV*physics_JetTauEtmiss*v9/
* group.phys-higgs.hhskim.data12_8TeV*physics_JetTauEtmiss*v9/

MC
--

* group.phys-higgs.hhskim.mc11_7TeV*v9/
* group.phys-higgs.hhskim.mc12_8TeV*v9/

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


Skimmer Source
--------------

.. literalinclude:: ../../hhskim.py
   :emphasize-lines: 127-259
