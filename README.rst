.. -*- mode: rst -*-

News
----
The master branch is now dedicated to analysis of the new atlas data format (xAOD).
For production of skims out of the old format (D3PD), please refer to the d3pd branch.

Dependencies
------------

* Latest 5.X `ROOT <http://root.cern.ch/drupal/>`_ with PyROOT enabled.

* `rootpy <https://github.com/rootpy/rootpy>`_::

   git clone git://github.com/rootpy/rootpy.git
   cd rootpy
   python setup.py install --user

* `goodruns <http://pypi.python.org/pypi/goodruns>`_::

   pip install --user goodruns

* `PyYAML <https://pypi.python.org/pypi/PyYAML>`_::

   pip install --user pyyaml

* `ConfigObj <http://www.voidspace.org.uk/python/configobj.html>`_::

   pip install --user configobj

* `externaltools <https://github.com/htautau/externaltools>`_::

   git clone git://github.com/htautau/externaltools.git

* `lumi <https://github.com/htautau/lumi>`_::

   git clone git://github.com/htautau/lumi.git

* `hhntup <https://github.com/htautau/hhntup>`_::

   git clone git://github.com/htautau/hhntup.git

* `TauSpinnerTool
  <https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauSpinnerTool>`_::

   svn co svn+ssh://${USER}@svn.cern.ch/reps/atlasoff/PhysicsAnalysis/TauID/TauSpinnerTool/trunk TauSpinnerTool


Place externaltools, lumi, and TauSpinnerTool in the same directory containing
hhntup to satisfy the symlinks in hhana. See the README in externaltools and
TauSpinnerTool for further instructions. Use at least Python version 2.6 (2.7
is preferred).

xAOD Migration
--------------
* Analysis release used::
  
   Base, 2.0.14

* Dataset used::

   mc14_8TeV.147808.PowhegPythia8_AU2CT10_Ztautau.merge.AOD.e2372_s1933_s1911_r5591_r5625

* To run the test::
  
   skim --local-test mc12_hadhad_xaod

Build and setup
---------------

Now build the C extension module for jet cleaning in the higgstautau package::

   make lib

Before running tests locally::

   source setup.sh


Skimming
--------

The skimming is performed by the ``hhskim.py`` script.

Run the skims on the grid (after setting up the panda client and your VOMS
proxy with the phys-higgs production role)::

    ./skim --yall mc11_hadhad mc12_hadhad \
                  data11_hadhad data12_hadhad \
                  embed11_hadhad embed12_hadhad


Running a local test of the skimming
------------------------------------

The samples are organized in several blocks defined in ``skims.cfg``.
Each block is written following the template::

   [block_name]
    student = hhskim.py
    dataset = dataset_block (defined in datasets.cfg) 
    version = version_number
    testinput = /path/to/the/input/files/for/test
    dest = SFU-LCG2_LOCALGROUPDISK,

For each block, modify the variable **testinput** according to your own setup.

Run the test::

    ./skim --yall block_name --local-test

The output will be created in the main directory as::

    hhskim_dataset_block.root

Creating ntuples
----------------

After the skims are finished and downloaded, update the paths in
``higgstautau/datasets_config.yml`` and update the datasets database::

    ./dsdb --reset hh

Then launch the batch jobs that create all the analysis ntuples (nominal and
systematics) with::

    ./run-all
