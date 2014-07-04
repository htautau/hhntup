.. -*- mode: rst -*-

Dependencies
------------

* Latest 5.X `ROOT <http://root.cern.ch/drupal/>`_ with PyROOT enabled.

* `rootpy <https://github.com/rootpy/rootpy>`_::

   git clone git://github.com/rootpy/rootpy.git
   cd rootpy
   python setup.py install --user

* `goodruns <http://pypi.python.org/pypi/goodruns>`_::

   pip install --user goodruns

* `PyYAML <https://pypi.python.org/pypi/PyYAML>`::

   pip install --user pyyaml

* `ConfigObj <http://www.voidspace.org.uk/python/configobj.html>`::

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


Creating ntuples
----------------

After the skims are finished and downloaded, update the paths in
``higgstautau/datasets_config.yml`` and update the datasets database::

    ./dsdb --reset hh

Then launch the batch jobs that create all the analysis ntuples (nominal and
systematics) with::

    ./run-all
