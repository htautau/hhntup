.. -*- mode: rst -*-

Dependencies
------------

You need to install/clone/checkout these packages:

* `rootpy <https://github.com/rootpy/rootpy>`_
* `goodruns <http://pypi.python.org/pypi/goodruns>`_
* `externaltools <https://github.com/htautau/externaltools>`_
* `lumi <https://github.com/htautau/lumi>`_
* `TauSpinnerTool
  <https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauSpinnerTool>`_

Place externaltools, lumi, and TauSpinnerTool in the same directory containing
hhntup. Install rootpy and goodruns into your ``PYTHONUSERBASE`` directory,
usually ``~/.local`` with ``python setup.py install --user`` from inside each
package.

Use at least Python version 2.6 (2.7 is preferred).


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
