.. -*- mode: rst -*-

Installation
============

You need to install these packages:

* `rootpy <https://github.com/rootpy/rootpy>`_
* `goodruns <http://pypi.python.org/pypi/goodruns/2.0>`_

Download each package with::

   git clone git://github.com/ndawe/${packagename}.git

or::

   svn checkout http://svn.github.com/ndawe/${packagename}.git ${packagename}

then follow the README in each package to install them
(use the ``--user`` option after ``setup.py``)::

   cd rootpy
   python setup.py install --user
   cd ..

   cd atlastools
   python setup.py install --user
   cd ..

   cd goodruns
   python setup.py install --user
   cd ..

and be sure to use at least Python version 2.6 (2.7 is preferred).
DO NOT PLACE THESE PACKAGES IN THIS DIRECTORY! Put them somewhere else,
such as ``~/python-modules`` or where you keep your python package sources.

Now build the C extension module for jet cleaning in the higgstautau package::

   cd higgstautau/jetcleaning
   ./setup.py build_ext --inplace

Before submitting grid jobs be sure to create a panda.[user].cfg and then::

   ln -s panda.[user].cfg panda.cfg

or use an existing panda.[user].cfg.

Before running tests locally remember to::

   source setup.sh


External Tools
==============

higgspy depends on various external tools to perform systematics etc. See
``EXTERNALTOOLS_README.txt`` for further instructions.

Also see the Htautau common ntuple production:
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauCommonNtupleProduction


Recommendations
===============

https://twiki.cern.ch/twiki/bin/viewauth/Atlas/DataPreparationCheckListForPhysicsAnalysis
