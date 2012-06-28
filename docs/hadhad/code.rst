Code
====

.. include:: links.rst

The main code is hosted by `bitbucket.org <http://bitbucket.org>`_ in private
`Git <http://git-scm.com/>`_ repositories and mirrored by the official ATLAS
Subversion repositories:

* Skimming/event selection code: `higgspy <https://bitbucket.org/ndawe/higgspy>`_
  (`atlasphys/Physics/Higgs/HSG4/software/common/higgspy_svn <https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/higgspy_svn/trunk>`_)

* BDT analysis code: `higgs-tau-mva <https://bitbucket.org/ndawe/higgs-tau-mva>`_
  (`atlasphys/Physics/Higgs/HSG4/software/hadhad/bdtmva_svn <https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/hadhad/bdtmva_svn/trunk>`_)

To contribute or view the code on `bitbucket.org`_ you need an account there and
permission granted by the author. The copy in the ATLAS Subversion repositories
should only be used as a reference for viewing the code and should not be used
for development and grid jobs since it may not have been updated recently.

These analysis packages depend heavily on:

* `rootpy`_: *a more feature-rich and pythonic
  interface with the ROOT libraries on top of the existing PyROOT bindings.*
* `goodruns`_: *an implementation of an
  ATLAS "good run list" (GRL) reader/writer in Python*
* `atlastools`_: various ATLAS utilities
  (grid job submission etc.)
* `yellowhiggs`_: *Interface for the
  CERN Yellow Report*


External Tools
--------------

The above software also depends on various external tools written for the ATLAS
community. A package called "externaltools" has been created to aid in keeping
track of which version of each package should be used in the analysis. Checkout
externaltools here::

   svn co svn+ssh://${USER}@svn.cern.ch/reps/atlasphys/Physics/Higgs/HSG4/software/hadhad/externaltools/trunk externaltools

Before submitting grid jobs checkout the externaltools package, follow the
README to checkout all packages and build them. Then copy the newly created
externaltools directory (inside the externaltools package) into the higgspy
package. "externaltools" is in the .gitignore file of higgspy so these files
will not be tracked.
