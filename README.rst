.. -*- mode: rst -*-

Dependencies
============

You need to install these packages:

* `rootpy <https://github.com/rootpy/rootpy>`_
* `goodruns <http://pypi.python.org/pypi/goodruns>`_

and be sure to use at least Python version 2.6 (2.7 is preferred).

Now build the C extension module for jet cleaning in the higgstautau package::

   make lib

Before running tests locally::

   source setup.sh
