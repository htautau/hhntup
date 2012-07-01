Background Estimation
=====================

QCD and Ztautau are normalized (fit to data) in each category separately
and they are normalized simultaneously (in the same fit). The fit is
performed in a low MMC mass bin :math:`80 < MMC < 110` GeV.

* Take::

   OS QCD = SQCD × (SS Data − SZtautau × SS Ztautau − SS Other Bkg)

* Two free parameters: `SQCD`, `SZtautau`
* Fit OS QCD + SZtautau × OS Ztautau + OS Other Bkg to OS data.
* Fit according to a weighted log likelihood maximized with MIGRAD.
* Tried fitting 1D and 2D track (tau1, tau2) distribution and 2D tau
  BDT distribution.
* Most success on 2D BDT distributions. Convergence in both
  categories and reasonable scale factors.

Fitting the Track Distribution
------------------------------

In the 0/1-jet category before (left) and after (right) the fit:

.. image:: images/bg_est/2d_track_fit_01jet.png
	:width: 450px

.. image:: images/bg_est/2d_track_fit_01jet_after.png
	:width: 450px


In the 2+ jet category before (left) and after (right) the fit:

.. image:: images/bg_est/2d_track_fit_2jet.png
	:width: 450px

.. image:: images/bg_est/2d_track_fit_2jet_after.png
	:width: 450px


Fitting the Tau BDT Distribution
--------------------------------

In the 0/1-jet category before (left) and after (right) the fit:

.. image:: images/bg_est/2d_bdt_fit_01jet.png
	:width: 450px

.. image:: images/bg_est/2d_bdt_fit_01jet_after.png
	:width: 450px


In the 2+ jet category before (left) and after (right) the fit:

.. image:: images/bg_est/2d_bdt_fit_2jet.png
	:width: 450px

.. image:: images/bg_est/2d_bdt_fit_2jet_after.png
	:width: 450px


