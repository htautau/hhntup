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


Fitting the Tau BDT Distribution
--------------------------------

In the VBF category before (left) and after (right) the fit:

.. image:: plots/background_estimation/2d_bdt_fit_vbf.png
	:width: 450px

.. image:: plots/background_estimation/2d_bdt_fit_vbf_after.png
	:width: 450px


In the boosted category before (left) and after (right) the fit:

.. image:: plots/background_estimation/2d_bdt_fit_boosted.png
	:width: 450px

.. image:: plots/background_estimation/2d_bdt_fit_boosted_after.png
	:width: 450px

In the non-boosted category before (left) and after (right) the fit:

.. image:: plots/background_estimation/2d_bdt_fit_ggf.png
	:width: 450px

.. image:: plots/background_estimation/2d_bdt_fit_ggf_after.png
	:width: 450px

