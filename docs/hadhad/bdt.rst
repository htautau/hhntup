BDT Training
============

A boosted decision tree (BDT) is trained in each category to distinguish signal
from background using a selection of discriminating features.

To avoid bias on the training sample, the events used for training are not
included in the final results. The fraction of the datasets used for training is
adjustable and is currently 50%. Only data in the same-sign control region is
used for training in the QCD background model so all data events in the
opposite-sign signal region are included in the final results. 50% of the MC
signal and background samples in the opposite-sign signal region are discarded
after training. To further reduce bias on the training sample, the BDT
hyperparameters are tuned according to a proceedure outlined below.

Model Selection
---------------

The optimal BDT training parameters (hyperparameters) are selected by a
5-fold cross-validation at each point in a grid search on the number of trees
and the minimum number of events required at a leaf node. The hyperparameter
values that yield a BDT with the lowest classification error on the testing
sample define the final BDT used in the analysis in each category.


Output Distributions
--------------------

The output distributions of the BDT in each category is only shown in a low MMC
mass control region: :math:`80 < MMC < 110` GeV.


.. image:: plots/classify/var_vbf_event_bdt_score_control.png
	:width: 450px

.. image:: plots/classify/var_boosted_event_bdt_score_control.png
	:width: 450px

.. image:: plots/classify/var_ggf_event_bdt_score_control.png
	:width: 450px

