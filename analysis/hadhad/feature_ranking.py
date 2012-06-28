#!/usr/bin/env python

"""
=========================================
Feature importances with forests of trees
=========================================

This examples shows the use of forests of trees to evaluate the importance of
features. The red plots are the feature importances of each individual tree,
and the blue plot is the feature importance of the whole forest.
"""
print __doc__

import numpy as np

from sklearn.ensemble import ExtraTreesClassifier

import samples
import features

signals, backgrounds = samples.get_samples('2jet', purpose='train')


X_train, X_test,\
w_train, w_test,\
y_train, y_test = samples.make_classification(
        *(samples.make_train_test(signals, backgrounds,
            branches=features.hh_2jet_vars)))


# Build a forest and compute the feature importances
forest = ExtraTreesClassifier(n_estimators=250,
                              n_jobs=-1,
                              min_samples_leaf=100,
                              compute_importances=True,
                              random_state=0)

print X_train
print y_train
print w_train


forest.fit(X_train, y_train, sample_weight=w_train)
importances = forest.feature_importances_
indices = np.argsort(importances)[::-1]

names = features.hh_2jet_vars
names = [names[i] for i in indices]
latex_names = [features.VARIABLES[name]['title'] for name in names]

# Print the feature ranking
print "Feature ranking:"

n_features = len(features.hh_2jet_vars)

for f in xrange(n_features):
    print "%d. %s (%f)" % (f + 1, names[f], importances[indices[f]])

# Plot the feature importances of the trees and of the forest
import pylab as pl
pl.figure()
pl.title("Feature importances")

for tree in forest.estimators_:
    pl.plot(xrange(n_features), tree.feature_importances_[indices], "r")

pl.plot(xrange(n_features), importances[indices], "b")
pl.xticks(range(len(latex_names)), latex_names, rotation=-30,
          rotation_mode='anchor', ha='left', va='top')
pl.savefig('ranking_random_forest.png', bbox_inches='tight')
