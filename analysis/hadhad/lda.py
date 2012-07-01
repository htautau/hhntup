#!/usr/bin/env python

"""
============================
Linear Discriminant Analysis
============================

Linear Discriminant Analysis (LDA) tries to identify attributes that account for
the most variance between classes. In particular, LDA, in constrast to PCA, is a
supervised method, using known class labels.
"""
print __doc__

import numpy as np
import pylab as pl
from matplotlib.ticker import NullFormatter

from sklearn.lda import LDA

import samples
import features

signals, backgrounds = samples.get_samples('2jet', purpose='train')

X_train, X_test,\
w_train, w_test,\
y_train, y_test = samples.make_classification(
        *(samples.make_train_test(signals, backgrounds,
            branches=features.hh_2jet_vars,
            train_fraction=.5,
            same_size_train=True,
            same_size_test=True)),
        standardize=True)

print X_train
print X_test

print w_train
print w_test

# standardize
X_train = standardize(X_train)
X_test = standardize(X_test)

lda = LDA(n_components=2)
lda.fit(X_train, y_train)
X_test_lda = lda.transform(X_test)

# plot LDA output
pl.figure(figsize=(8,10))

# definitions for the axes
left, width = 0.1, 0.8
bottom, height = 0.1, 0.6
bottom_h = bottom+height+0.02

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]

axScatter = pl.axes(rect_scatter)
pl.xlabel('Principal Component [arb. units]')
pl.ylabel('Secondary Component [arb. units]')

axHistx = pl.axes(rect_histx)

# no labels
nullfmt = NullFormatter()
axHistx.xaxis.set_major_formatter(nullfmt)

# the scatter plots:
for c, i, target_name in zip("rb", target_values, target_names):
    axScatter.scatter(X_test_lda[y_test == i, 0], X_test_lda[y_test == i, 1],
               c=c, label=target_name,
               s=w_test*10,
               alpha=0.9)

axScatter.legend()

axHistx.hist(X_test_lda[y_test == 1, 0], bins=30, range=axScatter.get_xlim(),
        color='r', label='Signal')
axHistx.hist(X_test_lda[y_test == 0, 0], bins=30, range=axScatter.get_xlim(),
        color='b', label='Background')

axHistx.set_xlim( axScatter.get_xlim() )

#axHistx.legend()
pl.title('Linear Discriminant Analysis')
pl.savefig('lda.png')
