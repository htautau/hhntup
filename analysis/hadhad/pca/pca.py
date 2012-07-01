#!/usr/bin/env python

"""
============================
Principal Component Analysis
============================

Principal Component Analysis (PCA) applied to this data identifies the
combination of attributes (principal components, or directions in the feature
space) that account for the most variance in the data.
"""
print __doc__

import numpy as np
import pylab as pl
from matplotlib.ticker import NullFormatter

from sklearn.decomposition import PCA
from sklearn import svm

import samples
import features


def perform_pca(channel):

    signals, backgrounds = samples.get_samples(channel, purpose='train')

    if channel == '01jet':
        branches = features.hh_01jet_vars
    else:
        branches = features.hh_2jet_vars

    X_train, X_test,\
    w_train, w_test,\
    y_train, y_test = samples.make_classification(
            *(samples.make_train_test(signals, backgrounds,
                branches=branches,
                train_fraction=.5,
                max_sig_train=2000,
                max_bkg_train=2000,
                max_sig_test=2000,
                max_bkg_test=2000,
                same_size_train=True,
                same_size_test=True,
                norm_sig_to_bkg_train=True,
                norm_sig_to_bkg_test=True)),
            standardize=True)

    print X_train
    print X_test

    print w_train
    print w_test

    print w_train.min(), w_train.max()

    pca = PCA(n_components=2)
    # fit only on background
    pca.fit(X_train[y_train == 0])
    X_train_pca = pca.transform(X_train)
    X_test_pca = pca.transform(X_test)

    xmin = X_test_pca[:, 0].min()
    xmax = X_test_pca[:, 0].max()
    ymin = X_test_pca[:, 1].min()
    ymax = X_test_pca[:, 1].max()

    width = xmax - xmin
    height = ymax - ymin

    xmin -= width*.1
    xmax += width*.1
    ymin -= height*.1
    ymax += height*.1

    # fit support vector machine on output of PCA
    clf = svm.SVC(C=100, gamma=.01, probability=True, scale_C=True)
    clf.fit(X_train_pca, y_train, sample_weight=w_train)

    # plot the decision function
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 500), np.linspace(ymin, ymax, 500))

    Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    channel_name = samples.CHANNEL_NAMES[channel]
    target_names = ['%s Signal' % channel_name,
                    '%s Background' % channel_name]
    target_values = [1, 0]

    # Percentage of variance explained for each components
    print 'explained variance ratio (first two components):', \
        pca.explained_variance_ratio_

    # plot PCA and SVM output
    pl.figure()
    # plot support vector machine decision function
    pl.set_cmap(pl.cm.jet)
    pl.contourf(xx, yy, Z, alpha=0.75)

    for c, i, target_name in zip("rb", target_values, target_names):
        pl.scatter(X_test_pca[y_test == i, 0], X_test_pca[y_test == i, 1],
                   c=c, label=target_name,
                   s=w_test[y_test == i]*10,
                   alpha=0.9)

    pl.xlim((xmin, xmax))
    pl.ylim((ymin, ymax))
    pl.legend()
    pl.xlabel('Principal Component [arb. units]')
    pl.ylabel('Secondary Component [arb. units]')
    pl.title('Principal Component Analysis\n'
             'and Support Vector Machine Decision Function')
    pl.savefig('pca_%s.png' % channel)


    # testing:
    signals, backgrounds = samples.get_samples(channel, purpose='test',
            mass=125)

    signal_train, signal_weight_train, \
    signal_test, signal_weight_test, \
    background_train, background_weight_train, \
    background_test, background_weight_test = samples.make_train_test(
                signals, backgrounds,
                branches=branches,
                train_fraction=.5,
                norm_sig_to_bkg_train=False,
                norm_sig_to_bkg_test=False)

    sample_test = np.concatenate((background_test, signal_test))
    sample_test = samples.std(sample_test)
    background_test, signal_test = sample_test[:len(background_test)], \
                                   sample_test[len(background_test):]

    signal_test = pca.transform(signal_test)
    background_test = pca.transform(background_test)

    pl.figure()
    pl.hist(clf.predict_proba(background_test)[:,-1],
            weights=background_weight_test, bins=30, range=(0, 1),
            label='Background', color='b')
    pl.hist(clf.predict_proba(signal_test)[:,-1],
            weights=signal_weight_test*10, bins=30, range=(0, 1),
            label='Signal x 10', color='r')
    pl.legend()
    pl.ylabel('Events')
    pl.xlabel('Support Vector Machine Signal Probability')
    pl.savefig('pca_svm_score_%s.png' % channel)


if __name__ == '__main__':

    perform_pca('2jet')
    #perform_pca('01jet')
