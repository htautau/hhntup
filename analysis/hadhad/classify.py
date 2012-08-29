#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--use-cache', action='store_true', default=False)
parser.add_argument('--unblind', action='store_true', default=False)
parser.add_argument('--train-frac', type=float, default=.5)
parser.add_argument('--cor', action='store_true', default=False)
args = parser.parse_args()


from pprint import pprint

from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

from features import *
import samples
from samples import *
from categories import CATEGORIES
import bkg_scales_cache
from systematics import iter_systematics

import numpy as np

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import IndexLocator, FuncFormatter

from rootpy.plotting import Hist
from rootpy.io import open as ropen

import pickle
import os

from tabulartext import PrettyTable
from utils import *


QUICK = False

# grid search params
if QUICK:
    # quick search for testing
    MIN_SAMPLES_LEAF = range(100, 120, 10)
    N_ESTIMATORS = range(10, 15, 2)
else:
    # full search
    MIN_SAMPLES_LEAF = range(10, 100, 10) + range(100, 2050, 50)
    N_ESTIMATORS = range(20, 1001, 30)

LIMITS_DIR = os.getenv('HIGGSTAUTAU_LIMITS_DIR')
if not LIMITS_DIR:
    sys.exit('You did not source setup.sh!')
LIMITS_DIR = os.path.join(LIMITS_DIR, 'hadhad', 'data')


def plot_grid_scores(grid_scores, best_point, params, name,
                     label_all_bins=False,
                     label_all_ticks=False,
                     format='png'):

    param_names = sorted(grid_scores[0][0].keys())
    param_values = dict([(pname, []) for pname in param_names])
    for pvalues, score, cv_scores in grid_scores:
        for pname in param_names:
            param_values[pname].append(pvalues[pname])

    for pname in param_names:
        param_values[pname] = np.unique(param_values[pname]).tolist()

    scores = np.empty(shape=[len(param_values[pname]) for pname in param_names])

    for pvalues, score, cv_scores in grid_scores:
        index = []
        for pname in param_names:
            index.append(param_values[pname].index(pvalues[pname]))
        scores.itemset(tuple(index), score)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmap = cm.get_cmap('jet', 100) # jet doesn't have white color
    #cmap.set_bad('w') # default value is 'k'
    ax.imshow(scores, interpolation="nearest", cmap=cmap)

    if label_all_ticks:
        plt.xticks(range(len(param_values[param_names[1]])),
                param_values[param_names[1]])
        plt.yticks(range(len(param_values[param_names[0]])),
                param_values[param_names[0]])
    else:
        trees = param_values[param_names[1]]
        def tree_formatter(x, pos):
            return str(trees[int(x)])

        leaves = param_values[param_names[0]]
        def leaf_formatter(x, pos):
            return str(leaves[int(x)])

        ax.xaxis.set_major_formatter(FuncFormatter(tree_formatter))
        ax.yaxis.set_major_formatter(FuncFormatter(leaf_formatter))
        ax.xaxis.set_major_locator(IndexLocator(2, 0))
        ax.yaxis.set_major_locator(IndexLocator(2, 0))
        xlabels = ax.get_xticklabels()
        for label in xlabels:
            label.set_rotation(45)

    ax.set_xlabel(params[param_names[1]], fontsize=20, position=(1., 0.), ha='right')
    ax.set_ylabel(params[param_names[0]], fontsize=20, position=(0., 1.), va='top')

    ax.set_frame_on(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    if label_all_bins:
        for row in range(scores.shape[0]):
            for col in range(scores.shape[1]):
                decor={}
                if ((param_values[param_names[0]].index(best_point[param_names[0]])
                     == row) and
                    (param_values[param_names[1]].index(best_point[param_names[1]])
                     == col)):
                    decor = dict(weight='bold',
                                 bbox=dict(boxstyle="round,pad=0.5",
                                           ec='black',
                                           fill=False))
                plt.text(col, row, "%.3f" % (scores[row][col]), ha='center',
                         va='center', **decor)
    else:
        # circle the best bin and label the parameters
        pass

    plt.savefig("grid_scores_%s.%s" % (name, format), bbox_inches='tight')
    plt.clf()


def apply_clf(clf,
              sample,
              category,
              region,
              branches,
              cuts=None,
              train_fraction=args.train_frac,
              systematic='NOMINAL'):

    if isinstance(sample, QCD):
        scores, weight = sample.scores(
                clf,
                category=category,
                region=region,
                branches=branches,
                train_fraction=train_fraction,
                cuts=cuts,
                systematic=systematic)
    elif isinstance(sample, Data):
        data_sample = data.ndarray(
                category=category,
                region=region,
                branches=branches,
                include_weight=False,
                cuts=cuts)
        scores = clf.predict_proba(data_sample)[:,-1]
        weight = np.ones(len(scores))
    else: # MC
        train, test = sample.train_test(
                category=category,
                region=region,
                branches=branches,
                train_fraction=train_fraction,
                cuts=cuts,
                systematic=systematic)
        weight = test['weight']
        input = np.vstack(test[branch] for branch in branches).T
        scores = clf.predict_proba(input)[:,-1]
    return scores, weight


def plot_clf(
        clf,
        backgrounds,
        category,
        region,
        branches,
        signals=None,
        signal_scale=1.,
        data=None,
        cuts=None,
        train_fraction=args.train_frac,
        name=None,
        draw_histograms=True,
        draw_data=args.unblind,
        save_histograms=False,
        bins=10,
        systematic='NOMINAL'):

    max_score = 0.
    min_score = 1.

    bkg_scores = []
    for bkg in backgrounds:
        scores, weight = apply_clf(
            clf,
            bkg,
            category=category,
            region=region,
            branches=branches,
            cuts=cuts,
            train_fraction=train_fraction,
            systematic=systematic)
        _min = scores.min()
        _max = scores.max()
        if _min < min_score:
            min_score = _min
        if _max > max_score:
            max_score = _max
        bkg_scores.append((bkg.label, scores, weight))

    if signals is not None:
        sig_scores = []
        for sig in signals:
            scores, weight = apply_clf(
                clf,
                sig,
                category=category,
                region=region,
                branches=branches,
                cuts=cuts,
                train_fraction=train_fraction,
                systematic=systematic)
            _min = scores.min()
            _max = scores.max()
            if _min < min_score:
                min_score = _min
            if _max > max_score:
                max_score = _max
            sig_scores.append((sig.label, scores, weight))
    else:
        sig_scores = None

    if data is not None and draw_data:
        data_scores, _ = apply_clf(
            clf,
            data,
            category=category,
            region=region,
            branches=branches,
            cuts=cuts,
            train_fraction=train_fraction)
        _min = data_scores.min()
        _max = data_scores.max()
        if _min < min_score:
            min_score = _min
        if _max > max_score:
            max_score = _max
    else:
        data_scores = None

    padding = (max_score - min_score) / (2 * bins)
    min_score -= padding
    max_score += padding
    hist_template = Hist(bins, min_score, max_score)

    bkg_hists = []
    for label, scores, weight in bkg_scores:
        hist = hist_template.Clone(title=label)
        for score, w in zip(scores, weight):
            hist.Fill(score, w)
        bkg_hists.append(hist)

    if signals is not None:
        sig_hists = []
        for label, scores, weight in sig_scores:
            hist = hist_template.Clone(title=label)
            for score, w in zip(scores, weight):
                hist.Fill(score, w)
            sig_hists.append(hist)
    else:
        sig_hists = None

    if data is not None and draw_data:
        data_hist = hist_template.Clone(title=data.label)
        map(data_hist.Fill, data_scores)
        print "Data events: %d" % sum(data_hist)
        print "Model events: %f" % sum(sum(bkg_hists))
        for hist in bkg_hists:
            print hist.GetTitle(), sum(hist)
        print "Data / Model: %f" % (sum(data_hist) / sum(sum(bkg_hists)))
    else:
        data_hist = None

    if draw_histograms:
        output_name = 'event_bdt_score'
        if name is not None:
            output_name += '_' + name
        draw(data=data_hist,
             model=bkg_hists,
             signal=sig_hists,
             signal_scale=signal_scale,
             category=category,
             category_name=info['name'],
             name="BDT Score",
             output_name=output_name,
             range=(min_score, max_score),
             show_ratio=data_hist is not None)


if __name__ == '__main__':

    use_cache = args.use_cache
    train_fraction = args.train_frac
    bins = 20

    #ztautau = MC_Ztautau()
    ztautau = Embedded_Ztautau()

    mc_ewk = MC_EWK()
    mc_top = MC_Top()
    mc_diboson = MC_Diboson()

    vbf_all = MC_VBF(mass=None)
    ggf_all = MC_ggF(mass=None)
    wh_all = MC_WH(mass=None)
    zh_all = MC_ZH(mass=None)

    backgrounds = [
        mc_top,
        mc_diboson,
        mc_ewk,
        ztautau,
    ]

    data = Data()

    qcd = QCD(data=data,
              mc=backgrounds[:])

    backgrounds.insert(0, qcd)

    signals = [
        vbf_all,
        ggf_all,
        wh_all,
        zh_all
    ]

    for category, info in sorted(CATEGORIES.items(), key=lambda item: item[0]):
        if category == 'preselection':
            continue
        print category

        # QCD shape region SS or !OS
        qcd.shape_region = info['qcd_shape_region']
        target_region = info['target_region']
        cuts = Cut()

        qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error = \
        bkg_scales_cache.get_scales(category)

        qcd.scale = qcd_scale
        ztautau.scale = ztautau_scale

        branches = info['features']

        if args.cor:
            branches = branches + ['mass_mmc_tau1_tau2']

        # split into testing and training samples
        sample_train, sample_test,\
        sample_weight_train, sample_weight_test,\
        labels_train, labels_test = make_classification(
            signals, backgrounds,
            category=category,
            region=target_region,
            branches=branches,
            train_fraction=None,
            max_sig_train=None,
            max_bkg_train=None,
            max_sig_test=None,
            max_bkg_test=None,
            norm_sig_to_bkg_train=True,
            norm_sig_to_bkg_test=False,
            same_size_train=True,
            same_size_test=False,
            standardize=False)

        if args.cor:
            # draw a linear correlation matrix
            samples.correlations(
                signal=sample_test[labels_test==1],
                signal_weight=sample_weight_test[labels_test==1],
                background=sample_test[labels_test==0],
                background_weight=sample_weight_test[labels_test==0],
                branches=branches,
                category=category)
            continue

        # train a classifier
        if use_cache and os.path.isfile('clf_%s.pickle' % category):
            # use a previously trained classifier
            with open('clf_%s.pickle' % category, 'r') as f:
                clf = pickle.load(f)
            print clf

        else:
            n_fold = 5
            cv = StratifiedKFold(labels_train, n_fold)
            # train a new BDT
            clf = AdaBoostClassifier(DecisionTreeClassifier(),
                    beta=.5,
                    compute_importances=True)
            # see top of file for grid search param constants
            grid_params = {
                'base_estimator__min_samples_leaf': MIN_SAMPLES_LEAF,
                'n_estimators': N_ESTIMATORS
            }
            #clf = SVC(probability=True, scale_C=True)
            # first grid search min_samples_leaf for the maximum n_estimators
            grid_clf = GridSearchCV(clf, grid_params, score_func=precision_score,
                                    n_jobs=-1)
            grid_clf.fit(sample_train, labels_train, sample_weight=sample_weight_train,
                         cv=StratifiedKFold(labels_train, n_fold))
            clf = grid_clf.best_estimator_
            grid_scores = grid_clf.grid_scores_
            """
            for
                for limit in xrange(1, clf.n_estimators):
                    for train, test in cv:
            """
            print "Classification report for the best estimator: "
            print clf
            y_true, y_pred = labels_test, clf.predict(sample_test)
            print "Tuned for 'precision' with optimal value: %0.3f" % precision_score(y_true, y_pred)
            print classification_report(y_true, y_pred)
            # plot a grid of the scores
            plot_grid_scores(
                grid_scores,
                best_point={
                    'base_estimator__min_samples_leaf':
                    clf.base_estimator.min_samples_leaf,
                    'n_estimators':
                    clf.n_estimators},
                params={
                    'base_estimator__min_samples_leaf':
                    'min leaf',
                    'n_estimators':
                    'trees'},
                name=category)

            if 'base_estimator__min_samples_leaf' in grid_params:
                # scale up the min-leaf and retrain on the whole set
                min_samples_leaf = grid_clf.best_estimator_.base_estimator.min_samples_leaf
                n_estimators=grid_clf.best_estimator_.n_estimators
                clf = AdaBoostClassifier(
                        DecisionTreeClassifier(
                            min_samples_leaf=int(min_samples_leaf * n_fold / float(n_fold - 1))),
                        n_estimators=n_estimators,
                        beta=.5, compute_importances=True)
                clf.fit(sample_train, labels_train,
                        sample_weight=sample_weight_train)
                print
                print "After scaling up min_leaf"
                print clf
                y_true, y_pred = labels_test, clf.predict(sample_test)
                print "Classification report: "
                print "Tuned for 'precision' with optimal value: %0.3f" % precision_score(y_true, y_pred)
                print classification_report(y_true, y_pred)

            if hasattr(clf, 'feature_importances_'):
                importances = clf.feature_importances_
                indices = np.argsort(importances)[::-1]
                print "Feature ranking:"
                print r"\begin{tabular}{c|c|c}"
                table = PrettyTable(["Rank", "Variable", "Importance"])
                print r"\hline\hline"
                print r"Rank & Variable & Importance\\"
                for f, feature in enumerate(branches):
                    table.add_row([f+1, feature, '%.3f' % importances[indices[f]]])
                    print r"%d & %s & %.3f\\" % (f + 1, VARIABLES[feature]['title'], importances[indices[f]])
                    #print "%d. %s (%f)" % (f + 1, feature, importances[indices[f]])
                print r"\end{tabular}"
                print
                print table.get_string(hrules=1)

            with open('clf_%s.pickle' % category, 'w') as f:
                pickle.dump(clf, f)

        # compare data and the model in a low mass control region
        cuts = Cut('mass_mmc_tau1_tau2 < 110')
        plot_clf(clf,
                 backgrounds,
                 category,
                 target_region,
                 branches,
                 signals=None,
                 data=data,
                 cuts=cuts,
                 train_fraction=train_fraction,
                 draw_data=True,
                 name='control')

        # show the background model and 125 GeV signal above mass control region
        cuts = Cut('mass_mmc_tau1_tau2 > 110')
        signals = [
            MC_VBF(mass=125),
            MC_ggF(mass=125),
            MC_WH(mass=125),
            MC_ZH(mass=125),
        ]
        plot_clf(
            clf,
            backgrounds,
            category,
            target_region,
            branches,
            signals=signals,
            signal_scale=20,
            cuts=cuts,
            train_fraction=train_fraction,
            name='ROI')

        # Create histograms for the limit setting with HistFactory
        # Include all systematic variations
        min_score = 1.
        max_score = 0.

        # apply on data
        data_scores, _ = apply_clf(
            clf,
            data,
            category=category,
            region=target_region,
            branches=branches,
            cuts=cuts)
        _min = data_scores.min()
        _max = data_scores.max()
        if _min < min_score:
            min_score = _min
        if _max > max_score:
            max_score = _max

        bkg_scores = {}
        sig_scores = {}

        for sys_variations in iter_systematics(
                channel='hadhad',
                include_nominal=True):

            if sys_variations == 'NOMINAL':
                sys_term = sys_variations
            else:
                sys_term = '_'.join(sys_variations)

            # apply on all backgrounds
            bkg_scores[sys_term] = []
            for bkg in backgrounds:
                scores, weight = apply_clf(
                    clf,
                    bkg,
                    category=category,
                    region=target_region,
                    branches=branches,
                    cuts=cuts,
                    systematic=sys_term)
                bkg_scores[sys_term].append((bkg.name, scores, weight))
                _min = scores.min()
                _max = scores.max()
                if _min < min_score:
                    min_score = _min
                if _max > max_score:
                    max_score = _max

            # apply on all signal masses
            sig_scores[sys_term] = []
            for mass in MC_Higgs.MASS_POINTS:
                for mode in (MC_VBF, MC_ggF, MC_WH, MC_ZH):
                    signal = mode(mass=mass)
                    scores, weight = apply_clf(clf,
                        signal,
                        category=category,
                        region=target_region,
                        branches=branches,
                        cuts=cuts)
                    _min = scores.min()
                    _max = scores.max()
                    if _min < min_score:
                        min_score = _min
                    if _max > max_score:
                        max_score = _max
                    sig_scores[sys_term].append(
                        ('Signal_%d_%s' %
                            (mass, signal.mode), scores, weight))

        padding = (max_score - min_score) / (2 * bins)
        min_score -= padding
        max_score += padding
        hist_template = Hist(bins, min_score, max_score)

        with ropen(os.path.join(LIMITS_DIR, '%s.root' % category),
                   'recreate') as f:

            data_hist = hist_template.Clone(name=data.name)
            map(data_hist.Fill, data_scores)
            f.cd()
            data_hist.Write()
            total_data = sum(data_hist)

            for sys_variations in iter_systematics(
                    channel='hadhad',
                    include_nominal=True):

                if sys_variations == 'NOMINAL':
                    sys_term = sys_variations
                    suffix = ''
                else:
                    sys_term = '_'.join(sys_variations)
                    suffix = '_' + sys_term

                bkg_hists = []
                for bkg_name, scores, weight in bkg_scores[sys_term]:
                    hist = hist_template.Clone(name=bkg_name + suffix)
                    for score, w in zip(scores, weight):
                        hist.Fill(score, w)
                    bkg_hists.append(hist)

                sig_hists = []
                for sig_name, scores, weight in sig_scores[sys_term]:
                    hist = hist_template.Clone(name=sig_name + suffix)
                    for score, w in zip(scores, weight):
                        hist.Fill(score, w)
                    sig_hists.append(hist)

                total_model = sum(sum(bkg_hists))
                print "Systematic: %s  Data / Model: %.5f" % (
                    sys_term, total_data / total_model)
                f.cd()
                for bkg in bkg_hists:
                    bkg.Write()
                for sig in sig_hists:
                    sig.Write()
