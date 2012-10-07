#!/usr/bin/env python

import argparse
from categories import CATEGORIES, CONTROLS, DEFAULT_LOW_MASS, DEFAULT_HIGH_MASS
from variables import VARIABLES


class formatter_class(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


parser = argparse.ArgumentParser(formatter_class=formatter_class)
"""
General Options
"""
parser.add_argument('--actions', nargs='*', choices=['plot', 'train'],
        default=[],
        help='only perform these actions')
parser.add_argument('--no-systematics', action='store_false',
        dest='systematics',
        help="turn off systematics",
        default=True)
parser.add_argument('--categories', nargs='*', default=CATEGORIES.keys(),
        help='which categories to draw plot or train in')
parser.add_argument('--controls', nargs='*', default=CONTROLS.keys(),
        help='which controls to draw plots in')
parser.add_argument('--only-controls', action='store_true', default=False,
        help='only draw control plots. no category plots.')
parser.add_argument('--no-controls', action='store_true', default=False,
        help='do not plot controls')
parser.add_argument('--unblind', action='store_true', default=False,
        help='plot the data in the signal region of the classifier output')
parser.add_argument('--embedding', action='store_true', default=False,
        help='use embedding instead of ALPGEN')

"""
Mass Regions Options
"""
parser.add_argument('--low-mass-cut', type=int,
        default=DEFAULT_LOW_MASS,
        help='the low mass window cut. '
        'Norms of Z and QCD are fit below this and '
        'the signal region of the classifier output is above this')
parser.add_argument('--high-mass-cut', type=int,
        default=DEFAULT_HIGH_MASS,
        help='the high mass window cut. '
        'Norms of Z and QCD are fit above this and '
        'the signal region of the classifier output is below this')
parser.add_argument('--no-sideband-in-control',
        dest='high_sideband_in_control',
        action='store_false',
        default=True,
        help='Exclude the high mass sideband in the mass control and include '
        'it in the signal region')

"""
Fitting Options
"""
parser.add_argument('--refit',
        action='store_false', dest='use_fit_cache',
        help="do not use cached background scale factors "
        "and instead recalculate them",
        default=True)
parser.add_argument('--fit-param', choices=('bdt', 'track', 'track1d'),
        default='track',
        help='parameters used to determine normalization of QCD and Z')
parser.add_argument('--draw-fit', action='store_true', default=False,
        help='draw the QCD/Z norm fit results')

"""
Training Options
"""
parser.add_argument('--retrain',
        action='store_false', dest='use_clf_cache',
        help="do not use cached classifier "
        "and instead train a new one",
        default=True)
parser.add_argument('--nfold', type=int, default=5,
        help='the number of folds in the cross-validation')
parser.add_argument('--clf-bins', dest='bins', type=int, default=10,
        help='the number of bins to use in the limit histograms and plots of '
        'the final classifier output')
parser.add_argument('--train-fraction', type=float, default=.5,
        help='the fraction of events used for training and excluded from the '
        'final limit histograms')
parser.add_argument('--train-categories', nargs='*', default=[],
        help='only train in these categories')
parser.add_argument('--quick-train', action='store_true', default=False,
        help='perform a very small grid search for testing purposes')
parser.add_argument('--grid-search', action='store_true', default=False,
        help='perform a grid-searched cross validation')
parser.add_argument('--forest-feature-ranking',
        action='store_true', default=False,
        help='Use a random forest to perform a feature ranking.')
parser.add_argument('--cor', action='store_true', default=False,
        help='draw correlation plots')

"""
Plotting Options
"""
parser.add_argument('--plots', nargs='*',
        help='only draw these plots. see the keys in variables.py')
parser.add_argument('--plot-cut', default=None, nargs='?',
        help='extra cut to be applied on the plots, but excluded from the '
        'QCD/Z normaliation and training and classifier output')
parser.add_argument('--plot-expr', default=None, nargs='?',
        help='expression to plot, instead of predefined ones in variables.py')
parser.add_argument('--plot-name', default=None, nargs='?',
        help='name of expr')
parser.add_argument('--plot-min', type=float, default=0, nargs='?',
        help='minimum of expr')
parser.add_argument('--plot-max', type=float, default=1, nargs='?',
        help='maximum of expr')
parser.add_argument('--plot-bins', type=int, default=20, nargs='?',
        help='number of bins to plot expr in')
parser.add_argument('--root', action='store_true', default=False,
        help='draw plots with ROOT. default is matplotlib')
parser.add_argument('--suffix', default=None, nargs='?',
        help='suffix to add to any output files or plots')
parser.add_argument('--output-formats', default=['png'], nargs='+',
        choices=('png', 'eps', 'pdf'),
        help='output formats')

args = parser.parse_args()

# root imports
import ROOT
ROOT.gROOT.SetBatch(True)

# local imports
from utils import *
from samples import *
import samples
from background_estimation import qcd_ztautau_norm
from classify import *
from config import plots_dir
from utils import *
import variables
from systematics import SYSTEMATICS
from categories import MassRegions

# stdlib imports
import pickle
import os
from pprint import pprint

# numpy imports
import numpy as np

# matplotlib imports
from matplotlib import pyplot as plt

# rootpy imports
from rootpy.plotting import Hist
from rootpy.io import open as ropen
from rootpy.extern.tabulartext import PrettyTable


def staged_score(self, X, y, sample_weight, n_estimators=-1):
    """
    calculate maximum signal significance
    """
    bins = 50
    for p in self.staged_predict_proba(X, n_estimators=n_estimators):

        scores = p[:,-1]

        # weighted mean accuracy
        y_pred = scores >= .5
        acc = np.average((y_pred == y), weights=sample_weight)

        min_score, max_score = scores.min(), scores.max()
        b_hist = Hist(bins, min_score, max_score + 0.0001)
        s_hist = b_hist.Clone()

        scores_s, w_s = scores[y==1], sample_weight[y==1]
        scores_b, w_b = scores[y==0], sample_weight[y==0]

        # fill the histograms
        for s, w in zip(scores_s, w_s):
            s_hist.Fill(s, w)
        for s, w in zip(scores_b, w_b):
            b_hist.Fill(s, w)

        # reverse cumsum
        #bins = list(b_hist.xedges())[:-1]
        s_counts = np.array(s_hist)
        b_counts = np.array(b_hist)
        S = s_counts[::-1].cumsum()[::-1]
        B = b_counts[::-1].cumsum()[::-1]

        # S / sqrt(S + B)
        s_sig = np.divide(list(S), np.sqrt(list(S + B)))

        #max_bin = np.argmax(np.ma.masked_invalid(significance)) #+ 1
        #max_sig = significance[max_bin]
        #max_cut = bins[max_bin]

        s_sig_max = np.max(np.ma.masked_invalid(s_sig))
        yield s_sig_max * acc


LIMITS_DIR = os.getenv('HIGGSTAUTAU_LIMITS_DIR')
if not LIMITS_DIR:
    sys.exit('You did not source setup.sh!')
LIMITS_DIR = os.path.join(LIMITS_DIR, 'hadhad', 'data')

PLOTS_DIR = plots_dir(__file__)

if args.embedding:
    ztautau = Embedded_Ztautau(systematics=args.systematics)
else:
    ztautau = MC_Ztautau(systematics=args.systematics)
others = Others(systematics=args.systematics)
data = Data(markersize=1.2)

higgs_125 = Higgs(
        mass=125,
        systematics=args.systematics,
        scale=50,
        linecolor='red',
        linestyle='dashed')

if 'train' in args.actions:

    # all modes, all masses
    signals_train = [
        Higgs(systematics=args.systematics),
    ]

    # all modes, 125GeV mass
    signal_eval = Higgs(
            mass=125,
            systematics=args.systematics,
            linecolor='red',
            linestyle='dashed')

figures = {}

categories_controls = (
        sorted(CATEGORIES.items(), key=lambda item: item[0]) +
        sorted(CONTROLS.items(), key=lambda item: item[0]))

output_suffix = '_%sfit' % args.fit_param
if args.embedding:
    output_suffix += '_embedding'
if args.suffix:
    output_suffix += '_%s' % args.suffix

mass_regions = MassRegions(
        low=args.low_mass_cut,
        high=args.high_mass_cut,
        high_sideband_in_control=args.high_sideband_in_control)

control_region = mass_regions.control_region
signal_region = mass_regions.signal_region
train_region = mass_regions.train_region

for category, cat_info in categories_controls:

    if category not in args.categories and category not in args.controls:
        continue

    if args.only_controls and category not in args.controls:
        continue

    if args.no_controls and category in args.controls:
        continue

    print
    print "=" * 40
    if category in args.categories:
        print "%s category" % category
    else:
        print "%s control region" % category
    print "=" * 40

    print "Cuts: %s" % cat_info['cuts']

    # QCD shape region SS or !OS
    qcd_shape_region = cat_info['qcd_shape_region']
    target_region = cat_info['target_region']

    qcd = QCD(data=data, mc=[others, ztautau],
          shape_region=qcd_shape_region)

    figures[category] = {}

    qcd.scale = 1.
    ztautau.scale = 1.

    # determine normalization of QCD and Ztautau
    # in each category separately
    qcd_ztautau_norm(
        ztautau=ztautau,
        others=others,
        qcd=qcd,
        data=data,
        category=category,
        target_region=target_region,
        mass_regions=mass_regions,
        bins=cat_info['fitbins'],
        draw=args.draw_fit,
        use_cache=args.use_fit_cache,
        param=args.fit_param,
        systematics=SYSTEMATICS if args.systematics else None,
        root=args.root)

    if 'plot' in args.actions:
        cuts = Cut(args.plot_cut)

        if args.plot_expr is not None:
            VARS = {tuple(args.plot_expr.split(',')):
                    {'title': args.plot_name,
                     'range': (args.plot_min, args.plot_max),
                     'bins': args.plot_bins,
                     'cats': ['GGF'],
                     'filename': 'expr_' + args.plot_name.replace(' ', '_')}}
        else:
            VARS = VARIABLES

        for expr, var_info in VARS.items():

            if category in args.controls and 'GGF' not in var_info['cats']:
                continue
            elif category in args.categories and category.upper() not in var_info['cats']:
                continue
            elif args.plots and expr not in args.plots:
                continue

            print
            print "plotting %s in %s category" % (expr, category)
            print cat_info['cuts'] & cuts

            bins = var_info['bins']
            min, max = var_info['range']

            if 'scale' in var_info:
                expr = "%s * %f" % (expr, var_info['scale'])

            other_hist = others.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            qcd_hist = qcd.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            ztautau_hist = ztautau.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            bkg_hists = [qcd_hist, other_hist, ztautau_hist]

            data_hist = data.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            signal_hist = higgs_125.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            print "Data events: %d" % sum(data_hist)
            print "Model events: %f" % sum(sum(bkg_hists))
            for hist in bkg_hists:
                print hist.GetTitle(), sum(hist)
            print "Data / Model: %f" % (sum(data_hist) / sum(sum(bkg_hists)))

            output_name = var_info['filename'] + output_suffix
            if cuts:
                output_name += '_' + cuts.safe()

            fig = draw(
                    data=data_hist,
                    model=bkg_hists,
                    signal=signal_hist,
                    name=var_info['root'] if args.root else var_info['title'],
                    output_name=output_name,
                    category_name=cat_info['name'],
                    category=category,
                    units=var_info.get('units', None),
                    range=var_info['range'],
                    show_ratio=True,
                    show_qq=False,
                    dir=PLOTS_DIR,
                    systematics=SYSTEMATICS if args.systematics else None,
                    root=args.root,
                    output_formats=args.output_formats)
            figures[category][expr] = fig

    if category not in args.categories:
        # don't train a classifier in a control region
        continue

    if 'train' in args.actions:

        # scikit-learn imports
        from sklearn.cross_validation import StratifiedKFold
        from sklearn.grid_search import GridSearchCV
        from sklearn.ensemble.grid_search import BoostGridSearchCV
        from sklearn.metrics import classification_report
        from sklearn.metrics import precision_score
        from sklearn.ensemble import AdaBoostClassifier, ExtraTreesClassifier
        from sklearn.tree import DecisionTreeClassifier

        backgrounds = [
            qcd,
            others,
            ztautau,
        ]

        # define training and test samples
        branches = cat_info['features']
        clf_filename = 'clf_%s%s.pickle' % (category, output_suffix)

        # train a classifier
        if (args.use_clf_cache
                and (category not in args.train_categories)
                and os.path.isfile(clf_filename)):
            # use a previously trained classifier
            print "using a previously trained classifier..."
            with open(clf_filename, 'r') as f:
                clf = pickle.load(f)
            print clf
        else:
            print "training a new classifier..."
            print "using these features:"
            print
            for branch in branches:
                print branch
            print

            if args.cor:
                branches = branches + ['mass_mmc_tau1_tau2']

            # split into testing and training samples
            sample_train, sample_test,\
            sample_weight_train, sample_weight_test,\
            labels_train, labels_test = make_classification(
                signals_train, backgrounds,
                category=category,
                region=target_region,
                branches=branches,
                train_fraction=args.train_fraction,
                max_sig_train=None,
                max_bkg_train=None,
                max_sig_test=None,
                max_bkg_test=None,
                norm_sig_to_bkg_train=True,
                norm_sig_to_bkg_test=False,
                same_size_train=True,
                same_size_test=False,
                remove_negative_train_weights=True,
                standardize=False,
                cuts=train_region)

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

            print
            print "plotting input variables as they are given to the BDT"
            # draw plots of the input variables
            for i, branch in enumerate(branches):
                print "plotting %s ..." % branch
                branch_data = sample_train[:,i]
                if 'scale' in variables.VARIABLES[branch]:
                    branch_data *= variables.VARIABLES[branch]['scale']
                _min, _max = branch_data.min(), branch_data.max()
                plt.figure()
                plt.hist(branch_data[labels_train==0],
                        bins=20, range=(_min, _max),
                        weights=sample_weight_train[labels_train==0],
                        label='Background', histtype='stepfilled',
                        alpha=.5)
                plt.hist(branch_data[labels_train==1],
                        bins=20, range=(_min, _max),
                        weights=sample_weight_train[labels_train==1],
                        label='Signal', histtype='stepfilled', alpha=.5)
                label = variables.VARIABLES[branch]['title']
                if 'units' in variables.VARIABLES[branch]:
                    label += ' [%s]' % variables.VARIABLES[branch]['units']
                plt.xlabel(label)
                plt.legend()
                plt.savefig('train_var_%s_%s.png' % (category, branch))

            print "plotting sample weights ..."
            _min, _max = sample_weight_train.min(), sample_weight_train.max()
            plt.figure()
            plt.hist(sample_weight_train[labels_train==0],
                    bins=20, range=(_min, _max),
                    label='Background', histtype='stepfilled',
                    alpha=.5)
            plt.hist(sample_weight_train[labels_train==1],
                    bins=20, range=(_min, _max),
                    label='Signal', histtype='stepfilled', alpha=.5)
            plt.xlabel('sample weight')
            plt.legend()
            plt.savefig('train_sample_weight_%s.png' % category)

            # train a new BDT

            # grid search params
            min_leaf_high = int((sample_train.shape[0] / 4.) *
                    (args.nfold - 1.) / args.nfold)
            min_leaf_low = max(10, int(min_leaf_high / 100.))

            if args.forest_feature_ranking:
                # perform a feature ranking with random forests
                # Build a forest and compute the feature importances
                forest = ExtraTreesClassifier(
                        n_estimators=200,
                        n_jobs=-1,
                        min_samples_leaf=(min_leaf_high + min_leaf_low) / 2,
                        compute_importances=True,
                        random_state=0)

                forest.fit(sample_train, labels_train,
                        sample_weight=sample_weight_train)
                importances = forest.feature_importances_
                indices = np.argsort(importances)[::-1]

                latex_names = [variables.VARIABLES[branch]['title']
                        for branch in branches]

                # Print the feature ranking
                print "random forest feature ranking:"

                n_features = len(branches)

                for i, idx in enumerate(indices):
                    print "%d. %s (%f)" % (
                            i + 1, branches[idx], importances[idx])

                # plot the feature importances of the trees and of the forest
                plt.figure()
                plt.title("Feature importances")

                for tree in forest.estimators_:
                    plt.plot(xrange(n_features),
                            tree.feature_importances_[indices], "r")

                plt.plot(xrange(n_features), importances[indices], "b")
                plt.xticks(range(len(latex_names)),
                        [latex_names[idx] for idx in indices], rotation=-30,
                          rotation_mode='anchor', ha='left', va='top')
                plt.savefig('ranking_random_forest_%s.png' % category,
                        bbox_inches='tight')
                continue

            if args.grid_search:
                if args.quick_train:
                    # quick search for testing
                    min_leaf_step = max((min_leaf_high - min_leaf_low) / 20, 1)
                    MAX_N_ESTIMATORS = 500
                    MIN_N_ESTIMATORS = 10

                else:
                    # larger search
                    min_leaf_step = max((min_leaf_high - min_leaf_low) / 100, 1)
                    MAX_N_ESTIMATORS = 1000
                    MIN_N_ESTIMATORS = 10

                MIN_SAMPLES_LEAF = range(
                        min_leaf_low, min_leaf_high, min_leaf_step)

                grid_params = {
                    'base_estimator__min_samples_leaf': MIN_SAMPLES_LEAF,
                }

                #AdaBoostClassifier.staged_score = staged_score

                clf = AdaBoostClassifier(
                        DecisionTreeClassifier(),
                        compute_importances=True,
                        learn_rate=.5)

                grid_clf = BoostGridSearchCV(
                        clf, grid_params,
                        max_n_estimators=MAX_N_ESTIMATORS,
                        min_n_estimators=MIN_N_ESTIMATORS,
                        # can use default ClassifierMixin score
                        #score_func=precision_score,
                        cv = StratifiedKFold(labels_train, args.nfold),
                        n_jobs=-1)

                print
                print "performing a grid search over these parameter values:"
                for param, values in grid_params.items():
                    print param.split('__')[-1], values
                    print '--'
                print "Minimum number of classifiers: %d" % MIN_N_ESTIMATORS
                print "Maximum number of classifiers: %d" % MAX_N_ESTIMATORS
                print
                print "training new classifiers ..."
                grid_clf.fit(
                        sample_train, labels_train,
                        sample_weight=sample_weight_train)
                clf = grid_clf.best_estimator_
                grid_scores = grid_clf.grid_scores_

                print "Best score: %f" % grid_clf.best_score_
                print "Best Parameters:"
                print grid_clf.best_params_

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

                """
                if 'base_estimator__min_samples_leaf' in grid_params:
                    # scale up the min-leaf and retrain on the whole set
                    min_samples_leaf = grid_clf.best_estimator_.base_estimator.min_samples_leaf
                    n_estimators=grid_clf.best_estimator_.n_estimators
                    clf = AdaBoostClassifier(
                            DecisionTreeClassifier(
                                min_samples_leaf=int(min_samples_leaf * args.nfold / float(args.nfold - 1))),
                            n_estimators=n_estimators,
                            compute_importances=True)
                    clf.fit(sample_train, labels_train,
                            sample_weight=sample_weight_train)
                    print
                    print "After scaling up min_leaf"
                    print clf
                """
            else:
                print "training a new classifier ..."

                if category == 'vbf':
                    min_samples_leaf=200
                    n_estimators=50
                else:
                    min_samples_leaf=150
                    n_estimators=20

                clf = AdaBoostClassifier(
                    DecisionTreeClassifier(
                        min_samples_leaf=min_samples_leaf),
                    compute_importances=True,
                    learn_rate=.5,
                    n_estimators=n_estimators)

                clf.fit(sample_train, labels_train,
                        sample_weight=sample_weight_train)

            with open(clf_filename, 'w') as f:
                pickle.dump(clf, f)

        if hasattr(clf, 'feature_importances_'):
            importances = clf.feature_importances_
            indices = np.argsort(importances)[::-1]
            print "Feature ranking:"
            print r"\begin{tabular}{c|c|c}"
            table = PrettyTable(["Rank", "Variable", "Importance"])
            print r"\hline\hline"
            print r"Rank & Variable & Importance\\"
            for f, idx in enumerate(indices):
                table.add_row([f + 1,
                    branches[idx],
                    '%.3f' % importances[idx]])
                print r"%d & %s & %.3f\\" % (f + 1,
                    VARIABLES[branches[idx]]['title'],
                    importances[idx])
            print r"\end{tabular}"
            print
            print table.get_string(hrules=1)

        # show the background model and data in the control region
        print "plotting classifier output in control region..."
        print control_region
        # data scores
        data_scores, _ = data.scores(clf,
                branches,
                train_fraction=args.train_fraction,
                category=category,
                region=target_region,
                cuts=control_region)

        # determine min and max scores
        min_score = 1.
        max_score = 0.
        _min = data_scores.min()
        _max = data_scores.max()
        if _min < min_score:
            min_score = _min
        if _max > max_score:
            max_score = _max

        # background model scores
        bkg_scores = []
        for bkg in backgrounds:
            scores_dict = bkg.scores(clf,
                    branches,
                    train_fraction=args.train_fraction,
                    category=category,
                    region=target_region,
                    cuts=control_region)

            for sys_term, (scores, weights) in scores_dict.items():
                assert len(scores) == len(weights)
                if len(scores) == 0:
                    continue
                _min = np.min(scores)
                _max = np.max(scores)
                if _min < min_score:
                    min_score = _min
                if _max > max_score:
                    max_score = _max

            bkg_scores.append((bkg, scores_dict))

        print "minimum score: %f" % min_score
        print "maximum score: %f" % max_score

        # prevent bin threshold effects
        min_score -= 0.00001
        max_score += 0.00001

        # add a bin above max score and below min score for extra beauty
        score_width = max_score - min_score
        bin_width = score_width / args.bins
        min_score -= bin_width
        max_score += bin_width

        # compare data and the model in a low mass control region
        plot_clf(
            background_scores=bkg_scores,
            category=category,
            category_name=cat_info['name'],
            signal_scores=None,
            data_scores=(data, data_scores),
            draw_data=True,
            name='control' + output_suffix,
            bins=args.bins + 2,
            min_score=min_score,
            max_score=max_score,
            systematics=SYSTEMATICS if args.systematics else None)

        # show the background model and 125 GeV signal in the signal region
        print "Plotting classifier output in signal region..."
        print signal_region
        # data scores
        data_scores, _ = data.scores(clf,
                branches,
                train_fraction=args.train_fraction,
                category=category,
                region=target_region,
                cuts=signal_region)

        # determine min and max scores
        min_score = 1.
        max_score = 0.
        _min = data_scores.min()
        _max = data_scores.max()
        if _min < min_score:
            min_score = _min
        if _max > max_score:
            max_score = _max

        # background model scores
        bkg_scores = []
        for bkg in backgrounds:
            scores_dict = bkg.scores(clf,
                    branches,
                    train_fraction=args.train_fraction,
                    category=category,
                    region=target_region,
                    cuts=signal_region)

            for sys_term, (scores, weights) in scores_dict.items():
                assert len(scores) == len(weights)
                if len(scores) == 0:
                    continue
                _min = np.min(scores)
                _max = np.max(scores)
                if _min < min_score:
                    min_score = _min
                if _max > max_score:
                    max_score = _max

            bkg_scores.append((bkg, scores_dict))

        # signal scores for M=125
        signal_scores_eval = signal_eval.scores(clf,
                branches,
                train_fraction=args.train_fraction,
                category=category,
                region=target_region,
                cuts=signal_region)

        min_score_signal = 1.
        max_score_signal = 0.

        for sys_term, (scores, weights) in signal_scores_eval.items():
            assert len(scores) == len(weights)
            if len(scores) == 0:
                continue
            _min = np.min(scores)
            _max = np.max(scores)
            if _min < min_score:
                min_score = _min
            if _max > max_score:
                max_score = _max
            if _min < min_score_signal:
                min_score_signal = _min
            if _max > max_score_signal:
                max_score_signal = _max

        print "minimum score: %f" % min_score
        print "maximum score: %f" % max_score
        print "minimum signal score: %f" % min_score_signal
        print "maximum signal score: %f" % max_score_signal

        # prevent bin threshold effects
        min_score -= 0.00001
        max_score += 0.00001
        min_score_signal -= 0.00001
        max_score_signal += 0.00001

        # add a bin above max score and below min score for extra beauty
        score_width_signal = max_score_signal - min_score_signal
        bin_width_signal = score_width_signal / args.bins
        min_score_signal -= bin_width_signal
        max_score_signal += bin_width_signal

        plot_clf(
            background_scores=bkg_scores,
            category=category,
            category_name=cat_info['name'],
            signal_scores=(signal_eval, signal_scores_eval),
            signal_scale=50,
            name='ROI' + output_suffix,
            bins=args.bins + 2,
            min_score=min_score_signal,
            max_score=max_score_signal,
            systematics=SYSTEMATICS if args.systematics else None)

        print "creating histograms for limits"
        bkg_scores = dict([(bkg.name, (bkg, scores_dict))
                for (bkg, scores_dict) in bkg_scores])

        for mass in Higgs.MASS_POINTS:
            print "%d GeV mass hypothesis" % mass
            # create separate signal. background and data histograms for each
            # mass hypothesis since the binning is optimized for each mass
            # individually.
            # The binning is determined by first locating the BDT cut value at
            # which the signal significance is maximized (S / sqrt(B)).
            # Everything above that cut is put in one bin. Everything below that
            # cut is put into N variable width bins such that the background is
            # flat.
            min_score_signal = 1.
            max_score_signal = 0.
            sig_scores = {}
            # signal scores
            for mode in Higgs.MODES.keys():
                sig = Higgs(mode=mode, mass=mass,
                        systematics=args.systematics)
                scores_dict = sig.scores(clf,
                        branches,
                        train_fraction=args.train_fraction,
                        category=category,
                        region=target_region,
                        cuts=signal_region)

                for sys_term, (scores, weights) in scores_dict.items():
                    assert len(scores) == len(weights)
                    if len(scores) == 0:
                        continue
                    _min = np.min(scores)
                    _max = np.max(scores)
                    if _min < min_score:
                        min_score = _min
                    if _max > max_score:
                        max_score = _max
                    if _min < min_score_signal:
                        min_score_signal = _min
                    if _max > max_score_signal:
                        max_score_signal = _max

                sig_scores['Signal_%s' % mode] = (sig, scores_dict)

            print "minimum signal score: %f" % min_score_signal
            print "maximum signal score: %f" % max_score_signal
            # prevent bin threshold effects
            min_score_signal -= 0.00001
            max_score_signal += 0.00001

            # using 50 bins determine location that maximizes signal
            # significance
            bkg_hist = Hist(50, min_score_signal, max_score_signal)
            sig_hist = bkg_hist.Clone()
            # fill background
            for bkg_sample, scores_dict in bkg_scores.values():
                for score, w in zip(*scores_dict['NOMINAL']):
                    bkg_hist.Fill(score, w)
            # fill signal
            for sig_sample, scores_dict in sig_scores.values():
                for score, w in zip(*scores_dict['NOMINAL']):
                    sig_hist.Fill(score, w)
            # determine maximum significance
            sig, max_sig, max_cut = significance(sig_hist, bkg_hist)
            print "maximum signal significance of %f at %f" % (max_sig, max_cut)
            # define one bin above max_cut and 5 below max_cut
            trans_bins = list(np.linspace(min_score_signal, max_cut))
            trans_bins.append(max_score_signal)

            hist_template = Hist(trans_bins)

            root_filename = '%s%s.root' % (category, output_suffix)
            f = ropen(os.path.join(LIMITS_DIR, root_filename), 'recreate')

            if args.unblind:
                data_hist = hist_template.Clone(name=data.name + '_%s' % mass)
                map(data_hist.Fill, data_scores)
                f.cd()
                data_hist.Write()

            for d in (bkg_scores, sig_scores):
                for name, (samp, scores_dict) in d.items():
                    for sys_term, (scores, weights) in scores_dict.items():
                        if sys_term == 'NOMINAL':
                            suffix = ''
                        else:
                            suffix = '_' + '_'.join(sys_term)
                        hist = hist_template.Clone(
                                name=name + ('_%d' % mass) + suffix)
                        for score, w in zip(scores, weights):
                            hist.Fill(score, w)
                        f.cd()
                        hist.Write()
            f.close()

# save all variable plots in one large multipage pdf
if 'plot' in args.actions and set(args.categories) == set(CATEGORIES.keys()) and not args.plots:
    # only create multipage pdf of all plots if we created all plots
    from matplotlib.backends.backend_pdf import PdfPages
    import datetime
    now = datetime.datetime.today()
    # put all plots in a multipage pdf
    for category, exprs in figures.items():
        pdf = PdfPages(os.path.join(
            PLOTS_DIR, 'variables_%s%s.pdf' % (category, output_suffix)))
        for expr, fig in sorted(exprs.items(), key=lambda x: x[0]):
            pdf.savefig(fig)
        d = pdf.infodict()
        # set pdf metadata
        d['Title'] = 'Features'
        d['Author'] = 'Noel Dawe'
        d['Subject'] = 'Higgs tautau hh features'
        d['Keywords'] = 'higgs tau'
        d['CreationDate'] = now
        d['ModDate'] = now
        pdf.close()
