#!/usr/bin/env python

import argparse
from categories import CATEGORIES, CONTROLS
from variables import VARIABLES


class formatter_class(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawTextHelpFormatter):
    pass


parser = argparse.ArgumentParser(formatter_class=formatter_class)
parser.add_argument('--no-fit-cache',
        action='store_false', dest='use_fit_cache',
        help="do not use cached background scale factors "
        "and instead recalculate them",
        default=True)
parser.add_argument('--no-clf-cache',
        action='store_false', dest='use_clf_cache',
        help="do not use cached classifier "
        "and instead train a new one",
        default=True)
parser.add_argument('--actions', nargs='*', choices=['plot', 'train'],
        default=[],
        help='only perform these actions')
parser.add_argument('--no-systematics', action='store_false',
        dest='systematics',
        help="turn off systematics",
        default=True)
parser.add_argument('--mass-cut', type=int, default=110,
        help='the mass window cut. norms of Z and QCD are fit below this and '
        'the signal region of the classifier output is above this')
parser.add_argument('--nfold', type=int, default=5,
        help='the number of folds in the cross-validation')
parser.add_argument('--clf-bins', dest='bins', type=int, default=10,
        help='the number of bins to use in the limit histograms and plots of '
        'the final classifier output')
parser.add_argument('--cor', action='store_true', default=False,
        help='draw correlation plots')
parser.add_argument('--unblind', action='store_true', default=False,
        help='plot the data in the signal region of the classifier output')
parser.add_argument('--embedding', action='store_true', default=False,
        help='use embedding instead of ALPGEN')
parser.add_argument('--train-fraction', type=float, default=.5,
        help='the fraction of events used for training and excluded from the '
        'final limit histograms')
parser.add_argument('--quick-train', action='store_true', default=False,
        help='perform a very small grid search for testing purposes')
parser.add_argument('--categories', nargs='*', default=CATEGORIES.keys(),
        help='which categories to draw plot or train in')
parser.add_argument('--controls', nargs='*', default=CONTROLS.keys(),
        help='which controls to draw plots in')
parser.add_argument('--only-controls', action='store_true', default=False,
        help='only draw control plots. no category plots.')
parser.add_argument('--forest-feature-ranking',
        action='store_true', default=False,
        help='Use a random forest to perform a feature ranking.')
parser.add_argument('--train-categories', nargs='*', default=[],
        help='only train in these categories')
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

for category, cat_info in categories_controls:

    if category not in args.categories and category not in args.controls:
        continue

    if args.only_controls and category not in args.controls:
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
    qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error = qcd_ztautau_norm(
        ztautau=ztautau,
        backgrounds=[others],
        data=data,
        category=category,
        target_region=target_region,
        qcd_shape_region=qcd_shape_region,
        mass_cut=args.mass_cut,
        bins=cat_info['fitbins'],
        use_cache=args.use_fit_cache)

    qcd.scale = qcd_scale
    qcd.scale_error = qcd_scale_error
    ztautau.scale = ztautau_scale
    ztautau.scale_error = ztautau_scale_error

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

            output_name = var_info['filename']
            if args.embedding:
                output_name += '_embedding'

            fig = draw(
                    data=data_hist,
                    model=bkg_hists,
                    signal=signal_hist,
                    name=var_info['title'],
                    output_name=output_name + '_' + cuts.safe() if cuts else output_name,
                    category_name=cat_info['name'],
                    category=category,
                    units=var_info.get('units', None),
                    range=var_info['range'],
                    show_ratio=True,
                    show_qq=False,
                    model_colour_map=None,
                    dir=PLOTS_DIR,
                    systematics=SYSTEMATICS if args.systematics else None,
                    root=args.root)
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

        if args.embedding:
            clf_filename = 'clf_%s_embedding.pickle' % category
        else:
            clf_filename = 'clf_%s.pickle' % category

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

            if args.forest_feature_ranking:
                # perform a feature ranking with random forests
                # Build a forest and compute the feature importances
                forest = ExtraTreesClassifier(
                        n_estimators=100,
                        n_jobs=-1,
                        max_depth=3,
                        min_samples_leaf=10,
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

                for f in xrange(n_features):
                    print "%d. %s (%f)" % (
                            f + 1, branches[f], importances[indices[f]])

                # plot the feature importances of the trees and of the forest
                plt.figure()
                plt.title("Feature importances")

                for tree in forest.estimators_:
                    plt.plot(xrange(n_features),
                            tree.feature_importances_[indices], "r")

                plt.plot(xrange(n_features), importances[indices], "b")
                plt.xticks(range(len(latex_names)), latex_names, rotation=-30,
                          rotation_mode='anchor', ha='left', va='top')
                plt.savefig('ranking_random_forest_%s.png' % category,
                        bbox_inches='tight')
                continue

            # train a new BDT
            clf = AdaBoostClassifier(
                    DecisionTreeClassifier(),
                    compute_importances=True,
                    learn_rate=.5)

            # grid search params
            if args.quick_train:
                # quick search for testing
                MIN_SAMPLES_LEAF = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
                N_ESTIMATORS = [
                        1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
                grid_params = {
                    'base_estimator__min_samples_leaf': MIN_SAMPLES_LEAF,
                    'n_estimators': N_ESTIMATORS
                }
                grid_clf = GridSearchCV(
                        clf, grid_params,
                        # can use default ClassifierMixin score
                        #score_func=precision_score,
                        cv = StratifiedKFold(labels_train, args.nfold),
                        n_jobs=-1)
            else:
                # full search
                min_leaf_high = int((sample_train.shape[0] / 2.) *
                        (args.nfold - 1.) / args.nfold)
                min_leaf_low = max(10, int(min_leaf_high / 30.))
                min_leaf_step = max((min_leaf_high - min_leaf_low) / 100, 1)
                MIN_SAMPLES_LEAF = range(
                        min_leaf_low, min_leaf_high, min_leaf_step)
                MAX_N_ESTIMATORS = 1000
                grid_params = {
                    'base_estimator__min_samples_leaf': MIN_SAMPLES_LEAF,
                }
                grid_clf = BoostGridSearchCV(
                        clf, grid_params,
                        max_n_estimators=MAX_N_ESTIMATORS,
                        # can use default ClassifierMixin score
                        #score_func=precision_score,
                        cv = StratifiedKFold(labels_train, args.nfold),
                        n_jobs=-1)

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

        # Create histograms for the limit setting with HistFactory
        # Include all systematic variations

        cuts = Cut('mass_mmc_tau1_tau2 < %d' % args.mass_cut)
        print "plotting classifier output in control region..."
        print cuts
        # data scores
        data_scores, _ = data.scores(clf,
                branches,
                train_fraction=args.train_fraction,
                category=category,
                region=target_region,
                cuts=cuts)

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
                    cuts=cuts)

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

        # compare data and the model in a low mass control region
        plot_clf(
            background_scores=bkg_scores,
            category=category,
            category_name=cat_info['name'],
            signal_scores=None,
            data_scores=(data, data_scores),
            draw_data=True,
            name='control_embedding' if args.embedding else 'control',
            bins=args.bins,
            min_score=min_score,
            max_score=max_score,
            systematics=SYSTEMATICS if args.systematics else None)

        # show the background model and 125 GeV signal in the signal region
        cuts = Cut('mass_mmc_tau1_tau2 > %d' % args.mass_cut)
        print "Plotting classifier output in signal region..."
        print cuts
        # data scores
        data_scores, _ = data.scores(clf,
                branches,
                train_fraction=args.train_fraction,
                category=category,
                region=target_region,
                cuts=cuts)

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
                    cuts=cuts)

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
                cuts=cuts)

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

        print "minimum score: %f" % min_score
        print "maximum score: %f" % max_score

        plot_clf(
            background_scores=bkg_scores,
            category=category,
            category_name=cat_info['name'],
            signal_scores=(signal_eval, signal_scores_eval),
            signal_scale=50,
            name='ROI_embedding' if args.embedding else 'ROI',
            bins=args.bins,
            min_score=min_score,
            max_score=max_score,
            systematics=SYSTEMATICS if args.systematics else None)

        print "creating histograms for limits"

        # signal scores for all masses and modes
        sig_scores = {}
        for mass in Higgs.MASS_POINTS:
            for mode in Higgs.MODES.keys():
                sig = Higgs(mode=mode, mass=mass,
                        systematics=args.systematics)
                scores_dict = sig.scores(clf,
                        branches,
                        train_fraction=args.train_fraction,
                        category=category,
                        region=target_region,
                        cuts=cuts)

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

                name = 'Signal_%d_%s' % (mass, mode)
                sig_scores[name] = (sig, scores_dict)

        print "minimum score: %f" % min_score
        print "maximum score: %f" % max_score

        #padding = (max_score - min_score) / (2 * args.bins)
        min_score -= 0.00001
        max_score += 0.00001
        hist_template = Hist(args.bins, min_score, max_score)

        if args.embedding:
            root_filename = '%s_embedding.root' % category
        else:
            root_filename = '%s.root' % category

        with ropen(os.path.join(LIMITS_DIR, root_filename), 'recreate') as f:

            data_hist = hist_template.Clone(name=data.name)
            map(data_hist.Fill, data_scores)
            f.cd()
            data_hist.Write()

            bkg_scores = dict([(bkg.name, (bkg, scores_dict))
                for (bkg, scores_dict) in bkg_scores])

            for d in (bkg_scores, sig_scores):
                for name, (samp, scores_dict) in d.items():
                    for sys_term, (scores, weights) in scores_dict.items():
                        if sys_term == 'NOMINAL':
                            suffix = ''
                        else:
                            sys_term = '_'.join(sys_term)
                            suffix = '_' + sys_term
                        hist = hist_template.Clone(name=name + suffix)
                        for score, w in zip(scores, weights):
                            hist.Fill(score, w)
                        f.cd()
                        hist.Write()

# save all variable plots in one large multipage pdf
if 'plot' in args.actions and set(args.categories) == set(CATEGORIES.keys()) and not args.plots:
    # only create multipage pdf of all plots if we created all plots
    from matplotlib.backends.backend_pdf import PdfPages
    import datetime
    now = datetime.datetime.today()
    # put all plots in a multipage pdf
    for category, exprs in figures.items():
        pdf = PdfPages(os.path.join(PLOTS_DIR, 'features_%s.pdf' % category))
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
