#!/usr/bin/env python

from argparse import ArgumentParser
from categories import CATEGORIES
from variables import VARIABLES

parser = ArgumentParser()
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
parser.add_argument('--actions', nargs='*', choices=['plot', 'train'])
parser.add_argument('--no-systematics', action='store_false',
        dest='systematics',
        help="turn off systematics",
        default=True)
parser.add_argument('--nfold', type=int, default=5)
parser.add_argument('--cor', action='store_true', default=False)
parser.add_argument('--unblind', action='store_true', default=False)
parser.add_argument('--train-frac', type=float, default=.5)
parser.add_argument('--categories', nargs='*', default=CATEGORIES.keys())
parser.add_argument('--plots', nargs='*')
args = parser.parse_args()

# root imports
import ROOT
ROOT.gROOT.SetBatch(True)

# local imports
from utils import *
from matplotlib import cm
from samples import *
import samples
import bkg_scales_cache
from background_estimation import qcd_ztautau_norm
from classify import *
from config import plots_dir
from utils import *

# stdlib imports
import pickle
import os
from pprint import pprint

# scikit-learn imports
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

# numpy imports
import numpy as np

# matplotlib imports
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import IndexLocator, FuncFormatter

# rootpy imports
from rootpy.plotting import Hist
from rootpy.io import open as ropen
from tabulartext import PrettyTable


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

PLOTS_DIR = plots_dir(__file__)

bins=20

#ztautau   = MC_Ztautau(systematics=args.systematics)
ztautau = Embedded_Ztautau(systematics=args.systematics)
others = Others(systematics=args.systematics)
data = Data(markersize=1.2)

higgs_125 = Higgs(
        masses=[125],
        systematics=args.systematics,
        scale=50,
        linecolor='red',
        linestyle='dashed')

if 'train' in args.actions:

    # all modes, all masses
    signals_train = [
        Higgs(),
    ]

    # all modes, 125GeV mass
    signals_test = [
        Higgs(mass=125),
    ]

figures = {}

for category, cat_info in sorted(CATEGORIES.items(), key=lambda item: item[0]):

    if category not in args.categories:
        continue

    # QCD shape region SS or !OS
    qcd_shape_region = cat_info['qcd_shape_region']
    target_region = cat_info['target_region']

    qcd = QCD(data=data, mc=[others, ztautau],
          shape_region=qcd_shape_region)

    figures[category] = {}

    #cuts = Cut('80 < mass_mmc_tau1_tau2 < 110')
    cuts = Cut()

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
        use_cache=args.use_fit_cache)

    qcd.scale = qcd_scale
    qcd.scale_error = qcd_scale_error
    ztautau.scale = ztautau_scale
    ztautau.scale_error = ztautau_scale_error

    if 'plot' in args.actions:
        for expr, var_info in VARIABLES.items():

            if category.upper() not in var_info['cats']:
                continue

            if args.plots and expr not in args.plots:
                continue

            print
            print "plotting %s in %s category" % (expr, category)

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

            fig = draw(
                    data=data_hist, model=bkg_hists,
                    signal=signal_hist,
                    name=var_info['title'],
                    output_name=var_info['filename'],
                    category_name=cat_info['name'],
                    category=category,
                    units=var_info.get('units', None),
                    range=var_info['range'],
                    show_ratio=True,
                    show_qq=False,
                    model_colour_map=None,
                    dir=PLOTS_DIR,
                    systematics=SYSTEMATICS if args.systematics else None)
            figures[category][expr] = fig

    if category == 'preselection':
        # don't train a classifier at preselection
        continue

    if 'train' in args.actions:

        backgrounds = [
            qcd,
            ztautau,
            others,
        ]

        # define training and test samples
        branches = cat_info['features']

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
        if args.use_clf_cache and os.path.isfile('clf_%s.pickle' % category):
            # use a previously trained classifier
            with open('clf_%s.pickle' % category, 'r') as f:
                clf = pickle.load(f)
            print clf
        else:
            cv = StratifiedKFold(labels_train, args.nfold)
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
                         cv=StratifiedKFold(labels_train, args.nfold))
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
                            min_samples_leaf=int(min_samples_leaf * args.nfold / float(args.nfold - 1))),
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
        plot_clf(
            clf,
            backgrounds,
            category,
            target_region,
            branches,
            signals=None,
            data=data,
            cuts=cuts,
            train_fraction=args.train_fraction,
            draw_data=True,
            name='control')

        # show the background model and 125 GeV signal above mass control region
        cuts = Cut('mass_mmc_tau1_tau2 > 110')
        plot_clf(
            clf,
            backgrounds,
            category,
            target_region,
            branches,
            signals=signals_test,
            signal_scale=20,
            cuts=cuts,
            train_fraction=args.train_fraction,
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
            for mass in Higgs.MASS_POINTS:
                for mode in Higgs.MODES.keys():
                    scores, weight = apply_clf(clf,
                        Higgs(mode=mode, mass=mass),
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
                            (mass, mode), scores, weight))

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
