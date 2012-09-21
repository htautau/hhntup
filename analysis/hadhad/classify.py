import numpy as np
# for reproducibilty
# especially for test/train set selection
np.random.seed(1987) # my birth year ;)

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.ticker import IndexLocator, FuncFormatter

from rootpy.plotting import Hist
from rootpy.io import open as ropen

from samples import *
from utils import draw

def plot_grid_scores(
        grid_scores, best_point, params, name,
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


def plot_clf(
        background_scores,
        category,
        category_name,
        signal_scores=None,
        signal_scale=1.,
        data_scores=None,
        name=None,
        draw_histograms=True,
        draw_data=False,
        save_histograms=False,
        bins=10,
        min_score=0,
        max_score=1,
        systematics=None):

    hist_template = Hist(bins, min_score - 0.00001, max_score + 0.00001)

    bkg_hists = []
    for bkg, scores_dict in background_scores:
        hist = hist_template.Clone(title=bkg.label)
        scores, weight = scores_dict['NOMINAL']
        for score, w in zip(scores, weight):
            hist.Fill(score, w)
        hist.decorate(**bkg.hist_decor)
        hist.systematics = {}
        for sys_term in scores_dict.keys():
            if sys_term == 'NOMINAL':
                continue
            sys_hist = hist_template.Clone()
            scores, weight = scores_dict[sys_term]
            for score, w in zip(scores, weight):
                sys_hist.Fill(score, w)
            hist.systematics[sys_term] = sys_hist
        bkg_hists.append(hist)

    if signal_scores is not None:
        signal, signal_scores = signal_scores
        sig_hist = hist_template.Clone(title=signal.label)
        scores, weight = signal_scores['NOMINAL']
        for score, w in zip(scores, weight):
            sig_hist.Fill(score, w)
        sig_hist.decorate(**signal.hist_decor)
        sig_hist.systematics = {}
        for sys_term in signal_scores.keys():
            if sys_term == 'NOMINAL':
                continue
            sys_hist = hist_template.Clone()
            scores, weight = signal_scores[sys_term]
            for score, w in zip(scores, weight):
                sys_hist.Fill(score, w)
            sig_hist.systematics[sys_term] = sys_hist
    else:
        sig_hist = None

    if data_scores is not None and draw_data:
        data, data_scores = data_scores
        data_hist = hist_template.Clone(title=data.label)
        data_hist.decorate(**data.hist_decor)
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
             signal=sig_hist,
             signal_scale=signal_scale,
             category=category,
             category_name=category_name,
             name="BDT Score",
             output_name=output_name,
             range=(min_score, max_score),
             show_ratio=data_hist is not None,
             model_colour_map=None,
             systematics=systematics)


def std(X):

    return (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)


def make_classification(
        signals,
        backgrounds,
        category,
        region,
        branches,
        train_fraction,
        cuts=None,
        max_sig_train=None,
        max_bkg_train=None,
        max_sig_test=None,
        max_bkg_test=None,
        norm_sig_to_bkg_train=True,
        norm_sig_to_bkg_test=False,
        same_size_train=True,
        same_size_test=False,
        standardize=False,
        remove_negative_train_weights=False,
        systematic='NOMINAL'):

    signal_train_arrs = []
    signal_weight_train_arrs = []
    signal_test_arrs = []
    signal_weight_test_arrs = []

    for signal in signals:
        train, test = signal.train_test(
            category=category,
            region=region,
            branches=branches,
            train_fraction=train_fraction,
            cuts=cuts,
            systematic=systematic)
        signal_weight_train_arrs.append(train['weight'])
        signal_weight_test_arrs.append(test['weight'])

        signal_train_arrs.append(
            np.vstack(train[branch] for branch in branches).T)
        signal_test_arrs.append(
            np.vstack(test[branch] for branch in branches).T)

    background_train_arrs = []
    background_weight_train_arrs = []
    background_test_arrs = []
    background_weight_test_arrs = []

    for background in backgrounds:
        train, test = background.train_test(
            category=category,
            region=region,
            branches=branches,
            train_fraction=train_fraction,
            cuts=cuts,
            systematic=systematic)
        background_weight_train_arrs.append(train['weight'])
        background_weight_test_arrs.append(test['weight'])

        background_train_arrs.append(
            np.vstack(train[branch] for branch in branches).T)
        background_test_arrs.append(
            np.vstack(test[branch] for branch in branches).T)

    signal_train = np.concatenate(signal_train_arrs)
    signal_weight_train = np.concatenate(signal_weight_train_arrs)
    signal_test = np.concatenate(signal_test_arrs)
    signal_weight_test = np.concatenate(signal_weight_test_arrs)

    background_train = np.concatenate(background_train_arrs)
    background_weight_train = np.concatenate(background_weight_train_arrs)
    background_test = np.concatenate(background_test_arrs)
    background_weight_test = np.concatenate(background_weight_test_arrs)

    if remove_negative_train_weights:
        # remove samples from the training sample with a negative weight
        signal_train = signal_train[signal_weight_train >= 0]
        background_train = background_train[background_weight_train >= 0]

        signal_weight_train = signal_weight_train[signal_weight_train >= 0]
        background_weight_train = background_weight_train[background_weight_train >= 0]

    if max_sig_train is not None and max_sig_train < len(signal_train):
        subsample = np.random.permutation(max_sig_train)[:len(signal_train)]
        signal_train = signal_train[subsample]
        signal_weight_train = signal_weight_train[subsample]

    if max_bkg_train is not None and max_bkg_train < len(background_train):
        subsample = np.random.permutation(max_bkg_train)[:len(background_train)]
        background_train = background_train[subsample]
        background_weight_train = background_weight_train[subsample]

    if max_sig_test is not None and max_sig_test < len(signal_test):
        subsample = np.random.permutation(max_sig_test)[:len(signal_test)]
        signal_test = signal_test[subsample]
        signal_weight_test = signal_weight_test[subsample]

    if max_bkg_test is not None and max_bkg_test < len(background_test):
        subsample = np.random.permutation(max_bkg_test)[:len(background_test)]
        background_test = background_test[subsample]
        background_weight_test = background_weight_test[subsample]

    if same_size_train:
        if len(background_train) > len(signal_train):
            # random subsample of background so it's the same size as signal
            subsample = np.random.permutation(
                len(background_train))[:len(signal_train)]
            background_train = background_train[subsample]
            background_weight_train = background_weight_train[subsample]
        elif len(background_train) < len(signal_train):
            # random subsample of signal so it's the same size as background
            subsample = np.random.permutation(
                len(signal_train))[:len(background_train)]
            signal_train = signal_train[subsample]
            signal_weight_train = signal_weight_train[subsample]

    if same_size_test:
        if len(background_test) > len(signal_test):
            # random subsample of background so it's the same size as signal
            subsample = np.random.permutation(
                len(background_test))[:len(signal_test)]
            background_test = background_test[subsample]
            background_weight_test = background_weight_test[subsample]
        elif len(background_test) < len(signal_test):
            # random subsample of signal so it's the same size as background
            subsample = np.random.permutation(
                len(signal_test))[:len(background_test)]
            signal_test = signal_test[subsample]
            signal_weight_test = signal_weight_test[subsample]

    if norm_sig_to_bkg_train:
        # normalize signal to background
        signal_weight_train *= (
            background_weight_train.sum() / signal_weight_train.sum())

    if norm_sig_to_bkg_test:
        # normalize signal to background
        signal_weight_test *= (
            background_weight_test.sum() / signal_weight_test.sum())

    print "Training Samples:"
    print "Signal: %d events, %s features" % signal_train.shape
    print "Sum(signal weights): %f" % signal_weight_train.sum()
    print "Background: %d events, %s features" % background_train.shape
    print "Sum(background weight): %f" % background_weight_train.sum()
    print
    print "Test Samples:"
    print "Signal: %d events, %s features" % signal_test.shape
    print "Sum(signal weights): %f" % signal_weight_test.sum()
    print "Background: %d events, %s features" % background_test.shape
    print "Sum(background weight): %f" % background_weight_test.sum()

    # create training/testing samples
    sample_train = np.concatenate((background_train, signal_train))
    sample_test = np.concatenate((background_test, signal_test))

    sample_weight_train = np.concatenate(
        (background_weight_train, signal_weight_train))
    sample_weight_test = np.concatenate(
        (background_weight_test, signal_weight_test))

    if standardize:
        sample_train = std(sample_train)
        sample_test = std(sample_test)

    labels_train = np.concatenate(
        (np.zeros(len(background_train)), np.ones(len(signal_train))))
    labels_test = np.concatenate(
        (np.zeros(len(background_test)), np.ones(len(signal_test))))

    # random permutation of training sample
    perm = np.random.permutation(len(labels_train))
    sample_train = sample_train[perm]
    sample_weight_train = sample_weight_train[perm]
    labels_train = labels_train[perm]

    # split the dataset in two equal parts respecting label proportions
    #train, test = iter(StratifiedKFold(labels, 2)).next()
    return sample_train, sample_test,\
        sample_weight_train, sample_weight_test,\
        labels_train, labels_test
