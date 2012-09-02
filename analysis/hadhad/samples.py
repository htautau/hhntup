import os
import sys
import atexit
from operator import add
import math

import numpy as np
from numpy.lib import recfunctions

# for reproducibilty
np.random.seed(1987) # my birth year ;)

from higgstautau.hadhad.periods import total_lumi
from higgstautau import datasets
from higgstautau.decorators import cached_property, memoize_method
from higgstautau import samples as samples_db

# Higgs cross sections
import yellowhiggs

from rootpy.plotting import Hist, Canvas, HistStack
from rootpy.io import open as ropen
from rootpy.tree import Tree, Cut
from rootpy.utils import asrootpy
from rootpy import root2array as r2a
from rootpy.math.stats.correlation import correlation_plot

import categories
import features


NTUPLE_PATH = os.getenv('HIGGSTAUTAU_NTUPLE_DIR')
if not NTUPLE_PATH:
    sys.exit("You did not source setup.sh")
NTUPLE_PATH = os.path.join(NTUPLE_PATH, 'hadhad')
DEFAULT_STUDENT = 'HHProcessor'
DEFAULT_TREENAME = 'higgstautauhh'
TOTAL_LUMI = total_lumi()
TAUTAUHADHADBR = 0.412997
VERBOSE = False
DB_HH = datasets.Database(name='datasets_hh', verbose=VERBOSE)
DB_TAUID = datasets.Database(name='datasets_tauid', verbose=VERBOSE)
FILES = {}
WORKING_POINT = 'Tight'
ID = Cut('tau1_JetBDTSig%s==1 && tau2_JetBDTSig%s==1' %
         (WORKING_POINT, WORKING_POINT))
NOID = Cut('tau1_JetBDTSig%s!=1 && tau2_JetBDTSig%s!=1' %
           (WORKING_POINT, WORKING_POINT))
OS = Cut('tau1_charge * tau2_charge == -1')
NOT_OS = Cut('tau1_charge * tau2_charge != -1')
SS = Cut('tau1_charge * tau2_charge == 1')
# mass_jet1_jet2 > 100000
TEMPFILE = ropen('tempfile.root', 'recreate')


#@atexit.register
def cleanup():

    TEMPFILE.Close()
    os.unlink(TEMPFILE.GetName())


def std(X):

    return (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)


def correlations(signal, signal_weight,
                 background, background_weight,
                 branches, channel):

    # draw correlation plots
    names = [features.VARIABLES[branch]['title'] for branch in branches]
    correlation_plot(signal, signal_weight, names,
                     "correlation_signal_%s" % channel)
    correlation_plot(background, background_weight, names,
                     "correlation_background_%s" % channel)


def get_samples(masses=None, modes=None, embedding=True):

    if embedding:
        ztautau = Embedded_Ztautau()
    else:
        ztautau = MC_Ztautau()
    ewk = EWK()
    others = Others()

    backgrounds = (
        ztautau,
        ewk,
        others,
    )

    signals = (
        Higgs(masses=masses, modes=modes)
    )
    return signals, backgrounds


def make_classification(
        signals,
        backgrounds,
        category,
        region,
        branches,
        cuts=None,
        train_fraction=None,
        max_sig_train=None,
        max_bkg_train=None,
        max_sig_test=None,
        max_bkg_test=None,
        norm_sig_to_bkg_train=True,
        norm_sig_to_bkg_test=False,
        same_size_train=True,
        same_size_test=False,
        standardize=False,
        systematic='NOMINAL'):

    signal_train_arrs = []
    signal_weight_train_arrs = []
    signal_test_arrs = []
    signal_weight_test_arrs = []

    for signal in signals:
        train, test = signal.train_test(
            category, region,
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
            category, region,
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


class Sample(object):

    REGIONS = {
        'ALL': Cut(),
        'OS': OS,
        '!OS': NOT_OS,
        'SS': SS,
        'OS-ID': OS & ID,
        '!OS-ID': NOT_OS & ID,
        'SS-ID': SS & ID,
        'OS-NOID': OS & NOID,
        '!OS-NOID': NOT_OS & NOID,
        'SS-NOID': SS & NOID}

    CATEGORIES = dict([
        (name, Cut('category==%d' % info['code']))
        if info['code'] is not None
        else (name, Cut(''))
        for name, info in categories.CATEGORIES.items()])

    WEIGHT_BRANCHES = [
        'mc_weight',
        'pileup_weight',
        'ggf_weight',
        # effic high and low already accounted for in TAUBDT_UP/DOWN
        'tau1_efficiency_scale_factor',
        'tau2_efficiency_scale_factor',
    ]

    WEIGHT_SYSTEMATICS = {
        'TRIGGER': {
            'UP': [
                'tau1_trigger_scale_factor_high',
                'tau2_trigger_scale_factor_high'],
            'DOWN': [
                'tau1_trigger_scale_factor_low',
                'tau2_trigger_scale_factor_low'],
            'NOMINAL': [
                'tau1_trigger_scale_factor',
                'tau2_trigger_scale_factor']},
        'FAKERATE': {
            'UP': [
                'tau1_fakerate_scale_factor_high',
                'tau2_fakerate_scale_factor_high'],
            'DOWN': [
                'tau1_fakerate_scale_factor_low',
                'tau2_fakerate_scale_factor_low'],
            'NOMINAL': [
            	'tau1_fakerate_scale_factor',
                'tau2_fakerate_scale_factor']},
    }


    def __init__(self, scale=1., cuts=None,
                 student=DEFAULT_STUDENT,
                 treename=DEFAULT_TREENAME,
                 **hist_decor):

        self.scale = scale
        if cuts is None:
            self._cuts = Cut()
        else:
            self._cuts = cuts
        self.student = student
        self.treename = treename
        self.hist_decor = hist_decor
        if isinstance(self, Higgs):
            self.hist_decor['fillstyle'] = 'hollow'
        else:
            self.hist_decor['fillstyle'] = 'solid'

    def check_systematic(self, systematic):

        if systematic != 'NOMINAL' and isinstance(self, Data):
            raise TypeError('Do not apply systematics on data!')

    def get_weight_branches(self, systematic):

        self.check_systematic(systematic)
        weight_branches = Sample.WEIGHT_BRANCHES[:]
        if systematic == 'NOMINAL':
            systerm = None
            variation = 'NOMINAL'
        else:
            systerm, variation = systematic.split('_')
        for term, variations in Sample.WEIGHT_SYSTEMATICS.items():
            if term == systerm:
                weight_branches += variations[variation]
            else:
                weight_branches += variations['NOMINAL']
        return weight_branches

    def iter_weight_branches(self):

        for type, variations in Sample.WEIGHT_SYSTEMATICS.items():
            for variation in variations:
                if variation == 'NOMINAL':
                    continue
                term = '%s_%s' % (type, variation)
                yield self.get_weight_branches(term), tuple([term])

    def cuts(self, category, region):

        return (Sample.CATEGORIES[category] &
                categories.CATEGORIES[category]['cuts'] &
                Sample.REGIONS[region] & self._cuts)

    def train_test(self,
                   category,
                   region,
                   branches,
                   train_fraction=None,
                   cuts=None,
                   systematic='NOMINAL'):
        """
        Return recarray for training and for testing
        """
        self.check_systematic(systematic)
        weight_branches = self.get_weight_branches(systematic)

        if train_fraction is not None:

            assert 0 < train_fraction < 1.

            if isinstance(self, MC):
                branches = branches + weight_branches

            train_arrs = []
            test_arrs = []

            for tree in self.trees(
                    category,
                    region,
                    cuts=cuts,
                    systematic=systematic):
                arr = r2a.tree_to_recarray(
                    tree,
                    branches=branches,
                    include_weight=True,
                    weight_name='weight')
                if isinstance(self, MC):
                    # merge the three weight columns
                    arr['weight'] *= reduce(np.multiply,
                            [arr[br] for br in weight_branches])
                    # remove the mc_weight and pileup_weight fields
                    arr = recfunctions.rec_drop_fields(
                        arr, weight_branches)
                split_idx = int(train_fraction * float(arr.shape[0]))
                arr_train, arr_test = arr[:split_idx], arr[split_idx:]
                # scale the weights to account for train_fraction
                arr_train['weight'] *= 1. / train_fraction
                arr_test['weight'] *= 1. / (1. - train_fraction)
                train_arrs.append(arr_train)
                test_arrs.append(arr_test)
            arr_train, arr_test = np.hstack(train_arrs), np.hstack(test_arrs)
        else:
            arr = self.recarray(
                category,
                region,
                branches,
                include_weight=True,
                cuts=cuts,
                systematic=systematic)
            arr_train, arr_test = arr, arr

        return arr_train, arr_test

    def recarray(self,
                 category,
                 region,
                 branches,
                 include_weight=True,
                 cuts=None,
                 systematic='NOMINAL'):

        self.check_systematic(systematic)

        weight_branches = self.get_weight_branches(systematic)
        if include_weight and isinstance(self, MC):
            branches = branches + weight_branches

        try:
            arr = r2a.tree_to_recarray(
                self.trees(
                    category,
                    region,
                    cuts=cuts,
                    systematic=systematic),
                branches=branches,
                include_weight=include_weight,
                weight_name='weight')
        except IOError, e:
            raise IOError("%s: %s" % (self.__class__.__name__, e))

        if include_weight and isinstance(self, MC):
            # merge the three weight columns
            arr['weight'] *= reduce(np.multiply,
                    [arr[br] for br in weight_branches])

            # remove the mc_weight and pileup_weight fields
            arr = recfunctions.rec_drop_fields(
                arr, weight_branches)
        return arr

    def ndarray(self,
                category,
                region,
                branches,
                include_weight=True,
                cuts=None,
                systematic='NOMINAL'):

        self.check_systematic(systematic)
        return r2a.recarray_to_ndarray(
                   self.recarray(
                       category,
                       region,
                       branches=branches,
                       include_weight=include_weight,
                       cuts=cuts,
                       systematic=systematic))


class Data(Sample):

    def __init__(self, **kwargs):

        super(Data, self).__init__(scale=1., **kwargs)
        self.DATA_FILE = ropen('.'.join([os.path.join(NTUPLE_PATH, self.student, self.student),
                                'data-JetTauEtmiss.root']))
        self.data = self.DATA_FILE.Get(self.treename)
        self.label = ('2011 Data $\sqrt{s} = 7$ TeV\n'
                      '$\int L dt = %.2f$ fb$^{-1}$' % (TOTAL_LUMI / 1e3))
        self.name = 'Data'

    def draw(self, expr, category, region, bins, min, max, cuts=None):

        hist = Hist(bins, min, max, title=self.label, name=self.name,
                **self.hist_decor)
        self.draw_into(hist, expr, category, region, cuts=cuts)
        return hist

    def draw_into(self, hist, expr, category, region, cuts=None):

        self.data.draw(expr, self.cuts(category, region) & cuts, hist=hist)

    def trees(self,
              category,
              region,
              cuts=None,
              systematic='NOMINAL'):

        self.check_sytematic(systematic)
        TEMPFILE.cd()
        return [asrootpy(self.data.CopyTree(self.cuts(category, region) & cuts))]


class Signal:
    pass


class Background:
    pass


class MC(Sample):

    SYSTEMATICS = [
        ('TAUBDT_UP',),
        ('TAUBDT_DOWN',),
        ('JES_UP', 'TES_UP'),
        ('JES_DOWN','TES_DOWN'),
        ('JER_UP',),
        ('MFS_UP',),
        ('MFS_DOWN',),
        ('ISOL_UP',),
        ('ISOL_DOWN',),
    ]

    SYSTEMATICS_BY_WEIGHT = [
        ('TRIGGER_UP',),
        ('TRIGGER_DOWN',),
        ('FAKERATE_UP',),
        ('FAKERATE_DOWN',),
    ]

    def __init__(self,
            systematics=True,
            systematics_terms=None,
            systematics_samples=None,
            db=DB_HH,
            **kwargs):

        super(MC, self).__init__(**kwargs)
        self.db = db
        self.datasets = []
        self.systematics = systematics

        for i, name in enumerate(self.samples):

            ds = self.db[name]
            trees = {}
            weighted_events = {}

            trees['NOMINAL'] = None
            weighted_events['NOMINAL'] = None

            if isinstance(self, Embedded_Ztautau):
                events_bin = 0
            else:
                events_bin = 1

            if ds.name in FILES and 'NOMINAL' in FILES[ds.name]:
                rfile = FILES[ds.name]['NOMINAL']
                trees['NOMINAL'] = rfile.Get(self.treename)
                weighted_events['NOMINAL'] = rfile.cutflow[events_bin]
            else:
                rfile = ropen('.'.join([
                    os.path.join(NTUPLE_PATH, self.student, self.student), ds.name, 'root']))
                trees['NOMINAL'] = rfile.Get(self.treename)
                weighted_events['NOMINAL'] = rfile.cutflow[events_bin]
                if ds.name not in FILES:
                    FILES[ds.name] = {}
                FILES[ds.name]['NOMINAL'] = rfile

            if self.systematics:

                unused_terms = MC.SYSTEMATICS[:]

                if systematics_terms:
                    for sys_term in systematics_terms:

                        # merge terms such as JES_UP,TES_UP (embedding) and TES_UP (MC)
                        actual_sys_term = sys_term
                        for term in unused_terms:
                            if set(term) & set(sys_term):
                                if len(sys_term) < len(term):
                                    print "merging %s and %s" % (term, sys_term)
                                    sys_term = term
                                break

                        trees[sys_term] = None
                        weighted_events[sys_term] = None

                        if ds.name in FILES and sys_term in FILES[ds.name]:
                            rfile = FILES[ds.name][sys_term]
                            trees[sys_term] = rfile.Get(self.treename)
                            weighted_events[sys_term] = rfile.cutflow[events_bin]
                        else:
                            rfile = ropen('.'.join([
                                os.path.join(NTUPLE_PATH, self.student, self.student),
                                '_'.join([ds.name, '_'.join(actual_sys_term)]), 'root']))
                            trees[sys_term] = rfile.Get(self.treename)
                            weighted_events[sys_term] = rfile.cutflow[events_bin]
                            if ds.name not in FILES:
                                FILES[ds.name] = {}
                            FILES[ds.name][sys_term] = rfile

                        unused_terms.remove(sys_term)

                if systematics_samples and name in systematics_samples:
                    for sample_name, sys_term in systematics_samples[name].items():

                        sys_term = tuple(sys_term.split(','))
                        sys_ds = self.db[sample_name]
                        trees[sys_term] = None
                        weighted_events[sys_term] = None

                        if sys_ds.name in FILES and sys_term in FILES[sys_ds.name]:
                            rfile = FILES[sys_ds.name][sys_term]
                            trees[sys_term] = rfile.Get(self.treename)
                            weighted_events[sys_term] = rfile.cutflow[events_bin]
                        else:
                            rfile = ropen('.'.join([
                                os.path.join(NTUPLE_PATH, self.student, self.student),
                                sys_ds.name, 'root']))
                            trees[sys_term] = rfile.Get(self.treename)
                            weighted_events[sys_term] = rfile.cutflow[events_bin]
                            if sys_ds.name not in FILES:
                                FILES[sys_ds.name] = {}
                            FILES[sys_ds.name][sys_term] = rfile

                        unused_terms.remove(sys_term)

                if unused_terms:
                    #print "UNUSED TERMS for %s:" % self.name
                    #print unused_terms

                    for term in unused_terms:
                        trees[term] = None # flag to use NOMINAL
                        weighted_events[term] = None # flag to use NOMINAL

            if isinstance(self, Higgs):
                # use yellowhiggs for cross sections
                xs = yellowhiggs.xsbr(
                        7, self.masses[i],
                        self.modes[i], 'tautau')[0] * TAUTAUHADHADBR
                kfact = 1.
                effic = 1.
            elif isinstance(self, Embedded_Ztautau):
                xs, kfact, effic = 1., 1., 1.
            else:
                xs, kfact, effic = ds.xsec_kfact_effic
            if VERBOSE:
                print ds.name, xs, kfact, effic,
            self.datasets.append((ds, trees, weighted_events, xs, kfact, effic))

    @property
    def label(self):

        l = self._label
        if self.scale != 1. and not isinstance(self,
                (MC_Ztautau, Embedded_Ztautau)):
            l += r' ($\sigma_{SM} \times %g$)' % self.scale
        return l

    def draw(self, expr, category, region, bins, min, max, cuts=None):

        hist = Hist(bins, min, max, title=self.label, name=self.name,
                **self.hist_decor)
        self.draw_into(hist, expr, category, region, cuts=cuts)
        return hist

    def draw_into(self, hist, expr, category, region, cuts=None):

        selection = self.cuts(category, region) & cuts
        if isinstance(expr, (list, tuple)):
            exprs = expr
        else:
            exprs = (expr,)

        if hasattr(hist, 'systematics'):
            sys_hists = hist.systematics
        else:
            sys_hists = {}

        for ds, sys_trees, sys_events, xs, kfact, effic in self.datasets:

            nominal_tree = sys_trees['NOMINAL']
            nominal_events = sys_events['NOMINAL']

            nominal_weight = (TOTAL_LUMI * self.scale *
                    xs * kfact * effic / nominal_events)

            nominal_weighted_selection = (
                '%f * %s * (%s)' %
                (nominal_weight,
                 ' * '.join(self.get_weight_branches('NOMINAL')),
                 selection))

            if VERBOSE:
                print nominal_weighted_selection

            # fill nominal histogram
            for expr in exprs:
                nominal_tree.Draw(expr, nominal_weighted_selection, hist=hist)

            if not self.systematics:
                continue

            # iterate over systematic variation trees
            for sys_term in sys_trees.iterkeys():

                # skip the nominal tree
                if sys_term == 'NOMINAL':
                    continue

                sys_hist = hist.Clone()

                sys_tree = sys_trees[sys_term]
                sys_event = sys_events[sys_term]

                if sys_tree is not None:

                    sys_hist.Reset()

                    sys_weight = (TOTAL_LUMI * self.scale *
                            xs * kfact * effic / sys_event)

                    sys_weighted_selection = (
                        '%f * %s * (%s)' %
                        (sys_weight,
                         ' * '.join(self.get_weight_branches('NOMINAL')),
                         selection))

                    for expr in exprs:
                        sys_tree.Draw(expr, sys_weighted_selection, hist=sys_hist)

                if sys_term not in sys_hists:
                    sys_hists[sys_term] = sys_hist
                else:
                    sys_hists[sys_term] += sys_hist

            # iterate over weight systematics on the nominal tree
            for weight_branches, sys_term in self.iter_weight_branches():

                sys_hist = hist.Clone()
                sys_hist.Reset()

                weighted_selection = (
                    '%f * %s * (%s)' %
                    (nominal_weight,
                     ' * '.join(weight_branches),
                     selection))

                for expr in exprs:
                    nominal_tree.Draw(expr, weighted_selection, hist=sys_hist)

                if sys_term not in sys_hists:
                    sys_hists[sys_term] = sys_hist
                else:
                    sys_hists[sys_term] += sys_hist

        # set the systematics
        hist.systematics = sys_hists

    def trees(self, category, region, cuts=None,
              systematic='NOMINAL'):

        TEMPFILE.cd()
        selection = self.cuts(category, region) & cuts
        trees = []
        for ds, sys_trees, sys_events, xs, kfact, effic in self.datasets:
            tree = sys_trees[systematic]
            events = sys_events[systematic]

            weight = TOTAL_LUMI * self.scale * xs * kfact * effic / events

            selected_tree = asrootpy(tree.CopyTree(selection))
            selected_tree.SetWeight(weight)
            trees.append(selected_tree)
        return trees

    def events(self, selection='', systematic='NOMINAL'):

        total = 0.
        for ds, sys_trees, sys_events, xs, kfact, effic in self.datasets:
            tree = sys_trees[systematic]
            events = sys_events[systematic]

            weight = TOTAL_LUMI * self.scale * xs * kfact * effic / events

            total += weight * tree.GetEntries(selection)
        return total

    def iter(self, selection='', systematic='NOMINAL'):

        TEMPFILE.cd()
        for ds, sys_trees, sys_events, xs, kfact, effic in self.datasets:
            tree = sys_trees[systematic]
            events = sys_events[systematic]

            weight = TOTAL_LUMI * self.scale * xs * kfact * effic / events

            if selection:
                selected_tree = asrootpy(tree.CopyTree(selection))
            else:
                selected_tree = tree
            for event in selected_tree:
                yield weight, event


class MC_Ztautau(MC, Background):

    def __init__(self, color='#00a4ff', **kwargs):
        """
        Instead of setting the k factor here
        the normalization is determined by a fit to the data
        """
        yml = samples_db.BACKGROUNDS['hadhad']['ztautau']
        self.name = 'Ztautau'
        self._label = yml['latex']
        self.samples = yml['samples']
        syst = samples_db.SYSTEMATICS['hadhad'][yml['systematics']]
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(MC_Ztautau, self).__init__(
                color=color,
                systematics_terms=systematics_terms,
                **kwargs)


class Embedded_Ztautau(MC, Background):

    def __init__(self, color='#00a4ff', **kwargs):
        """
        Instead of setting the k factor here
        the normalization is determined by a fit to the data
        """
        yml = samples_db.BACKGROUNDS['hadhad']['embedded_ztautau']
        self.name = 'Ztautau'
        self._label = yml['latex']
        self.samples = yml['samples']
        systematics_samples = yml['systematics_samples']
        syst = samples_db.SYSTEMATICS['hadhad'][yml['systematics']]
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(Embedded_Ztautau, self).__init__(
                color=color,
                systematics_samples=systematics_samples,
                systematics_terms=systematics_terms,
                **kwargs)


class EWK(MC, Background):

    def __init__(self, color='#ff9f71', **kwargs):

        yml = samples_db.BACKGROUNDS['hadhad']['ewk']
        self.name = 'EWK'
        self._label = yml['latex']
        self.samples = yml['samples']
        syst = samples_db.SYSTEMATICS['hadhad'][yml['systematics']]
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(EWK, self).__init__(
                color=color,
                systematics_terms=systematics_terms,
                **kwargs)


class Top(MC, Background):

    def __init__(self, color='#0000ff', **kwargs):

        yml = samples_db.BACKGROUNDS['hadhad']['top']
        self.name = 'Top'
        self._label = yml['latex']
        self.samples = yml['samples']
        syst = samples_db.SYSTEMATICS['hadhad'][yml['systematics']]
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(Top, self).__init__(
                color=color,
                systematics_terms=systematics_terms,
                **kwargs)


class Diboson(MC, Background):

    def __init__(self, color='#ffd075', **kwargs):

        yml = samples_db.BACKGROUNDS['hadhad']['diboson']
        self.name = 'Diboson'
        self._label = yml['latex']
        self.samples = yml['samples']
        syst = samples_db.SYSTEMATICS['hadhad'][yml['systematics']]
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(Diboson, self).__init__(
                color=color,
                systematics_terms=systematics_terms,
                **kwargs)


class Others(MC, Background):

    def __init__(self, color='#ff7700', **kwargs):

        yml_diboson = samples_db.BACKGROUNDS['hadhad']['diboson']
        yml_top = samples_db.BACKGROUNDS['hadhad']['top']
        yml_ewk = samples_db.BACKGROUNDS['hadhad']['ewk']
        self.samples = (yml_diboson['samples'] +
                        yml_top['samples'] +
                        yml_ewk['samples'])
        self._label = 'Others'
        self.name = 'Others'
        syst = samples_db.SYSTEMATICS['hadhad']['mc']
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(Others, self).__init__(
                color=color,
                systematics_terms=systematics_terms,
                **kwargs)


class Higgs(MC, Signal):

    MASS_POINTS = range(100, 155, 5)

    MODES = {
        'ggf': ('ggH', 'PowHegPythia_'),
        'vbf': ('VBFH', 'PowHegPythia_'),
        'zh': ('ZH', 'Pythia'),
        'wh': ('WH', 'Pythia'),
    }

    def __init__(self, modes=None, masses=None, **kwargs):

        if masses is None:
            masses = Higgs.MASS_POINTS
        else:
            assert len(masses) > 0
            for mass in masses:
                assert mass in Higgs.MASS_POINTS
            assert len(set(masses)) == len(masses)

        if modes is None:
            modes = Higgs.MODES.keys()
        else:
            assert len(modes) > 0
            for mode in modes:
                assert mode in Higgs.MODES
            assert len(set(modes)) == len(modes)

        str_mass = ''
        if len(masses) == 1:
            str_mass = str(masses[0])

        str_mode = ''
        if len(modes) == 1:
            str_mode = str(modes[0]) + ' '

        self._label = r'%s$H%s\rightarrow\tau_{h}\tau_{h}$' % (
                str_mode, str_mass)

        self.name = '{mode}Signal{mass}'.format(
                mass=str_mass,
                mode=str_mode.strip())

        self.samples = []
        self.masses = []
        self.modes = []
        for mode in modes:
            mode_str = Higgs.MODES[mode][0]
            generator = Higgs.MODES[mode][1]
            for mass in masses:
                self.samples.append('%s%s%d_tautauhh.mc11c' % (
                    generator, mode_str, mass))
                self.masses.append(mass)
                self.modes.append(mode)

        syst = samples_db.SYSTEMATICS['hadhad']['mc']
        systematics_terms = [tuple(term.split(',')) for term in syst]
        super(Higgs, self).__init__(
                systematics_terms=systematics_terms,
                **kwargs)


class QCD(Sample):

    def __init__(self, data, mc,
                 scale=1.,
                 shape_region='SS',
                 cuts=None,
                 color='#59d454'):

        super(QCD, self).__init__(scale=scale, color=color)
        self.data = data
        self.mc = mc
        self.name = 'QCD'
        self.label = 'QCD Multi-jet'
        self.scale = 1.
        self.shape_region = shape_region

    def draw(self, expr, category, region, bins, min, max, cuts=None):

        hist = Hist(bins, min, max, title=self.label, name=self.name,
                **self.hist_decor)
        self.draw_into(hist, expr, category, region, cuts=cuts)
        return hist

    def draw_into(self, hist, expr, category, region, cuts=None):

        MC_bkg_notOS = hist.Clone()
        for mc in self.mc:
            mc.draw_into(MC_bkg_notOS, expr, category, self.shape_region,
                         cuts=cuts)

        data_hist = hist.Clone()
        self.data.draw_into(data_hist, expr,
                            category, self.shape_region, cuts=cuts)

        hist += (data_hist - MC_bkg_notOS) * self.scale
        if hasattr(MC_bkg_notOS, 'systematics'):
            hist.systematics = {}
            for sys_term, sys_hist in MC_bkg_notOS.systematics.items():
                hist.systematics[sys_term] = (data_hist - sys_hist) * self.scale

        hist.SetTitle(self.label)

    def scores(self,
               clf,
               category,
               region,
               branches,
               train_fraction,
               cuts=None,
               systematic='NOMINAL'):

        # SS data
        train, test = self.data.train_test(category=category,
                                           region=self.shape_region,
                                           branches=branches,
                                           train_fraction=train_fraction,
                                           cuts=cuts)
        weight = test['weight']
        sample = np.vstack(test[branch] for branch in branches).T
        scores = clf.predict_proba(sample)[:,-1]

        # subtract SS MC
        for mc in self.mc:
            # didn't train on MC here if using SS or !OS
            train, test = mc.train_test(category=category,
                                        region=self.shape_region,
                                        branches=branches,
                                        train_fraction=train_fraction,
                                        cuts=cuts,
                                        systematic=systematic)
            sample = np.vstack(test[branch] for branch in branches).T
            scores = np.concatenate((scores, clf.predict_proba(sample)[:,-1]))
            weight = np.concatenate((weight, test['weight'] * -1))

        weight *= self.scale
        return scores, weight

    def trees(self, category, region, cuts=None,
              systematic='NOMINAL'):

        TEMPFILE.cd()
        trees = [asrootpy(self.data.data.CopyTree(
                    self.data.cuts(category,
                                   region=self.shape_region) & cuts))]
        for mc in self.mc:
            _trees = mc.trees(
                    category,
                    region=self.shape_region,
                    cuts=cuts,
                    systematic=systematic)
            for tree in _trees:
                tree.Scale(-1)
            trees += _trees

        for tree in trees:
            tree.Scale(self.scale)
        return trees


class MC_TauID(MC):

    def __init__(self, **kwargs):

        self.name = 'TauID'
        self._label = 'TauID'
        self.samples = ['PythiaWtaunu_incl.mc11c']
        super(MC_TauID, self).__init__(student='TauIDProcessor',
                db=DB_TAUID, **kwargs)


if __name__ == '__main__':

    from pyAMI.query import print_table
    signals, backgrounds = get_samples('2jet', purpose='train')

    table = [['dataset',
              'mean sigma [pb]', 'min sigma [pb]', 'max sigma [pb]',
              'sigma factor',
              'filter effic', 'K factor']]
    format = "%s %.8f %.8f %.8f %.6f %.6f %.6f"
    for sample in signals + backgrounds:
        if isinstance(sample, QCD):
            continue
        for datasets in sample.datasets:
            for ds, tree, events, (xsec, kfact, effic) in datasets:
                row = format % (ds.ds, xsec*1E3, kfact, ds.xsec_factor, effic)
                table.append(row.split())
    print
    print
    table.sort(key=lambda row: row[0])
    print_table(table)
