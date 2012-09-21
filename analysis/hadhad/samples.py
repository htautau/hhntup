# std lib imports
import os
import sys
import atexit
from operator import add
import math

# numpy imports
import numpy as np
from numpy.lib import recfunctions

# for reproducibilty
# especially for test/train set selection
np.random.seed(1987) # my birth year ;)

# pytables imports
import tables

# higgstautau imports
from higgstautau.hadhad.periods import total_lumi
from higgstautau import datasets
from higgstautau.decorators import cached_property, memoize_method
from higgstautau import samples as samples_db

# rootpy imports
from rootpy.plotting import Hist, Canvas, HistStack
from rootpy.io import open as ropen
from rootpy.tree import Tree, Cut
from rootpy.utils import asrootpy
from rootpy import root2array as r2a
from rootpy.math.stats.correlation import correlation_plot

# local imports
import categories
import variables
import systematics

# Higgs cross sections
import yellowhiggs


NTUPLE_PATH = os.getenv('HIGGSTAUTAU_NTUPLE_DIR')
if not NTUPLE_PATH:
    sys.exit("You did not source setup.sh")
NTUPLE_PATH = os.path.join(NTUPLE_PATH, 'hadhad')
DEFAULT_STUDENT = 'HHProcessor'
TOTAL_LUMI = total_lumi()
TAUTAUHADHADBR = 0.412997
VERBOSE = False
DB_HH = datasets.Database(name='datasets_hh', verbose=VERBOSE)
DB_TAUID = datasets.Database(name='datasets_tauid', verbose=VERBOSE)
FILES = {}
WORKING_POINT = 'Tight'
ID = Cut('(tau1_JetBDTSig%s==1) && (tau2_JetBDTSig%s==1)' %
         (WORKING_POINT, WORKING_POINT))
NOID = Cut('(tau1_JetBDTSig%s!=1) && (tau2_JetBDTSig%s!=1)' %
           (WORKING_POINT, WORKING_POINT))
OS = Cut('tau1_charge * tau2_charge == -1')
NOT_OS = Cut('tau1_charge * tau2_charge != -1')
SS = Cut('tau1_charge * tau2_charge == 1')
# mass_jet1_jet2 > 100000
TEMPFILE = ropen('tempfile.root', 'recreate')


def get_file(student, hdf=False):

    if hdf:
        ext = '.h5'
    else:
        ext = '.root'
    filename = student + ext
    if filename in FILES:
        return FILES[filename]
    if hdf:
        student_file = tables.openFile(
                os.path.join(NTUPLE_PATH, student, filename))
    else:
        student_file = ropen(
                os.path.join(NTUPLE_PATH, student, filename), 'READ')
    FILES[filename] = student_file
    return student_file


@atexit.register
def cleanup():

    TEMPFILE.Close()
    os.unlink(TEMPFILE.GetName())
    for filehandle in FILES.values():
        filehandle.close()


def correlations(signal, signal_weight,
                 background, background_weight,
                 branches, category):

    # draw correlation plots
    names = [variables.VARIABLES[branch]['title'] for branch in branches]
    correlation_plot(signal, signal_weight, names,
                     "correlation_signal_%s" % category)
    correlation_plot(background, background_weight, names,
                     "correlation_background_%s" % category)


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

    CATEGORIES = dict(
        [(name, Cut('category==%d' % info['code']))
         for name, info in categories.CATEGORIES.items()] +
        [(name, Cut(''))
         for name, info in categories.CONTROLS.items()])

    WEIGHT_BRANCHES = [
        'mc_weight',
        'pileup_weight',
        #'ggf_weight', # <= affects limits a lot!
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
                 **hist_decor):

        self.scale = scale
        if cuts is None:
            self._cuts = Cut()
        else:
            self._cuts = cuts
        self.student = student
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
        elif len(systematic) > 1:
            # no support for this yet...
            systerm = None
            variation = 'NOMINAL'
        else:
            systerm, variation = systematic[0].split('_')
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
                term = ('%s_%s' % (type, variation),)
                yield self.get_weight_branches(term), term

    def cuts(self, category, region):

        if category in categories.CATEGORIES:
            return (Sample.CATEGORIES[category] &
                    categories.CATEGORIES[category]['cuts'] &
                    Sample.REGIONS[region] & self._cuts)
        elif category in categories.CONTROLS:
            return (Sample.CATEGORIES[category] &
                    categories.CONTROLS[category]['cuts'] &
                    Sample.REGIONS[region] & self._cuts)
        else:
            raise ValueError(
                    'no such category or control region: %s' % category)

    def train_test(self,
                   branches,
                   category,
                   region,
                   train_fraction,
                   cuts=None,
                   systematic='NOMINAL'):
        """
        Return recarray for training and for testing
        """
        self.check_systematic(systematic)
        assert 0 < train_fraction < 1, "Train fraction must be between 0 and 1"
        arrays = self.tables(
                branches,
                category,
                region,
                include_weight=True,
                cuts=cuts,
                systematic=systematic)
        test_arrs = []
        train_arrs = []
        for arr in arrays:
            # split into test and train samples
            split_idx = int(train_fraction * arr.shape[0])
            arr_train, arr_test = arr[:split_idx], arr[split_idx:]
            test_arrs.append(arr_test)
            train_arrs.append(arr_train)
        arr_train = np.concatenate(train_arrs)
        arr_test = np.concatenate(test_arrs)
        # scale the weights to account for train_fraction
        arr_train['weight'] *= (1. / train_fraction)
        arr_test['weight'] *= (1. / (1. - train_fraction))
        return arr_train, arr_test

    def recarray(self,
                 branches,
                 category,
                 region,
                 include_weight=True,
                 cuts=None,
                 systematic='NOMINAL'):

        self.check_systematic(systematic)
        arrays = self.tables(
                branches,
                category,
                region,
                include_weight=include_weight,
                cuts=cuts,
                systematic=systematic)
        return np.concatenate(arrays)


class Data(Sample):

    def __init__(self, **kwargs):

        super(Data, self).__init__(scale=1., **kwargs)

        rfile = get_file(self.student)
        h5file = get_file(self.student, hdf=True)
        self.data = rfile.data_JetTauEtmiss
        self.h5data = h5file.root.data_JetTauEtmiss
        self.label = ('2011 Data $\sqrt{s} = 7$ TeV\n'
                      '$\int L dt = %.2f$ fb$^{-1}$' % (TOTAL_LUMI / 1e3))
        self.name = 'Data'

    def draw(self, expr, category, region, bins, min, max, cuts=None):

        hist = Hist(bins, min, max, title=self.label, **self.hist_decor)
        self.draw_into(hist, expr, category, region, cuts=cuts)
        return hist

    def draw_into(self, hist, expr, category, region, cuts=None):

        self.data.draw(expr, self.cuts(category, region) & cuts, hist=hist)

    def scores(self, clf, branches, train_fraction,
            category, region, cuts=None):

        if region == 'OS':
            # target region. Not trained on, so use all of it
            data_recarray = self.recarray(
                    branches=branches,
                    category=category,
                    region=region,
                    include_weight=True,
                    cuts=cuts,
                    systematic='NOMINAL')
        else:
            # ignore training sample
            _, data_recarray = self.train_test(
                    branches=branches,
                    category=category,
                    region=region,
                    train_fraction=train_fraction,
                    cuts=cuts)
        data_weights = data_recarray['weight']
        data_sample = np.vstack(data_recarray[branch] for branch in branches).T
        data_scores = clf.predict_proba(data_sample)[:,-1]
        return data_scores, data_weights

    def trees(self,
              category,
              region,
              cuts=None,
              systematic='NOMINAL'):

        self.check_systematic(systematic)
        TEMPFILE.cd()
        tree = asrootpy(self.data.CopyTree(self.cuts(category, region) & cuts))
        tree.userdata.weight_branches = []
        return [tree]

    def tables(self,
              branches,
              category,
              region,
              cuts=None,
              include_weight=True,
              systematic='NOMINAL'):

        self.check_systematic(systematic)
        selection = self.cuts(category, region) & cuts
        # read the table with a selection
        table = self.h5data.readWhere(selection.where())
        # add weight field
        if include_weight:
            # data is not weighted
            weights = np.ones(table.shape[0], dtype='f4')
            table = recfunctions.rec_append_fields(table, names='weight',
                    data=weights,
                    dtypes='f4')
            branches = branches + ['weight']
        # drop all branches except branches (+ weight)
        table = table[branches]
        return [table]


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
        rfile = get_file(self.student)
        h5file = get_file(self.student, hdf=True)

        for i, name in enumerate(self.samples):

            ds = self.db[name]
            treename = name.replace('.', '_')
            treename = treename.replace('-', '_')

            trees = {}
            tables = {}
            weighted_events = {}

            if isinstance(self, Embedded_Ztautau):
                events_bin = 0
                events_hist_suffix = '_cutflow_event'
            else:
                events_bin = 1
                events_hist_suffix = '_cutflow'

            trees['NOMINAL'] = rfile.Get(treename)
            tables['NOMINAL'] = getattr(h5file.root, treename)
            weighted_events['NOMINAL'] = rfile.Get(treename + events_hist_suffix)[events_bin]

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

                        sys_name = treename + '_' + '_'.join(actual_sys_term)
                        trees[sys_term] = rfile.Get(sys_name)
                        tables[sys_term] = getattr(h5file.root, sys_name)
                        weighted_events[sys_term] = rfile.Get(sys_name + events_hist_suffix)[events_bin]

                        unused_terms.remove(sys_term)

                if systematics_samples and name in systematics_samples:
                    for sample_name, sys_term in systematics_samples[name].items():

                        print "%s -> %s %s" % (name, sample_name, sys_term)

                        sys_term = tuple(sys_term.split(','))
                        sys_ds = self.db[sample_name]
                        sample_name = sample_name.replace('.', '_')
                        sample_name = sample_name.replace('-', '_')

                        trees[sys_term] = rfile.Get(sample_name)
                        tables[sys_term] = getattr(h5file.root, sample_name)
                        weighted_events[sys_term] = getattr(rfile,
                                sample_name + events_hist_suffix)[events_bin]

                        unused_terms.remove(sys_term)

                if unused_terms:
                    #print "UNUSED TERMS for %s:" % self.name
                    #print unused_terms

                    for term in unused_terms:
                        trees[term] = None # flag to use NOMINAL
                        tables[term] = None
                        weighted_events[term] = None # flag to use NOMINAL

            if isinstance(self, Higgs):
                # use yellowhiggs for cross sections
                xs = yellowhiggs.xsbr(
                        7, self.masses[i],
                        self.modes[i], 'tautau')[0] * TAUTAUHADHADBR
                kfact = 1.
                effic = 1.
            elif isinstance(self, Embedded_Ztautau):
                xs, kfact, effic = 100., 1., 1.
            else:
                xs, kfact, effic = ds.xsec_kfact_effic
            if VERBOSE:
                print ds.name, xs, kfact, effic,
            self.datasets.append((ds, trees, tables, weighted_events, xs, kfact, effic))

    @property
    def label(self):

        l = self._label
        if self.scale != 1. and not isinstance(self,
                (MC_Ztautau, Embedded_Ztautau)):
            l += r' ($\sigma_{SM} \times %g$)' % self.scale
        return l

    def draw(self, expr, category, region, bins, min, max, cuts=None):

        hist = Hist(bins, min, max, title=self.label, **self.hist_decor)
        self.draw_into(hist, expr, category, region, cuts=cuts)
        return hist

    def draw_into(self, hist, expr, category, region, cuts=None):

        if isinstance(expr, (list, tuple)):
            exprs = expr
        else:
            exprs = (expr,)

        if hasattr(hist, 'systematics'):
            sys_hists = hist.systematics
        else:
            sys_hists = {}

        selection = self.cuts(category, region) & cuts

        for ds, sys_trees, sys_tables, sys_events, xs, kfact, effic in self.datasets:

            nominal_tree = sys_trees['NOMINAL']
            nominal_events = sys_events['NOMINAL']


            nominal_weight = (TOTAL_LUMI * self.scale *
                    xs * kfact * effic / nominal_events)

            #print nominal_tree.GetEntries(selection), nominal_weight, nominal_tree.GetWeight()

            nominal_weighted_selection = (
                '%f * %s * (%s)' %
                (nominal_weight,
                 ' * '.join(self.get_weight_branches('NOMINAL')),
                 selection))

            if VERBOSE:
                print nominal_weighted_selection

            current_hist = hist.Clone()
            current_hist.Reset()

            # fill nominal histogram
            for expr in exprs:
                nominal_tree.Draw(expr, nominal_weighted_selection,
                        hist=current_hist)

            hist += current_hist

            if not self.systematics:
                continue

            # iterate over systematic variation trees
            for sys_term in sys_trees.iterkeys():

                # skip the nominal tree
                if sys_term == 'NOMINAL':
                    continue

                sys_hist = current_hist.Clone()

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

                sys_hist = current_hist.Clone()
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

            # QCD + Ztautau fit error
            if isinstance(self, Ztautau):
                up_fit = current_hist.Clone()
                up_fit *= ((self.scale + self.scale_error) / self.scale)
                down_fit = current_hist.Clone()
                down_fit *= ((self.scale - self.scale_error) / self.scale)
                if ('ZFIT_UP',) not in sys_hists:
                    sys_hists[('ZFIT_UP',)] = up_fit
                    sys_hists[('ZFIT_DOWN',)] = down_fit
                else:
                    sys_hists[('ZFIT_UP',)] += up_fit
                    sys_hists[('ZFIT_DOWN',)] += down_fit
            else:
                for _term in [('ZFIT_UP',), ('ZFIT_DOWN',)]:
                    if _term not in sys_hists:
                        sys_hists[_term] = current_hist.Clone()
                    else:
                        sys_hists[_term] += current_hist.Clone()

            for _term in [('QCDFIT_UP',), ('QCDFIT_DOWN',)]:
                if _term not in sys_hists:
                    sys_hists[_term] = current_hist.Clone()
                else:
                    sys_hists[_term] += current_hist.Clone()

        # set the systematics
        hist.systematics = sys_hists

    def scores(self, clf, branches, train_fraction,
            category, region, cuts=None, scores_dict=None):

        if scores_dict is None:
            scores_dict = {}

        for systematic in systematics.iter_systematics(True):

            if not self.systematics and systematic != 'NOMINAL':
                continue

            arr_train, arr_test = self.train_test(
                    branches=branches,
                    category=category,
                    region=region,
                    train_fraction=train_fraction,
                    cuts=cuts,
                    systematic=systematic)

            sample = np.vstack(arr_test[branch] for branch in branches).T
            scores = clf.predict_proba(sample)[:,-1]
            weights = arr_test['weight']

            if systematic not in scores_dict:
                scores_dict[systematic] = (scores, weights)
            else:
                prev_scores, prev_weights = scores_dict[systematic]
                scores_dict[systematic] = (
                        np.concatenate((prev_scores, scores)),
                        np.concatenate((prev_weights, weights)))
        return scores_dict

    def trees(self, category, region, cuts=None, systematic='NOMINAL'):
        """
        This is where all the magic happens...
        """
        TEMPFILE.cd()
        selection = self.cuts(category, region) & cuts
        weight_branches = self.get_weight_branches(systematic)
        if systematic in MC.SYSTEMATICS_BY_WEIGHT:
            systematic = 'NOMINAL'

        trees = []
        for ds, sys_trees, sys_tables, sys_events, xs, kfact, effic in self.datasets:

            if systematic in (('ZFIT_UP',), ('ZFIT_DOWN',),
                              ('QCDFIT_UP',), ('QCDFIT_DOWN',)):
                tree = sys_trees['NOMINAL']
                events = sys_events['NOMINAL']
            else:
                tree = sys_trees[systematic]
                events = sys_events[systematic]

                if tree is None:
                    tree = sys_trees['NOMINAL']
                    events = sys_events['NOMINAL']

            scale = self.scale
            if isinstance(self, Ztautau):
                if systematic == ('ZFIT_UP',):
                    scale = self.scale + self.scale_error
                elif systematic == ('ZFIT_DOWN',):
                    scale = self.scale - self.scale_error
            weight = scale * TOTAL_LUMI * xs * kfact * effic / events

            selected_tree = asrootpy(tree.CopyTree(selection))
            #print selected_tree.GetEntries(), weight
            selected_tree.SetWeight(weight)
            selected_tree.userdata.weight_branches = weight_branches
            #print self.name, selected_tree.GetEntries(), selected_tree.GetWeight()
            trees.append(selected_tree)
        return trees

    def tables(self, branches, category, region, cuts=None,
            include_weight=True, systematic='NOMINAL'):
        """
        This is where all the magic happens...
        """
        TEMPFILE.cd()
        selection = self.cuts(category, region) & cuts
        weight_branches = self.get_weight_branches(systematic)
        if systematic in MC.SYSTEMATICS_BY_WEIGHT:
            systematic = 'NOMINAL'
        if include_weight:
            branches = branches + ['weight']
        tables = []
        for ds, sys_trees, sys_tables, sys_events, xs, kfact, effic in self.datasets:

            if systematic in (('ZFIT_UP',), ('ZFIT_DOWN',),
                              ('QCDFIT_UP',), ('QCDFIT_DOWN',)):
                table = sys_tables['NOMINAL']
                events = sys_events['NOMINAL']
            else:
                table = sys_tables[systematic]
                events = sys_events[systematic]

                if table is None:
                    table = sys_tables['NOMINAL']
                    events = sys_events['NOMINAL']

            scale = self.scale
            if isinstance(self, Ztautau):
                if systematic == ('ZFIT_UP',):
                    scale = self.scale + self.scale_error
                elif systematic == ('ZFIT_DOWN',):
                    scale = self.scale - self.scale_error
            weight = scale * TOTAL_LUMI * xs * kfact * effic / events

            # read the table with a selection
            table = table.readWhere(selection.where())
            # add weight field
            if include_weight:
                weights = np.ones(table.shape[0], dtype='f4') * weight
                table = recfunctions.rec_append_fields(table,
                        names='weight',
                        data=weights,
                        dtypes='f4')
                # merge the weight columns
                table['weight'] *= reduce(np.multiply,
                        [table[br] for br in weight_branches])
            # drop all branches except branches (+ weight)
            table = table[branches]
            tables.append(table)
        return tables

    def events(self, selection='', systematic='NOMINAL'):

        total = 0.
        for ds, sys_trees, sys_tables, sys_events, xs, kfact, effic in self.datasets:
            tree = sys_trees[systematic]
            events = sys_events[systematic]

            weight = TOTAL_LUMI * self.scale * xs * kfact * effic / events

            total += weight * tree.GetEntries(selection)
        return total

    def iter(self, selection='', systematic='NOMINAL'):

        TEMPFILE.cd()
        for ds, sys_trees, sys_tables, sys_events, xs, kfact, effic in self.datasets:
            tree = sys_trees[systematic]
            events = sys_events[systematic]

            weight = TOTAL_LUMI * self.scale * xs * kfact * effic / events

            if selection:
                selected_tree = asrootpy(tree.CopyTree(selection))
            else:
                selected_tree = tree
            for event in selected_tree:
                yield weight, event


class Ztautau:
    pass


class MC_Ztautau(MC, Ztautau, Background):

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
        self.scale_error = 0.
        super(MC_Ztautau, self).__init__(
                color=color,
                systematics_terms=systematics_terms,
                **kwargs)


class Embedded_Ztautau(MC, Ztautau, Background):

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
        self.scale_error = 0.
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

    def __init__(self,
            mode=None, modes=None,
            mass=None, masses=None, **kwargs):

        if masses is None:
            if mass is not None:
                assert mass in Higgs.MASS_POINTS
                masses = [mass]
            else:
                masses = Higgs.MASS_POINTS
        else:
            assert len(masses) > 0
            for mass in masses:
                assert mass in Higgs.MASS_POINTS
            assert len(set(masses)) == len(masses)

        if modes is None:
            if mode is not None:
                assert mode in Higgs.MODES
                modes = [mode]
            else:
                modes = Higgs.MODES.keys()
        else:
            assert len(modes) > 0
            for mode in modes:
                assert mode in Higgs.MODES
            assert len(set(modes)) == len(modes)

        str_mass = ''
        if len(masses) == 1:
            str_mass = '(%s)' % str(masses[0])

        str_mode = ''
        if len(modes) == 1:
            str_mode = str(modes[0]) + ' '

        self._label = r'%s$H%s\rightarrow\tau_{\mathrm{had}}\tau_{\mathrm{had}}$' % (
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
        self.scale_error = 0.
        self.shape_region = shape_region

    def draw(self, expr, category, region, bins, min, max, cuts=None):

        hist = Hist(bins, min, max, title=self.label, **self.hist_decor)
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
            if not hasattr(hist, 'systematics'):
                hist.systematics = {}
            for sys_term, sys_hist in MC_bkg_notOS.systematics.items():
                scale = self.scale
                if sys_term == ('FIT_UP',):
                    scale = self.scale + self.scale_error
                elif sys_term == ('FIT_DOWN',):
                    scale = self.scale - self.scale_error
                qcd_hist = (data_hist - sys_hist) * scale
                if sys_term not in hist.systematics:
                    hist.systematics[sys_term] = qcd_hist
                else:
                    hist.systematics[sys_term] += qcd_hist

        hist.SetTitle(self.label)

    def scores(self,
               clf,
               branches,
               train_fraction,
               category,
               region,
               cuts=None):

        # SS data
        data_scores, data_weights = self.data.scores(
                clf,
                branches,
                train_fraction,
                category,
                region=self.shape_region,
                cuts=cuts)

        scores_dict = {}
        # subtract SS MC
        for mc in self.mc:
            mc.scores(
                    clf,
                    branches,
                    train_fraction,
                    category,
                    region=self.shape_region,
                    cuts=cuts,
                    scores_dict=scores_dict)

        for sys_term in scores_dict.keys():
            sys_scores, sys_weights = scores_dict[sys_term]
            scale = self.scale
            if sys_term == ('QCDFIT_UP',):
                scale += self.scale_error
            elif sys_term == ('QCDFIT_DOWN',):
                scale -= self.scale_error
            # subtract SS MC
            sys_weights *= -1 * scale
            # add SS data
            sys_scores = np.concatenate((sys_scores, np.copy(data_scores)))
            sys_weights = np.concatenate((sys_weights, data_weights * scale))
            scores_dict[sys_term] = (sys_scores, sys_weights)
        return scores_dict

    def trees(self, category, region, cuts=None,
              systematic='NOMINAL'):

        TEMPFILE.cd()
        data_tree = asrootpy(
                self.data.data.CopyTree(
                    self.data.cuts(
                        category,
                        region=self.shape_region) & cuts))
        data_tree.userdata.weight_branches = []
        trees = [data_tree]
        for mc in self.mc:
            _trees = mc.trees(
                    category,
                    region=self.shape_region,
                    cuts=cuts,
                    systematic=systematic)
            for tree in _trees:
                tree.Scale(-1)
            trees += _trees

        scale = self.scale
        if systematic == ('QCDFIT_UP',):
            scale += self.scale_error
        elif systematic == ('QCDFIT_DOWN',):
            scale -= self.scale_error

        for tree in trees:
            tree.Scale(scale)
        return trees

    def tables(self, branches, category, region, cuts=None,
            include_weight=True, systematic='NOMINAL'):

        assert include_weight == True

        data_tables = self.data.tables(
                branches=branches,
                category=category,
                region=self.shape_region,
                cuts=cuts,
                include_weight=include_weight,
                systematic='NOMINAL')
        arrays = data_tables

        for mc in self.mc:
            _arrays = mc.tables(
                    branches=branches,
                    category=category,
                    region=self.shape_region,
                    cuts=cuts,
                    include_weight=include_weight,
                    systematic=systematic)
            for array in _arrays:
                array['weight'] *= -1
            arrays.extend(_arrays)

        scale = self.scale
        if systematic == ('QCDFIT_UP',):
            scale += self.scale_error
        elif systematic == ('QCDFIT_DOWN',):
            scale -= self.scale_error

        for array in arrays:
            array['weight'] *= scale
        return arrays


class MC_TauID(MC):

    def __init__(self, **kwargs):

        self.name = 'TauID'
        self._label = 'TauID'
        self.samples = ['PythiaWtaunu_incl.mc11c']
        super(MC_TauID, self).__init__(student='TauIDProcessor',
                db=DB_TAUID, **kwargs)


if __name__ == '__main__':

    from background_estimation import qcd_ztautau_norm

    # tests
    category='boosted'
    shape_region = '!OS'
    target_region = 'OS'

    ztautau = MC_Ztautau(systematics=False)
    others = Others(systematics=False)
    data = Data()
    qcd = QCD(data=data, mc=[others, ztautau],
          shape_region=shape_region)

    qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error = qcd_ztautau_norm(
        ztautau=ztautau,
        backgrounds=[others],
        data=data,
        category=category,
        target_region=target_region,
        qcd_shape_region=shape_region,
        use_cache=True)

    qcd.scale = qcd_scale
    qcd.scale_error = qcd_scale_error
    ztautau.scale = ztautau_scale
    ztautau.scale_error = ztautau_scale_error

    expr = 'tau1_BDTJetScore'
    cuts = None
    bins = 20
    min, max = -0.5, 1.5

    print '--- Others'
    other_hist = others.draw(
            expr,
            category, target_region,
            bins, min, max,
            cuts=cuts)
    print "--- QCD"
    qcd_hist = qcd.draw(
            expr,
            category, target_region,
            bins, min, max,
            cuts=cuts)
    print '--- Z'
    ztautau_hist = ztautau.draw(
            expr,
            category, target_region,
            bins, min, max,
            cuts=cuts)
    print '--- Data'
    data_hist = data.draw(
            expr,
            category, target_region,
            bins, min, max,
            cuts=cuts)

    print "Data: %f" % sum(data_hist)
    print "QCD: %f" % sum(qcd_hist)
    print "Z: %f" % sum(ztautau_hist)
    print "Others: %f" % sum(other_hist)
    print "Data / Model: %f" % (sum(data_hist) / (sum(qcd_hist) +
        sum(ztautau_hist) + sum(other_hist)))

    # test scores
    from categories import CATEGORIES
    import pickle
    import numpy as np

    branches = CATEGORIES[category]['features']

    train_frac = .5

    with open('clf_%s.pickle' % category, 'r') as f:
        clf = pickle.load(f)
        print clf
    print '--- Others'
    other_scores, other_weights = others.scores(
            clf, branches,
            train_frac,
            category, target_region,
            cuts=cuts)['NOMINAL']
    print '--- QCD'
    qcd_scores, qcd_weights = qcd.scores(
            clf, branches,
            train_frac,
            category, target_region,
            cuts=cuts)['NOMINAL']
    print '--- Z'
    ztautau_scores, ztautau_weights = ztautau.scores(
            clf, branches,
            train_frac,
            category, target_region,
            cuts=cuts)['NOMINAL']
    print '--- Data'
    data_scores, data_weights = data.scores(
            clf, branches,
            train_frac,
            category, target_region,
            cuts=cuts)

    print "Data: %d" % (len(data_scores))
    print "QCD: %f" % np.sum(qcd_weights)
    print "Z: %f" % np.sum(ztautau_weights)
    print "Others: %f" % np.sum(other_weights)
