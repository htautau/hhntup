#!/usr/bin/env python

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

"""
This module generates a database of all MC and data datasets
"""
USE_PYAMI = True
try:
    from pyAMI.client import AMIClient
    from pyAMI.query import get_dataset_xsec_min_max_effic, \
                            get_dataset_info, \
                            get_provenance, \
                            get_periods, \
                            get_runs
    from pyAMI import query
    from pyAMI.auth import AMI_CONFIG, create_auth_config
except ImportError:
    USE_PYAMI = False
    print "Warning: pyAMI is not installed."
    print "Cross-section retrieval will be disabled."


import sys
from rootpy.io import open as ropen, DoesNotExist
from rootpy.plotting import Hist
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from operator import itemgetter
import logging

import re
import glob
import os
import cPickle as pickle
import atexit
import fnmatch

from decorators import cached_property
import yaml

from atlastools.datasets import DATA, MC

YEAR = 11
GRL = 'grl/data11_7TeV.periodAllYear_DetStatus-v36-pro10_CoolRunQuery-00-04-08_Higgs_tautau_lh.xml'

"""
LepHad constants
"""
MC_LEPHAD_PATH = '/global/mtm/data/SKIMS/MC'
MC_LEPHAD_PREFIX = 'user.mtm.muLHSkim'
MC_LEPHAD_FILE_PATTERN = '*.root*'

DATA_LEPHAD_PATH = '/global/mtm/data/SKIMS/DATA'
DATA_LEPHAD_PREFIX = 'user.mtm.muLHSkim'
DATA_LEPHAD_MUON_STREAM = 'Muons'
DATA_LEPHAD_EGAMMA_STREAM = 'Egamma'
DATA_LEPHAD_FILE_PATTERN = '*.root*'

"""
HadHad constants
"""
MC_HADHAD_PATH = '/global/endw/mc11_7TeV/higgs_tautau_hh_reskim_p851'
MC_HADHAD_PREFIX = 'user.NoelDawe.HHSkim'
MC_HADHAD_FILE_PATTERN = '*HHSkim*.root*'

DATA_HADHAD_PATH = '/global/endw/data11_7TeV/higgs_tautau_hh_skim_p851'
DATA_HADHAD_PREFIX = 'user.NoelDawe.HTauSkim'
DATA_HADHAD_STREAM = 'JetTauEtmiss'
DATA_HADHAD_FILE_PATTERN = '*HTauSkim*.root*'

"""
Common constants
"""
MC_TREENAME = 'tau'
DATA_TREENAME = 'tau'

DATA_PATTERN = re.compile('^(?P<prefix>\S+.)?data11_7TeV.'
                          '(?P<run>\d+).physics_'
                          '(?P<stream>\S+)?'
                          '.merge.NTUP_TAUMEDIUM.(?P<tag>\w+)'
                          '(.(?P<suffix>\S+))?$')

MC_TAG_PATTERN1 = re.compile('^e(?P<evnt>\d+)_'
                             's(?P<digi>\d+)_'
                             's(?P<digimerge>\d+)_'
                             'r(?P<reco>\d+)_'
                             'r(?P<recomerge>\d+)_'
                             'p(?P<ntup>\d+)$')

# not all valid samples have a recomerge tag:
MC_TAG_PATTERN2 = re.compile('^e(?P<evnt>\d+)_'
                             's(?P<digi>\d+)_'
                             's(?P<digimerge>\d+)_'
                             'r(?P<reco>\d+)_'
                             'p(?P<ntup>\d+)$')

"""
MC11a/b/c categories are defined here
Each MC dataset is automatically classified
acccording to these categories by matching the reco
and merge tags of the dataset name
"""
# order by decreasing preference
MC_CATEGORIES = {
    'mc11a': {'reco':  (2730, 2731),
              'merge': (2780, 2700)},
    'mc11b': {'reco':  (2920, 2923),
              'merge': (3063, 2993, 2900)},
    'mc11c': {'reco':  (3043, 3060),
              'merge': (3109, 3063, 2993)},
}

HERE = os.path.dirname(os.path.abspath(__file__))


"""
Any datasets which don't have the provenance stored properly in AMI
should be hardcoded here (it happens)
"""
DS_NOPROV = {}

"""
Cross-sections are cached so that we don't need to keep asking AMI
for them over and over
"""
XSEC_CACHE_FILE = os.path.join(HERE, 'xsec_cache')
XSEC_CACHE_MODIFIED = False
XSEC_CACHE = {}


if USE_PYAMI:
    amiclient = AMIClient()
    if not os.path.exists(AMI_CONFIG):
        create_auth_config()
    amiclient.read_config(AMI_CONFIG)


class Database(dict):

    def __init__(self, name='datasets'):

        super(Database, self).__init__()
        self.name = name
        self.filepath = os.path.join(HERE, '%s.yml' % self.name)
        if os.path.isfile(self.filepath):
            with open(self.filepath) as db:
                print "Loading '%s' dataset database..." % self.name
                d = yaml.load(db)
                if d:
                    self.update(d)
        self.modified = False

    def write(self):

        if self.modified:
            with open(self.filepath, 'w') as db:
                print "Saving '%s' dataset database to disk..." % self.name
                yaml.dump(dict(self), db)

    def reset(self):

        return self.clear()

    def clear(self):

        # erase all datasets in database
        print "Resetting '%s' dataset database..." % self.name
        super(Database, self).clear()
        self.modified = True

    def validate(self,
                 pattern=None,
                 datatype=MC):

        ds = {}
        for name, info in self.items():
            if info.datatype == datatype:
                if pattern is None or fnmatch.fnmatch(name, pattern):
                    ds[name] = info
        incomplete = []
        for name, info in ds.items():
            print "Validating %s..." % name
            complete = validate_single((name, info), child=False)
            print "Complete: %s" % complete
            print '-'*50
            if not complete:
                incomplete.append(info.ds)
        """
        pool = Pool(processes=cpu_count())
        for result, complete in pool.map(validate_single, sorted(ds.items(), key=itemgetter(0))):
            print result
            print "Complete: %s" % complete
            print '-'*50
            if not complete:
                all_complete = False
        """
        if not incomplete:
            print "ALL DATASETS ARE COMPLETE"
        else:
            print "SOME DATASETS ARE NOT COMPLETE"
            print "INCOMPLETE DATASETS:"
            for ds in incomplete:
                print ds

    def scan(self,
             mc_path=None,
             mc_prefix=None,
             mc_pattern=None,
             data_path=None,
             data_prefix=None,
             data_pattern=None,
             versioned=False,
             deep=False):
        """
        Update the dataset database
        """
        print "Updating '%s' dataset database..." % self.name
        self.modified = True

        if mc_path is not None:
            if versioned:
                pattern1 = ('^mc11_7TeV.(?P<id>\d+)'
                            '.(?P<name>\w+)(.merge.NTUP_TAUMEDIUM)?'
                            '.(?P<tag>e\d+_s\d+_s\d+_r\d+_r\d+_p\d+)'
                            '.v(?P<version>\d+).(?P<suffix>\S+)$')
                pattern2 = ('^mc11_7TeV.(?P<id>\d+)'
                            '.(?P<name>\w+)(.merge.NTUP_TAUMEDIUM)?'
                            '.(?P<tag>e\d+_s\d+_s\d+_r\d+_p\d+)'
                            '.v(?P<version>\d+).(?P<suffix>\S+)$')
            else:
                pattern1 = ('^mc11_7TeV.(?P<id>\d+)'
                            '.(?P<name>\w+)(.merge.NTUP_TAUMEDIUM)?'
                            '.(?P<tag>e\d+_s\d+_s\d+_r\d+_r\d+_p\d+)'
                            '_(?P<suffix>\S+)$')
                pattern2 = ('^mc11_7TeV.(?P<id>\d+)'
                            '.(?P<name>\w+)(.merge.NTUP_TAUMEDIUM)?'
                            '.(?P<tag>e\d+_s\d+_s\d+_r\d+_p\d+)'
                            '_(?P<suffix>\S+)$')

            if mc_prefix:
                pattern1 = ('^%s.' % mc_prefix) + pattern1[1:]
                pattern2 = ('^%s.' % mc_prefix) + pattern2[1:]

            MC_PATTERN1 = re.compile(pattern1)
            MC_PATTERN2 = re.compile(pattern2)

            if deep:
                mc_dirs = get_all_dirs_under(mc_path)
            else:
                if mc_prefix:
                    mc_dirs = glob.glob(os.path.join(mc_path, mc_prefix) + '*')
                else:
                    mc_dirs = glob.glob(os.path.join(mc_path, '*'))

            for dir in mc_dirs:
                dirname, basename = os.path.split(dir)
                match  = re.match(MC_PATTERN1, basename)
                match2 = re.match(MC_PATTERN2, basename)

                if (match2 and not match): match = match2

                if match:
                    name = match.group('name')
                    tag = match.group('tag')
                    try:
                        version = int(match.group('version'))
                    except IndexError:
                        version = 0
                    tag_match = re.match(MC_TAG_PATTERN1, tag)
                    tag_match2 = re.match(MC_TAG_PATTERN2, tag)
                    MC_TAG_PATTERN = MC_TAG_PATTERN1

                    if (tag_match2 and not tag_match) :
                        tag_match = tag_match2
                        MC_TAG_PATTERN = MC_TAG_PATTERN2

                    if not tag_match:
                        print "Dataset not tag-matched: %s" % basename
                        continue
                    cat = None
                    for cat_name, cat_params in MC_CATEGORIES.items():
                        if int(tag_match.group('reco')) in cat_params['reco']:
                            cat = cat_name
                            break
                    if cat is None:
                        print "Dataset does not match a category: %s" % basename
                        continue
                    name += '.' + cat
                    dataset = self.get(name, None)
                    if dataset is not None and version == dataset.version:
                        if tag != dataset.tag:
                            this_reco = int(tag_match.group('reco'))
                            other_reco = int(re.match(dataset.tag_pattern,
                                                      dataset.tag).group('reco'))
                            use_mergetag = True
                            try:
                                this_merge = int(tag_match.group('recomerge'))
                                other_merge = int(re.match(dataset.tag_pattern,
                                                           dataset.tag).group('recomerge'))
                            except IndexError:
                                use_mergetag = False

                            cat_params = MC_CATEGORIES[cat]
                            reco_tags = list(cat_params['reco'])
                            merge_tags = list(cat_params['merge'])
                            assert(this_reco in reco_tags and other_reco in reco_tags)
                            take_this = False
                            if reco_tags.index(this_reco) < reco_tags.index(other_reco):
                                take_this = True
                            elif use_mergetag and this_reco == other_reco and \
                                    merge_tags.index(this_merge) < merge_tags.index(other_merge):
                                take_this = True

                            if take_this:
                                print "taking %s over %s" % (basename, dataset.ds)
                                DATASETS[name] = Dataset(name=name,
                                                datatype=MC,
                                                treename=MC_TREENAME,
                                                ds='mc11_7TeV.' +
                                                match.group('id') + '.' +
                                                match.group('name') +
                                                '.merge.NTUP_TAUMEDIUM.' +
                                                match.group('tag'),
                                                id=int(match.group('id')),
                                                category=cat,
                                                version=version,
                                                tag_pattern=MC_TAG_PATTERN.pattern,
                                                tag=tag,
                                                dirs=[dir],
                                                file_pattern=mc_pattern)
                        else:
                            dataset.dirs.append(dir)
                    elif dataset is None or (dataset is not None and version > dataset.version):
                        self[name] = Dataset(name=name,
                                            datatype=MC,
                                            treename=MC_TREENAME,
                                            ds='mc11_7TeV.' +
                                            match.group('id') + '.' +
                                            match.group('name') +
                                            '.merge.NTUP_TAUMEDIUM.' +
                                            match.group('tag'),
                                            id=int(match.group('id')),
                                            category=cat,
                                            version=version,
                                            tag_pattern=MC_TAG_PATTERN.pattern,
                                            tag=tag,
                                            dirs=[dir],
                                            file_pattern=mc_pattern)
                else:
                    print "Dataset not matched: %s" % basename

        if data_path is not None:
            if data_prefix:
                data_dirs = glob.glob(os.path.join(data_path, data_prefix) + '*')
            else:
                data_dirs = glob.glob(os.path.join(data_path, '*'))

            self['data'] = Dataset(name='data',
                                       datatype=DATA,
                                       treename=DATA_TREENAME,
                                       ds='data',
                                       id=1,
                                       # The GRL is the same for both lephad and hadhad analyses
                                       grl=GRL,
                                       dirs=data_dirs,
                                       file_pattern=data_pattern)
            # TODO create datasets for each run and each period
            runs = {}
            for dir in data_dirs:
                if os.path.isdir(dir):
                    match = re.match(DATA_PATTERN, dir)
                    if match:
                        run = int(match.group('run'))
                        tag = match.group('tag')
                        if run not in runs:
                            runs[run] = {'tag': tag, 'dirs': [dir]}
                        else:
                            runs[run]['dirs'].append(dir)
                            if tag != runs[run]['tag']:
                                print 'multiple copies of run with different tags: %s' % runs[run]['dirs']
                    else:
                        print "this dir does not match valid ds name: %s" % dir
                else:
                    print "this is not a dir: %s" % dir
            for run, info in runs.items():
                name = 'data-%d' % run
                self[name] = Dataset(name=name,
                                         datatype=DATA,
                                         treename=DATA_TREENAME,
                                         ds=name,
                                         id=1,
                                         grl=GRL,
                                         dirs=info['dirs'],
                                         file_pattern=data_pattern)
            if USE_PYAMI:
                run_periods = get_periods(amiclient, year=YEAR, level=2)
                run_periods = [p.name for p in run_periods]
                period_runs = {}
                for period in run_periods:
                    if period == 'VdM':
                        continue
                    _runs = get_runs(amiclient, periods=period, year=YEAR)
                    for run in _runs:
                        period_runs[run] = period
                periods = {}
                for run, info in runs.items():
                    _period = period_runs[run]
                    if _period in periods:
                        periods[_period] += info['dirs']
                    else:
                        periods[_period] = info['dirs'][:]
                for period, dirs in periods.items():
                    name = 'data-%s' % period
                    self[name] = Dataset(name=name,
                                             datatype=DATA,
                                             treename=DATA_TREENAME,
                                             ds=name,
                                             id=1,
                                             grl=GRL,
                                             dirs=dirs,
                                             file_pattern=data_pattern)

    def search(self, pattern):

        data = []
        patterns = pattern
        if not isinstance(pattern, (list, tuple)):
            patterns = [pattern]
        for name, ds in self.items():
            for pattern in patterns:
                if re.match(pattern, name) or fnmatch.fnmatch(name, pattern):
                    data.append(ds)
                    continue
        return data


class Dataset(yaml.YAMLObject):

    yaml_tag = u'!Dataset'

    def __init__(self, name, datatype, treename,
                 ds, dirs,
                 file_pattern='*.root*',
                 id=None,
                 category=None,
                 version=None,
                 tag_pattern=None,
                 tag=None,
                 grl=None):

        self.name = name
        self.datatype = datatype
        self.treename = treename
        self.id = id
        self.ds = ds
        self.category = category
        self.version = version
        self.tag_pattern = tag_pattern
        self.tag = tag
        self.dirs = dirs
        self.file_pattern = file_pattern
        self.grl = grl

    @cached_property
    def xsec_effic(self):

        global XSEC_CACHE_MODIFIED

        if self.datatype == DATA:
            return 1., 1., 1., 1.
        if self.name in XSEC_CACHE:
            return XSEC_CACHE[self.name]
        elif USE_PYAMI:
            if self.ds in DS_NOPROV:
                xsec, xsec_min, xsec_max, effic = get_dataset_xsec_min_max_effic(amiclient, DS_NOPROV[self.ds])
            else:
                xsec, xsec_min, xsec_max, effic = get_dataset_xsec_min_max_effic(amiclient, self.ds)
            #if 'tautauhh' in self.name:
            #    xsec *= 0.412997
            XSEC_CACHE[self.name] = (xsec, xsec_min, xsec_max, effic)
            XSEC_CACHE_MODIFIED = True
            return (xsec, xsec_min, xsec_max, effic)
        else:
            return (None, None, None, None)

    @cached_property
    def xsec_factor(self):

        if 'tautauhh' in self.name:
            return 0.412997
        return 1.

    @cached_property
    def files(self):

        _files = []
        for dir in self.dirs:
            for path, dirs, files in os.walk(dir):
                _files += [os.path.join(path, f) for f in
                           fnmatch.filter(files, self.file_pattern)]
        return _files


def dataset_constructor(loader, node):

    return Dataset(**loader.construct_mapping(node))

yaml.add_constructor(u'!Dataset', dataset_constructor)



if os.path.isfile(XSEC_CACHE_FILE):
    with open(XSEC_CACHE_FILE) as cache:
        print "Loading cross-section cache..."
        XSEC_CACHE = pickle.load(cache)


@atexit.register
def write_cache():

    if XSEC_CACHE_MODIFIED:
        with open(XSEC_CACHE_FILE, 'w') as cache:
            print "Saving cross-section cache to disk..."
            pickle.dump(XSEC_CACHE, cache)


def validate_single(args, child=True):

    if child:
        from cStringIO import StringIO

        sys.stdout = out = StringIO()
        sys.stderr = out

    name = args[0]
    info = args[1]
    complete = True

    try:
        dirs = info.dirs
        root_files = []
        for dir in dirs:
            root_files += glob.glob(os.path.join(dir, '*.root*'))
        events = 0
        for fname in root_files:
            try:
                with ropen(fname) as rfile:
                    tree = rfile.tau
                    try: # skimmed dataset
                        cutflow = rfile.cutflow
                        events += int(cutflow[0])
                    except DoesNotExist: # unskimmed dataset
                        events += tree.GetEntries()
            except IOError:
                print "Currupt file: %s" % fname
                pass
        # determine events in original ntuples
        # use first dir
        ds_name = info.ds
        print 'NTUP: ' + ds_name
        ds_info = get_dataset_info(amiclient, ds_name)
        ntuple_events = int(ds_info.info['totalEvents'])
        try:
            # determine events in AODs
            prov = get_provenance(amiclient, ds_name, type='AOD')
            AOD_ds = prov.values()[0][0]
            print 'AOD: ' + AOD_ds
            AOD_ds_info = get_dataset_info(amiclient, AOD_ds)
            AOD_events = int(AOD_ds_info.info['totalEvents'])
        except IndexError:
            print 'AOD: UNKNOWN'
            AOD_events = ntuple_events
        print name
        print "\tevts\tNTUP\tAOD"
        print "\t%i\t%i\t%i" % (events, ntuple_events, AOD_events)
        if events != ntuple_events:
            print "NTUP MISMATCH"
        if events != AOD_events:
            print "AOD MISMATCH"
        if events != ntuple_events and (events != AOD_events or AOD_events == 0):
            print "MISSING EVENTS"
            complete = False
        if child:
            return out.getvalue(), complete
        return complete
    except Exception, e:
        import traceback
        print "dataset %s exception" % name
        traceback.print_exception(*sys.exc_info())
        if child:
            return out.getvalue(), False
        return False


def get_all_dirs_under(path):
    """
    Get list of all directories under path
    """
    dirs = []
    for dirpath, dirnames, filenames in os.walk(path):

        _dirnames = []
        for dirname in dirnames:
            fullpath = os.path.join(dirpath, dirname)
            # check if this dir contains other dirs
            subdirs_exist = False
            subdirs = os.listdir(fullpath)
            for subdir in subdirs:
                if os.path.isdir(os.path.join(fullpath, subdir)):
                    subdirs_exist = True
                    break
            if subdirs_exist:
                _dirnames.append(dirname)
            else:
                # this must be a dataset, don't walk into this dir
                dirs.append(fullpath)
        # only recurse on directories containing subdirectories
        dirnames = _dirnames

    return dirs


if __name__ == '__main__':

    """
    Update the database
    """
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--deep', action='store_true', default=False)
    parser.add_argument('--reset', action='store_true', default=False)
    parser.add_argument('--versioned', action='store_true', default=False)
    parser.add_argument('--validate', action='store_true', default=False)
    parser.add_argument('--validate-pattern', default=None)
    parser.add_argument('--info', action='store_true', default=False)

    parser.add_argument('--year', type=int, default=None)
    parser.add_argument('--grl', default=None)

    parser.add_argument('--mc-path', default=None)
    parser.add_argument('--mc-prefix', default=None)
    parser.add_argument('--mc-pattern', default='*.root*')

    parser.add_argument('--data-path', default=None)
    parser.add_argument('--data-prefix', default=None)
    parser.add_argument('--data-pattern', default='*.root*')

    parser.add_argument('--name', default='datasets')

    parser.add_argument('analysis',
                        choices=('mulh', 'elh', 'hh', 'custom'),
                        default='custom', nargs='?')
    args = parser.parse_args()

    if args.year is not None:
        YEAR = args.year
    if args.grl is not None:
        GRL = args.grl

    if args.analysis == 'hh':
        mc_path = MC_HADHAD_PATH
        mc_prefix = MC_HADHAD_PREFIX
        mc_pattern = MC_HADHAD_FILE_PATTERN
        data_path = DATA_HADHAD_PATH
        data_prefix = DATA_HADHAD_PREFIX
        data_pattern = DATA_HADHAD_FILE_PATTERN
        args.versioned = True
        args.name = 'datasets_hh'
    elif args.analysis == 'mulh':
        mc_path = MC_LEPHAD_PATH
        mc_prefix = MC_LEPHAD_PREFIX
        mc_pattern = MC_LEPHAD_FILE_PATTERN
        data_path = DATA_LEPHAD_PATH
        data_prefix = DATA_LEPHAD_PREFIX
        data_pattern = DATA_LEPHAD_FILE_PATTERN
        args.versioned = True
        args.name = 'datasets_mulh'
    elif args.analysis == 'elh':
        mc_path = MC_LEPHAD_PATH
        mc_prefix = MC_LEPHAD_PREFIX
        mc_pattern = MC_LEPHAD_FILE_PATTERN
        data_path = DATA_LEPHAD_PATH
        data_prefix = DATA_LEPHAD_PREFIX
        data_pattern = DATA_LEPHAD_FILE_PATTERN
        args.versioned = True
        args.name = 'datasets_elh'
    else: # custom
        mc_path = args.mc_path
        mc_prefix = args.mc_prefix
        mc_pattern = args.mc_pattern
        data_path = args.data_path
        data_prefix = args.data_prefix
        data_pattern = args.data_pattern

    db = Database(name=args.name)

    if args.reset:
        db.clear()

    if mc_path is not None or data_path is not None:
        db.scan(mc_path=mc_path,
                mc_prefix=mc_prefix,
                mc_pattern=mc_pattern,
                data_path=data_path,
                data_prefix=data_prefix,
                data_pattern=data_pattern,
                deep=args.deep,
                versioned=args.versioned)

    if args.validate or args.validate_pattern is not None:
        # check for missing events etc...
        db.validate(pattern=args.validate_pattern,
                    datatype=MC)
    else:
        if args.info:
            print "%i datasets in database" % len(db)
        else:
            for name in sorted(db.keys()):
                print "%s => %s" % (name, db[name].ds)
    db.write()
