#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--grl', default=None)
parser.add_argument('dirs', nargs='+')
args = parser.parse_args()

import re
import os
import sys
from pyAMI import query
from rootpy.io import open as ropen
from rootpy.plotting import Hist
import glob
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from operator import itemgetter
import logging
import ROOT
ROOT.gROOT.SetBatch(True)


#logger=mp.log_to_stderr(logging.DEBUG)


def validate_run(args):

    import sys
    from cStringIO import StringIO

    sys.stdout = out = StringIO()
    sys.stderr = out

    run = args[0]
    info = args[1]
    skim_complete = True

    try:
        dirs = info['dirs']
        root_files = []
        for dir in dirs:
            root_files += glob.glob(os.path.join(dir, '*.root*'))
        events = 0
        for fname in root_files:
            try:
                with ropen(fname) as rfile:
                    tree = rfile.tau
                    cutflow = rfile.cutflow
                    events += int(cutflow[0])
                    #events += tree.GetEntries()
            except IOError:
                print "Currupt file: %s" % fname
                pass
        # determine events in original ntuples
        # use first dir
        ds_name = dirs[0]
        match = re.match(pattern, ds_name)
        ds_name = ds_name[len(match.group('prefix')) + 1 : -1 * (len(match.group('suffix')) + 1)]
        print 'NTUP: ' + ds_name
        ds_info = get_dataset_info(client, ds_name)
        ntuple_events = int(ds_info.info['totalEvents'])
        try:
            # determine events in AODs
            prov = get_provenance(client, ds_name, type='AOD')
            AOD_ds = prov.values()[0][0]
            print 'AOD: ' + AOD_ds
            AOD_ds_info = get_dataset_info(client, AOD_ds)
            AOD_events = int(AOD_ds_info.info['totalEvents'])
        except IndexError:
            print 'AOD: UNKNOWN'
            AOD_events = ntuple_events
        print "run\tevts\tNTUP\tAOD"
        print '%i\t%i\t%i\t%i' % (run, events, ntuple_events, AOD_events)
        if events != ntuple_events:
            print "NTUP MISMATCH"
        if events != AOD_events:
            print "AOD MISMATCH"
        if events != ntuple_events and events != AOD_events:
            print "MISSING EVENTS"
            skim_complete = False
        return out.getvalue(), skim_complete
    except Exception, e:
        import traceback
        print "run %i exception" % run
        traceback.print_exception(*sys.exc_info())
        return out.getvalue(), False

pattern = re.compile('^(?P<prefix>\S+).data11_7TeV.(?P<run>\d+).physics_JetTauEtmiss.merge.NTUP_TAUMEDIUM.(?P<tag>\w+).(?P<suffix>\S+)$')

skim_complete = True

runs = {}
for dir in args.dirs:
    if os.path.isdir(dir):
        match = re.match(pattern, dir)
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

# check that all runs in GRL are present
if args.grl:
    from goodruns import GRL
    grl = GRL(args.grl)
    print "validating against runs in GRL..."
    for run in grl:
        if run not in runs:
            print "WARNING: run %i in GRL but not in skim" % run
            skim_complete = False

from pyAMI.client import AMIClient
from pyAMI.auth import AMI_CONFIG, create_auth_config
from pyAMI.query import get_dataset_info, get_provenance

client = AMIClient()
if not os.path.exists(AMI_CONFIG):
    create_auth_config()
client.read_config(AMI_CONFIG)

pool = Pool(processes=cpu_count())
for result, complete in pool.map(validate_run, sorted(runs.items(), key=itemgetter(0))):
    print result
    print "Complete: %s" % complete
    print '-'*50
    if not complete:
        skim_complete = False

if skim_complete:
    print "SKIM IS COMPLETE"
else:
    print "SKIM IS NOT COMPLETE"
