from argparse import ArgumentParser

from higgstautau import datasets
from higgstautau.hadhad import periods
from higgstautau.hadhad.periods import total_lumi
from tabulartext import PrettyTable
from rootpy.io import open as ropen
import math
import os


def get_parser():

    parser = ArgumentParser()
    parser.add_argument('--format', choices=('latex', 'text'), default='text')
    parser.add_argument('--short', action='store_true', default=False)
    parser.add_argument('--mass', type=int, default=125)
    parser.add_argument('--proc', default='HHProcessor')
    parser.add_argument('--db')
    parser.add_argument('--noweight', action='store_true', default=False)
    parser.add_argument('--verbose', action='store_true', default=False)
    #parser.add_argument('--dir', default='samples/midpt')
    return parser


def make_cutflow(samples,
                 args,
                 num_format="%.1f"):

    lumi = total_lumi()

    data_num_format = "%d"

    filters = None
    db = datasets.Database(args.db)

    # build table
    cutflow_table = {}
    for i, (latex_name, text_name, sample) in enumerate(samples):
        matched_samples = db.search(sample)
        total_cutflow = None
        total = 0
        for ds in matched_samples:
            with ropen(os.path.join(args.dir, '.'.join([args.proc, ds.name, 'root']))) as rfile:
                cutflow = rfile.cutflow_event
                cutflow.SetDirectory(0)
                ds_total = rfile.cutflow[0]

                if filters is None:
                    filters = [cutflow.GetXaxis().GetBinLabel(j + 1) for j in xrange(len(cutflow))]

                # scale MC by lumi and xsec
                if ds.datatype != datasets.DATA and not args.noweight:
                    events = rfile.cutflow[0]
                    xsec, xsec_min, xsec_max, effic = ds.xsec_effic
                    weight = 1E3 * lumi * xsec / (effic * events)
                    if args.verbose:
                        print '-' * 30
                        print ds.name
                        print "xsec: %f [nb]" % xsec
                        print "effic: %f" % effic
                        print "events: %d" % events
                        print "lumi: %f [1/pb]" % lumi
                        print "weight (1E3 * lumi * xsec / (effic * events)): %f" % weight
                    cutflow *= weight
                    ds_total *= weight

                total += ds_total
                if total_cutflow is None:
                    total_cutflow = cutflow
                else:
                    total_cutflow += cutflow
        cutflow_table[i] = list(zip(total_cutflow, total_cutflow.yerrh()))
        cutflow_table[i].insert(0, (total, math.sqrt(total)))

    filters[0] = 'Skim'
    filters.insert(0, 'Total')

    cutflows = [cutflow_table[i] for i in xrange(len(samples))]

    print
    if args.format == 'text':
        print "Higgs mass of %d GeV" % args.mass
        print "Integrated luminosity of %.3f fb^-1" % (lumi/1000.)
        sample_names = [sample[1] for sample in samples]
        table = PrettyTable(['Filter'] + sample_names)
        if args.short:
            for i, row in enumerate(zip(*cutflows)):
                table.add_row([filters[i]] + [num_format % passing for passing, error in row])
        else:
            for i, row in enumerate(zip(*cutflows)):
                table.add_row([filters[i]] + ["%s(%s)" % (num_format % passing,
                                                          num_format % error)
                                                          for passing, error in row])
        print table.get_string(hrules=1)
    else:
        print "Higgs mass of %d GeV\\\\" % args.mass
        print "Integrated luminosity of %.3f fb$^{-1}$\\\\" % (lumi/1000.)
        sample_names = [sample[0] for sample in samples]
        print r'\begin{center}'
        print r'\begin{scriptsize}'
        print r'\begin{tabular}{%s}' % ('|'.join(['c' for i in xrange(len(samples) + 1)]))
        print " & ".join(['Filter'] + sample_names) + '\\\\'
        print r'\hline\hline'
        if args.short:
            for i, row in enumerate(zip(*cutflows)):
                print " & ".join([filters[i]] + [num_format % passing
                                  for passing, error in row]) + '\\\\'
        else:
            for i, row in enumerate(zip(*cutflows)):
                print " & ".join([filters[i]] + ["$%s\pm%s$" % (num_format % passing,
                                                                num_format % error)
                                  for passing, error in row]) + '\\\\'
        print r'\end{tabular}'
        print r'\end{scriptsize}'
        print r'\end{center}'
    print
