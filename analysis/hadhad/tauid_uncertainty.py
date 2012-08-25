#!/usr/bin/env python

from rootpy.tree import Cut
from higgstautau.tauid import uncertainty


def efficiency(sample, selection, prong, category):

    category = category.replace(
            'tau_numberOfVertices', 'number_of_good_vertices')

    #total = sample.events(Cut('trueTau1_nProng==%d' % prong) & category)
    total = sample.events(Cut('tau1_numTrack==%d' % prong) & category)

    passing = 0.
    cut = Cut('tau1_numTrack==%d' % prong) & category
    for weight, event in sample.iter(cut):
        if (event.tau1_BDTJetScore >
            selection.Eval(event.tau1_fourvect.Pt())):
            passing += weight
    return passing / total


def efficiency_uncertainty(sample, selection, prong, category):

    category = category.replace(
            'tau_numberOfVertices', 'number_of_good_vertices')

    #total = sample.events(Cut('trueTau1_nProng==%d' % prong) & category)
    total = sample.events(Cut('tau1_numTrack==%d' % prong) & category)

    passing_high = 0.
    passing_low = 0.
    cut = Cut('tau1_numTrack==%d' % prong) & category
    for weight, event in sample.iter(cut):
        pt = event.tau1_fourvect.Pt()
        high_score, low_score = uncertainty(
                event.tau1_BDTJetScore,
                pt,
                event.tau1_numTrack,
                event.number_of_good_vertices)
        selection_cut = selection.Eval(pt)
        if high_score > selection_cut:
            passing_high += weight
        if low_score > selection_cut:
            passing_low += weight
    return passing_high / total, passing_low / total


def efficiency_validation(sample, selection, prong, category):

    category = category.replace(
            'tau_numberOfVertices', 'number_of_good_vertices')

    #total = sample.events(Cut('trueTau1_nProng==%d' % prong) & category)
    total = sample.events(Cut('tau1_numTrack==%d' % prong) & category)

    passing = 0.
    cut = Cut('tau1_numTrack==%d' % prong) & category
    for weight, event in sample.iter(cut):
        if validate == 'loose':
            if event.tau1_JetBDTSigLoose == 1:
                passing += weight
        elif validate == 'medium':
            if event.tau1_JetBDTSigMedium == 1:
                passing += weight
        elif validate == 'tight':
            if event.tau1_JetBDTSigTight == 1:
                passing += weight
        else:
            raise ValueError("invalid working point: %s" % selection)
    return passing / total


if __name__ == '__main__':

    from samples import MC_TauID

    from higgstautau.tauid.p851 import selection, nvtx_to_category, CATEGORIES
    from higgstautau.tauid.common import LEVELS, PRONGS
    from higgstautau.tauid import EFFIC_UNCERT_2011 as EFFIC_UNCERT

    from rootpy.io import open as ropen

    with ropen('bdt_uncertainty.root', 'recreate') as f:
        ztautau = MC_TauID(systematics=False)
        for prong in PRONGS:
            for cat_str, category in CATEGORIES.items():

                def binary_search(target, selection, shift, min_error=0.0001,
                        min_change=0.000001,
                        reverse=False):

                    print "target:", target
                    if reverse:
                        a = 1.
                        b = 0.
                    else:
                        a = 0.
                        b = 1.
                    curr_effic = 1E100
                    prev_effic = None
                    curr_shift = None
                    while abs(curr_effic - target) > min_error:
                        mid = (a + b) / 2.
                        curr_shift = mid * shift
                        curr_effic = efficiency(ztautau,
                                selection + curr_shift,
                                prong, category)
                        print curr_effic
                        if prev_effic is not None and abs(curr_effic - prev_effic) < min_change:
                            break
                        if curr_effic > target:
                            a = mid
                        else:
                            b = mid
                        prev_effic = curr_effic
                    print efficiency(ztautau, selection + curr_shift, prong, category)
                    print list(curr_shift.y())
                    print list((selection + curr_shift).y())
                    print "=" * 20
                    return curr_shift

                loose = selection('loose', prong, cat_str)
                medium = selection('medium', prong, cat_str)
                tight = selection('tight', prong, cat_str)

                shift_loose_low = medium - loose
                shift_loose_high = 0. - loose

                shift_medium_low = tight - medium
                shift_medium_high = loose - medium

                shift_tight_low = 1. - tight
                shift_tight_high = medium - tight

                uncert_loose = EFFIC_UNCERT['loose'][prong]
                uncert_medium = EFFIC_UNCERT['medium'][prong]
                uncert_tight = EFFIC_UNCERT['tight'][prong]

                target_loose = efficiency(ztautau, loose, prong, category)
                target_medium = efficiency(ztautau, medium, prong, category)
                target_tight = efficiency(ztautau, tight, prong, category)

                print target_loose
                print efficiency_validation(ztautau, 'loose', prong, category)

                print target_medium
                print efficiency_validation(ztautau, 'medium', prong, category)

                print target_tight
                print efficiency_validation(ztautau, 'tight', prong, category)

                target_loose_high = target_loose + uncert_loose
                target_loose_low = target_loose - uncert_loose

                target_medium_high = target_medium + uncert_medium
                target_medium_low = target_medium - uncert_medium

                target_tight_high = target_tight + uncert_tight
                target_tight_low = target_tight - uncert_tight

                shift_loose_low = binary_search(target_loose_low, loose,
                        shift_loose_low)
                shift_loose_high = binary_search(target_loose_high, loose,
                        shift_loose_high, reverse=True)

                shift_medium_low = binary_search(target_medium_low, medium,
                        shift_medium_low)
                shift_medium_high = binary_search(target_medium_high, medium,
                        shift_medium_high, reverse=True)

                shift_tight_low = binary_search(target_tight_low, tight,
                        shift_tight_low)
                shift_tight_high = binary_search(target_tight_high, tight,
                        shift_tight_high, reverse=True)

                f.cd()
                shift_loose_high.name = 'loose_high_%dp_%s' % (prong, cat_str)
                shift_loose_high.Write()
                shift_loose_low.name = 'loose_low_%dp_%s' % (prong, cat_str)
                shift_loose_low.Write()

                shift_medium_high.name = 'medium_high_%dp_%s' % (prong, cat_str)
                shift_medium_high.Write()
                shift_medium_low.name = 'medium_low_%dp_%s' % (prong, cat_str)
                shift_medium_low.Write()

                shift_tight_high.name = 'tight_high_%dp_%s' % (prong, cat_str)
                shift_tight_high.Write()
                shift_tight_low.name = 'tight_low_%dp_%s' % (prong, cat_str)
                shift_tight_low.Write()
