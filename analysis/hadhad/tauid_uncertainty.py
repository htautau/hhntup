#!/usr/bin/env python

from samples import MC_Ztautau
from tauid.p851.selection import selection, nvtx_to_category, LEVELS, \
    CATEGORIES, PRONGS
from rootpy.tree import Cut
from rootpy.io import open as ropen

"""
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012
Values are percents / 100
"""
EFFIC_UNCERT = {
    'loose': {
        1: 0.04,
        3: 0.08,
    },
    'medium': {
        1: 0.05,
        3: 0.08,
    },
    'tight': {
        1: 0.04,
        3: 0.07,
    },
}


def uncertainty(score, pt, prong, nvtx):

    loose = selection('loose', prong, nvtx).Eval(pt)
    medium = selection('medium', prong, nvtx).Eval(pt)
    tight = selection('tight', prong, nvtx).Eval(pt)

    if score < loose:
        raise ValueError(
            'No uncertainties defined for scores lower than loose')

    if score < medium:
        return selection_uncertainty('loose', pt, prong, nvtx)
    elif score < tight:
        return selection_uncertainty('medium', pt, prong, nvtx)
    else:
        return selection_uncertainty('tight', pt, prong, nvtx)


def selection_uncertainty(level, pt, prong, nvtx):

    return UNCERT[level][prong][nvtx_to_category(nvtx)].Eval(pt)


if __name__ == '__main__':

    def efficiency(sample, selection, prong, category, validate=False):

        category = category.replace(
                'tau_numberOfVertices', 'number_of_good_vertices')

        total = (sample.events(Cut('trueTau1_nProng==%d' % prong) & category) +
                 sample.events(Cut('trueTau2_nProng==%d' % prong) & category))
        passing = 0.
        cut = (Cut('tau1_numTrack==%d' % prong) |
               Cut('tau2_numTrack==%d' % prong)) & category
        for weight, event in sample.iter(cut):
            if event.tau1_numTrack == prong:
                if validate == 'loose':
                    if event.tau1_JetBDTSigLoose == 1:
                        passing += weight
                elif validate == 'medium':
                    if event.tau1_JetBDTSigMedium == 1:
                        passing += weight
                elif validate == 'tight':
                    if event.tau1_JetBDTSigTight == 1:
                        passing += weight
                elif (event.tau1_BDTJetScore >
                    selection.Eval(event.tau1_fourvect.Pt())):
                    passing += weight
            if event.tau2_numTrack == prong:
                if validate == 'loose':
                    if event.tau2_JetBDTSigLoose == 1:
                        passing += weight
                elif validate == 'medium':
                    if event.tau2_JetBDTSigMedium == 1:
                        passing += weight
                elif validate == 'tight':
                    if event.tau2_JetBDTSigTight == 1:
                        passing += weight
                elif (event.tau2_BDTJetScore >
                    selection.Eval(event.tau2_fourvect.Pt())):
                    passing += weight
        return passing / total

    min_error = 0.005

    with ropen('bdt_uncertainty.root', 'recreate') as f:
        ztautau = MC_Ztautau(systematics=False, student='TauIDProcessor')
        for prong in PRONGS:
            for cat_str, category in CATEGORIES.items():

                loose = selection('loose', prong, cat_str)
                medium = selection('medium', prong, cat_str)
                tight = selection('tight', prong, cat_str)

                # binary search alpha x (medium - loose)
                #shift = medium - loose
                #uncert = EFFIC_UNCERT['loose'][prong]
                #print efficiency(ztautau, loose, prong, category)
                #print efficiency(ztautau, loose, prong, category,
                #        validate='loose')
                #shift.name = 'loose_%dp_%s' % (prong, cat_str)
                #shift.Write()

                # binary search alpha x (tight - medium)
                shift_medium_low = tight - medium
                shift_medium_high = medium - tight

                shift_tight_low = 1. - tight
                shift_tight_high = tight - 1.

                uncert_medium = EFFIC_UNCERT['medium'][prong]
                uncert_tight = EFFIC_UNCERT['tight'][prong]

                target_medium = efficiency(ztautau, medium, prong, category)
                target_tight = efficiency(ztautau, tight, prong, category)

                print target_medium
                print efficiency(ztautau, medium, prong, category,
                        validate='medium')

                target_medium_high = target_medium + uncert_medium
                target_medium_low = target_medium - uncert_medium

                print "target:", target_medium_low
                a = 0.
                b = 1.
                curr_effic = 1E100
                while abs(curr_effic - target_medium_low) > min_error:
                    mid = (a + b) / 2.
                    curr_shift = mid * shift_medium_low
                    curr_effic = efficiency(ztautau,
                            medium + curr_shift,
                            prong, category)
                    print curr_effic
                    if curr_effic > target_medium_low:
                        a = mid
                    else:
                        b = mid
                print efficiency(ztautau, medium + curr_shift, prong, category)

                print "=" * 20

                shift_medium_high.name = 'medium_high_%dp_%s' % (prong, cat_str)
                shift_medium_high.Write()
                shift_medium_low.name = 'medium_low_%dp_%s' % (prong, cat_str)
                shift_medium_low.Write()

                # binary search alpha x (1. - tight)
                print target_tight
                print efficiency(ztautau, tight, prong, category,
                        validate='tight')
                shift_tight_high.name = 'tight_high_%dp_%s' % (prong, cat_str)
                shift_tight_high.Write()
                shift_tight_low.name = 'tight_low_%dp_%s' % (prong, cat_str)
                shift_tight_low.Write()

else:

    UNCERT = {}

