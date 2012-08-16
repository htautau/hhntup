#!/usr/bin/env python

from samples import MC_Ztautau
from tauid.p851.selection import selection, nvtx_to_category, LEVELS, \
    CATEGORIES, PRONGS
from rootpy.tree import Cut
from rootpy.io import open as ropen

"""
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012
Values are percents
"""
EFFIC_UNCERT = {
    'loose': {
        1: 4.,
        3: 8.,
    },
    'medium': {
        1: 5.,
        3: 8.,
    },
    'tight': {
        1: 4.,
        3: 7.,
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
        cut = (Cut('tau1_numTrack==%d && tau1_matched' % prong) |
               Cut('tau2_numTrack==%d && tau2_matched' % prong)) & category
        for weight, event in sample.iter(cut):
            if event.tau1_numTrack == prong and event.tau1_matched:
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
            if event.tau2_numTrack == prong and event.tau2_matched:
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


    with ropen('bdt_uncertainty.root', 'recreate') as f:
        ztautau = MC_Ztautau(systematics=False)
        for prong in PRONGS:
            for cat_str, category in CATEGORIES.items():

                loose = selection('loose', prong, cat_str)
                medium = selection('medium', prong, cat_str)
                tight = selection('tight', prong, cat_str)

                '''
                # binary search alpha x (medium - loose)
                shift = medium - loose
                uncert = EFFIC_UNCERT['loose'][prong]
                print efficiency(ztautau, loose, prong, category)
                shift.name = 'loose_%dp_%s' % (prong, cat_str)
                shift.Write()
                '''

                # binary search alpha x (tight - medium)
                shift = tight - medium
                uncert = EFFIC_UNCERT['medium'][prong]
                print efficiency(ztautau, medium, prong, category)
                print efficiency(ztautau, medium, prong, category,
                        validate='medium')
                shift.name = 'medium_%dp_%s' % (prong, cat_str)
                shift.Write()

                # binary search alpha x (1. - tight)
                shift = 1. - tight
                uncert = EFFIC_UNCERT['tight'][prong]
                print efficiency(ztautau, tight, prong, category)
                print efficiency(ztautau, tight, prong, category,
                        validate='tight')
                shift.name = 'tight_%dp_%s' % (prong, cat_str)
                shift.Write()

else:

    UNCERT = {}

