from rootpy.tree.filtering import EventFilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
from .models import TrueTauBlock
from . import track_counting


class TauLeadSublead(EventFilter):

    def __init__(self, lead=35*GeV, sublead=25*GeV, **kwargs):

        super(TauLeadSublead, self).__init__(**kwargs)
        """
        Leading and subleading tau pT thresholds
        """
        self.lead = lead
        self.sublead = sublead

    def passes(self, event):
        # sort in descending order by pT
        event.taus.sort(key=lambda tau: tau.pt, reverse=True)
        # only keep leading two taus
        event.taus.slice(0, 2)
        # Event passes if the highest pT tau is above the leading
        # pT threshold and the next subleading tau pT is above the subleading pT theshold
        return event.taus[0].pt > self.lead and event.taus[1].pt > self.sublead


class Triggers(EventFilter):
    """
    See lowest unprescaled triggers here:
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled#Taus_electron_muon_MET
    """
    def __init__(self, year, old_skim=False, **kwargs):

        if year == 2011:
            if old_skim:
                self.passes = self.passes_11_old
            else:
                self.passes = self.passes_11
        elif year == 2012:
            self.passes = self.passes_12
        else:
            raise ValueError("No triggers defined for year %d" % year)
        super(Triggers, self).__init__(**kwargs)

    def passes_11(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)

    def passes_11_old(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1_EMULATED
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1_EMULATED
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)

    def passes_12(self, event):
        try:
            return event.EF_tau29Ti_medium1_tau20Ti_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e

        # TODO use tau27Ti_m1_tau18Ti_m1_L2loose for period E
        # need emulaion, SFs for this


class ElectronVeto(EventFilter):

    def passes(self, event):

        for el in event.electrons:
            pt = el.cl_E / cosh(el.tracketa)
            if pt <= 15 * GeV: continue
            if not ((abs(el.tracketa) < 1.37) or (1.52 < abs(el.tracketa) < 2.47)): continue
            if el.author not in (1, 3): continue
            if not abs(el.charge) == 1: continue
            if el.mediumPP != 1: continue
            if (el.OQ & 1446) != 0: continue
            return False
        return True


from ..filters import muon_has_good_track

class MuonVeto(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(MuonVeto, self).__init__(**kwargs)

    def passes(self, event):

       for muon in event.muons:
           if muon.pt <= 10 * GeV:
               continue
           if abs(muon.eta) >= 2.5:
               continue
           if muon.loose != 1:
               continue
           if not muon_has_good_track(muon, self.year):
               continue
           return False
       return True


class TaudR(EventFilter):

    def __init__(self, dr=3.2, **kwargs):

        super(TaudR, self).__init__(**kwargs)
        self.dr = dr

    def passes(self, event):

        assert len(event.taus) == 2
        tau1, tau2 = event.taus
        return utils.dR(tau1.eta, tau1.phi, tau2.eta, tau2.phi) < self.dr


class TruthMatching(EventFilter):

    def __init__(self, tree, **kwargs):

        self.tree = tree
        super(TruthMatching, self).__init__(**kwargs)

    def passes(self, event):

        # Truth-matching
        # match only with visible true taus
        event.truetaus.select(lambda tau: tau.vis_Et > 10 * GeV and abs(tau.vis_eta) < 2.5)

        if len(event.truetaus) > 2:
            print "ERROR: too many true taus: %i" % len(event.truetaus)
            for truetau in event.truetaus:
                print "truth (pT: %.4f, eta: %.4f, phi: %.4f)" % (truetau.pt, truetau.eta, truetau.phi),
                if truetau.tauAssoc_index >= 0:
                    matched_tau = event.taus.getitem(truetau.tauAssoc_index)
                    print " ==> reco (pT: %.4f, eta: %.4f, phi: %.4f)" % (matched_tau.pt, matched_tau.eta, matched_tau.phi),
                    print "dR = %.4f" % truetau.tauAssoc_dr
                else:
                    print ""
            tree.error = True

        tau1, tau2 = event.taus
        unmatched_reco = range(2)
        unmatched_truth = range(event.truetaus.len())
        matched_truth = []
        for i, tau in enumerate((tau1, tau2)):
            matching_truth_index = tau.trueTauAssoc_index
            if matching_truth_index >= 0:
                unmatched_reco.remove(i)
                # check that this tau / true tau was not previously matched
                if matching_truth_index not in unmatched_truth or \
                   matching_truth_index in matched_truth:
                    print "ERROR: match collision!"
                    tau1.matched_collision = True
                    tau2.matched_collision = True
                    self.tree.trueTau1_matched_collision = True
                    self.tree.trueTau2_matched_collision = True
                    self.tree.error = True
                else:
                    unmatched_truth.remove(matching_truth_index)
                    matched_truth.append(matching_truth_index)
                    tau.matched = True
                    tau.matched_dR = tau.trueTauAssoc_dr
                    setattr(self.tree, "trueTau%i_matched" % (i+1), 1)
                    setattr(self.tree, "trueTau%i_matched_dR" % (i+1), event.truetaus.getitem(matching_truth_index).tauAssoc_dr)
                    TrueTauBlock.set(self.tree, i+1, event.truetaus.getitem(matching_truth_index))

        for i, j in zip(unmatched_reco, unmatched_truth):
            TrueTauBlock.set(self.tree, i+1, event.truetaus.getitem(j))

        self.tree.mass_vis_true_tau1_tau2 = (self.tree.trueTau1_fourvect_vis + self.tree.trueTau2_fourvect_vis).M()
        return True


class TauTrackRecounting(EventFilter):

    def __init__(self, year, **kwargs):
        self.year = year
        super(TauTrackRecounting, self).__init__(**kwargs)

   def passes(self, event):
        for tau in event.taus:
            tau.ntrack_full = track_counting.count_tracks(
                    tau, event, self.year)
        return True
