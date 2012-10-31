from rootpy.tree.filtering import EventFilter
from .models import TrueTauBlock


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
