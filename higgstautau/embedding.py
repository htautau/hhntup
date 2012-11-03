from rootpy.tree.filtering import EventFilter


class EmbeddingPileupPatch(EventFilter):

    def passes(self, event):

        # fix averageIntPerXing
        # HACK: stored in pz of mc particle with pdg ID 39
        # https://twiki.cern.ch/twiki/bin/viewauth/Atlas/EmbeddingTools
        averageIntPerXing = None
        for p in event.mc:
            if p.pdgId == 39:
                averageIntPerXing = p.fourvect.Pz()
                break
        if averageIntPerXing is not None:

            #print 'Current value %f' % event.averageIntPerXing
            event.averageIntPerXing = averageIntPerXing
            print 'Proposed value %f' % averageIntPerXing
            print 'Accept value %f' % event.averageIntPerXing

        else:
            print "pdgID 39 not found! Skipping event..."
            return None
        return True
