# --
# November 6th, 2014: Not converted to XAOD yet
# --
from rootpy.tree.filtering import EventFilter
from . import log; log = log[__name__]

import ROOT

__all__ = [
    'EmbeddingPileupPatch',
    'EmbeddingIsolation',
    'EmbeddingCorrections'
]

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
            event.averageIntPerXing = averageIntPerXing
        else:
            log.warning("pdgID 39 not found! Skipping event...")
            # ignore event
            return None
        return True


class EmbeddingIsolation(EventFilter):
    """
    2012 embedding isolation systematics
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/EmbeddingTools
    """
    def __init__(self, tree, **kwargs):
        self.tree = tree
        super(EmbeddingIsolation, self).__init__(**kwargs)

    def passes(self, event):
        # no isolation
        isolation = 0
        found = False
        # find particle with pdgid = 82
        for p in event.mc:
            if p.pdgId == 82:
                found = True
                if p.fourvect.Px() > 0:
                    # default isolation
                    isolation = 1
                    if p.fourvect.Py() > 0:
                        # tight isolation
                        isolation = 2
                break
        if not found:
            log.warning("pdgID 82 not found! Skipping event...")
            # ignore event
            return None
        self.tree.embedding_isolation = isolation
        return True


class EmbeddingCorrections(EventFilter):
    # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauEmbeddedCorrections2013
    def __init__(self, tree, year, passthrough=False, **kwargs):
        super(EmbeddingCorrections, self).__init__(passthrough=passthrough, **kwargs)
        self.year = year
        self.tree = tree
        if not passthrough:
            from externaltools import EmbeddedCorrections
            from externaltools import TrigMuonEfficiency
            from externaltools import ElectronEfficiencyCorrection
            from externaltools import HSG4LepLepTriggerSF
            from externaltools import MuonEfficiencyCorrections
            if year == 2011:
                self.tool = ROOT.EmbeddedCorrections.Embedded7TeV(
                    EmbeddedCorrections.get_resource('2DMaps_7TeV.root'),
                    TrigMuonEfficiency.RESOURCE_PATH,
                    ElectronEfficiencyCorrection.RESOURCE_PATH,
                    'muon_trigger_sf_mc11c.root',
                    'rel17p0.v01',
                    EmbeddedCorrections.get_resource('TriggerEventNumberDimuons7TeV.root'),
                    EmbeddedCorrections.get_resource('trigeff_2011.root'))
            else:
                self.tool = ROOT.EmbeddedCorrections.Embedded(
                    EmbeddedCorrections.get_resource('2DMaps.root'),
                    TrigMuonEfficiency.RESOURCE_PATH,
                    'muon_trigger_sf_2012_AtoL.p1328.root',
                    ElectronEfficiencyCorrection.RESOURCE_PATH,
                    HSG4LepLepTriggerSF.RESOURCE_PATH,
                    MuonEfficiencyCorrections.RESOURCE_PATH,
                    'STACO_CB_plus_ST_2012_SF.txt.gz')

    def passes(self, event):
        if self.year == 2011:
            self.tool.SetupEmbeddedEvent(
                event.mc_pt,
                event.mc_eta,
                event.mc_phi,
                event.mc_m,
                event.mc_pdgId,
                # possibly a random run number (for MC embedding)
                self.tree.RunNumber,
                event.EventNumber,
                self.tree.lbn)
        else:
            self.tool.SetupEmbeddedEvent(
                event.mc_pt,
                event.mc_eta,
                event.mc_phi,
                event.mc_m,
                event.mc_pdgId,
                # possibly a random run number (for MC embedding)
                self.tree.RunNumber)
        # Retrieve the unfolding weight
        self.tree.embedding_reco_unfold = self.tool.GetEmbeddingRecoUnfolding()
        # Access the trigger unfolding weight
        self.tree.embedding_trigger_weight = self.tool.GetEmbeddingTriggerWeight()
        # Get the original mass of the dimuon event
        self.tree.embedding_dimuon_mass = self.tool.GetOriginalZ().M()
        return True
