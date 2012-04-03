"""
==================
Skimming procedure
==================

A) For DATA skimming:

    1) Trigger

        The lowest un-prescaled trigger can be found in
        https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled

        - Period A-I  :  tau29_medium1_tau20_medium1 || tau100_medium || xe60_noMu
        - Period J    :  tau29_medium1_tau20_medium1 || tau100_medium || xe60_tight_noMu
        - Period K    :  tau29_medium1_tau20_medium1 || tau125_medium1 || xe60_tight_noMu
        - Period L-M  :  tau29T_medium1_tau20T_medium1 || tau125_medium1 || xe60_verytight_noMu


    2) Two LOOSE taus

        - tau_author!=2 && tau_pT > 18 GeV && tau_numTrack > 0
        - tau_JetBDTLoose==1 || tau_tauLlhLoose==1

        * Note pT>18GeV cut was chosen to be able to fluctuate the TES.

        With these two selection (1 & 2), we expect factor 50 reduction which
        results in ~400-500 GB disk space at 5 fb-1.


B) For MC skimming:

    1) Trigger: OR of above.

==========
Histograms
==========

Some information has to be kept during the skimming.

Just after the "Trigger" requirement, the histograms/trees should be made according
to each trigger decision.  (three histograms for each variable below)

 * number of events within GRL (to check consistency of LumiCalc)
 * mu-distribution
 * # of vertices
 * # of LOOSE taus (pt>18GeV, LLH||BDT-loose)
 * # of EF tau trigger object  (no matching is necessary)
 * tau pT spectrum
 * EF tau pT spectrum  (no matching is necessary)
 * MET
 * anything else???
"""


import ROOT
from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from rootpy.tree.filtering import EventFilter, EventFilterList
from atlastools.batch import ATLASStudent
from rootpy.tree import Tree, TreeChain
from mixins import TauFourMomentum


ROOT.gErrorIgnoreLevel = ROOT.kFatal


class Triggers(EventFilter):

    def passes(self, event):
        """
        - Period A-I  :  tau29_medium1_tau20_medium1 || tau100_medium || xe60_noMu
        - Period J    :  tau29_medium1_tau20_medium1 || tau100_medium || xe60_tight_noMu
        - Period K    :  tau29_medium1_tau20_medium1 || tau125_medium1 || xe60_tight_noMu
        - Period L-M  :  tau29T_medium1_tau20T_medium1 || tau125_medium1 || xe60_verytight_noMu
        """
        if 177531 <= event.RunNumber <= 186493: # Period A-I
            return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau100_medium or event.EF_xe60_noMu
        elif 186516 <= event.RunNumber <= 186755: # Period J
            return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau100_medium or event.EF_xe60_tight_noMu 
        elif 186873 <= event.RunNumber <= 187815: # Period K
            return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau125_medium1 or event.EF_xe60_tight_noMu
        elif 188902 <= event.RunNumber <= 191933: # Period L-M
            return event.EF_tau29T_medium1_tau20T_medium1 or event.EF_tau125_medium1 or event.EF_xe60_verytight_noMu
        else:
            raise ValueError("No trigger condition defined for run %s" % event.RunNumber)
        

class TwoGoodLooseTaus(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.author != 2 and tau.numTrack > 0 and
                                      tau.pt > 18*GeV and (tau.tauLlhLoose == 1 or tau.JetBDTSigLoose == 1))
        return len(event.taus) > 1


class HTauSkim(ATLASStudent):

    def work(self):
        
        # initialize the TreeChain of all input files
        intree = TreeChain(self.fileset.treename,
                          files=self.fileset.files,
                          events=self.events)
        
        outtree = Tree(name=self.fileset.treename, file=self.output)
        outtree.set_buffer(intree.buffer, create_branches=True, visible=False)

        # copy TrigConfTree from first file in input list
        # make sure input files are all from the same run
        # use the prun --useContElementBoundary option if the
        # input container consists of run datasets
        metadirname = '%sMeta' % self.fileset.treename
        trigconf = intree.file['%s/TrigConfTree' % metadirname]
        metadir = self.output.mkdir(metadirname)
        metadir.cd()
        newtrigconf = trigconf.CloneTree(-1,'fast')
        newtrigconf.Write()
        self.output.cd()
        
        # set the event filters
        self.event_filters = EventFilterList([
            Triggers(),
            TwoGoodLooseTaus(passthrough = self.fileset.datatype == datasets.MC)
        ])
        intree.filters += self.event_filters

        # define tau collection
        intree.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)

        # entering the main event loop...
        for event in intree:
            
            outtree.Fill()

        # flush any remaining baskets to disk
        outtree.FlushBaskets()
        outtree.Write()
