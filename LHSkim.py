"""
==================
Skimming procedure
==================

For DATA skimming:
------------------

    1) Triggers

        See muTriggers and muMCTriggers in lephad/filters.py

    2) One Loose Tau

        - tau_author!=2 && tau_pT > 18 GeV && tau_numTrack > 0
        - tau_JetBDTLoose==1 || tau_tauLlhLoose==1

        * Note pT>18GeV cut was chosen to be able to fluctuate the TES.


For MC skimming:
----------------

    1) Triggers: OR of above.


======================
Extra Trees/Histograms
======================

Some information has to be kept during the skimming.

Just after the "Trigger" requirement, the histograms/trees should be made according
to each trigger decision.  (three histograms for each variable below)

 * number of events within GRL (to check consistency of LumiCalc)
 * mu-distribution
 * # of vertices
 * # of LOOSE taus (pt>18GeV, LLH loose || BDT loose)
 * # of EF tau trigger object  (no matching is necessary)
 * tau pT spectrum
 * EF tau pT spectrum  (no matching is necessary)
 * MET
 * anything else???
"""


import ROOT
import math

from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from atlastools.batch import ATLASStudent

from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel, TreeBuffer
from rootpy.types import *
from rootpy.io import open as ropen
from rootpy.plotting import Hist
from rootpy.registry import lookup_demotion

from higgstautau.mixins import TauFourMomentum
from higgstautau.mixins import FourMomentum
from higgstautau.lephad.filters import muTriggers, eTriggers, MCTriggers, \
                                       data_triggers, mc_triggers, \
                                       tau_skimselection, muon_skimselection, electron_skimselection, vertex_selection
import goodruns


ROOT.gErrorIgnoreLevel = ROOT.kFatal


class SkimExtraModel(TreeModel):

    number_of_good_vertices = IntCol()
    number_of_good_taus = IntCol()
    number_of_good_muons = IntCol()
    number_of_good_electrons = IntCol()


class SkimExtraTauPtModel(TreeModel):

    tau_pt = FloatCol()
    muon_pt = FloatCol()


class LHSkim(ATLASStudent):

    def work(self):

        # initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                          files=self.files,
                          events=self.events)

        outtree = Tree(name=self.metadata.treename,
                       file=self.output,
                       model=SkimExtraModel)
        outtree.set_buffer(intree.buffer, create_branches=True, visible=False)

        # Define additional trimming
        blocks_to_remove= ['ph_*',
                           'cl_*',
                           'trk_*',
                           'tau_otherTrk_*',
                           'jet_AntiKt6LCTopo_*',
                           'jet_AntiKt4TopoEM_*',
                           'mu_muid*',
                           'mu_calo*'
                           ]
        
        variables_to_keep = outtree.glob('*', prune=blocks_to_remove) + data_triggers
        outtree.activate(variables_to_keep, exclusive=True)

        if self.metadata.datatype == datasets.DATA:
            # outtree_extra holds info for events not included in the skim
            outtree_extra = Tree(name=self.metadata.treename + '_failed_skim_after_trigger',
                                 file=self.output,
                                 model=SkimExtraModel + SkimExtraTauPtModel)

            extra_variables = [
                'trig_EF_tau_pt',
                'actualIntPerXing',
                'averageIntPerXing',
                'MET_RefFinal_phi',
                'MET_RefFinal_et'
                'MET_RefFinal_sumet',
                'MET_LocHadTopo_phi',
                'MET_LocHadTopo_et',
                'MET_LocHadTopo_sumet',
                'EventNumber',
                'RunNumber',
                'lbn'
            ] + data_triggers

            outtree_extra.set_buffer(intree.buffer, variables=extra_variables, create_branches=True, visible=False)

        # merge TrigConfTrees
        metadirname = '%sMeta' % self.metadata.treename
        trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
        map(trigconfchain.Add, self.files)
        metadir = self.output.mkdir(metadirname)
        metadir.cd()
        trigconfchain.Merge(self.output, -1, 'fast keep')
        self.output.cd()

        if self.metadata.datatype == datasets.DATA:
            # merge GRL XML strings
            grls = []
            merged_grl = goodruns.GRL()
            for fname in self.files:
                merged_grl |= goodruns.GRL('%s:/Lumi/%s' % (fname, self.metadata.treename))
            lumi_dir = self.output.mkdir('Lumi')
            lumi_dir.cd()
            xml_string= ROOT.TObjString(merged_grl.str())
            xml_string.Write(self.metadata.treename)
            self.output.cd()

        # set the event filters
        trigger_filter = None
        if self.metadata.datatype == datasets.DATA:
            if self.metadata.title == datasets.MUON:
                trigger_filter = muTriggers()
            if self.metadata.title == datasets.ELEC:
                trigger_filter = eTriggers()
        else:
            trigger_filter = MCTriggers()

        # define collections for preselection
        intree.define_collection(name='taus', prefix='tau_', size='tau_n', mix=TauFourMomentum)
        intree.define_collection(name='muons', prefix='mu_staco_', size='mu_staco_n')
        intree.define_collection(name='electrons', prefix='el_', size='el_n')
        intree.define_collection(name='vertices', prefix='vxp_', size='vxp_n')

        #Cut Flow counters
        nevents = 0
        nevents_passing_trigger = 0
        nevents_with_good_vertices = 0
        nevents_with_good_muons = 0
        nevents_with_good_electrons = 0
        nevents_with_good_taus = 0
        nevents_with_good_lephad = 0
        
        # entering the main event loop...
        for event in intree:
            nevents += 1
            if trigger_filter(event):
                nevents_passing_trigger +=1
                #Vertex requirements
                event.vertices.select(lambda vxp: vertex_selection(vxp))
                number_of_good_vertices = len(event.vertices)
                if number_of_good_vertices > 0:
                    nevents_with_good_vertices +=1

                # Tau Preselection
                # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Taus
                event.taus.select(lambda tau : tau_skimselection(tau))
                number_of_good_taus = len(event.taus)
                if number_of_good_taus > 0:
                    nevents_with_good_taus +=1

                # Muon Preselection
                # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Muons
                event.muons.select(lambda mu : muon_skimselection(mu)) 
                number_of_good_muons = len(event.muons)
                if number_of_good_muons > 0:
                    nevents_with_good_muons +=1

                # Electron Preselection
                # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Electrons
                event.electrons.select(lambda e : electron_skimselection(e)) 
                number_of_good_electrons = len(event.electrons)
                if number_of_good_electrons > 0:
                    nevents_with_good_electrons +=1

                if ((number_of_good_taus > 0 and self.metadata.datatype == datasets.DATA) and
                    ((number_of_good_muons > 0 and self.metadata.title == datasets.MUON) or \
                    (number_of_good_electrons > 0 and self.metadata.title == datasets.ELEC))) or \
                    (self.metadata.datatype == datasets.MC):
                    nevents_with_good_lephad +=1
                    outtree.number_of_good_vertices = number_of_good_vertices
                    outtree.number_of_good_taus = number_of_good_taus
                    outtree.number_of_good_muons = number_of_good_muons
                    outtree.number_of_good_electrons = number_of_good_electrons
                    outtree.Fill()
                else:
                    outtree_extra.number_of_good_vertices = number_of_good_vertices
                    outtree_extra.number_of_good_taus = number_of_good_taus
                    outtree_extra.number_of_good_muons = number_of_good_muons
                    outtree_extra.number_of_good_electrons = number_of_good_electrons
                    if event.taus:
                        # There can be at most one good tau if this event failed the skim
                        outtree_extra.tau_pt = event.taus[0].pt
                    else:
                        outtree_extra.tau_pt = -1111.
                    outtree_extra.Fill()

                    if event.muons:
                        # There can be at most one good muon if this event failed the skim
                        outtree_extra.muon_pt = event.muons[0].pt
                    else:
                        outtree_extra.muon_pt = -1111.
                    outtree_extra.Fill()

        self.output.cd()

        # store the original number of events
        cutflow = Hist(7, 0, 7, name='cutflow', type='D')
        cutflow[0] = nevents
        cutflow[1] = nevents_passing_trigger
        cutflow[2] = nevents_with_good_vertices
        cutflow[3] = nevents_with_good_muons
        cutflow[4] = nevents_with_good_electrons
        cutflow[5] = nevents_with_good_taus
        cutflow[6] = nevents_with_good_lephad
        cutflow.Write()

        # flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()
        if self.metadata.datatype == datasets.DATA:
            outtree_extra.FlushBaskets()
            outtree_extra.Write()
