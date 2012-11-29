"""
Add instructions to run the skim here
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

from higgstautau.mixins import TauFourMomentum, ElectronFourMomentum
from higgstautau.lephad.filters import tau_skimselection, muon_skimselection, electron_skimselection, \
                                       OverlapCheck, PrepareInputTree

from higgstautau.patches import ElectronIDpatch, TauIDpatch
import goodruns


ROOT.gErrorIgnoreLevel = ROOT.kFatal

YEAR = 2011

REMOVE = [
    'cl_*',
    'trk_*',
    'tau_otherTrk_*',
    'jet_AntiKt6LCTopo_*',
    'mu_muid*',
    'mu_calo*'
]

# override globs above
KEEP = [
    'trk_atTJVA_*'
    ]

class SkimExtraModel(TreeModel):

    number_of_good_taus = IntCol()
    number_of_good_muons = IntCol()
    number_of_good_electrons = IntCol()
    el_trk_d3pdindex = ROOT.vector('int')
    mu_staco_trk_d3pdindex = ROOT.vector('int')

class LHSkim(ATLASStudent):

    def work(self):

        # Determine year
        YEAR = 2012
        if self.metadata.year == 2011:
            YEAR = 2011

        # Initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                          files=self.files,
                          events=self.events)

        self.output.cd()
        outtree = Tree(name=self.metadata.treename,
                       model=SkimExtraModel)
#        outtree.set_buffer(intree.buffer, create_branches=True, visible=False)

        outtree.set_buffer(
                intree.buffer,
                ignore_branches=intree.glob(REMOVE,exclude=KEEP),
                create_branches=True,
                transfer_objects=True,
            visible=False)

        # Merge TrigConfTrees
        metadirname = '%sMeta' % self.metadata.treename
        trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
        map(trigconfchain.Add, self.files)
        metadir = self.output.mkdir(metadirname)
        metadir.cd()
        trigconfchain.Merge(self.output, -1, 'fast keep')
        self.output.cd()


        # merge GRL XML strings
        if self.metadata.datatype == datasets.DATA:
            grls = []
            merged_grl = goodruns.GRL()
            for fname in self.files:
                merged_grl |= goodruns.GRL('%s:/Lumi/%s' % (fname, self.metadata.treename))
            lumi_dir = self.output.mkdir('Lumi')
            lumi_dir.cd()
            xml_string= ROOT.TObjString(merged_grl.str())
            xml_string.Write(self.metadata.treename)
            self.output.cd()


        # Define collections for preselection
        intree.define_collection(name='taus', prefix='tau_', size='tau_n', mix=TauFourMomentum)
        intree.define_collection(name='muons', prefix='mu_staco_', size='mu_staco_n')
        intree.define_collection(name='electrons', prefix='el_', size='el_n', mix=ElectronFourMomentum)
        intree.define_object(name='isLTT', prefix='')
        intree.define_object(name='leptonType', prefix='')

        # set the event filters

        eventFilterList = [PrepareInputTree()]

        if YEAR == 2012:
            eventFilterList.append(ElectronIDpatch())
            eventFilterList.append(TauIDpatch(YEAR))

        event_filters = EventFilterList(eventFilterList)

        self.filters['event'] = event_filters
        intree.filters += event_filters


        # Cut Flow counters
        nevents = 0
        nevents_mc_weights = 0
        nevents_passing_trigger = 0
        nevents_with_good_muons = 0
        nevents_with_good_electrons = 0
        nevents_with_good_taus = 0
        nevents_with_good_lephad = 0


        # Entering the main event loop...
        for event in intree:
            nevents += 1

            ## Findin the track index for the electron
            outtree.el_trk_d3pdindex.clear()
            for i_el in xrange(intree.el_n):
                _el_trk_d3pdindex=-1
                _el_trk_mindist = 999;
                for i_trk in xrange(intree.trk_n):
                    _trk_mindist  = math.pow(intree.trk_theta[i_trk] - intree.el_tracktheta[i_el],2)
                    _trk_mindist += math.pow(intree.trk_phi[i_trk] - intree.el_trackphi[i_el],2)
                    if _trk_mindist < _el_trk_mindist :
                       _el_trk_mindist = _trk_mindist
                       _el_trk_d3pdindex = i_trk
                outtree.el_trk_d3pdindex.push_back(_el_trk_d3pdindex)


            ## Finding the track index for the muon
            outtree.mu_staco_trk_d3pdindex.clear()
            for i_mu in xrange(intree.mu_staco_n):
                _mu_staco_trk_d3pdindex=-1
                _mu_staco_trk_mindist = 999;
                for i_trk in xrange(intree.trk_n):
                    _trk_mindist  = math.pow(intree.trk_theta[i_trk] - intree.mu_staco_tracktheta[i_mu],2)
                    _trk_mindist += math.pow(intree.trk_phi[i_trk] - intree.mu_staco_trackphi[i_mu],2)
                    if _trk_mindist < _mu_staco_trk_mindist :
                       _mu_staco_trk_mindist = _trk_mindist
                       _mu_staco_trk_d3pdindex = i_trk
                outtree.mu_staco_trk_d3pdindex.push_back(_mu_staco_trk_d3pdindex)

            if self.metadata.datatype == datasets.MC:
                nevents_mc_weights += event.mc_event_weight

            if True: #trigger_filter(event):
                nevents_passing_trigger +=1

                # Tau Skim Selection
                event.taus.select(lambda tau : tau_skimselection(tau))
                number_of_good_taus = len(event.taus)
                if number_of_good_taus > 0:
                    nevents_with_good_taus +=1

                # Muon Skim Selection
                event.muons.select(lambda mu : muon_skimselection(mu))
                number_of_good_muons = len(event.muons)
                if number_of_good_muons > 0:
                    nevents_with_good_muons +=1

                # Electron Skim Selection
                event.electrons.select(lambda e : electron_skimselection(e))
                number_of_good_electrons = len(event.electrons)
                if number_of_good_electrons > 0:
                    nevents_with_good_electrons +=1

                # Steer object selection for data and MC
                hasTau = number_of_good_taus > 0
                hasLeptonData = (number_of_good_muons > 0 and self.metadata.title == datasets.MUON and OverlapCheck(event, DoMuonCheck = True)) or \
                                (number_of_good_electrons > 0 and self.metadata.title == datasets.ELEC and OverlapCheck(event, DoElectronCheck = True))
                hasLeptonMC = (number_of_good_muons > 0 and OverlapCheck(event, DoMuonCheck = True)) or \
                              (number_of_good_electrons > 0 and OverlapCheck(event, DoElectronCheck = True))
                isLepTau = self.metadata.title == datasets.TAU
                isData = self.metadata.datatype == datasets.DATA
                isMC = self.metadata.datatype == datasets.MC or self.metadata.datatype == datasets.EMBED

                # Fill the event
                if (isData and hasTau and hasLeptonData) or (isMC and hasTau and hasLeptonMC) or (isData and hasTau and hasLeptonMC and isLepTau):
                    nevents_with_good_lephad +=1
                    outtree.number_of_good_taus = number_of_good_taus
                    outtree.number_of_good_muons = number_of_good_muons
                    outtree.number_of_good_electrons = number_of_good_electrons
                    outtree.Fill()


        self.output.cd()

        # Store the original number of events
        cutflow = Hist(7, 0, 7, name='cutflow', type='D')
        cutflow[0] = nevents
        cutflow[1] = nevents_mc_weights
        cutflow[2] = nevents_passing_trigger
        cutflow[3] = nevents_with_good_muons
        cutflow[4] = nevents_with_good_electrons
        cutflow[5] = nevents_with_good_taus
        cutflow[6] = nevents_with_good_lephad
        cutflow.Write()

        # Flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()
