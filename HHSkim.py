"""
==================
Skimming procedure
==================

For DATA skimming:
------------------

    1) Triggers

        See Triggers in filters.py

    2) Two LOOSE taus

        - tau_author!=2 && tau_pT > 18 GeV && tau_numTrack > 0
        - tau_JetBDTLoose==1 || tau_tauLlhLoose==1

        * Note pT>18GeV cut was chosen to be able to fluctuate the TES.

        With these two selection (1 & 2), we expect factor 50 reduction which
        results in ~400-500 GB disk space at 5 fb-1.


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
from rootpy.tree import Tree, TreeChain, TreeModel
from rootpy.types import *
from rootpy.io import open as ropen
from rootpy.plotting import Hist

from higgstautau.mixins import TauFourMomentum
from higgstautau.hadhad.filters import SkimmingDataTriggers, SkimmingMCTriggers, \
                                       data_triggers, mc_triggers
import goodruns


ROOT.gErrorIgnoreLevel = ROOT.kFatal


class SkimExtraModel(TreeModel):

    number_of_good_vertices = IntCol()
    number_of_good_taus = IntCol()


class SkimExtraTauPtModel(TreeModel):

    tau_pt = FloatCol()

#TODO also create pileup reweighting files here
#TODO store mc_event_weight for events failing skim in MC

class HHSkim(ATLASStudent):

    def work(self):

        # initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                          files=self.files,
                          events=self.events)

        outtree = Tree(name=self.metadata.treename,
                       file=self.output,
                       model=SkimExtraModel)
        outtree.set_buffer(intree.buffer, create_branches=True, visible=False)

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
        if self.metadata.datatype == datasets.DATA:
            trigger_filter = SkimmingDataTriggers()
        else:
            trigger_filter = SkimmingMCTriggers()

        # define tau collection
        intree.define_collection(name='taus', prefix='tau_', size='tau_n', mix=TauFourMomentum)
        intree.define_collection(name='vertices', prefix='vxp_', size='vxp_n')

        nevents = 0
        # entering the main event loop...
        for event in intree:
            nevents += 1
            if trigger_filter(event):
                event.vertices.select(lambda vxp: (vxp.type == 1 and vxp.nTracks >= 4) or (vxp.type == 3 and vxp.nTracks >= 2))
                number_of_good_vertices = len(event.vertices)
                event.taus.select(lambda tau: tau.author != 2 and tau.numTrack > 0 and
                                              tau.pt > 18*GeV and
                                              (tau.tauLlhLoose == 1 or tau.JetBDTSigLoose == 1))
                number_of_good_taus = len(event.taus)
                if (number_of_good_taus > 1 and self.metadata.datatype == datasets.DATA) or \
                    self.metadata.datatype == datasets.MC:
                    outtree.number_of_good_vertices = number_of_good_vertices
                    outtree.number_of_good_taus = number_of_good_taus
                    outtree.Fill()
                else:
                    outtree_extra.number_of_good_vertices = number_of_good_vertices
                    outtree_extra.number_of_good_taus = number_of_good_taus
                    if event.taus:
                        # There can be at most one good tau if this event failed the skim
                        outtree_extra.tau_pt = event.taus[0].pt
                    else:
                        outtree_extra.tau_pt = -1111.
                    outtree_extra.Fill()

        self.output.cd()

        # store the original number of events
        cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        cutflow[0] = nevents
        cutflow.Write()

        # flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()
        if self.metadata.datatype == datasets.DATA:
            outtree_extra.FlushBaskets()
            outtree_extra.Write()


"""
From Zinonas

to obtain a reduction factor of ~ 10.4 :

ch.SetBranchStatus("cl_*",    0)
ch.SetBranchStatus("ph_*",    0)
ch.SetBranchStatus("jet_AntiKt4TopoEM_*",    0)
ch.SetBranchStatus("jet_AntiKt4LCTopo_*",    0)
ch.SetBranchStatus("jet_AntiKt6LCTopo_*",    0)
ch.SetBranchStatus("jet_flavor_*",    0)
ch.SetBranchStatus("jet_*Assoc*",    0)
ch.SetBranchStatus("tau_jet_*",    0)
ch.SetBranchStatus("tau_seedCalo_*",    0)
ch.SetBranchStatus("tau_track_atTJVA_*",    0)
ch.SetBranchStatus("tau_secvtx_*",    0)
ch.SetBranchStatus("tau_privtx_*",    0)
ch.SetBranchStatus("tau_calcVars_*",    0)
ch.SetBranchStatus("tau_seedTrk_*",    0)
ch.SetBranchStatus("tau_otherTrk_*",    0)
ch.SetBranchStatus("tau_Pi0Cluster_*",    0)
ch.SetBranchStatus("tau_cell_*",    0)
ch.SetBranchStatus("tau_cluster_*",    0)
ch.SetBranchStatus("tau_privtx_*",    0)
ch.SetBranchStatus("tau_*Assoc*",    0)
ch.SetBranchStatus("tau_MET_*",    0)
ch.SetBranchStatus("tau_EF_*",    0)
ch.SetBranchStatus("tau_L2_*",    0)
ch.SetBranchStatus("tau_L1_*",    0)
ch.SetBranchStatus("EF_2e*",    0)
ch.SetBranchStatus("EF_2mu*",    0)
ch.SetBranchStatus("EF_2j*",    0)
ch.SetBranchStatus("EF_xe*",    0)
ch.SetBranchStatus("EF_xs*",    0)
ch.SetBranchStatus("EF_e*",    0)
ch.SetBranchStatus("EF_mu*",    0)
ch.SetBranchStatus("EF_MU*",    0)
ch.SetBranchStatus("EF_g*",    0)
ch.SetBranchStatus("EF_j*",    0)
ch.SetBranchStatus("EF_g*",    0)
ch.SetBranchStatus("L1_*",    0)
ch.SetBranchStatus("L2_*",    0)
ch.SetBranchStatus("trig_L1_emtau_*",    0)
ch.SetBranchStatus("trig_L1_esum_*",    0)
ch.SetBranchStatus("trig_L2_met_*",    0)
ch.SetBranchStatus("trig_EF_met_*",    0)
ch.SetBranchStatus("trig_L2_feb_*",    0)
ch.SetBranchStatus("trig_L1_jet_*",    0)
ch.SetBranchStatus("trig_L2_jet_*",    0)
ch.SetBranchStatus("trig_EF_jet_*",    0)
ch.SetBranchStatus("trig_EF_trigmuonef_*",    0)
ch.SetBranchStatus("trig_EF_el_*",    0)
ch.SetBranchStatus("trig_RoI_EF_mu_*",    0)
ch.SetBranchStatus("trig_RoI_EF_e_*",    0)
ch.SetBranchStatus("trig_RoI_L2_tau_*",    0)
ch.SetBranchStatus("trig_EF_tau_seedCalo_*",    0)
ch.SetBranchStatus("trig_EF_tau_calcVars*",    0)
ch.SetBranchStatus("trig_L2_trk_*",    0)
ch.SetBranchStatus("trig_L2_tau_*",    0)
ch.SetBranchStatus("mc_*",    0)
ch.SetBranchStatus("mcevt_*",    0)
ch.SetBranchStatus("true*",    0)
ch.SetBranchStatus("muonTruth*",    0)
ch.SetBranchStatus("jet_antikt4truth_*",    0)
ch.SetBranchStatus("trk_*",    0)
ch.SetBranchStatus("collcand_*",    0)

ch.SetBranchStatus("el_*",    0)
ch.SetBranchStatus("el_cl_E",    1)
ch.SetBranchStatus("el_tracketa",    1)
ch.SetBranchStatus("el_trackphi",    1)
ch.SetBranchStatus("el_author",    1)
ch.SetBranchStatus("el_charge",    1)
ch.SetBranchStatus("el_loosePP",    1)
ch.SetBranchStatus("el_mediumPP",    1)
ch.SetBranchStatus("el_tightPP",    1)
ch.SetBranchStatus("el_OQ",    1)

ch.SetBranchStatus("mu_*",    0)
ch.SetBranchStatus("mu_staco_E",    1)
ch.SetBranchStatus("mu_staco_pt",    1)
ch.SetBranchStatus("mu_staco_eta",    1)
ch.SetBranchStatus("mu_staco_phi",    1)
ch.SetBranchStatus("mu_staco_loose",    1)
ch.SetBranchStatus("mu_staco_medium",    1)
ch.SetBranchStatus("mu_staco_tight",    1)
ch.SetBranchStatus("mu_staco_isSegmentTaggedMuon",    1)
ch.SetBranchStatus("mu_staco_expectBLayerHit",    1)
ch.SetBranchStatus("mu_staco_nBLHits",    1)
ch.SetBranchStatus("mu_staco_nPixHits",    1)
ch.SetBranchStatus("mu_staco_nPixelDeadSensors",    1)
ch.SetBranchStatus("mu_staco_nSCTHits",    1)
ch.SetBranchStatus("mu_staco_nSCTDeadSensors",    1)
ch.SetBranchStatus("mu_staco_nPixHoles",    1)
ch.SetBranchStatus("mu_staco_nSCTHoles",    1)
ch.SetBranchStatus("mu_staco_nTRTHits",    1)
ch.SetBranchStatus("mu_staco_nTRTOutliers",    1)
"""
