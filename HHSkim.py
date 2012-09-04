"""
==================
Skimming procedure
==================

For DATA skimming:
------------------

    1) Triggers

        if 177986 <= event.RunNumber <= 187815: # Periods B-K
            return event.EF_tau29_medium1_tau20_medium1
        elif 188902 <= event.RunNumber <= 191933: # Periods L-M
            return event.EF_tau29T_medium1_tau20T_medium1

        See Triggers in filters.py

    2) Two LOOSE taus

        - tau_author!=2 && tau_pT > 18 GeV && tau_numTrack > 0
        - tau_JetBDTLoose==1 || tau_tauLlhLoose==1

        * Note pT>18GeV cut was chosen to be able to fluctuate the TES.

        With these two selection (1 & 2), we expect factor 50 reduction which
        results in ~400-500 GB disk space at 5 fb-1.


For MC skimming:
----------------

    1) Triggers: same as above.


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


=================
Trigger Emulation
=================

Emulation is performed using the CoEPPTrigTool package

180164: EF_tau29_medium1_tau20_medium1_Hypo_00_02_42
183003: EF_tau29_medium1_tau20_medium1_Hypo_00_03_02
186169: EF_tau29_medium1_tau20_medium1_Hypo_00_03_02
189751: EF_tau29T_medium1_tau20T_medium1_Hypo_00_03_02


==============
Branch removal
==============

See branches_remove and branches_keep below
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

from higgstautau.mixins import TauFourMomentum, MCParticle
from higgstautau.hadhad.filters import Triggers
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.filters import vertex_selection
import higgstautau.skimming.hadhad as hhskimming
from hhskimming import branches as hhbranches
from hhskimming.models import *
from higgtautau.pileup import PileupTemplates, PileupReweight

import goodruns


#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class HHSkim(ATLASStudent):

    def work(self):

        if self.metadata.datatype != datasets.EMBED:
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
                merged_grl |= goodruns.GRL(
                        '%s:/Lumi/%s' % (fname, self.metadata.treename))
            lumi_dir = self.output.mkdir('Lumi')
            lumi_dir.cd()
            xml_string= ROOT.TObjString(merged_grl.str())
            xml_string.Write(self.metadata.treename)
            self.output.cd()

        Model = SkimExtraModel
        if self.metadata.datatype == datasets.MC:
            Model += TriggerEmulation

        outtree = Tree(
                name=self.metadata.treename,
                file=self.output,
                model=Model)

        event_filters = EventFilterList([
            EmbeddingPileupPatch(
                passthrough=self.metadata.datatype != datasets.EMBED),
            PileupTemplates(
                passthrough=self.metadata.datatype != datasets.MC),
            ExtraInfoTree(
                )
            TauTriggerEmulation(
                year=self.metadata.year,
                outtree=outtree,
                passthrough=self.metadata.datatype != datasets.MC),
            Triggers(
                datatype=self.metadata.datatype,
                year=self.metadata.year,
                skim=True),
            # the BDT bits are broken in the p1130 production, correct them
            # DON'T FORGET TO REMOVE THIS WHEN SWITCHING TO A NEWER
            # PRODUCTION TAG!!!
            TauIDpatch(
                year=self.metadata.year,
                passthrough=self.metadata.year != 2012),
            # patch electron ID for 2012
            ElectronIDpatch(
                passthrough=self.metadata.year != 2012),
            PileupReweight(
                passthrough=self.metadata.datatype != datasets.MC),
        ])

        self.filters['event'] = event_filters

        onfilechange = []

        if self.metadata.datatype == datasets.MC:
            # TODO get trigger emul tool here
            onfilechange.append(
                        (update_trigger_trees, (self, trigger_tool_wrapper,)))

        # initialize the TreeChain of all input files
        chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                events=self.events,
                onfilechange=onfilechange,
                filters=event_filters)

        outtree.set_buffer(
                chain.buffer,
                ignore_branches=chain.glob(
                    hhbranches.REMOVE,
                    prune=hhbranches.KEEP),
                create_branches=True,
                visible=False)

        """
        if self.metadata.datatype == datasets.DATA:
            # outtree_extra holds info for events not included in the skim
            outtree_extra = Tree(
                    name=self.metadata.treename + '_failed_skim_after_trigger',
                    file=self.output,
                    model=SkimExtraModel + SkimExtraTauPtModel)

            extra_variables = [
                'trig_EF_tau_pt',
                'actualIntPerXing',
                'averageIntPerXing',
                'MET_RefFinal_BDTMedium_phi',
                'MET_RefFinal_BDTMedium_et'
                'MET_RefFinal_BDTMedium_sumet',
                'EventNumber',
                'RunNumber',
                'lbn'
            ]

            if self.metadata.year % 1000 == 11:
                extra_variables += Triggers.triggers_11
            elif self.metadata.year % 1000 == 12:
                extra_variables += Triggers.triggers_12
            else:
                raise ValueError("No triggers defined for year %d" % year)

            outtree_extra.set_buffer(
                    chain.buffer,
                    branches=extra_variables,
                    create_branches=True,
                    visible=False)
        else:
            # write out some branches for all events
            # that failed the skim before trigger only for MC
            # Used to get the total pileup reweighting sum
            if self.metadata.datatype == datasets.MC:
                outtree_extra = Tree(
                        name=self.metadata.treename + '_failed_skim_before_trigger',
                        file=self.output)
            else: #embedding
                outtree_extra = Tree(
                        name=self.metadata.treename + '_failed_skim_before_selection',
                        file=self.output)
            extra_variables = [
                'actualIntPerXing',
                'averageIntPerXing',
                'RunNumber',
            ]
            outtree_extra.set_buffer(
                    chain.buffer,
                    branches=extra_variables,
                    create_branches=True, visible=False)
        """

        # define tau collection
        chain.define_collection(
                name='taus',
                prefix='tau_',
                size='tau_n',
                mix=TauFourMomentum)
        chain.define_collection(
                name="electrons",
                prefix="el_",
                size="el_n")
        chain.define_collection(
                name='vertices',
                prefix='vxp_',
                size='vxp_n')
        chain.define_collection(
                name="mc",
                prefix="mc_",
                size="mc_n",
                mix=MCParticle)

        nevents = 0
        nevents_mc_weight = 0

        # entering the main event loop...
        for event in chain:

            if self.metadata.datatype == datasets.MC:
                nevents_mc_weight += event.mc_event_weight

            event.vertices.select(vertex_selection)
            number_of_good_vertices = len(event.vertices)
            event.taus.select(lambda tau:
                    tau.author != 2 and
                    tau.numTrack > 0 and
                    tau.pt > 18*GeV and
                    (tau.tauLlhLoose == 1 or tau.JetBDTSigLoose == 1))
            number_of_good_taus = len(event.taus)

            if (number_of_good_taus > 1 and
                (self.metadata.datatype in (datasets.DATA, datasets.EMBED))) \
               or self.metadata.datatype == datasets.MC:
                outtree.number_of_good_vertices = number_of_good_vertices
                outtree.Fill()
            """
            elif self.metadata.datatype == datasets.DATA:
                outtree_extra.number_of_good_vertices = number_of_good_vertices
                if event.taus:
                    # There can be at most one good tau if this event failed the skim
                    outtree_extra.tau_pt = event.taus[0].pt
                else:
                    outtree_extra.tau_pt = -1111.
                outtree_extra.Fill()
            elif self.metadata.datatype == datasets.EMBED:
                outtree_extra.Fill()
            """

        self.output.cd()

        if self.metadata.datatype == datasets.MC:
            # store the original weighted number of events
            cutflow = Hist(2, 0, 2, name='cutflow', type='D')
            cutflow[1] = nevents_mc_weight
        else:
            cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        # store the original number of events
        cutflow[0] = nevents
        cutflow.Write()

        # flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()
        """
        outtree_extra.FlushBaskets()
        outtree_extra.Write()
        """
