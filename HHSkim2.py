import ROOT
import math

from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from atlastools.batch import ATLASStudent
from atlastools.filtering import GRLFilter

from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel
from rootpy.types import *
from rootpy.io import open as ropen

from higgstautau.mixins import *
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau.trigger import update_trigger_config, get_trigger_config

import goodruns

#ROOT.gErrorIgnoreLevel = ROOT.kFatal
VALIDATE = True
YEAR = 2011

class TriggerMatching(TreeModel):

    tau_trigger_match_index = ROOT.vector('int')
    tau_trigger_match_thresh = ROOT.vector('int')


class HHSkim2(ATLASStudent):
    """
    Apply cleaning
    """

    def work(self):

        # merge TrigConfTrees
        metadirname = '%sMeta' % self.metadata.treename
        trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
        map(trigconfchain.Add, self.files)
        metadir = self.output.mkdir(metadirname)
        metadir.cd()
        trigconfchain.Merge(self.output, -1, 'fast keep')
        self.output.cd()

        # merge the cutflow hists from the first skim
        cutflow = None
        for fname in self.files:
            with ropen(fname) as f:
                if cutflow is None:
                    cutflow = f.cutflow.Clone(name='cutflow')
                    cutflow.SetDirectory(0)
                else:
                    cutflow += f.cutflow
        self.output.cd()
        cutflow.Write()

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

            # merge outtree_extras from the first skim
            extra_trees = ROOT.TChain(self.metadata.treename +
                                      '_failed_skim_after_trigger')
            map(extra_trees.Add, self.files)
            extra_trees.Merge(self.output, -1, 'fast keep')
            self.output.cd()

        onfilechange = []

        # trigger config tool to read trigger info in the ntuples
        trigger_config = get_trigger_config()

        # update the trigger config maps on every file change
        onfilechange.append((update_trigger_config, (trigger_config,)))

        # initialize the TreeChain of all input files
        chain = TreeChain(self.metadata.treename,
                          files=self.files,
                          events=self.events,
                          onfilechange=onfilechange)

        Model = None
        if self.metadata.datatype == datasets.DATA:
            Model = TriggerMatching

        tree = Tree(name=self.metadata.treename,
                    file=self.output,
                    model=Model)

        tree.set_buffer(chain.buffer,
                        create_branches=True)

        # set the event filters
        event_filters = EventFilterList([
            GRLFilter(self.grl, passthrough=self.metadata.datatype != datasets.DATA),
            Triggers(datatype=self.metadata.datatype,
                     year=YEAR,
                     skim=False),
            PriVertex(),
            LArError(),
            LArHole(datatype=self.metadata.datatype),
            JetCleaning(),
            ElectronVeto(),
            MuonVeto(),
            TauAuthor(),
            TauHasTrack(),
            TauMuonVeto(),
            TauElectronVeto(),
            TauPT(),
            TauEta(),
            TauCrack(),
            TauLArHole(),
            TauIDMedium(),
            TauTriggerMatch(config=trigger_config,
                            year=YEAR,
                            datatype=self.metadata.datatype,
                            skim=True,
                            tree=tree),
            TauLeadSublead(lead=35*GeV,
                           sublead=25*GeV),
        ])

        self.filters['event'] = event_filters

        chain.filters += event_filters

        # define tree collections
        chain.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        chain.define_collection(name="taus_EF", prefix="trig_EF_tau_",
                                size="trig_EF_tau_n", mix=TauFourMomentum)

        # jet_* etc. is AntiKt4LCTopo_* in tau-perf D3PDs
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        chain.define_collection(name="truetaus", prefix="trueTau_", size="trueTau_n", mix=MCTauFourMomentum)
        chain.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        chain.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n")
        chain.define_collection(name="electrons", prefix="el_", size="el_n")
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        # define tree objects
        tree.define_object(name='tau1', prefix='tau1_')
        tree.define_object(name='tau2', prefix='tau2_')
        tree.define_object(name='jet1', prefix='jet1_')
        tree.define_object(name='jet2', prefix='jet2_')

        if VALIDATE:
            validate_log = open(self.metadata.name + '_validate.log', 'w')

        # entering the main event loop...
        for event in chain:
            if VALIDATE:
                print >> validate_log, event.EventNumber,
                for tau in event.taus:
                    print >> validate_log, tau.index,
                print >> validate_log
            tree.Fill()

        if VALIDATE:
            validate_log.close()

        self.output.cd()

        # flush any baskets remaining in memory to disk
        tree.FlushBaskets()
        tree.Write()
