import warnings
import numpy as np
warnings.filterwarnings('error', category=np.ComplexWarning)

# jet calibration sometimes gives jets with pT=0
#warnings.filterwarnings('error', 'transvers momentum = 0!', RuntimeWarning)

import ROOT
import math

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.math.physics.vector import Vector2, LorentzVector
from rootpy.plotting import Hist
from rootpy.extern.argparse import ArgumentParser

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.filtering import GRLFilter
from atlastools.batch import ATLASStudent

from higgstautau.mixins import *
from higgstautau.models import *
from higgstautau.hadhad.models import *
from higgstautau import eventview
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau import mass
from higgstautau.embedding import EmbeddingPileupPatch, EmbeddingIsolation
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.emulation import TauTriggerEmulation, update_trigger_trees
from higgstautau.trigger.matching import TauTriggerMatchIndex, TauTriggerMatchThreshold
from higgstautau.trigger.efficiency import TauTriggerEfficiency
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.overlap import TauJetOverlapRemoval
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.pileup import (PileupReweight, PileupDataScale,
                                averageIntPerXingPatch,
                                get_pileup_reweighting_tool)
from higgstautau.hadhad.objects import define_objects
from higgstautau import log; log = log[__name__]

from goodruns import GRL
import subprocess

#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class HHProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def __init__(self, options, **kwargs):

        super(HHProcessor, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-terms', default=None)
        parser.add_argument('--redo-mmc', default=False, action='store_true')
        #parser.add_argument('--use-numjets25', default=False, action='store_true')
        parser.add_argument('--student-verbose', default=False, action='store_true')
        self.args = parser.parse_args(options)
        if self.args.syst_terms is not None:
            self.args.syst_terms = set([
                eval('Systematics.%s' % term) for term in
                self.args.syst_terms.split(',')])

    @staticmethod
    def merge(inputs, output, metadata):

        # merge output trees
        root_output = output + '.root'
        subprocess.call(['hadd', root_output] + inputs)

        if metadata.datatype == datasets.DATA:
            # merge GRLs
            grl = GRL()
            for input in inputs:
                grl |= GRL('%s:/lumi' % input)
            grl.save('%s:/lumi' % root_output)

    def work(self):
        """
        This is the one function that all "ATLASStudent"s must implement.
        """
        datatype = self.metadata.datatype
        year = self.metadata.year
        redo_mmc = self.args.redo_mmc
        #use_numjets25 = self.args.use_numjets25
        verbose = self.args.student_verbose

        # get pileup reweighting tool
        pileup_tool = get_pileup_reweighting_tool(
            year=year,
            use_defaults=True)

        OutputModel = RecoTauBlock + RecoJetBlock + EventVariables

        if datatype == datasets.MC:
            # only create truth branches for MC
            OutputModel += TrueTauBlock

            # add branches for VBF Higgs associated partons
            if 'VBF' in self.metadata.name:
                OutputModel += PartonBlock

        if datatype == datasets.EMBED:
            # add embedding systematics branches
            OutputModel += EmbeddingBlock

        onfilechange = []
        count_funcs = {}

        if datatype in (datasets.MC, datasets.EMBED):

            def mc_weight_count(event):
                return event.mc_event_weight

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        trigger_config = None

        if datatype != datasets.EMBED:
            # trigger config tool to read trigger info in the ntuples
            trigger_config = get_trigger_config()

            # update the trigger config maps on every file change
            onfilechange.append((update_trigger_config, (trigger_config,)))

        if datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file, tree):

                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))

        if datatype == datasets.DATA:
            merged_cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        else:
            merged_cutflow = Hist(2, 0, 2, name='cutflow', type='D')

        def update_cutflow(student, cutflow, name, file, tree):

            year = student.metadata.year
            datatype = student.metadata.datatype
            if datatype == datasets.MC:
                cutflow[0] += file.cutflow_event[0]
                cutflow[1] += file.cutflow_event_mc_weight[0]
            else:
                cutflow[0] += file.cutflow_event[0]

        onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                events=self.events,
                read_branches_on_demand=True,
                cache=True,
                onfilechange=onfilechange)

        # create output tree
        self.output.cd()
        tree = Tree(name='higgstautauhh', model=OutputModel)

        copied_variables = [
                'actualIntPerXing',
                'averageIntPerXing',
                'number_of_good_vertices',
                'RunNumber',
                'EventNumber',
                'lbn']

        if 'ggH' in self.metadata.name or 'VBFH' in self.metadata.name:
            # branched needed for theory uncertainties
            # only for VBF and ggH
            copied_variables += [
                'jet_antikt4truth_n',
                'jet_antikt4truth_pt',
                'jet_antikt4truth_m',
                'jet_antikt4truth_eta',
                'jet_antikt4truth_phi',
                'mc_n',
                'mc_pdgId',
                'mc_child_index',
                'mc_parent_index',
                'mc_pt',
                'mc_eta',
                'mc_phi',
                'mc_m',
                'mc_status',
            ]

        tree.set_buffer(
                chain._buffer,
                branches=copied_variables,
                create_branches=True,
                visible=False)
        chain.always_read(copied_variables)

        # set the event filters
        event_filters = EventFilterList([
        ])

        self.filters['event'] = event_filters

        chain._filters += event_filters
