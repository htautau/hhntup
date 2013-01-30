import warnings
import numpy as np
warnings.filterwarnings('error', category=np.ComplexWarning)

# jet calibration sometimes gives jets with pT=0
#warnings.filterwarnings('error', 'transvers momentum = 0!', RuntimeWarning)

import ROOT
import math

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.math.physics.vector import Vector2
from rootpy.plotting import Hist
from rootpy.io import open as ropen
from rootpy.extern.argparse import ArgumentParser
from rootpy.extern.hep import pdg

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.filtering import GRLFilter
from atlastools.batch import ATLASStudent

from higgstautau.mixins import *
from higgstautau import hepmc
from higgstautau import tautools
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau import mass
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.emulation import TauTriggerEmulation, update_trigger_trees
from higgstautau.trigger.matching import TauTriggerMatchIndex, TauTriggerMatchThreshold
from higgstautau.trigger.efficiency import TauTriggerEfficiency

from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.pileup import PileupReweight
from higgstautau.hadhad.objects import define_objects
from higgstautau.ditaumass.models import *

from goodruns import GRL
import subprocess

#ROOT.gErrorIgnoreLevel = ROOT.kFatal


def get_taus(event):

    # attempt to get the resonance Z (23) or Higgs (25)
    resonance = tautools.get_particles(event, [25, 23],
            num_expected=1)

    if not resonance:
        return None, None

    # get the resonance just before the decay
    resonance = resonance[0].last_self

    # collect decay products (taus and photons)
    tau_decays = []
    for child in resonance.iter_children():
        if abs(child.pdgId) == pdg.tau_minus:
            # ignore status 3 taus in 2012 (something strange in the
            # MC record...)
            if child.status == 3:
                continue
            tau_decays.append(tautools.TauDecay(child))

    # check for incomplete tau decays
    for decay in tau_decays:
        if not decay.valid:
            print "invalid tau decay:"
            print decay
            decay.init.export_graphvis(
                    'decay_invalid_%d.dot' %
                    event.EventNumber)

    return resonance, tau_decays



class C3POProcessor(ATLASStudent):

    def __init__(self, options, **kwargs):

        super(C3POProcessor, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--student-verbose',
            dest='verbose',
            action='store_true', default=False)
        parser.add_argument('--syst-terms', default=None)
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
        verbose = self.args.verbose

        OutputModel = C3POEvent

        if datatype == datasets.MC:
            # only create truth branches for MC
            OutputModel += (
                    FourVectModel.prefix('resonance_') +
                    TrueTau.prefix('truetau1_') +
                    TrueTau.prefix('truetau2_'))

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

        # initialize the TreeChain of all input files
        # (each containing one tree named self.metadata.treename)
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
                'RunNumber',
                'EventNumber',
                'lbn']

        tree.set_buffer(
                chain._buffer,
                branches=copied_variables,
                create_branches=True,
                visible=False)

        chain.always_read(copied_variables)

        # set the event filters
        event_filters = EventFilterList([
            CoreFlags(
                count_funcs=count_funcs),
            TauSelected(2,
                count_funcs=count_funcs),
            TruthMatching(
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            MCWeight(
                datatype=datatype,
                tree=tree,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs)
        ])

        self.filters['event'] = event_filters

        chain._filters += event_filters

        define_objects(chain, year, skim=False)

        # define tree objects
        taus = [
            tree.define_object(name='tau1', prefix='tau1_'),
            tree.define_object(name='tau2', prefix='tau2_')]

        if datatype == datasets.MC:
            truetaus = [
                tree.define_object(name='truetau1', prefix='truetau1_'),
                tree.define_object(name='truetau2', prefix='truetau2_')]

            tree.define_object(name='resonance', prefix='resonance_')

        # entering the main event loop...
        for event in chain:

            # sort taus and jets in decreasing order by pT
            event.taus.sort(key=lambda tau: tau.pt, reverse=True)

            tau1, tau2 = event.taus

            # MET
            METx = event.MET.etx
            METy = event.MET.ety
            MET_vect = Vector2(METx, METy)
            MET = event.MET.et

            tree.MET = MET
            tree.MET_x = METx
            tree.MET_y = METy
            tree.MET_phi = event.MET.phi

            sumET = event.MET.sumet
            tree.sumET = sumET
            if sumET != 0:
                tree.MET_sig = ((2. * MET / GeV) /
                        (utils.sign(sumET) * sqrt(abs(sumET / GeV))))
            else:
                tree.MET_sig = -1.

            # use MMC values from skim
            mmc_mass = event.tau_MMC_mass
            mmc_resonance = event.tau_MMC_resonance
            mmc_met = Vector2(event.tau_MMC_MET_x, event.tau_MMC_MET_y)

            tree.mass_mmc_tau1_tau2 = mmc_mass
            tree.mmc_resonance.set_from(mmc_resonance)
            if mmc_mass > 0:
                tree.mmc_resonance_pt = mmc_resonance.Pt()
            tree.MET_mmc = mmc_met.Mod()
            tree.MET_mmc_x = mmc_met.X()
            tree.MET_mmc_y = mmc_met.Y()
            tree.MET_mmc_phi = math.pi - mmc_met.Phi()

            # truth matching
            if datatype == datasets.MC:

                resonance, tau_decays = get_taus(event)

                if resonance is not None:

                    FourVectModel.set(tree.resonance, resonance)

                    matched_taus = []
                    decays = tau_decays[:]
                    for itau, tau in enumerate(event.taus):
                        for idecay, tau_decay in enumerate(decays):
                            if tau.matches_vect(tau_decay.fourvect_visible):
                                tau_decay.matched = True
                                tau_decay.matched_object = tau
                                tau.matched = True
                                tau.matched_object = tau_decay
                                TrueTau.set(truetaus[itau], tau_decay,
                                        verbose=verbose)
                                decays.pop(idecay)
                                matched_taus.append(itau)
                                break

                    if len(decays) > 0:
                        for idecay, decay in enumerate(decays):
                            reco_idx = -1
                            remaining_idx = range(2)
                            for imatched in remaining_idx:
                                if imatched not in matched_taus:
                                    reco_idx = imatched
                                    remaining_idx.remove(imatched)
                                    break
                            TrueTau.set(truetaus[reco_idx], tau_decay,
                                    verbose=verbose)

                    if len(tau_decays) == 2:
                        # write truth met
                        fourvect_missing = (tau_decays[0].fourvect_missing +
                                            tau_decays[1].fourvect_missing)

                        tree.MET_true = fourvect_missing.Pt()
                        tree.MET_phi_true = fourvect_missing.Phi()
                        tree.MET_x_true = tree.MET_true * math.cos(tree.MET_phi_true)
                        tree.MET_y_true = tree.MET_true * math.sin(tree.MET_phi_true)

            # tau - vertex association
            tree.tau_same_vertex = (
                    tau1.privtx_x == tau2.privtx_x and
                    tau1.privtx_y == tau2.privtx_y and
                    tau1.privtx_z == tau2.privtx_z)

            # fill tau block
            for outtau, intau in zip(taus, event.taus):
                RecoTau.set(outtau, intau, verbose=verbose)

            # fill output ntuple
            tree.Fill(reset=True)

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

        if datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        merged_cutflow.Write()
