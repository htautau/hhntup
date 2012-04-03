import ROOT
import math
from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.batch import ATLASStudent
from rootpy.tree import Tree, TreeBuffer, TreeChain
from mixins import MCParticle
import tautools
from models import *

ROOT.gErrorIgnoreLevel = ROOT.kFatal


class MCTauProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def work(self):
         
        # initialize the TreeChain of all input files (each containing one tree named self.fileset.treename)
        tree = TreeChain(self.fileset.treename,
                         files=self.fileset.files,
                         events=self.events,
                         usecache=True,
                         cache_size=10000000,
                         learn_entries=30)
        
        # create output tree
        self.output.cd()
        
        # this tree will contain info pertaining to true tau decays
        # for possible use in the optimization of a missing mass calculator
        mc_tree = Tree(name = "tau_mc", model=TrueTau_MCBlock)

        tree.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
         
        for event in tree:
            
            """
            Need to get all MC tau final states to build ntuple for missing mass calculator pdfs
            """ 
            if self.fileset.datatype == datasets.MC:
                tau_decays = tautools.get_tau_decays(event)
                for decay in tau_decays:
                    hadronic = decay.hadronic
                    if hadronic:
                        mc_tree.hadronic = True
                        mc_tree.nprong = decay.nprong
                        mc_tree.npi0 = decay.npi0
                        mc_tree.nneutrals = decay.nneutrals
                        mc_tree.fourvect.set_from(decay.fourvect)
                        mc_tree.fourvect_vis.set_from(decay.fourvect_visible)
                        mc_tree.fourvect_miss.set_from(decay.fourvect_missing)
                        mc_tree.dR_tau_nu = decay.dR_tau_nu
                        mc_tree.dTheta3d_tau_nu = decay.dTheta3d_tau_nu
                        mc_tree.Fill(reset=True)
