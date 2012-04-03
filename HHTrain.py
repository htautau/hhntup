import ROOT
import math
from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from atlastools.batch import ATLASStudent
from rootpy.tree import Tree, TreeChain, TreeModel
from rootpy.types import *
from rootpy.io import open as ropen
import numpy as np


ROOT.gErrorIgnoreLevel = ROOT.kFatal


class ExtraVariables(TreeModel):

    sphericity = FloatCol()
    aplanarity = FloatCol()


class HHTrain(ATLASStudent):

    def work(self):

        # initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                          files=self.files,
                          events=self.events)

        outtree = Tree(name=self.metadata.treename,
                       file=self.output,
                       model=ExtraVariables)
        outtree.set_buffer(intree.buffer, create_branches=True, visible=False)

        # define tree objects
        intree.define_object(name='tau1', prefix='tau1_')
        intree.define_object(name='tau2', prefix='tau2_')
        intree.define_object(name='jet1', prefix='jet1_')
        intree.define_object(name='jet2', prefix='jet2_')

        # entering the main event loop...
        for event in intree:

            # require two additional jets
            if event.jet2.fourvect.P() < 20 * GeV:
                continue

            # sphericity tensor
            norm = event.tau1.fourvect.P()**2 + \
                   event.tau2.fourvect.P()**2 + \
                   event.jet1.fourvect.P()**2 + \
                   event.jet2.fourvect.P()**2
            S = np.zeros(shape=(3,3))
            for i in xrange(3):
                for j in xrange(3):
                    S[i][j] = event.tau1.fourvect.Vect()[i] * event.tau1.fourvect.Vect()[j] + \
                              event.tau2.fourvect.Vect()[i] * event.tau2.fourvect.Vect()[j] + \
                              event.jet1.fourvect.Vect()[i] * event.jet1.fourvect.Vect()[j] + \
                              event.jet2.fourvect.Vect()[i] * event.jet2.fourvect.Vect()[j]
            S /= norm
            eigvals, eigvects = np.linalg.eig(S)
            eigvals = sorted(eigvals)
            print eigvals

            # sphericity
            outtree.sphericity = (eigvals[0] + eigvals[1]) * 1.5

            # aplanarity
            outtree.aplanarity = eigvals[0] * 1.5

            outtree.Fill()

        self.output.cd()

        # flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()
