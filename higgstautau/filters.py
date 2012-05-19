"""
Event filters common to both hadhad and lephad go here
"""

from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *

from . import jetcleaning


class PriVertex(EventFilter):

    def passes(self, event):

        return any(ifilter(lambda vxp: vxp.type == 1 and vxp.nTracks >= 4, event.vertices))


class JetCleaning(EventFilter):

    def __init__(self, level=jetcleaning.LOOSER,
                 pt_thresh=20*GeV,
                 eta_max=4.5,
                 **kwargs):

        super(JetCleaning, self).__init__(**kwargs)
        self.level = level
        self.pt_thresh = pt_thresh
        self.eta_max = eta_max

    def passes(self, event):

        for jet in event.jets:

            if jet.pt <= self.pt_thresh or abs(jet.eta) >= self.eta_max: continue
            LArQmean = jet.AverageLArQF / 65535.0
            chf = jet.sumPtTrk / jet.pt
            if jetcleaning.is_bad(level=self.level,
                                  quality=jet.LArQuality,
                                  NegE=jet.NegativeE,
                                  emf=jet.emfrac,
                                  hecf=jet.hecf,
                                  time=jet.Timing,
                                  fmax=jet.fracSamplingMax,
                                  eta=jet.emscale_eta,
                                  chf=chf,
                                  HecQ=jet.HECQuality,
                                  LArQmean=LArQmean
                                 ): return False
        return True


class LArError(EventFilter):

    def passes(self, event):

        return event.larError <= 1


def in_lar_hole(eta, phi):

    return (-0.2 < eta < 1.6) and (-0.988 < phi < -0.392)


class LArHole(EventFilter):

    def __init__(self, datatype, **kwargs):

        super(LArHole, self).__init__(**kwargs)
        if datatype == datasets.DATA:
            self.passes = self.passes_data
        else:
            self.passes = self.passes_mc

    def passes_data(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        for jet in event.jets:
            if not jet.pt > 20 * GeV * (1 - jet.BCH_CORR_JET) / (1 - jet.BCH_CORR_CELL): continue
            if in_lar_hole(jet.eta, jet.phi): return False
        return True

    def passes_mc(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        for jet in event.jets:
            if not jet.pt > 20 * GeV: continue
            if in_lar_hole(jet.eta, jet.phi): return False
        return True
