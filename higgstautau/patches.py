from rootpy.tree.filtering import *
from rootpy.io import root_open as ropen
from math import *
import ROOT

#################################################
# Electron ID recalculation (p1130 tau D3PDs)
#################################################

from externaltools import egammaAnalysisUtils
from ROOT import isTightPlusPlus, isMediumPlusPlus, isLoosePlusPlus, egammaMenu#, egammaPID


class egammaPID(object):

    ConversionMatch_Electron = 1


class ElectronIDpatch(EventFilter):
    """
    Recalculates the electron ID based on electron variables
    """
    def passes(self, event):

        #iterate over all electron candidates
        for e in event.electrons:

            #Retrieve relevant branches
            eta          = e.etas2
            if eta == -999:
                continue
            eT           = e.cl_E/cosh(eta)
            rHad         = e.Ethad/eT
            rHad1        = e.Ethad1/eT
            Reta         = e.reta
            w2           = e.weta2
            f1           = e.f1
            f3           = e.f3
            wstot        = e.wstot
            DEmaxs1      = 0
            if (e.emaxs1 + e.Emax2) != 0:
                DEmaxs1      = (e.emaxs1 - e.Emax2)/(e.emaxs1 + e.Emax2)
            deltaEta     = e.deltaeta1
            deltaPhi     = e.deltaphi2
            eOverp       = e.cl_E*abs(e.trackqoverp)
            d0           = e.trackd0_physics
            TRratio      = e.TRTHighTOutliersRatio
            nTRT         = e.nTRTHits
            nTRTOutliers = e.nTRTOutliers
            nSi          = e.nSiHits
            nSiOutliers  = e.nSCTOutliers + e.nPixelOutliers
            nPix         = e.nPixHits
            nPixOutliers = e.nPixelOutliers
            nBlayer      = e.nBLHits
            nBlayerOutliers = e.nBLayerOutliers
            expectBlayer = bool(e.expectHitInBLayer)
            ConvBit      = e.isEM and (1 << egammaPID.ConversionMatch_Electron)

            #Correct the loosePP flag
            e.loosePP = isLoosePlusPlus(
                    eta, eT, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1,
                    deltaEta, nSi, nSiOutliers, nPix, nPixOutliers)

            #Correct the mediumPP flag
            e.mediumPP  = isMediumPlusPlus(
                    eta, eT, f3, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1,
                    deltaEta, d0, TRratio, nTRT, nTRTOutliers,
                    nSi, nSiOutliers, nPix, nPixOutliers, nBlayer,
                    nBlayerOutliers, expectBlayer)

            #Correct the tightPP flag
            e.tightPP = isTightPlusPlus(
                    eta, eT, f3, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1,
                    deltaEta, d0, TRratio, nTRT, nTRTOutliers,
                    nSi, nSiOutliers, nPix, nPixOutliers, nBlayer,
                    nBlayerOutliers, expectBlayer, eOverp, deltaPhi, ConvBit)

        return True



#################################################
# Tau ID recalculation
#################################################
from .tauid.p851.selection import selection as selection_2011
from .tauid.p1130.selection import selection as selection_2012

class TauIDpatch(EventFilter):
    """
    Recalculates the tau ID
    """
    def __init__(self, year, **kwargs):

        super(TauIDpatch, self).__init__(**kwargs)

        if year == 2011:
            self.passes = self.passes_2011
        elif year == 2012:
            self.loose_1p   = selection_2012('loose', 1)
            self.medium_1p  = selection_2012('medium', 1)
            self.tight_1p   = selection_2012('tight', 1)
            self.loose_3p   = selection_2012('loose', 3)
            self.medium_3p  = selection_2012('medium', 3)
            self.tight_3p   = selection_2012('tight', 3)
            self.passes = self.passes_2012
        else:
            raise ValueError("No tauid patch defined for year %d" % year)

    def passes_2011(self, event):

        #nvtx = event.number_of_good_vertices
        # assume vertex selection already applied!
        nvtx = len(event.vertices)

        for tau in event.taus:

            pt = tau.pt
            ntrack = tau.numTrack

            loose = selection_2011('loose', ntrack, nvtx).Eval(pt)
            medium = selection_2011('medium', ntrack, nvtx).Eval(pt)
            tight = selection_2011('tight', ntrack, nvtx).Eval(pt)

            tau.JetBDTSigLoose  = (tau.BDTJetScore > loose)
            tau.JetBDTSigMedium = (tau.BDTJetScore > medium)
            tau.JetBDTSigTight  = (tau.BDTJetScore > tight)

        return True

    def passes_2012(self, event):

        for tau in event.taus:
            pt = tau.pt

            if tau.numTrack <= 1:
                cut_loose  = self.loose_1p.Eval(pt)
                cut_medium = self.medium_1p.Eval(pt)
                cut_tight  = self.tight_1p.Eval(pt)
            else:
                cut_loose  = self.loose_3p.Eval(pt)
                cut_medium = self.medium_3p.Eval(pt)
                cut_tight  = self.tight_3p.Eval(pt)

            tau.JetBDTSigLoose  = (tau.BDTJetScore > cut_loose)
            tau.JetBDTSigMedium = (tau.BDTJetScore > cut_medium)
            tau.JetBDTSigTight  = (tau.BDTJetScore > cut_tight)

        return True
