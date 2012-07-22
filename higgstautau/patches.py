from rootpy.tree.filtering import *
from math import *
import ROOT

#################################################
# Electron ID recalculation (p1130 tau D3PDs)
#################################################

from externaltools import egammaAnalysisUtils
from ROOT import isTightPlusPlus, isMediumPlusPlus, egammaMenu

class ElectronIDpatch(EventFilter):
    """
    Recalculates the electron ID based on electron variables
    """

    def passes(self, event):

        #iterate over all electron candidates
        for e in event.electrons:

            #Retrieve relevant branches
            eta          = e.etas2
            if eta == -999: continue
            eT           = e.cl_E/cosh(eta)
            rHad         = e.Ethad/eT
            rHad1        = e.Ethad1/eT
            Reta         = e.reta
            w2           = e.weta2
            f1           = e.f1
            f3           = e.f3
            wstot        = e.wstot
            DEmaxs1      = -1
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
            ConvBit      = e.isEM

            #Correct the mediumPP flag
            e.mediumPP  = isMediumPlusPlus(eta, eT, f3, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1, deltaEta, d0, TRratio, nTRT, nTRTOutliers,
                                           nSi, nSiOutliers, nPix, nPixOutliers, nBlayer, nBlayerOutliers, expectBlayer)


            #Correct the tightPP flag
            e.tightPP = isTightPlusPlus(eta, eT, f3, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1, deltaEta, d0, TRratio, nTRT, nTRTOutliers,
                                        nSi, nSiOutliers, nPix, nPixOutliers, nBlayer, nBlayerOutliers, expectBlayer, eOverp, deltaPhi, ConvBit)
            
        return True
    


#################################################
# Tau ID recalculation (p1130 tau D3PDs)
#################################################

class TauIDpatch(EventFilter):
    """
    Recalculates the tau ID
    """

    def __init__(self, graph, **kwargs):

        super(TauIDpatch, self).__init__(**kwargs)

        #Load TGraphs

        f = ROOT.TFile(graph)

        self.loose_1p   = f.Get('loose_1p')
        self.medium_1p  = f.Get('medium_1p')
        self.tight_1p   = f.Get('tight_1p')
        self.loose_3p   = f.Get('loose_3p')
        self.medium_3p  = f.Get('medium_3p')
        self.tight_3p   = f.Get('tight_3p')


    def passes(self, event):

        for tau in event.taus:
            pt = tau.pt

            cut_loose  = 0
            cut_medium = 0
            cut_tight  = 0
            
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
