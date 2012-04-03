from rootpy.tree import TreeModel
from rootpy.types import *


class CurrentList(TreeModel):

    # tau_numTrack
    NUMTRACK = IntCol()
    # evt_calcVars_numGoodVertices
    NUM_PILEUP_AND_PRIMARY_VERTICES = IntCol()
    # tau_seedCalo_centFrac
    CENTFRAC = FloatCol()
    # tau_etOverPtLeadTrk
    ETOVERPTLEADTRK = FloatCol()
    # tau_seedCalo_trkAvgDist
    TRKAVGDIST = FloatCol()
    # tau_effTopoInvMass
    EFFTOPOINVMASS = FloatCol()
    # tau_ipSigLeadTrk
    IPSIGLEADTRK = FloatCol()
    # tau_seedCalo_lead3ClusterEOverAllClusterE
    LEAD3CLUSTEREOVERALLCLUSTERE = FloatCol()
    # tau_seedCalo_wideTrk_n
    NUMWIDETRACK = IntCol()
    # tau_calcVars_calRadius
    CALRADIUS = FloatCol()
    # tau_calcVars_drMax
    DRMAX = FloatCol()
    # tau_massTrkSys
    MASSTRKSYS = FloatCol()
    # tau_trFlightPathSig
    TRFLIGHTPATHSIG = FloatCol()
