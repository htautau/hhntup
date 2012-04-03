from rootpy.tree import TreeModel
from rootpy.types import *

class TruthTau(TreeModel):

    hadronic = BoolCol(default=False)
    nprong = IntCol(default=-1111)
    npi0 = IntCol(default=-1111)
    nneutrals = IntCol(default=-1111)
    pt = FloatCol(default=0)
    eta = FloatCol(default=-1111)
    phi = FloatCol(default=-1111)
    pt_vis = FloatCol(default=0)
    eta_vis = FloatCol(default=-1111)
    phi_vis = FloatCol(default=-1111)
    nu_pt = FloatCol(default=0)
    nu_eta = FloatCol(default=-1111)
    nu_phi = FloatCol(default=-1111)
    dR_tau_nu = FloatCol(default=-1111)
    dTheta3d_tau_nu = FloatCol(default=-1111)

reco_variables = (
    ("BDTJetScore", "F"),
    ("BDTEleScore", "F"),
    ("nPi0", "I"),
    ("pt", "F"),
    ("seedCalo_eta", "F"),
    ("seedCalo_phi", "F"),
    ("seedCalo_numTrack", "I"),
    ("charge", "I"),
)

truth_variables = (
    ("pt", "F"),
    ("m", "F"),
    ("eta", "F"),
    ("phi", "F"),
    ("vis_m", "F"),
    ("vis_eta", "F"),
    ("vis_phi", "F"),
    ("vis_Et", "F"),
    ("nProng", "I"),
    ("nPi0", "I"),
    ("charge", "I"),
)

common_variables = (
    ("matched", "I"),
    ("matched_dR", "F"),
    ("matched_collision", "B"),
)

# define branches for output ntuple
variables = [
    ("Mvis_tau1_tau2","F"),
    ("numJets","I"),
    ("jetDeltaEta", "F"),
    ("jetDeltaEta_matched", "F"),
    ("jet_AntiKt4TopoEM_matched", "VI"),
    ("jet_AntiKt4TopoEM_matched_dR", "VF"),
    ("numVertices","I"),
    ("MET","F"),
    ("MET_phi","F"),
    ("HT","F"),
    ("MMC_mass","F"),
    ("error", "B"),
    ("mu", "I"),
    ("selected", "B")
]

jet_variables = (
    ("E", "F"),
    ("pt", "F"),
    ("m", "F"),
    ("eta", "F"),
    ("phi", "F"),
    ("jvtxf", "F")
)

jet_extra_variables = (
    ("Et", "F"),
)

jet_matched_variables = (
    ("matched", "B"),
    ("matched_dR", "F"),
    ("matched_pdgId", "I")
)

parton_variables = (
    ("pt", "F"),
    ("eta", "F"),
    ("phi", "F"),
    ("pdgId", "I"),
)

parton_extra_variables = (
    ("Et", "F"),
)
