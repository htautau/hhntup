import math

VARIABLES = {
    'averageIntPerXing': {
        'title': r'$\langle\mu\rangle|_{LB,BCID}$',
        'filename': 'averageIntPerXing',
        'bins': 20,
        'range': (1, 21),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'actualIntPerXing': {
        'title': r'$\langle\mu\rangle|_{LB}(BCID)$',
        'filename': 'actualIntPerXing',
        'bins': 20,
        'range': (1, 21),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'sum_pt': {
        'title': r'$\sum p_T$',
        'filename': 'sum_pt',
        'bins': 20,
        'range': (50, 450),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'sum_pt_full': {
        'title': r'$\sum p_T$ (all)',
        'filename': 'sum_pt_full',
        'bins': 20,
        'range': (50, 450),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'mmc_resonance_pt': {
        'title': r'MMC Resonance $p_T$',
        'filename': 'mmc_resonance_pt',
        'bins': 20,
        'range': (0, 200),
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'numJets': {
        'title': r'Number of Jets with $p_T>25$ GeV',
        'filename': 'numjets',
        'bins': 7,
        'range': (-.5, 6.5),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET': {
        'title': r'$E^{miss}_{T}$',
        'filename': 'MET',
        'bins': 20,
        'range': (0, 100),
        'scale': 1./1000,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_x': {
        'title': r'$E^{miss}_{T_{x}}$',
        'filename': 'MET_x',
        'bins': 20,
        'range': (-75, 75),
        'scale': 1./1000,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_y': {
        'title': r'$E^{miss}_{T_{y}}$',
        'filename': 'MET_y',
        'bins': 20,
        'range': (-75, 75),
        'scale': 1./1000,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_phi': {
        'title': r'$E^{miss}_{T} \phi$',
        'filename': 'MET_phi',
        'bins': 20,
        'range': (-math.pi, math.pi),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_mmc': {
        'title': r'$E^{miss}_{T}$ MMC',
        'filename': 'MET_mmc',
        'bins': 20,
        'range': (0, 100),
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_mmc_x': {
        'title': r'$E^{miss}_{T_{x}}$ MMC',
        'filename': 'MET_mmc_x',
        'bins': 20,
        'range': (-75, 75),
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_mmc_y': {
        'title': r'$E^{miss}_{T_{y}}$ MMC',
        'filename': 'MET_mmc_y',
        'bins': 20,
        'range': (-75, 75),
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'MET_mmc_vec.Phi()': {
        'title': r'$E^{miss}_{T} \phi$ MMC',
        'filename': 'MET_mmc_phi',
        'bins': 20,
        'range': (-math.pi, math.pi),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'sphericity': {
        'title': r'sphericity',
        'filename': 'sphericity',
        'bins': 20,
        'range': (0, 1),
        'cats': ['VBF', 'BOOSTED']
    },
    'sphericity_full': {
        'title': r'sphericity (all)',
        'filename': 'sphericity_full',
        'bins': 20,
        'range': (0, 1),
        'cats': ['VBF', 'BOOSTED']
    },
    'aplanarity': {
        'title': r'aplanarity',
        'filename': 'aplanarity',
        'bins': 20,
        'range': (0, .15),
        'cats': ['VBF', 'BOOSTED']
    },
    'aplanarity_full': {
        'title': r'aplanarity (all)',
        'filename': 'aplanarity_full',
        'bins': 20,
        'range': (0, .15),
        'cats': ['VBF', 'BOOSTED']
    },
    'MET_centrality': {
        'title': r'$E^{miss}_{T}$ Centrality',
        'filename': 'met_centrality',
        'bins': 20,
        'range': (-math.sqrt(2), math.sqrt(2)),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'mass2_vis_tau1_tau2': {
        'title': r'$M^{vis}_{\tau_{1},\/\tau_{2}}$',
        'filename': 'mass_vis',
        'bins': 20,
        'range': (0, 250),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'mass_mmc_tau1_tau2': {
        'title': r'$M^{MMC}_{\tau_{1},\/\tau_{2}}$',
        'filename': 'mass_MMC',
        'bins': 20,
        'range': (50, 250),
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'mass_collinear_tau1_tau2': {
        'title': r'$M^{col}_{\tau_{1},\/\tau_{2}}$',
        'filename': 'mass_collinear',
        'bins': 20,
        'range': (0, 250),
        'units': 'GeV',
        'scale': 0.001,
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_fourvect.Pt()': {
        'title': r'$p_{T_{\tau_{1}}}$',
        'filename': 'tau1_pt',
        'bins': 20,
        'range': (20, 100),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_fourvect.Pt()': {
        'title': r'$p_{T_{\tau_{2}}}$',
        'filename': 'tau2_pt',
        'bins': 20,
        'range': (20, 100),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_numTrack': {
        'title': r'$\tau_{1}$ Number of Tracks',
        'filename': 'tau1_numTrack',
        'bins': 5,
        'range': (-.5, 4.5),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_numTrack': {
        'title': r'$\tau_{2}$ Number of Tracks',
        'filename': 'tau2_numTrack',
        'bins': 5,
        'range': (-.5, 4.5),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_x': {
        'title': r'$\tau_{1_{x}}$',
        'filename': 'tau1_x',
        'bins': 20,
        'range': (-3, 4),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_x': {
        'title': r'$\tau_{2_{x}}$',
        'filename': 'tau2_x',
        'bins': 20,
        'range': (-3, 4),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_jvtxf': {
        'title': r'$\tau_{1}$ JVF',
        'filename': 'tau1_jvf',
        'bins': 20,
        'range': (0, 1),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_jvtxf': {
        'title': r'$\tau_{2}$ JVF',
        'filename': 'tau2_jvf',
        'bins': 20,
        'range': (0, 1),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_BDTJetScore': {
        'title': r'$\tau_{1}$ BDT Score',
        'filename': 'tau1_BDTJetScore',
        'bins': 20,
        'range': (.55, 1.0001),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_BDTJetScore': {
        'title': r'$\tau_{2}$ BDT Score',
        'filename': 'tau2_BDTJetScore',
        'bins': 20,
        'range': (.55, 1.0001),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'cos_theta_tau1_tau2': {
        'title': r'$\cos(\alpha_{\tau_{1},\/\tau_{2}})$',
        'filename': 'cos_theta_tau1_tau2',
        'bins': 20,
        'range': (-1, 1),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'theta_tau1_tau2': {
        'title': r'$\alpha_{\tau_{1},\/\tau_{2}}$',
        'filename': 'theta_tau1_tau2',
        'bins': 20,
        'range': (0, math.pi),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'dR_tau1_tau2': {
        'title': r'$\Delta R_{\tau_{1},\/\tau_{2}}$',
        'filename': 'dr_tau1_tau2',
        'bins': 20,
        'range': (0., 6.),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'fabs(dPhi_tau1_tau2)': {
        'title': r'$\Delta \phi_{\tau_{1},\/\tau_{2}}$',
        'filename': 'dphi_tau1_tau2',
        'bins': 20,
        'range': (0., math.pi),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    ('tau1_fourvect.Eta()', 'tau2_fourvect.Eta()'): {
        'title': r'$\eta_{\tau_{1,\/2}}$',
        'filename': 'tau_eta',
        'bins': 20,
        'range': (-3, 3),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_charge': {
        'title': r'$\tau_1$ Charge',
        'filename': 'tau1_charge',
        'bins': 5,
        'range': (-2.5, 2.5),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_charge': {
        'title': r'$\tau_2$ Charge',
        'filename': 'tau2_charge',
        'bins': 5,
        'range': (-2.5, 2.5),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau1_centrality': {
        'title': r'$\tau_1$ Centrality',
        'filename': 'tau1_centrality',
        'bins': 20,
        'range': (0, 1),
        'cats': ['VBF']
    },
    'tau2_centrality': {
        'title': r'$\tau_2$ Centrality',
        'filename': 'tau2_centrality',
        'bins': 20,
        'range': (0, 1),
        'cats': ['VBF']
    },
    'tau1_fourvect.Eta()': {
        'title': r'$\eta_{\tau_{1}}$',
        'filename': 'tau1_eta',
        'bins': 20,
        'range': (-3, 3),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_fourvect.Eta()': {
        'title': r'$\eta_{\tau_{2}}$',
        'filename': 'tau2_eta',
        'bins': 20,
        'range': (-3, 3),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'jet1_fourvect.Eta()': {
        'title': r'$\eta_{jet_{1}}$',
        'filename': 'jet1_eta',
        'bins': 20,
        'range': (-5, 5),
        'cats': ['VBF', 'BOOSTED']
    },
    'jet2_fourvect.Eta()': {
        'title': r'$\eta_{jet_{2}}$',
        'filename': 'jet2_eta',
        'bins': 20,
        'range': (-5, 5),
        'cats': ['VBF']
    },
    'jet1_fourvect.Pt()': {
        'title': r'$p_{T_{jet_{1}}}$',
        'filename': 'jet1_pt',
        'bins': 20,
        'range': (20, 160),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF','BOOSTED']
    },
    'jet2_fourvect.Pt()': {
        'title': r'$p_{T_{jet_{2}}}$',
        'filename': 'jet2_pt',
        'bins': 20,
        'range': (20, 160),
        'scale': 0.001,
        'units': 'GeV',
        'cats': ['VBF']
    },
    'jet1_fourvect_boosted.Eta()': {
        'title': r'boosted $\eta_{jet_{1}}$',
        'filename': 'jet1_eta_boosted',
        'bins': 20,
        'range': (-5, 5),
        'cats': ['VBF']
    },
    'jet2_fourvect_boosted.Eta()': {
        'title': r'boosted $\eta_{jet_{2}}$',
        'filename': 'jet2_eta_boosted',
        'bins': 20,
        'range': (-5, 5),
        'cats': ['VBF']
    },
    ('jet1_fourvect.Eta()', 'jet2_fourvect.Eta()'): {
        'title': r'$\eta_{jet_{1,\/2}}$',
        'filename': 'jet_eta',
        'bins': 20,
        'range': (-5, 5),
        'cats': ['VBF']
    },
    ('jet1_fourvect_boosted.Eta()', 'jet2_fourvect_boosted.Eta()'): {
        'title': r'boosted $\eta_{jet_{1,\/2}}$',
        'filename': 'jet_eta_boosted',
        'bins': 20,
        'range': (-5, 5),
        'cats': ['VBF']
    },
    'dEta_jets': {
        'title': r'$|\Delta\eta_{jet_{1},\/jet_{2}}|$',
        'filename': 'dEta_jets',
        'bins': 20,
        'range': (0, 6),
        'cats': ['VBF']
    },
    'dEta_jets_boosted': {
        'title': r'boosted $|\Delta\eta_{jet_{1},\/jet_{2}}|$',
        'filename': 'dEta_jets_boosted',
        'bins': 20,
        'range': (0, 6),
        'cats': ['VBF']
    },
    'eta_product_jets': {
        'title': r'$\eta_{jet_{1}} \times \/ \eta_{jet_{2}}$',
        'filename': 'eta_product_jets',
        'bins': 20,
        'range': (-10, 10),
        'cats': ['VBF']
    },
    'mass_jet1_jet2': {
        'title': r'$M_{jet_{1},\/jet_{2}}$',
        'filename': 'M_jet1_jet2',
        'bins': 21,
        'range': (0, 600),
        'scale': 1./1000,
        'units': 'GeV',
        'cats': ['VBF']
    },
}


