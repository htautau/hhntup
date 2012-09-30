from rootpy.tree import Cut
from higgstautau.hadhad import categories

TAU1_MEDIUM = Cut('tau1_JetBDTSigMedium==1')
TAU2_MEDIUM = Cut('tau2_JetBDTSigMedium==1')
TAU1_TIGHT = Cut('tau1_JetBDTSigTight==1')
TAU2_TIGHT = Cut('tau2_JetBDTSigTight==1')

TAU1_CENTRAL = Cut('fabs(tau1_fourvect.Eta()) < 1.5')
TAU1_FORWARD = Cut('fabs(tau1_fourvect.Eta()) > 1.5')
TAU2_CENTRAL = Cut('fabs(tau2_fourvect.Eta()) < 1.5')
TAU2_FORWARD = Cut('fabs(tau2_fourvect.Eta()) > 1.5')

ID_MEDIUM = TAU1_MEDIUM & TAU2_MEDIUM
ID_TIGHT = TAU1_TIGHT & TAU2_TIGHT
ID_MEDIUM_TIGHT = (TAU1_MEDIUM & TAU2_TIGHT) | (TAU1_TIGHT & TAU2_MEDIUM)

Z_PEAK = Cut('80 < mass_mmc_tau1_tau2 < 120')
ID_MEDIUM_FORWARD_TIGHT_CENTRAL = (
        (TAU1_MEDIUM & TAU1_FORWARD & TAU2_TIGHT & TAU2_CENTRAL) |
        (TAU1_TIGHT & TAU1_CENTRAL & TAU2_MEDIUM & TAU2_FORWARD))

# low cut fixes mass, high cut removes QCD
DR_FIX = Cut('1.0 < dR_tau1_tau2 < 3.0')
MASS_FIX = Cut('mass_mmc_tau1_tau2 > 80')
MAX_NJET = Cut('numJets <= 3')
MET = Cut('MET > 20000')

COMMON_CUTS = MET & MASS_FIX & MAX_NJET # & DR_FIX

CATEGORIES = {
    'vbf': {
        'name': r'$\tau_{had}\tau_{had}$: VBF Category',
        'code': categories.CATEGORY_VBF,
        'cuts': COMMON_CUTS & ID_MEDIUM,
        'fitbins': 5,
        'qcd_shape_region': 'SS',
        'target_region': 'OS',
        'features': [
            'dEta_jets',
            #'dEta_jets_boosted',
            'eta_product_jets',
            #'eta_product_jets_boosted',
            'mass_jet1_jet2',
            #'sphericity',
            #'aplanarity',
            'tau1_centrality',
            'tau2_centrality',
            #'sphericity_boosted',
            #'aplanarity_boosted',
            #'tau1_centrality_boosted',
            #'tau2_centrality_boosted',
            #'cos_theta_tau1_tau2',
            'dR_tau1_tau2',
            'tau1_BDTJetScore',
            'tau2_BDTJetScore',
            #'tau1_x',
            #'tau2_x',
            'MET_centrality',
            'mmc_resonance_pt',
            #'sum_pt',
        ]
    },
    'boosted': {
        'name': r'$\tau_{had}\tau_{had}$: Boosted Category',
        'code': categories.CATEGORY_BOOSTED,
        'cuts': COMMON_CUTS & ID_MEDIUM,
        'fitbins': 5,
        'qcd_shape_region': 'SS',
        'target_region': 'OS',
        'features': [
            'sphericity',
            #'aplanarity',
            #'cos_theta_tau1_tau2',
            'dR_tau1_tau2',
            'tau1_BDTJetScore',
            'tau2_BDTJetScore',
            'tau1_x',
            'tau2_x',
            'MET_centrality',
            #'sum_pt',
            'mmc_resonance_pt',
        ]
    },
    'ggf': {
        'name': r'$\tau_{had}\tau_{had}$: Non-Boosted Category',
        'code': categories.CATEGORY_GGF,
        'cuts': COMMON_CUTS & ID_MEDIUM,
        'fitbins': 8,
        'qcd_shape_region': 'SS',
        'target_region': 'OS',
        'features': [
            #'cos_theta_tau1_tau2',
            'dR_tau1_tau2',
            'tau1_BDTJetScore',
            'tau2_BDTJetScore',
            'tau1_x',
            'tau2_x',
            'MET_centrality',
            #'mmc_resonance_pt',
        ]
    },
}

CONTROLS = {
    'preselection': {
        'name': r'$\tau_{had}\tau_{had}$: At Preselection',
        'code': None,
        'cuts': COMMON_CUTS & ID_MEDIUM,
        'fitbins': 10,
        'qcd_shape_region': 'SS',
        'target_region': 'OS',
    },
    'z': {
        'name': r'$\tau_{had}\tau_{had}$: Z Control Region',
        'code': None,
        'cuts': MET & Cut('dR_tau1_tau2<2.8') & Z_PEAK & ID_MEDIUM_TIGHT,
        'fitbins': 8,
        'qcd_shape_region': 'SS',
        'target_region': 'OS',
    }
}
