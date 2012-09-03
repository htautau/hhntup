#!/usr/bin/env python

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
        'range': (.55, 1),
        'cats': ['VBF', 'GGF', 'BOOSTED', 'PRESELECTION']
    },
    'tau2_BDTJetScore': {
        'title': r'$\tau_{2}$ BDT Score',
        'filename': 'tau2_BDTJetScore',
        'bins': 20,
        'range': (.55, 1),
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
        'range': (0.8, 3.0),
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

SYSTEMATICS = [
    ('TAUBDT_UP', 'TAUBDT_DOWN'),
    ('JES_UP,TES_UP', 'JES_DOWN,TES_DOWN'),
    ('JER_UP',),
    #('MFS_UP', 'MFS_DOWN'),
    #('ISOL_UP', 'ISOL_DOWN'),
    ('TRIGGER_UP', 'TRIGGER_DOWN'),
    ('FAKERATE_UP', 'FAKERATE_DOWN'),
]


if __name__ == '__main__':

    from argparse import ArgumentParser
    from categories import CATEGORIES

    parser = ArgumentParser()
    parser.add_argument('--no-cache', action='store_false', dest='use_cache',
            help="do not use cached background scale factors "
            "and instead recalculate them",
            default=True)
    parser.add_argument('--only-fit', action='store_true',
            help="only fit the Ztautau and QCD background, don't make plots",
            default=False)
    parser.add_argument('--no-systematics', action='store_false',
            dest='systematics',
            help="turn off systematics",
            default=True)

    parser.add_argument('--categories', nargs='*', default=CATEGORIES.keys())
    parser.add_argument('--plots', nargs='*')
    args = parser.parse_args()

    import ROOT
    ROOT.gROOT.SetBatch(True)
    from utils import *
    from matplotlib import cm
    from samples import *
    from background_estimation import qcd_ztautau_norm
    from config import plots_dir
    import os

    PLOTS_DIR = plots_dir(__file__)

    #ztautau   = MC_Ztautau(systematics=args.systematics)
    ztautau = Embedded_Ztautau(systematics=args.systematics)
    others = Others(systematics=args.systematics)

    higgs_125 = Higgs(
            masses=[125],
            systematics=args.systematics,
            scale=50,
            linecolor='red',
            linestyle='dashed')

    data = Data(markersize=2)

    figures = {}

    for category, cat_info in sorted(CATEGORIES.items(), key=lambda item: item[0]):

        if category not in args.categories:
            continue

        # QCD shape region SS or !OS
        qcd_shape_region = cat_info['qcd_shape_region']
        target_region = cat_info['target_region']

        qcd = QCD(data=data, mc=[others, ztautau],
              shape_region=qcd_shape_region)

        figures[category] = {}

        #cuts = Cut('80 < mass_mmc_tau1_tau2 < 110')
        cuts = Cut()

        qcd.scale = 1.
        ztautau.scale = 1.

        # determine normalization of QCD and Ztautau
        # in each category separately
        qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error = qcd_ztautau_norm(
            ztautau=ztautau,
            backgrounds=[others],
            data=data,
            category=category,
            target_region=target_region,
            qcd_shape_region=qcd_shape_region,
            use_cache=args.use_cache)

        if args.only_fit:
            continue

        qcd.scale = qcd_scale
        ztautau.scale = ztautau_scale

        for expr, var_info in VARIABLES.items():

            if category.upper() not in var_info['cats']:
                continue

            if args.plots and expr not in args.plots:
                continue

            print
            print "plotting %s in %s category" % (expr, category)

            bins = var_info['bins']
            min, max = var_info['range']

            if 'scale' in var_info:
                expr = "%s * %f" % (expr, var_info['scale'])

            other_hist = others.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            qcd_hist = qcd.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            print qcd_hist.systematics

            ztautau_hist = ztautau.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            bkg_hists = [qcd_hist, other_hist, ztautau_hist]

            data_hist = data.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            signal_hist = higgs_125.draw(
                    expr,
                    category, target_region,
                    bins, min, max,
                    cuts=cuts)

            print "Data events: %d" % sum(data_hist)
            print "Model events: %f" % sum(sum(bkg_hists))
            for hist in bkg_hists:
                print hist.GetTitle(), sum(hist)
            print "Data / Model: %f" % (sum(data_hist) / sum(sum(bkg_hists)))

            fig = draw(
                    data=data_hist, model=bkg_hists,
                    signal=signal_hist,
                    name=var_info['title'],
                    output_name=var_info['filename'],
                    category_name=cat_info['name'],
                    category=category,
                    units=var_info.get('units', None),
                    range=var_info['range'],
                    show_ratio=True,
                    show_qq=False,
                    model_colour_map=None,
                    dir=PLOTS_DIR,
                    systematics=SYSTEMATICS if args.systematics else None)
            figures[category][expr] = fig

    if set(args.categories) == set(CATEGORIES.keys()) and not args.plots:
        # only create multipage pdf of all plots if we created all plots
        from matplotlib.backends.backend_pdf import PdfPages
        import datetime
        now = datetime.datetime.today()
        # put all plots in a multipage pdf
        for category, exprs in figures.items():
            pdf = PdfPages(os.path.join(PLOTS_DIR, 'features_%s.pdf' % category))
            for expr, fig in sorted(exprs.items(), key=lambda x: x[0]):
                pdf.savefig(fig)
            d = pdf.infodict()
            d['Title'] = 'Features'
            d['Author'] = 'Noel Dawe'
            d['Subject'] = 'Higgs tautau hh features'
            d['Keywords'] = 'higgs tau'
            d['CreationDate'] = now
            d['ModDate'] = now
            pdf.close()
