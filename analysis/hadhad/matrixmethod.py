#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--normed', action='store_true', default=False)
parser.add_argument('--format', default='png')
args = parser.parse_args()

from rootpy.plotting import Hist, Canvas, HistStack
from rootpy.common import draw
from rootpy.io import open as ropen
from rootpy.tree import Tree, Cut
from rootpy.common import logon
import rootpy.root2matplotlib as rplt

import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rc
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42

from samples import *
import samples
from features import *

logon()

OUT_PREFIX = 'plots/matrixmethod'
FORMAT = args.format
NORMED = args.normed

for channel in samples.CHANNELS:
    print channel

    data = Data(channel=channel)
    mc_ztautau = MC_Ztautau(channel=channel)
    mc_ewk = MC_EWK(channel=channel)
    mc_ttbar = MC_TTbar(channel=channel)
    mc_singletop = MC_SingleTop(channel=channel)
    mc_diboson = MC_Diboson(channel=channel)

    #mc_dy = MC_DY(channel=channel)

    vbf = MC_VBF(channel=channel, mass=125, scale=10. if not NORMED else 1.)
    ggf = MC_ggF(channel=channel, mass=125, scale=10. if not NORMED else 1.)
    zh =  MC_ZH(channel=channel, mass=125, scale=10. if not NORMED else 1.)
    wh =  MC_WH(channel=channel, mass=125, scale=10. if not NORMED else 1.)

    qcd = QCD(data=data, mc=(mc_ztautau,
                             mc_ewk,
                             mc_ttbar,
                             mc_singletop,
                             mc_diboson))

    variables = samples.CHANNEL_VARIABLES[channel]

    variables += [
        'mass_mmc_tau1_tau2',
        'mass2_vis_tau1_tau2',
        'mass_collinear_tau1_tau2',
        'averageIntPerXing',
        'actualIntPerXing',
    ]

    for expr in variables:
        info = features.VARIABLES[expr]
        print expr
        bins = info['bins']
        min, max = info['range']

        if 'scale' in info:
            expr = "%s * %f" % (expr, info['scale'])

        if NORMED:
            fig = plt.figure(figsize=(8, 6), dpi=100)
            regions = ['OS-ID']
            base_index = 111
        else:
            fig = plt.figure(figsize=(16, 12), dpi=100)
            regions = ['OS-NOID',
                       'OS-ID',
                       'SS-NOID',
                       'SS-ID']
            base_index = 221
        qcd_hists = qcd.draw(expr, bins, min, max)

        for index, region in enumerate(regions):
            print region

            qcd_hist = qcd_hists[region]
            mc_ztautau_hist = sum(mc_ztautau.draw(expr, region, bins, min, max))
            #mc_dy_hist = sum(mc_dy.draw(expr, region, bins, min, max))
            mc_ewk_hist = sum(mc_ewk.draw(expr, region, bins, min, max))
            mc_ttbar_hist = sum(mc_ttbar.draw(expr, region, bins, min, max))
            mc_singletop_hist = sum(mc_singletop.draw(expr, region, bins, min, max))
            mc_diboson_hist = sum(mc_diboson.draw(expr, region, bins, min, max))

            data_hist = data.draw(expr, region, bins, min, max)

            vbf_hist = sum(vbf.draw(expr, region, bins, min, max))
            ggf_hist = sum(ggf.draw(expr, region, bins, min, max))
            zh_hist = sum(zh.draw(expr, region, bins, min, max))
            wh_hist = sum(wh.draw(expr, region, bins, min, max))

            background_hists = [
                    qcd_hist,
                    mc_ttbar_hist,
                    mc_singletop_hist,
                    mc_diboson_hist,
                    mc_ewk_hist,
                    mc_ztautau_hist,
                    #mc_dy_hist
            ]

            signal_hists = [
                    vbf_hist,
                    ggf_hist,
                    zh_hist,
                    wh_hist
            ]

            # set colour schemes
            for i, h in enumerate(background_hists):
                colour = cm.jet(1.*i/(len(background_hists)-1))
                h.SetFillColor(colour)
                h.SetLineColor(colour)

            for i, h in enumerate(signal_hists):
                colour = cm.summer(1.*i/(len(signal_hists)-1))
                h.SetFillColor(colour)
                h.SetLineColor(colour)

            if NORMED:

                for hatch, hist in zip(['\\','/'], signal_hists):
                    hist.SetFillStyle(hatch)
                    hist.SetLineColor(hist.GetFillColor())
                    hist.SetLineWidth(2)

                total_bkg = sum(sum(background_hists))
                for bg in background_hists:
                    bg /= total_bkg

                total_sig = sum(sum(signal_hists))
                for sig in signal_hists:
                    sig /= total_sig

                data_hist /= sum(data_hist)

            if NORMED:
                ax = fig.add_axes((0.16, 0.16, 0.8, 0.8))
            else:
                ax = plt.subplot(base_index + index)
            ax.get_frame().set_linewidth(2)
            rplt.bar(background_hists, linewidth=0, stacked=True, yerr='quadratic')
            if NORMED:
                rplt.bar(signal_hists, stacked=True, yerr='quadratic', fill=False, alpha=0.75)
            else:
                rplt.bar(signal_hists, stacked=True, yerr='quadratic', linewidth=0)
            rplt.errorbar(data_hist, fmt='o')
            if region == 'SS-ID': # legend
                label = '2jet selection' if channel == '2jet' else '01jet selection'
                if not NORMED:
                    l = plt.legend(loc='upper center', bbox_to_anchor=(-.1, 1.4))
                else:
                    l = plt.legend(loc='upper left')
                frame = l.get_frame()
                frame.set_alpha(.75)
                frame.set_linewidth(0)
                label = info['title']
                if 'units' in info:
                    label += ' [%s]' % info['units']
                ax.set_xlabel(label, fontsize=20, position=(1.0, 0.), ha='right')
            elif region == 'OS-ID':
                if not NORMED:
                    plt.title('%s ID' % samples.WORKING_POINT, fontsize=20)
                else:
                    label = info['title']
                    if 'units' in info:
                        label += ' [%s]' % info['units']
                    ax.set_xlabel(label, fontsize=20, position=(1.0, 0.), ha='right')
                    ax.set_ylabel('Normalized OS ID Events', fontsize=20, position=(0., 1.), va='top')
            elif region == 'SS-ID':
                ax.set_xlabel(info['title'], fontsize=20, position=(1.0, 0.), ha='right')
            elif region == 'OS-NOID':
                ax.set_ylabel('OS Events', fontsize=20, position=(0., 1.), va='top')
                plt.title('Failed %s ID' % samples.WORKING_POINT, fontsize=20)
            elif region == 'SS-NOID':
                ax.set_ylabel('SS Events', fontsize=20, position=(0., 1.), va='top')
            ax.yaxis.set_label_coords(-0.1, 1.)
            ax.xaxis.set_label_coords(1., -0.1)

        filename = '%s_%s' % (os.path.join(OUT_PREFIX, info['filename']), channel)
        if NORMED:
            filename += '_normed'
        fig.savefig('%s_%s.%s' % (filename, samples.WORKING_POINT, FORMAT),
                    bbox_inches='tight')
