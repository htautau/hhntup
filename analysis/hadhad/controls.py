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
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt

from matplotlib import rc
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42

import rootpy.root2matplotlib as rplt
from samples import *
import samples
from features import *
from rootpy.common import logon

logon()

OUT_PREFIX = 'plots/controls'
FORMAT = args.format
NORMED = args.normed

channel = '2jet'

data = Data(channel=channel)
mc_ztautau = MC_Ztautau(channel=channel)
mc_ewk = MC_EWK(channel=channel)
#mc_dy = MC_DY(channel=channel)

vbf = MC_VBF_125(channel=channel, scale=10. if not NORMED else 1.)
ggf = MC_ggF_125(channel=channel, scale=10. if not NORMED else 1.)

variables = samples.CHANNEL_VARIABLES[channel]

variables += [
    'mass_mmc_tau1_tau2',
    'mass2_vis_tau1_tau2',
    'mass_collinear_tau1_tau2',
    'averageIntPerXing',
    'actualIntPerXing',
    'MET'
]


CONTROLS = {
    'QCD': Cut('MET<10000') & samples.ID,
    'Z':   Cut('MET>30000 && mass_mmc_tau1_tau2 < 100') & samples.ID,
}


for control, cuts in CONTROLS.items():
    for expr in variables:
        info = features.VARIABLES[expr]
        print expr
        bins = info['bins']
        bins = 10
        min, max = info['range']

        if 'scale' in info:
            expr = "%s * %f" % (expr, info['scale'])

        fig = plt.figure(figsize=(11, 6), dpi=100)
        base_index = 111

        mc_ztautau_hist = sum(mc_ztautau.draw(expr, cuts & samples.OS, bins, min, max))
        #mc_dy_hist = sum(mc_dy.draw(expr, region, bins, min, max))
        mc_ewk_hist = sum(mc_ewk.draw(expr, cuts & samples.OS, bins, min, max))
        data_hist = data.draw(expr, cuts & samples.OS, bins, min, max)

        # get OS / SS factor
        qcd_OS = (data_hist - mc_ewk_hist) - mc_ztautau_hist

        qcd_hist = data.draw(expr, cuts & samples.SS, bins, min, max)
        qcd_hist.SetTitle('QCD')
        qcd_hist.SetFillColor('blue')
        qcd_hist.SetLineColor('blue')
        # MC subtraction
        mc_ztautau_hist_control = sum(mc_ztautau.draw(expr, cuts & samples.SS, bins, min, max))
        mc_ewk_hist_control = sum(mc_ewk.draw(expr, cuts & samples.SS, bins, min, max))
        qcd_hist -= mc_ztautau_hist_control
        qcd_hist -= mc_ewk_hist_control
        qcd_hist *= (sum(qcd_OS) / sum(qcd_hist))

        vbf_hist = sum(vbf.draw(expr, cuts & samples.OS, bins, min, max))
        ggf_hist = sum(ggf.draw(expr, cuts & samples.OS, bins, min, max))

        background_hists = [qcd_hist, mc_ztautau_hist, mc_ewk_hist,] #mc_dy_hist]
        signal_hists = [vbf_hist, ggf_hist]

        ax = fig.add_axes((0.08, 0.16, 0.6, 0.8))
        # Shink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        ax.get_frame().set_linewidth(2)
        rplt.bar(background_hists, linewidth=0, stacked=True, yerr='quadratic')
        rplt.bar(signal_hists, stacked=True, yerr='quadratic', linewidth=0)
        rplt.errorbar(data_hist, fmt='o')

        label = info['title']
        if 'units' in info:
            label += ' [%s]' % info['units']
        ax.set_xlabel(label, fontsize=20, position=(1.0, 0.), ha='right')
        ax.set_ylabel('Events', fontsize=20, position=(0., 1.), va='top')
        l = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        frame = l.get_frame()
        frame.set_alpha(.75)
        frame.set_linewidth(0)

        ax.yaxis.set_label_coords(-0.1, 1.)
        ax.xaxis.set_label_coords(1., -0.1)

        filename = '%s_%s_%s' % (os.path.join(OUT_PREFIX, info['filename']),
                                 channel, control)
        fig.savefig('%s.%s' % (filename, FORMAT))

