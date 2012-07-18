import os

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.font_manager as fm
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator, NullFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from rootpy.plotting import Hist
import rootpy.root2matplotlib as rplt
from rootpy.math.stats.qqplot import qqplot


def package_path(name):

    return os.path.splitext(os.path.abspath('latex/%s.sty' % name))[0]


LATEX_PREAMBLE = '''
\usepackage[EULERGREEK]{%s}
\sansmath
''' % package_path('sansmath')

"""
LATEX_PREAMBLE = '''
\usepackage[math-style=upright]{%s}
''' % package_path('unicode-math')
"""

plt.rcParams['ps.useafm'] = True
rc('text', usetex=True)
rc('font', family='sans-serif')
rc('text.latex', preamble=LATEX_PREAMBLE)
plt.rcParams['pdf.fonttype'] = 42


def set_colours(hists, colour_map=cm.jet):

    for i, h in enumerate(hists):
        colour = colour_map(1.*i/(len(hists)-1))
        h.SetColor(colour)


def root_axes(ax, no_xlabels=False, vscale=1.):

    ax.get_frame().set_linewidth(2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    if no_xlabels:
        ax.xaxis.set_major_formatter(NullFormatter())

    for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(2)

    ax.yaxis.set_label_coords(-0.1, 1.)
    ax.xaxis.set_label_coords(1., -0.15 / vscale)

    ax.tick_params(which='major', labelsize=15, length=8)
    ax.tick_params(which='minor', length=4)


def draw(model,
         name,
         category_name,
         output_name,
         category,
         data=None,
         signal=None,
         signal_scale=1.,
         units=None,
         range=None,
         model_colour_map=cm.winter,
         signal_colour_map=cm.autumn,
         show_ratio=False,
         show_qq=False,
         output_formats=None,
         dir='.'):

    if output_formats is None:
        output_formats = ('png', 'eps', 'pdf')
    if data is None:
        show_ratio=False
        show_qq=False
    vscale = 1.
    left_margin = 0.16
    bottom_margin = 0.16
    top_margin = 0.05
    right_margin = 0.05

    width = 1. - right_margin - left_margin
    height = 1. - top_margin - bottom_margin

    figheight = baseheight = 6.
    figwidth = basewidth = 8.

    ratio_abs_height = 1.5
    qq_abs_height = 6.
    hist_abs_height = 6.

    if show_ratio and not show_qq:
        figheight += ratio_abs_height

        vscale = figheight / baseheight
        bottom_margin /= vscale
        top_margin /= vscale

        ratio_height = ratio_abs_height / figheight
        hist_height = hist_abs_height / figheight - top_margin - bottom_margin

        rect_hist = [left_margin, bottom_margin + ratio_height, width, hist_height]
        rect_ratio = [left_margin, bottom_margin, width, ratio_height]

    elif show_qq and not show_ratio:
        figheight += qq_abs_height

        vscale = figheight / baseheight
        bottom_margin /= vscale
        top_margin /= vscale

        gap = bottom_margin
        qq_height = qq_abs_height / figheight - gap - top_margin
        hist_height = hist_abs_height / figheight - bottom_margin

        rect_qq = [left_margin, bottom_margin + hist_height + gap, width, qq_height]
        rect_hist = [left_margin, bottom_margin, width, hist_height]

    elif show_ratio and show_qq:
        figheight += ratio_abs_height + qq_abs_height

        vscale = figheight / baseheight
        bottom_margin /= vscale
        top_margin /= vscale

        ratio_height = ratio_abs_height / figheight

        gap = bottom_margin
        qq_height = qq_abs_height / figheight - gap - top_margin
        hist_height = hist_abs_height / figheight - bottom_margin

        rect_qq = [left_margin, bottom_margin + ratio_height + hist_height + gap, width, qq_height]
        rect_hist = [left_margin, bottom_margin + ratio_height, width, hist_height]
        rect_ratio = [left_margin, bottom_margin, width, ratio_height]

    else:
        rect_hist = [left_margin, bottom_margin, width, height]

    fig = plt.figure(figsize=(figwidth, figheight), dpi=100)
    prop = fm.FontProperties(size=14)

    hist_ax = plt.axes(rect_hist)

    set_colours(model, model_colour_map)
    if signal is not None:
        set_colours(signal, signal_colour_map)

    rplt.bar(model, linewidth=0, stacked=True, yerr='quadratic', axes=hist_ax,
             ypadding=(.5, .1))
    if signal is not None:
        if signal_scale != 1.:
            for sig in signal:
                sig *= signal_scale
                sig.SetTitle(r'%s $\times\/%d$' % (sig.GetTitle(), signal_scale))
        rplt.bar(signal, linewidth=0, stacked=True, yerr='quadratic',
                 axes=hist_ax, alpha=.8, ypadding=(.5, .1))
    if data is not None:
        rplt.errorbar(data, fmt='o', axes=hist_ax, ypadding=(.5, .1))

    if show_ratio:
        ratio_ax = plt.axes(rect_ratio)
        ratio_ax.axhline(y=1)
        rplt.errorbar(Hist.divide(data, sum(model), option='B'), fmt='o', axes=ratio_ax)
        ratio_ax.set_ylim((0, 2.))
        ratio_ax.yaxis.tick_right()

    if show_qq:
        qq_ax = plt.axes(rect_qq)
        gg_graph = qqplot(data, sum(model))
        gg_graph.SetTitle('QQ plot')
        y = np.array(list(gg_graph.y()))
        y_up = y + np.array(list(gg_graph.yerrh()))
        y_low = y - np.array(list(gg_graph.yerrl()))
        f = qq_ax.fill_between(list(gg_graph.x()),
                               y_low,
                               y_up,
                               interpolate=True,
                               facecolor='green',
                               linewidth=0,
                               label='68% CL band')
        #l = qq_ax.plot(xrange(-10, 10), xrange(-10, 10), 'b--')[0]
        diag = [max(gg_graph.xedgesl(0), min(y)),
                max(gg_graph.xedgesh(-1), max(y))]
        l = Line2D(diag, diag, color='b', linestyle='--')
        qq_ax.add_line(l)
        p, _, _ = rplt.errorbar(gg_graph, axes=qq_ax, snap_zero=False,
                                xerr=False, yerr=False)
        qq_ax.set_ylabel('Model', fontsize=20, position=(0., 1.), va='top')
        qq_ax.set_xlabel('Data', fontsize=20, position=(1., 0.), ha='right')
        leg = qq_ax.legend([p, Patch(facecolor='green', linewidth=0), l],
                           ['QQ plot', '68% CL band', 'Diagonal'],
                           loc='lower right', prop=prop)
        frame = leg.get_frame()
        frame.set_linewidth(0)
        qq_ax.set_xlim((gg_graph.xedgesl(0), gg_graph.xedgesh(-1)))
        qq_ax.set_ylim((min(y_low), max(y_up)))

    l = hist_ax.legend(prop=prop, title=category_name,
            loc='upper left',
            numpoints=1)
    frame = l.get_frame()
    #frame.set_alpha(.8)
    frame.set_fill(False) # eps does not support alpha values
    frame.set_linewidth(0)

    hist_ax.set_ylabel('Events', fontsize=20, position=(0., 1.), va='top')

    label = name
    if units is not None:
        label = '%s [%s]' % (label, units)

    base_ax = hist_ax
    if show_ratio:
        ratio_ax.set_ylabel('Data / Model', fontsize=20, position=(0., 1.), va='top')
        base_ax = ratio_ax
    base_ax.set_xlabel(label, fontsize=20, position=(1., 0.), ha='right')

    root_axes(hist_ax, no_xlabels=show_ratio)
    if show_ratio:
        root_axes(ratio_ax, vscale=1 if show_qq else vscale * .4)
    if show_qq:
        root_axes(qq_ax, vscale=vscale)

    if range is not None:
        hist_ax.set_xlim(range)
        if show_ratio:
            ratio_ax.set_xlim(range)

    for format in output_formats:
        plt.savefig(
            os.path.join(dir,
            'var_%s_%s.%s' %
            (category,
             output_name.lower().replace(' ', '_'),
             format)))
    plt.close(fig)
    return fig


def significance():
    """
    # plot the signal significance on the same axis
    sig_ax = ax.twinx()
    # reverse cumsum
    bins = list(bkg_hists[0].xedges())[:-1]
    sig_counts = np.array(list(sum(sig_hists)))
    bkg_counts = np.array(list(sum(bkg_hists)))

    S = (sig_counts/SIG_SCALE)[::-1].cumsum()[::-1]
    B = bkg_counts[::-1].cumsum()[::-1]
    # S / sqrt(B)
    significance = np.divide(list(S), np.sqrt(list(B)))

    max_bin = np.argmax(np.ma.masked_invalid(significance)) #+ 1
    max_sig = significance[max_bin]
    max_cut = bins[max_bin]

    print "Max signal significance %.2f at %.2f" % (max_sig, max_cut)

    sig_ax.plot(bins, significance, 'k--', label='Signal Significance')
    sig_ax.set_ylabel(r'Significance: $S (\sigma=\sigma_{SM}) / \sqrt{B}$',
            color='black', fontsize=20, position=(0., 1.), va='top')
    sig_ax.tick_params(axis='y', colors='red')

    plt.annotate('(%.2f, %.2f)' % (max_cut, max_sig), xy=(max_cut, max_sig),
            xytext=(max_cut + 0.1 * 1., max_sig),
                 arrowprops=dict(color='black', shrink=0.15),
                 ha='left', va='center', color='black')
    """
    pass
