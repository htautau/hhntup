import os
import sys
import math

import numpy as np

import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm
import matplotlib.font_manager as fm
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator, NullFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from rootpy.plotting import Canvas, Pad, Legend, Hist, HistStack
import rootpy.root2matplotlib as rplt
from rootpy.math.stats.qqplot import qqplot


class Tee(object):
    """
    http://stackoverflow.com/questions/616645/
    how-do-i-duplicate-sys-stdout-to-a-log-file-in-python/3423392#3423392
    """
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        sys.stdout = self.stdout
        self.file.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)


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

#plt.rcParams['ps.useafm'] = True
#rc('text', usetex=True)
#rc('font', family='sans-serif')
rc('text.latex', preamble=LATEX_PREAMBLE)
#plt.rcParams['pdf.fonttype'] = 42


def set_colours(hists, colour_map=cm.jet):

    for i, h in enumerate(hists):
        colour = colour_map((i + 1) / float(len(hists) + 1))
        h.SetColor(colour)


def format_legend(l):

    frame = l.get_frame()
    #frame.set_alpha(.8)
    frame.set_fill(False) # eps does not support alpha values
    frame.set_linewidth(0)


def root_axes(ax, no_xlabels=False, vscale=1.):

    ax.patch.set_linewidth(2)
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
         dir='.',
         systematics=None,
         root=False,
         width=8.):

    if output_formats is None:
        output_formats = ('png', 'eps', 'pdf')
    elif isinstance(output_formats, str):
        output_formats = output_formats.split(',')
    if data is None:
        show_ratio=False
        show_qq=False

    vscale = 1.
    left_margin = 0.16
    bottom_margin = 0.16
    top_margin = 0.05
    right_margin = 0.05
    ratio_sep_margin = 0.025
    ypadding = (.55, .1)

    width = 1. - right_margin - left_margin
    height = 1. - top_margin - bottom_margin

    figheight = baseheight = 6.
    figwidth = basewidth = width

    ratio_abs_height = 1.975
    qq_abs_height = 6.
    hist_abs_height = 6.

    if show_ratio and not show_qq:
        figheight += ratio_abs_height + ratio_sep_margin

        vscale = figheight / baseheight
        bottom_margin /= vscale
        top_margin /= vscale

        ratio_height = ratio_abs_height / figheight
        hist_height = (hist_abs_height / figheight
                       - top_margin - bottom_margin)

        rect_hist = [left_margin,
                     bottom_margin + ratio_height + ratio_sep_margin,
                     width, hist_height]
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

    if root:
        fig = Canvas(width=int(figwidth * 100), height=int(figheight * 100))
        fig.SetLeftMargin(0)
        fig.SetBottomMargin(0)
        fig.SetRightMargin(0)
        fig.SetTopMargin(0)
        hist_pad = Pad("top", "top", 0, rect_hist[1], 1, 1)
        hist_pad.SetBottomMargin(0)
        hist_pad.SetLeftMargin(rect_hist[0])
        hist_pad.SetRightMargin(1. - rect_hist[2] - rect_hist[0])
        hist_pad.SetTopMargin(1. - rect_hist[3] - rect_hist[1])
        hist_pad.Draw()
    else:
        fig = plt.figure(figsize=(figwidth, figheight), dpi=100)
        prop = fm.FontProperties(size=14)
        hist_ax = plt.axes(rect_hist)

    if model_colour_map is not None:
        set_colours(model, model_colour_map)

    if isinstance(signal, (list, tuple)) and signal_colour_map is not None:
        set_colours(signal, signal_colour_map)

    if root:
        # plot model stack with ROOT
        hist_pad.cd()
        model_stack = HistStack()
        for hist in model:
            hist.SetLineWidth(0)
            hist.format = 'hist'
            model_stack.Add(hist)
        model_stack.Draw()
    else:
        model_bars = rplt.bar(model, linewidth=0,
                stacked=True, axes=hist_ax,
                ypadding=ypadding)

    if signal is not None:
        if signal_scale != 1.:
            if isinstance(signal, (list, tuple)):
                for sig in signal:
                    sig *= signal_scale
                    sig.SetTitle(r'%s $\times\/%d$' % (sig.GetTitle(), signal_scale))
            else:
                signal *= signal_scale
                signal.SetTitle(r'%s $\times\/%d$' % (signal.GetTitle(),
                    signal_scale))

        if isinstance(signal, (list, tuple)):
            if root:
                pass
            else:
                signal_bars = rplt.bar(signal, linewidth=0,
                        stacked=True, yerr='quadratic',
                        axes=hist_ax, alpha=.8, ypadding=ypadding)
        else:
            if root:
                pass
            else:
                _, _, signal_bars = rplt.hist(signal,
                        histtype='stepfilled',
                        axes=hist_ax, ypadding=ypadding)

    if show_qq:
        qq_ax = plt.axes(rect_qq)
        gg_graph = qqplot(data, sum(model))
        gg_graph.SetTitle('QQ plot')
        y = np.array(list(gg_graph.y()))
        y_up = y + np.array(list(gg_graph.yerrh()))
        y_low = y - np.array(list(gg_graph.yerrl()))
        f = qq_ax.fill_between(
                list(gg_graph.x()),
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

    if systematics is not None:
        # draw systematics band
        # add separate variations in quadrature
        # also include stat error in quadrature
        total_model = sum(model)
        var_high = []
        var_low = []
        for variations in systematics:
            if len(variations) == 2:
                high, low = variations
            elif len(variations) == 1:
                high = variations[0]
                low = 'NOMINAL'
            else:
                raise ValueError(
                        "only one or two variations per term are allowed")
            total_high = model[0].Clone()
            total_high.Reset()
            total_low = total_high.Clone()
            total_max = total_high.Clone()
            total_min = total_high.Clone()
            for m in model:
                total_high += m.systematics[high]
                if low == 'NOMINAL':
                    total_low += m.Clone()
                else:
                    total_low += m.systematics[low]
            for i in xrange(len(total_high)):
                total_max[i] = max(total_high[i], total_low[i], total_model[i])
                total_min[i] = min(total_high[i], total_low[i], total_model[i])
            var_high.append(total_max)
            var_low.append(total_min)

        # include stat error variation
        total_model_stat_high = total_model.Clone()
        total_model_stat_low = total_model.Clone()
        for i in xrange(len(total_model)):
            total_model_stat_high[i] += total_model.yerrh(i)
            total_model_stat_low[i] -= total_model.yerrl(i)
        var_high.append(total_model_stat_high)
        var_low.append(total_model_stat_low)

        # sum variations in quadrature bin-by-bin
        high_band = total_model.Clone()
        high_band.Reset()
        low_band = high_band.Clone()
        for i in xrange(len(high_band)):
            sum_high = math.sqrt(
                    sum([(v[i] - total_model[i])**2 for v in var_high]))
            sum_low = math.sqrt(
                    sum([(v[i] - total_model[i])**2 for v in var_low]))
            high_band[i] = sum_high
            low_band[i] = sum_low

        if root:
            pass
        else:
            # draw band as hatched histogram with base of model - low_band
            # and height of high_band + low_band
            rplt.fill_between(total_model + high_band,
                        total_model - low_band,
                        edgecolor='yellow',
                        linewidth=0,
                        facecolor=(0,0,0,0),
                        hatch='////',
                        axes=hist_ax,
                        zorder=100)

    if data is not None:
        # draw data
        if root:
            hist_pad.cd()
            data.Draw('same E1')
        else:
            data_bars = rplt.errorbar(data,
                    fmt='o', axes=hist_ax,
                    ypadding=ypadding,
                    emptybins=False,
                    barsabove=True,
                    zorder=1000)
        # draw ratio plot
        if show_ratio:
            total_model = sum(model)
            numerator = data - total_model
            error_hist = Hist.divide(numerator, total_model, option='B')
            error_hist.linecolor = 'black'
            error_hist.linewidth = 1
            error_hist.fillstyle = 'hollow'
            error_hist *= 100
            if root:
                fig.cd()
                ratio_pad = Pad("ratio", "ratio",
                    0, 0, 1, rect_ratio[1] + rect_ratio[3])
                ratio_pad.SetBottomMargin(rect_ratio[1])
                ratio_pad.SetLeftMargin(rect_ratio[0])
                ratio_pad.SetRightMargin(1. - rect_ratio[2] - rect_ratio[0])
                ratio_pad.SetTopMargin(0)
                ratio_pad.Draw()
                ratio_pad.cd()
                error_hist.Draw('hist')
                error_hist.yaxis.SetLimits(-100, 100)
                error_hist.yaxis.SetRangeUser(-100, 100)
                xmin = model_stack.xaxis.GetXmin()
                xmax = model_stack.xaxis.GetXmax()
                error_hist.xaxis.SetLimits(xmin, xmax)
                error_hist.xaxis.SetRangeUser(xmin, xmax)
            else:
                ratio_ax = plt.axes(rect_ratio)
                ratio_ax.axhline(y=0, color='black')
                ratio_ax.axhline(y=50, color='black', linestyle='--')
                ratio_ax.axhline(y=-50, color='black', linestyle='--')
                rplt.hist(
                        error_hist,
                        axes=ratio_ax,
                        histtype='stepfilled')
                ratio_ax.set_ylim((-100., 100.))
                ratio_ax.set_xlim(hist_ax.get_xlim())
                #ratio_ax.yaxis.tick_right()
                ratio_ax.set_ylabel(r'$\frac{\rm{Data - Model}}{\rm{Model}}$ [\%]',
                        fontsize=20, position=(0., 1.), va='top')
            if systematics is not None:
                # plot band on ratio plot
                # uncertainty on top is data + model
                high_band_top = high_band.Clone()
                low_band_top = low_band.Clone()
                # quadrature sum of model uncert + data stat uncert in numerator
                for i in xrange(len(high_band_top)):
                    high_band_top[i] = math.sqrt(
                            high_band[i]**2 +
                            data.yerrh(i)**2)
                    low_band_top[i] = math.sqrt(
                            low_band[i]**2 +
                            data.yerrl(i)**2)
                # full uncert
                high_band_full = high_band.Clone()
                low_band_full = low_band.Clone()
                # quadrature sum of numerator and denominator
                for i in xrange(len(high_band_full)):
                    if numerator[i] == 0 or total_model[i] == 0:
                        high_band_full[i] = 0.
                        low_band_full[i] = 0.
                    else:
                        high_band_full[i] = abs(error_hist[i]) * math.sqrt(
                                (high_band_top[i] / numerator[i])**2 +
                                (high_band[i] / total_model[i])**2)
                        low_band_full[i] = abs(error_hist[i]) * math.sqrt(
                                (low_band_top[i] / numerator[i])**2 +
                                (low_band[i] / total_model[i])**2)
                if root:
                    pass
                else:
                    rplt.fill_between(
                        error_hist + high_band_full,
                        error_hist - low_band_full,
                        edgecolor='black',
                        linewidth=0,
                        facecolor=(0,0,0,0),
                        hatch='\\\\\\\\',
                        axes=ratio_ax)
    if root:
        hist_pad.cd()
        model_legend = Legend(len(model), pad=hist_pad,
                leftmargin=0.05, rightmargin=0.5)
        for hist in model:
            model_legend.AddEntry(hist, 'F')
        model_legend.SetHeader(category_name)
        model_legend.Draw()
    else:
        model_legend = hist_ax.legend(
                reversed(model_bars), [h.title for h in reversed(model)],
                prop=prop, title=category_name,
                loc='upper left',
                numpoints=1)
        format_legend(model_legend)

    if root:
        hist_pad.cd()
        right_legend = Legend(2 if signal is not None else 1, pad=hist_pad)
        right_legend.AddEntry(data, 'lep')
        if signal is not None:
            # TODO support list of signal
            right_legend.AddEntry(signal, 'F')
        right_legend.Draw()
    else:
        right_legend_bars = []
        right_legend_titles =[]

        if data is not None:
            right_legend_bars.append(data_bars)
            right_legend_titles.append(data.title)
        if signal is not None:
            if isinstance(signal, (list, tuple)):
                right_legend_bars += signal_bars
                right_legend_titles += [s.title for s in signal]
            else:
                right_legend_bars.append(signal_bars[0])
                right_legend_titles.append(signal.title)

        if right_legend_bars:
            right_legend = hist_ax.legend(
                    right_legend_bars,
                    right_legend_titles,
                    prop=prop,
                    loc='upper right',
                    numpoints=1)
            format_legend(right_legend)
            hist_ax.add_artist(model_legend)

    if units is not None:
        label = '%s [%s]' % (name, units)
        binwidths = list(set(['%.3g' % w for w in list(model[0].xwidth())]))
        if len(binwidths) == 1:
            # constant width bins
            ylabel = 'Events / %s [%s]' % (binwidths[0], units)
        else:
            ylabel = 'Events'
    else:
        label = name
        ylabel = 'Events'

    if root:
        model_stack.yaxis.SetTitle('Events')
        base_hist = model_stack
        if show_ratio:
            base_hist = error_hist
        base_hist.xaxis.SetTitle(label)
    else:
        hist_ax.set_ylabel(ylabel, fontsize=20, position=(0., 1.), va='top')
        base_ax = hist_ax
        if show_ratio:
            base_ax = ratio_ax
        base_ax.set_xlabel(label, fontsize=20, position=(1., 0.), ha='right')
        root_axes(hist_ax, no_xlabels=show_ratio)
        if show_ratio:
            root_axes(ratio_ax, vscale=1 if show_qq else vscale * .4)
        if show_qq:
            root_axes(qq_ax, vscale=vscale)

    if range is not None:
        if root:
            model_stack.xaxis.SetLimits(*range)
            model_stack.xaxis.SetRangeUser(*range)
        else:
            hist_ax.set_xlim(range)
        if show_ratio:
            if root:
                pass
            else:
                ratio_ax.set_xlim(range)

    filename = os.path.join(dir,
            'var_%s_%s' %
            (category,
             output_name.lower().replace(' ', '_')))
    if root:
        filename += '_root'
    for format in output_formats:
        if root:
            fig.SaveAs('%s.%s' % (filename, format))
        else:
            plt.savefig('%s.%s' % (filename, format))
    if not root:
        plt.close(fig)
    else:
        fig.OwnMembers()
    return fig

"""
def significance():

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
