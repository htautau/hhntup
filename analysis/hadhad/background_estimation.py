import os
from array import array

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import TF1, TF2, TLatex

from rootpy.tree import Cut
from rootpy.plotting import Hist, Hist2D, HistStack, Legend, Canvas

import numpy as np
from scipy.optimize import fmin_l_bfgs_b

from utils import set_colours, draw
import categories
import bkg_scales_cache
from config import plots_dir
import samples


class FitError(Exception):
    pass


def draw_fit(
        expr, bins,
        xmin, xmax,
        ymin, ymax,
        model,
        data,
        category,
        region,
        name,
        output_name,
        output_formats=('png', 'eps', 'pdf'),
        root=False,
        systematics=None,
        cuts=None,
        after=False):

    PLOTS_DIR = plots_dir(__file__)

    """
    bkg_hist = bkg_hist.ravel()
    bkg_control_hist = bkg_control_hist.ravel()

    ztautau_hist = ztautau_hist.ravel()
    ztautau_hist *= ztautau_scale

    ztautau_control_hist = ztautau_control_hist.ravel()
    ztautau_control_hist *= ztautau_scale

    data_hist = data_hist.ravel()
    data_control_hist = data_control_hist.ravel()

    qcd_hist = (data_control_hist
                - ztautau_control_hist
                - bkg_control_hist) * qcd_scale

    qcd_hist.title = 'QCD Multi-jet'
    """

    model_hists = []
    for sample in model:
        hist2d = sample.draw2d(expr, category, region,
                bins, xmin, xmax, bins, ymin, ymax, cuts)
        hist = hist2d.ravel()
        if hasattr(hist2d, 'systematics'):
            hist.systematics = {}
            for term, _hist in hist2d.systematics.items():
                hist.systematics[term] = _hist.ravel()
        model_hists.append(hist)

    data_hist2d = data.draw2d(expr, category, region,
            bins, xmin, xmax, bins, ymin, ymax, cuts)
    data_hist = data_hist2d.ravel()
    if hasattr(data_hist2d, 'systematics'):
        data_hist.systematics = {}
        for term, hist in data_hist2d.systematics.items():
            data_hist.systematics[term] = hist.ravel()

    if after:
        output_name += '_after'

    draw(model=model_hists,
        data=data_hist,
        name=name,
        category_name=category,
        category=category,
        show_ratio=True,
        systematics=systematics,
        root=root,
        dir=PLOTS_DIR,
        output_formats=output_formats,
        output_name=output_name)


    """
    c = Canvas()
    hists = [qcd_hist, bkg_hist, ztautau_hist]
    set_colours(hists)
    for h in hists:
        h.SetFillStyle('solid')

    stack = HistStack()
    stack.Add(qcd_hist)
    stack.Add(bkg_hist)
    stack.Add(ztautau_hist)
    legend = Legend(4, leftmargin=.6)
    #legend.SetHeader(categories.CATEGORIES[category]['name'])
    legend.SetBorderSize(0)
    legend.AddEntry(ztautau_hist, 'F')
    legend.AddEntry(bkg_hist, 'F')
    legend.AddEntry(qcd_hist, 'F')

    if False:
        legend.AddEntry(data_hist, 'F')
        stack.Draw('lego1 0')
        data_hist.Draw('lego 0 same')

        axis_max = max(stack.GetMaximum(), data_hist.GetMaximum())

        stack[0].GetZaxis().SetLimits(0, axis_max * 1.1)
        stack[0].GetZaxis().SetRangeUser(0, axis_max * 1.1)
        if xlabel is not None:
            stack.GetXaxis().SetTitle(xlabel)
            stack.GetXaxis().SetTitleOffset(2.5)
        if ylabel is not None:
            stack.GetYaxis().SetTitle(ylabel)
            stack.GetYaxis().SetTitleOffset(2.5)
    else:
        legend.AddEntry(data_hist, 'LEP')
        stack.Draw('hist E1')
        data_hist.Draw('same E1')

        axis_max = max(stack.maximum(), data_hist.maximum())

        stack.GetYaxis().SetLimits(0, axis_max * 1.1)
        stack.GetYaxis().SetRangeUser(0, axis_max * 1.1)
        if xlabel is not None:
            stack[0].SetXTitle(xlabel)

    stack.SetMinimum(0)
    stack.SetMaximum(axis_max * 1.1)
    data_hist.SetMinimum(0)
    data_hist.SetMaximum(axis_max * 1.1)
    for f in data_hist.GetListOfFunctions():
        f.SetMinimum(0)
        f.SetMaximum(axis_max * 1.1)

    legend.Draw()

    if model_func is not None:
        chi2 = model_func.GetChisquare()
        ndf = model_func.GetNDF()
        qcd_scale = model_func.GetParameter('QCD_scale')
        qcd_scale_error = model_func.GetParError(0)
        ztautau_scale = model_func.GetParameter('Ztautau_scale')
        ztautau_scale_error = model_func.GetParError(1)
        fit_stats1 = TLatex(.2, .8, '#frac{#chi^{2}}{dof} = %.3f' % (chi2 / ndf))
        fit_stats2 = TLatex(.2, .6, 'QCD Scale = %.3f #pm %.3f' % (qcd_scale, qcd_scale_error))
        fit_stats3 = TLatex(.2, .4, 'Ztautau Scale = %.3f #pm %.3f' % (ztautau_scale, ztautau_scale_error))
        fit_stats1.Draw()
        fit_stats2.Draw()
        fit_stats3.Draw()
        name += '_after'

    c.Modified()
    c.Update()
    c.OwnMembers()
    c.Draw()
    for format in formats:
        c.SaveAs(os.path.join(PLOTS_DIR, "%s.%s" % (name, format)))
    """

def qcd_ztautau_norm(
        ztautau,
        others,
        qcd,
        data,
        category,
        target_region,
        mass_regions,
        cuts=None,
        bins=10,
        draw=False,
        use_cache=True,
        param='BDT',
        systematics=None,
        root=False):

    is_embedded = isinstance(ztautau, samples.Embedded_Ztautau)
    param = param.upper()

    if use_cache and bkg_scales_cache.has_category(
            category, is_embedded, param):
        qcd_scale, qcd_scale_error, ztautau_scale, ztautau_scale_error = \
                 bkg_scales_cache.get_scales(category, is_embedded, param)
        qcd.scale = qcd_scale
        qcd.scale_error = qcd_scale_error
        ztautau.scale = ztautau_scale
        ztautau.scale_error = ztautau_scale_error
        return

    qcd_shape_region = qcd.shape_region

    print "fitting scale factors for embedding: %s" % str(is_embedded)
    print "fitting scale factors for %s category" % category

    if param == 'BDT':
        xmin, xmax = .6, 1
        ymin, ymax = .55, 1
        expr = 'tau2_BDTJetScore:tau1_BDTJetScore'
        xlabel = '#tau_{1} BDT Score'
        ylabel = '#tau_{2} BDT Score'
        name = 'BDT Score Grid'
    elif param == 'TRACK':
        xmin, xmax = 1, 6
        ymin, ymax = 1, 6
        bins = 5 # ignore bins args above
        expr = 'tau2_ntrack_full:tau1_ntrack_full'
        xlabel = '#tau_{1} Number of Tracks'
        ylabel = '#tau_{2} Number of Tracks'
        name = 'Number of Tracks Grid'
    else:
        raise ValueError('No fit defined for %s parameters.' % param)

    ndim = 2

    output_name = "%dd_%s_fit_%s" % (ndim, param, category)
    if is_embedded:
        output_name += '_embedding'

    print "performing %d-dimensional fit using %s" % (ndim, expr)
    print "using %d bins on each axis" % bins

    assert(ndim in (1, 2))
    control = mass_regions.control_region
    control &= cuts

    print "fitting scale factors in control region: %s" % control

    if ndim == 1:
        hist = Hist(bins, min, max, name='fit_%s' % category)
    else:
        hist = Hist2D(bins, xmin, xmax, bins, ymin, ymax, name='fit_%s' % category)

    ztautau_hist = hist.Clone(title=ztautau.label)
    ztautau_hist_control = hist.Clone(title=ztautau.label)

    bkg_hist = hist.Clone(title=others.label)
    bkg_hist_control = hist.Clone(title=others.label)

    data_hist = hist.Clone(title=data.label)
    data_hist_control = hist.Clone(title=data.label)

    ztautau.draw_into(
            ztautau_hist,
            expr,
            category, target_region,
            cuts=control)

    ztautau.draw_into(
            ztautau_hist_control,
            expr,
            category, qcd_shape_region,
            cuts=control)

    others.draw_into(
            bkg_hist,
            expr,
            category, target_region,
            cuts=control)

    others.draw_into(
            bkg_hist_control,
            expr,
            category, qcd_shape_region,
            cuts=control)

    data.draw_into(
            data_hist,
            expr,
            category, target_region,
            cuts=control)

    data.draw_into(
            data_hist_control,
            expr,
            category, qcd_shape_region,
            cuts=control)

    if draw:
        draw_fit(
                expr, bins,
                xmin, xmax,
                ymin, ymax,
                model=[
                    qcd,
                    others,
                    ztautau],
                data=data,
                category=category,
                region=target_region,
                output_name=output_name,
                name=name,
                after=False,
                systematics=systematics,
                cuts=control,
                root=root)

    class Model(object):

        def __init__(self, ndim=1):

            self.ndim=ndim
            if ndim == 1:
                self.func = TF1('model_%s' % category, self, 0, 1, 2)
            else:
                self.func = TF2('model_%s' % category, self, 0, 1, 0, 1, 2)

            self.func.SetParName(0, 'QCD_scale')
            self.func.SetParName(1, 'Ztautau_scale')
            self.func.SetParameter(0, 1.)
            self.func.SetParameter(1, 1.2)

        def __call__(self, args, p):

            model = ( (data_hist_control
                       - ztautau_hist_control * p[1]
                       - bkg_hist_control) * p[0]
                    + ztautau_hist * p[1] + bkg_hist)
            bin = model.FindBin(*(list(args)[:self.ndim]))
            return model.GetBinContent(bin)

    model = Model(ndim=ndim)
    model_func = model.func
    model_func.SetLineWidth(0)
    fit_result = data_hist.Fit(model_func, 'WL')

    qcd_scale = model_func.GetParameter('QCD_scale')
    qcd_scale_error = model_func.GetParError(0)
    ztautau_scale = model_func.GetParameter('Ztautau_scale')
    ztautau_scale_error = model_func.GetParError(1)
    #data_hist.GetFunction('model_%s' % category).Delete()

    # check norm in control region
    qcd_hist = (data_hist_control
                - ztautau_hist_control * ztautau_scale
                - bkg_hist_control) * qcd_scale

    factor = data_hist.Integral() / (qcd_hist + ztautau_hist * ztautau_scale +
            bkg_hist).Integral()

    # check norm overall
    ztautau_hist_overall = hist.Clone(title=ztautau.label)
    ztautau_hist_control_overall = hist.Clone(title=ztautau.label)

    bkg_hist_overall = hist.Clone(title=others.label)
    bkg_hist_control_overall = hist.Clone(title=others.label)

    data_hist_overall = hist.Clone(title=data.label)
    data_hist_control_overall = hist.Clone(title=data.label)

    ztautau.draw_into(
            ztautau_hist_overall,
            expr,
            category, target_region)

    ztautau.draw_into(
            ztautau_hist_control_overall,
            expr,
            category, qcd_shape_region)

    others.draw_into(
            bkg_hist_overall,
            expr,
            category, target_region)

    others.draw_into(
            bkg_hist_control_overall,
            expr,
            category, qcd_shape_region)

    data.draw_into(
            data_hist_overall,
            expr,
            category, target_region)

    data.draw_into(
            data_hist_control_overall,
            expr,
            category, qcd_shape_region)

    qcd_hist_overall = (data_hist_control_overall
                - ztautau_hist_control_overall * ztautau_scale
                - bkg_hist_control_overall) * qcd_scale

    overall_factor = data_hist_overall.Integral() / (qcd_hist_overall +
            ztautau_hist_overall * ztautau_scale + bkg_hist_overall).Integral()

    print
    print "data / model in this control region: %.3f" % factor
    print "data / model overall: %.3f" % overall_factor
    print

    #qcd_scale *= factor
    #ztautau_scale *= factor

    qcd.scale = qcd_scale
    qcd.scale_error = qcd_scale_error
    ztautau.scale = ztautau_scale
    ztautau.scale_error = ztautau_scale_error

    if draw:
        draw_fit(
                expr, bins,
                xmin, xmax,
                ymin, ymax,
                model=[
                    qcd,
                    others,
                    ztautau],
                data=data,
                category=category,
                region=target_region,
                name=name,
                output_name=output_name,
                after=True,
                systematics=systematics,
                cuts=control,
                root=root)

    bkg_scales_cache.set_scales(
            category, is_embedded, param,
            qcd_scale, qcd_scale_error,
            ztautau_scale, ztautau_scale_error)
