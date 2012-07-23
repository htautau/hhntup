from rootpy.tree import Cut
from rootpy.plotting import Hist, Hist2D, HistStack, Legend, Canvas
import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from utils import set_colours
import categories
import ROOT
ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetOptFit(1111)

from ROOT import TF1, TF2, TLatex
from array import array
import bkg_scales_cache


class FitError(Exception):
    pass


def draw_fit(
        data_hist,
        data_control_hist,
        ztautau_hist,
        ztautau_control_hist,
        bkg_hist,
        bkg_control_hist,
        category,
        name,
        ndim=1,
        format='png',
        xlabel=None,
        ylabel=None,
        model_func=None,
        qcd_scale=1,
        ztautau_scale=1):

    ztautau_hist = ztautau_hist.Clone()
    ztautau_hist *= ztautau_scale

    ztautau_control_hist = ztautau_control_hist.Clone()
    ztautau_control_hist *= ztautau_scale

    qcd_hist = (data_control_hist
                - ztautau_control_hist
                - bkg_control_hist) * qcd_scale
    qcd_hist.SetTitle('QCD')

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
    legend.SetHeader(categories.CATEGORIES[category]['name'])
    legend.SetBorderSize(0)
    legend.AddEntry(qcd_hist, 'F')
    legend.AddEntry(bkg_hist, 'F')
    legend.AddEntry(ztautau_hist, 'F')
    legend.AddEntry(data_hist, 'F')

    if ndim == 2:
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

        stack.Draw('hist E1')
        data_hist.Draw('E1 same')

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
        fit_stats1 = TLatex(-.9, .8, '#frac{#chi^{2}}{dof} = %.3f' % (chi2 / ndf))
        fit_stats2 = TLatex(-.9, .6, 'QCD Scale = %.3f #pm %.3f' % (qcd_scale, qcd_scale_error))
        fit_stats3 = TLatex(-.9, .4, 'Ztautau Scale = %.3f #pm %.3f' % (ztautau_scale, ztautau_scale_error))
        fit_stats1.Draw()
        fit_stats2.Draw()
        fit_stats3.Draw()
        name += '_after'

    c.Modified()
    c.Update()
    c.OwnMembers()
    c.Draw()
    c.SaveAs("%s.%s" % (name, format))


def qcd_ztautau_norm(qcd,
                     ztautau,
                     backgrounds,
                     data,
                     category='preselection',
                     target_region='OS',
                     control_region='SS',
                     cuts=None,
                     use_cache=True):

    if use_cache and bkg_scales_cache.has_category(category):
        return bkg_scales_cache.get_scales(category)

    print "fitting scale factors for category %s" % category
    min, max = .55, 1
    bins = categories.CATEGORIES[category]['fitbins']
    expr = 'tau2_BDTJetScore:tau1_BDTJetScore'
    fit_name = 'bdt'
    xlabel = '#tau_{1} BDT Score'
    ylabel = '#tau_{2} BDT Score'

    #bins, min, max = 6, -0.5, 5.5
    #expr = 'tau2_numTrack:tau1_numTrack'
    #fit_name = 'track'
    #xlabel = '#tau_{1} Number of Tracks'
    #ylabel = '#tau_{2} Number of Tracks'

    ndim = 2

    name = "%dd_%s_fit_%s" % (ndim, fit_name, category)

    assert(ndim in (1, 2))
    control = Cut('mass_mmc_tau1_tau2 < 100')
    control &= cuts

    if ndim == 1:
        hist = Hist(bins, min, max, name='fit_%s' % category)
    else:
        hist = Hist2D(bins, min, max, bins, min, max, name='fit_%s' % category)

    ztautau_hist = hist.Clone(title='Ztautau')
    ztautau_hist_control = hist.Clone(title='Ztautau')

    bkg_hist = hist.Clone(title='Other Bkg')
    bkg_hist_control = hist.Clone(title='Other Bkg')

    data_hist = hist.Clone(title='Data')
    data_hist_control = hist.Clone(title='Data')

    ztautau.draw_into(
            ztautau_hist,
            expr,
            category, target_region,
            cuts=control)

    ztautau.draw_into(
            ztautau_hist_control,
            expr,
            category, control_region,
            cuts=control)

    for b in backgrounds:
        b.draw_into(
                bkg_hist, expr,
                category, target_region,
                cuts=control)
        b.draw_into(
                bkg_hist_control, expr,
                category, control_region,
                cuts=control)

    data.draw_into(
            data_hist,
            expr,
            category, target_region,
            cuts=control)

    data.draw_into(
            data_hist_control,
            expr,
            category, control_region,
            cuts=control)

    draw_fit(data_hist,
             data_hist_control,
             ztautau_hist,
             ztautau_hist_control,
             bkg_hist,
             bkg_hist_control,
             category,
             name=name,
             xlabel=xlabel,
             ylabel=ylabel,
             ndim=ndim)

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

    # correct for scale
    qcd_hist = (data_hist_control
                - ztautau_hist_control * ztautau_scale
                - bkg_hist_control) * qcd_scale

    factor = data_hist.Integral() / (qcd_hist + ztautau_hist * ztautau_scale +
            bkg_hist).Integral()
    print factor
    print qcd_scale, ztautau_scale
    qcd_scale *= factor
    ztautau_scale *= factor

    draw_fit(data_hist,
             data_hist_control,
             ztautau_hist,
             ztautau_hist_control,
             bkg_hist,
             bkg_hist_control,
             category,
             name=name,
             xlabel=xlabel,
             ylabel=ylabel,
             model_func=model_func,
             ndim=ndim,
             qcd_scale=qcd_scale,
             ztautau_scale=ztautau_scale)

    bkg_scales_cache.set_scales(category, qcd_scale, ztautau_scale)
    return qcd_scale, ztautau_scale
