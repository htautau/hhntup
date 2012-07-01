from ROOT import TGraph, TGraphAsymmErrors
from ROOT import TCanvas, TLegend
from ROOT import TLatex
from ROOT import gPad


class LimitPlot(object):

    def __init__(self, masses, lumi, format='png'):

        self.masses = masses
        self.lumi = lumi
        if not isinstance(format, (list, tuple)):
            format = [format]
        self.format = format

    def draw_exclusion(self, limit, name):

        Obs  = TGraph(len(self.masses))
        Exp  = TGraph(len(self.masses))
        Exp1 = TGraphAsymmErrors(len(self.masses))
        Exp2 = TGraphAsymmErrors(len(self.masses))

        for i, mH in enumerate(self.masses):
            m = float(mH.split("G")[0])
            Obs.SetPoint(i,m,limit[mH]['Obs'])
            Exp.SetPoint(i,m,limit[mH]['Exp'])
            Exp1.SetPoint(i,m,limit[mH]['Exp'])
            Exp2.SetPoint(i,m,limit[mH]['Exp'])
            Exp1.SetPointError(i,0.,0.,(limit[mH]['Exp']-limit[mH]['sig1-']),(limit[mH]['sig1+']-limit[mH]['Exp']))
            Exp2.SetPointError(i,0.,0.,(limit[mH]['Exp']-limit[mH]['sig2-']),(limit[mH]['sig2+']-limit[mH]['Exp']))

        c = TCanvas("Test", "Test",600,600)
        c.SetBorderMode(0)
        leg = TLegend(0.59,0.75,0.92,0.94)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetLineColor(10)
        leg.SetShadowColor(kWhite)
        leg.SetTextSize(0.04)
        leg.SetBorderSize(0)
        Obs.SetMarkerStyle(20)
        Exp1.SetLineColor(kBlue)
        Exp2.SetLineColor(kRed)
        Exp2.GetXaxis().SetTitle('m_{H} [GeV]')
        Exp2.GetYaxis().SetTitle('95% CL. limit on #sigma/#sigma_{H}^{SM}')
        Exp2.GetYaxis().SetTitleOffset(1.05)
        Exp2.GetYaxis().SetTitleSize(0.05)
        Exp2.GetXaxis().SetTitleOffset(1.05)
        Exp2.GetXaxis().SetTitleSize(0.05)

        Exp2.GetXaxis().SetRangeUser(90.,160.)
        Exp2.GetYaxis().SetRangeUser(0.,40.)
        Exp2.SetFillColor(kYellow)
        Exp2.SetLineColor(kBlack)

        Exp1.SetFillColor(kGreen)
        Exp.SetLineWidth(2)
        Exp.SetLineStyle(2)

        Exp2.Draw('ACE3')
        Exp1.Draw('sameCE3')
        Exp.Draw('sameC')
        Obs.Draw('sameLP')
        gPad.RedrawAxis()

        leg.AddEntry(Obs, 'Observed CLs',"LP")
        leg.AddEntry(Exp, 'Expected CLs',"LP")
        leg.AddEntry(Exp1, '#pm 1#sigma',"F")
        leg.AddEntry(Exp2, '#pm 2#sigma',"F")
        leg.Draw()

        lumilatex = TLatex()
        lumilatex.SetNDC(true)
        lumilatex.SetTextAlign(12)
        lumilatex.SetTextSize(0.033)
        lumilatex.DrawLatex(0.19,0.85,"#sqrt{s} = 7 TeV,  #intLdt = %.1f fb^{-1}" % (self.lumi))
        lumilatex.DrawLatex(0.69,0.21, "#font[72]{ATLAS} Internal")
        lumilatex2 = TLatex()
        lumilatex2.SetNDC(true)
        lumilatex2.SetTextSize(0.04)
        lumilatex2.DrawLatex(0.19,0.9, "#tau_{had}#tau_{had} channel")

        c.Update()

        for f in self.format:
            c.SaveAs('%s.%s' % (name, f))
