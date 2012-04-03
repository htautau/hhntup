from rootpy.plotting import Hist2D, Canvas
from rootpy.plotting.shapes import Ellipse
import math
import ROOT

ROOT.gStyle.SetOptStat(0)


def draw(event, vects, radii=None, colors=None, etamin=-5, etamax=5, phimin=-math.pi, phimax=math.pi):
    """
    This function creates a 2D view of an event over eta and phi
    """
    c = Canvas()
    display = Hist2D(100, etamin, etamax, 100, phimin, phimax)
    display.SetXTitle('#eta')
    display.SetYTitle('#phi')
    display.Draw()
    things = []
    for i, vect in enumerate(vects):
        radius = .4
        if radii:
            radius = radii[i]
        thing = Ellipse(vect.Eta(), vect.Phi(), radius, radius)
        thing.SetLineWidth(3)
        if colors:
            thing.SetLineColor(colors[i])
            thing.SetFillStyle(0)
        thing.Draw()
        things.append(thing)
    c.SaveAs('eventviews/event_%i.png' % event.EventNumber)
