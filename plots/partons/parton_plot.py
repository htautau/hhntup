from rootpy.io import open
from rootpy.plotting import Hist, Canvas
from rootpy.tree import Tree
from atlastools import style as  atlastyle
from rootpy.common import logon, set_style
logon(batch=True)
set_style(atlastyle.get_style())



f = open("VBFH130hh.root")
t = f.Get("VBFH130hh")


h1 = Hist(100, -6, 6)
h1.SetXTitle("#eta")
h2 = h1.Clone()
h2.SetXTitle("#eta")

t.Draw("parton1_eta","selected",hist=h1)
t.Draw("parton2_eta","selected",hist=h2)


h1.SetFillColor('red')
h1.SetFillStyle('\\')
h2.SetFillColor('blue')
h2.SetFillStyle('/')

c = Canvas()
h1.Draw("hist")
h2.Draw("hist same")

c.SaveAs("parton_eta.eps")

h3 = Hist(100, 0, 1)
h3.SetXTitle("Jet Vertex Fraction")
t.Draw("jet_AntiKt4TopoEM_jvtxf","selected && jet_AntiKt4TopoEM_jvtxf > 0",hist=h3)
c.Clear()
h3.Draw("hist")
c.SaveAs("jvf.eps")
