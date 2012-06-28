from ROOT import *
from TreeLoader import TreeLoader
from Sample import Sample
from sampleGroups import *

groups = [
    AlpgenJimmyWmunu,
    AlpgenJimmyWtaunu,
    AlpgenJimmyZtautau_pt20,
    AlpgenJimmyZmumu_pt20,
    TTBar,
    McAtNlo_JIMMY_WZ,
    McAtNlo_JIMMY_WW,
    McAtNlo_JIMMY_ZZ,
    [PowHegPythia_VBFH125_tautaulh],
    DataGroup
    ]

for group in groups:
    totalNEvents = 0
    totalWeightedNEvents = 0
    print '====================================================='
    for sample in group:
        loader = TreeLoader(sample.treeName)
        loader.add(sample)
        tree = loader.getTree()
        nEvents = tree.GetEntries()
        print sample.name, ':', nEvents

        totalNEvents += nEvents
        totalWeightedNEvents += nEvents* sample.weight

    print '------------'
    print 'Total :', totalNEvents
    print 'Total weighted :', totalWeightedNEvents 

        
