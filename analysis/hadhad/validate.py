#!/usr/bin/env python

from rootpy.tree import Cut
from samples import *

tau1_matched = Cut('tau1_matched')
tau2_matched = Cut('tau2_matched')
both_taus_matched = tau1_matched & tau2_matched


def matching_test(sample):

    print sample.label
    total_events = sample.events()
    print "total events: ", total_events
    print "tau1 matched: ", sample.events(tau1_matched) / total_events
    print "tau2 matched: ", sample.events(tau2_matched) / total_events
    print "both matched: ", sample.events(both_taus_matched) / total_events


def fakerate_test(sample):

    assert sample.events(-tau1_matched) == sample.events('tau1_fakerate_scale_factor < 1')
    assert sample.events(tau1_matched) == sample.events('tau1_fakerate_scale_factor == 1')


ztautau = MC_Ztautau(systematics=False)
matching_test(ztautau)
fakerate_test(ztautau)

signals = []
for mass in Higgs.MASS_POINTS:
    for mode in Higgs.MODES:
        signal = Higgs(mass=mass, mode=mode,
                systematics=False)
        matching_test(signal)
        fakerate_test(signal)

