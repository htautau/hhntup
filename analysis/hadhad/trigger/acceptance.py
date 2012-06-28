#!/usr/bin/env python

from rootpy.io import open as ropen
from rootpy.tree import Tree, Cut
from prettytable import PrettyTable
import operator

f = ropen("${data}/higgs/VBFH130hh/VBFH130hh.root")
tree = f.Get("VBFH130hh")

total = float(tree.GetEntries())

triggers = tree.glob("L1_*")

trigger_table = []

for trigger in triggers:
    accept = 100. * tree.GetEntries(trigger) / float(total)
    trigger_table.append((trigger, accept))

trigger_table.sort(reverse=True, key=lambda e: e[1])

table = PrettyTable(["Trigger", "Acceptance"])
for trigger, accept in trigger_table:
    table.add_row([trigger, "%.1f%%" % accept])
print table

print tree.GetEntries("EF_tau29_medium1_tau20_medium1 || EF_e20_medium || EF_e60_loose || EF_xe60_noMu") / total
print tree.GetEntries("selected && (EF_tau29_medium1_tau20_medium1 || EF_e20_medium || EF_e60_loose || EF_xe60_noMu)") / total

print tree.GetEntries("EF_tau29_medium1_tau20_medium1") / total
print tree.GetEntries("selected && EF_tau29_medium1_tau20_medium1") / total

with open('triggers.txt') as f:
    triggers = [l.strip() for l in f.readlines()]

triggers_existing = []
# first check that all triggers are available in tree
for trigger in triggers:
    if not trigger in tree.iterbranchnames():
        print "missing %s" % trigger 
    else:
        triggers_existing.append(trigger)
triggers = triggers_existing

def removed(triggers, trigger):

    return [t for t in triggers if t != trigger]

def or_triggers(triggers):

    return reduce(operator.__or__, [Cut(t) for t in triggers]) #snap!

while True:
    cut = or_triggers(triggers)
    accept = tree.GetEntries(cut)
    print cut
    print "Acceptance: %.2f%%" % (100 * accept / total)
    print "After selection: %.2f%%" % (100 * tree.GetEntries(cut & "selected") / total)
    # remove worst trigger
    if len(triggers) > 1:
        accept_reduce = [(trigger, tree.GetEntries(or_triggers(removed(triggers, trigger)))) for trigger in triggers]
        worst_trigger = max(accept_reduce, key=lambda t: t[1])[0]
        print "will remove %s in next round..." % worst_trigger
        triggers.remove(worst_trigger)
    else:
        break
    print "-"*20
