#!/usr/bin/env python

categories = ['ggf', 'boosted', 'vbf']
fitmethod = 'trackfit'
masses = range(100, 155, 5)
lumi = '1.'
lumi_rel_err = '0.039'

template_channel = ''.join(open('template_channel.xml', 'r').readlines())
template_combination = ''.join(open('template_combination.xml', 'r').readlines())
template_combination_all = ''.join(open('template_combination_all.xml', 'r').readlines())

for mass in masses:
    for category in categories:
        with open('hh_channel_%(category)s_%(mass)d.xml' % locals(), 'w') as f:
            f.write(template_channel % locals())
        with open('hh_combination_%(category)s_%(mass)d.xml' % locals(), 'w') as f:
            f.write(template_combination % locals())
    with open('hh_combination_%(mass)d.xml' % locals(), 'w') as f:
        f.write(template_combination_all % locals())
