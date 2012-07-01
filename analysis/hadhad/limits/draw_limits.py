#!/usr/bin/env python
from ROOT import *
from tools.plotting.Graph1D import Graph1D
from tools.plotting import palette


for category in ['ggf', 'boosted', 'vbf', 'combined']:
    #Load data points
    mass      = range(100, 155, 5)
    sigminus2 = []
    sigminus1 = []
    median    = []
    sigplus1  = []
    sigplus2  = []

    with open('bandPlotData_%s.txt' % category) as inF:
        for line in inF.readlines():
             if line.find('-2sigma') != -1:
                  sigminus2.append( float(line.split()[-1]) )
             if line.find('-1sigma') != -1:
                  sigminus1.append( float(line.split()[-1]) )
             if line.find('Median') != -1:
                  median.append( float(line.split()[-1]) )
             if line.find('+1sigma') != -1:
                  sigplus1.append( float(line.split()[-1]) )
             if line.find('+2sigma') != -1:
                  sigplus2.append( float(line.split()[-1]) )

    print mass
    print sigminus2
    print sigminus1
    print median
    print sigplus1
    print sigplus2

    # Instantiate the Graphs
    g = Graph1D('LimitBand_%s' % category, 'm_{H}[GeV]', '95% CL. #sigma/#sigma_{H}^{SM}')
    g.add(palette.bandyellow, 'band', '#pm2#sigma')
    g.add(palette.bandgreen, 'band', '#pm1#sigma')
    g.add('#000000', 'dashline', 'Expected CLs')
    g.addLine(97.7 , 1, 152.3, 1)
    g.addLabel('ATLASPreliminary', 0.20, 0.87)
    g.addLabel('#tau_{h}#tau_{h} channel', 0.6, 0.87)
    g.addLabel('%s category' % category, 0.6, 0.79)
    g.addLabel('#sqrt{s} = 7 TeV', 0.20, 0.54)
    g.addLabel('#intL dt = ' + str(round(4661.06/1000, 2  )) + 'fb^{-1}', 0.20, 0.44)

    for i in range(len(mass)):

        sigminus2error = abs(median[i] - sigminus2[i])
        sigplus2error  = abs(median[i] - sigplus2[i])

        sigminus1error = abs(median[i] - sigminus1[i])
        sigplus1error  = abs(median[i] - sigplus1[i])

        g.fill(2, mass[i], median[i])
        g.fill(1, mass[i], median[i], 0, 0, sigminus1error, sigplus1error)
        g.fill(0, mass[i], median[i], 0, 0, sigminus2error, sigplus2error)

    g.draw('L')
