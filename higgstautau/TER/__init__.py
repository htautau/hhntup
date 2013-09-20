from rootpy.utils.lock import lock
import os
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))

with lock(HERE):
    ROOT.gSystem.CompileMacro(os.path.join(HERE, 'TER.C'),
        'k',
        'TER',
        '/tmp')

from ROOT import TER

MyCrystallBall = TER.MyCrystallBall
TERafla = TER.TERafla
TERSigma = TER.TERSigma
EtaBin = TER.EtaBin
