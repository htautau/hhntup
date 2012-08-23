from rootpy.io import open as ropen
from ..common import nprong, PRONGS, LEVELS
import os


HERE = os.path.dirname(os.path.abspath(__file__))


P1130_FILE = ropen(os.path.join(HERE, 'ParametrizedBDTSelection.root'))


def selection(level, prong):

    prong = nprong(prong)
    return P1130_FILE.Get('%s_%dp' % (level, prong))
