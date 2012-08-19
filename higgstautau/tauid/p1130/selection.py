from rootpy.io import open as ropen
from ..common import nprong
import os


HERE = os.path.dirname(os.path.abspath(__file__))


f = ropen(os.path.join(HERE, 'ParametrizedBDTSelection.root'))


def selection(level, prong):

    prong = nprong(prong)
    return f.Get('%s_%dp' % (level, prong))
