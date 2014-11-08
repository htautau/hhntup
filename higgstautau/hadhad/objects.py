from ..mixins import *
from . import log; log = log[__name__]
from .. import datasets

def define_objects(tree, datatype=None):

    tree.define_collection('taus', 'TauRecContainer', mix=TauFourMomentum)    
    tree.define_collection('electrons', 'ElectronCollection')    
    tree.define_collection('vertices', 'PrimaryVertices')
    tree.define_collection('jets', 'AntiKt4LCTopoJets')
    tree.define_collection('jets_EM', 'AntiKt4EMTopoJets')
    tree.define_collection('MET', 'MET_RefFinal')

    if datatype != datasets.DATA:
        tree.define_collection('mc', 'TruthParticle')
        tree.define_collection('truejets', 'AntiKt4TruthJets')

        # collection filtered in hhskim
        tree.define_collection('truetaus', 'TruthParticle')

