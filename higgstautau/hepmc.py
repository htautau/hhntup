"""
This module contains various utility functions to extract
information from the MC event record.
"""

def get_VBF_partons(event):
    """
    Find and return the two outgoing VBF partons in a VBF event

    Barcodes:
    3,4,5,6 are the partons before VBF Higgs production
    7 is the Higgs
    8, 9 are the associated quark/gluons after Higgs production

    Does not work for 2012 samples...

    mc_status: 
    mc_status==21 is incoming partons
    mc_status==23 are outgoing partons 
    mc_status==22 are higgles

    This only works though for PowhegPythia we think...
    """

    partons = [p for p in event.mc if p.status==23]
    return partons

def get_VBF_partins(event):
    """
    Find and return the two incoming VBF partons in a VBF event
    """
    partins = [p for p in event.mc if p.status==21]
    return partins

def get_MC_higgles(event):
    """
    Find and return the higgs boson

    uses the optimum possible discriminator.
    """
    return [h for h in event.mc if h.status==22]
