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
    """
    return [p for p in event.mc if p.barcode in (8, 9)]
