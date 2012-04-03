from atlastools import pdg

def get_VBF_partons(event):
    """
    Barcodes:
    3,4,5,6 are the partons before VBF Higgs production
    7 is the Higgs
    8, 9 are the associated quark/gluons after Higgs production
    """ 
    return [p for p in event.mc if p.barcode in (8, 9)]
