from rootpy.hep import pdg

def get_VBF_partons(event):
    """
    Barcodes:
    3,4,5,6 are the partons before VBF Higgs production
    7 is the Higgs
    8, 9 are the associated quark/gluons after Higgs production
    
    Return partons sorted by eta (increasing) 
    """ 
    return sorted([p for p in event.mc if p.barcode in (8, 9)], key=lambda p: p.eta)
