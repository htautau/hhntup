def TTweight(OS, sample):
    """
    Return the TT background estimation normalization weight
    """

    if sample.name.find('T1') == -1 and sample.name.find('TTbar') == -1: return 1.0

    if OS: return 0.677986917525
    else: return 1.11204781386

    return 1.0