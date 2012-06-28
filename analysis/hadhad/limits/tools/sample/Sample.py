"""
A class to define samples
"""

class Sample():
    """ A class to define samples characteristics """

    def __init__(self,
                 path,
                 treeNameTrain,
                 treeNameTest,
                 sampleType,
		 Xsection = 1.0,
                 FilterEff = 1.0,
                 kFactor = 1.0,
                 Lumi = 4661.06*1000,
		 D3PDName = '',
		 D3PDnEvents = 1,
		 SkimName = '',
		 name = ''):
        """Define data members"""
        self.path = path
        self.treeNameTrain = treeNameTrain
        self.treeNameTest = treeNameTest
        self.sampleType = sampleType
        self.Xsection = Xsection
        self.FilterEff = FilterEff
        self.kFactor = kFactor
        self.Lumi = Lumi
        self.D3PDName = D3PDName
        self.D3PDnEvents = D3PDnEvents
        self.SkimName = SkimName
        self.name = name
        self.weight = (self.Xsection*self.Lumi*self.FilterEff*self.kFactor)/self.D3PDnEvents
