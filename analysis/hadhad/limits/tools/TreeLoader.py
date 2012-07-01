#TreeLoader Class
from ROOT import TFile, TChain
from sample.Sample import Sample

class TreeLoader:
    """A class to handle the loading of D3PDs"""

    #-----------------------------------------------------------
    def __init__(self, treeName = 'tau'):
        """Constructor"""
        self._tree = TChain(treeName)


    #-----------------------------------------------------------
    def add(self, data, Test = True, Train = True):
        """Adding samples to the internal TChain"""

        # Add single root files in a text file
        if isinstance(data, str):
            f = open(data)
            lines = f.readlines()
            for line in lines:
                self._tree.Add(line.replace('\n',''))


        # Add objects of the sample class
        elif isinstance(data, Sample):
            if Test:
                tcTest = TChain(data.treeNameTest)
                tcTest.Add(data.path)
                self._tree.Add(tcTest)
            if Train:
                tcTrain = TChain(data.treeNameTrain)
                tcTrain.Add(data.path)
                self._tree.Add(tcTrain)

        else:
            raise IOError('File must be a list of .root files ending with .txt or a Sample')


    #-----------------------------------------------------------
    def getTree(self):
        """Obtain the TTree"""

        return self._tree


    #-----------------------------------------------------------
    #Eventually add options to operate a selection on the tree
