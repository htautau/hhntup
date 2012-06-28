import array
from ROOT import std

class VariableContainer():
    """
    Define list of variables to be used in BDTs
    """

    def __init__(self):
        self.variables = {}
        self.orderedNameList = []

    def Add(self, varName, varType):

        if varType == 'F':
            setattr(self, varName, array.array('f', [0]))
        elif varType == 'I':
            setattr(self, varName, array.array('i', [0]))
        else:
            print 'Configuration.py : Please provide a valid variable type (\'I\' for integer or \'F\' for float)'
            return
        
        self.variables[varName] = varType
        self.orderedNameList.append(varName)

        return

    def SetValue(self, varName, value):
        if not varName in self.variables.iterkeys():
            print 'Variable ', varName, ' is not defined'
            return

        getattr(self, varName)[0] = value

        return

    def GetValues(self, tree):
        values = std.vector('double')()
        for variable in self.orderedNameList:
            values.push_back(getattr(tree, variable))

        return values
        
        


class PlotInfo():
    """
    Grouping plotting parameters together
    """
    def __init__(self, VarName, Nbins, Xlo, Xhi, Label, Factor):
        self.varName = VarName
        self.nbins = Nbins
        self.xlo = Xlo
        self.xhi = Xhi
        self.label = Label
        self.factor = Factor
    


class PlotInfoContainer():
    """
    Define a list of variable distributions to plot with plotting parameters
    """

    def __init__(self):
        self.plotInfos = []

    def Add(self, VarName, Nbins, Xlo, Xhi, Label, Factor = 1.0):
        self.plotInfos.append(PlotInfo(VarName, Nbins, Xlo, Xhi, Label, Factor))
