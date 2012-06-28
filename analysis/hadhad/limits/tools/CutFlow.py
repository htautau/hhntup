"""
A class to define cutflows
"""

class CutFlow():
    """A class to keep track of cuts"""
    def __init__(self):
        """ Constructor """
        self.CutFlowDict = {}
        self.numberOfCuts = 0


    def add(self, cutName):
        """Add a new counter"""
        self.CutFlowDict[cutName] = [0, self.numberOfCuts]
        self.numberOfCuts += 1


    def increment(self, cutName, weight = 1):
        """Increment one counter"""
        self.CutFlowDict[cutName][0] += weight


    def printCounters(self):
        """Print the counters in an orderly fashion"""

        numberOfSpaces = 0

        for key in self.CutFlowDict.iterkeys():
            nChar = len(key)
            if nChar > numberOfSpaces: numberOfSpaces = nChar

        print '  === Cut Flow ===  '

        for i in range(0, self.numberOfCuts):
            for key, value in self.CutFlowDict.iteritems():
                if value[1] == i:
                    print key.ljust(numberOfSpaces, ' ') + ' : ' + str(value[0])


    def reset(self):
        """Reset all counters to 0"""
        for key, value in self.CutFlowDict.iteritems():
            self.CutFlowDict[key][0] = 0
        
