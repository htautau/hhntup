#Plot2D class

from ROOT import *
import math
from palette import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()

styleDict = {'small'  : [20, 0.3, ''],
             'medium' : [20, 0.5, ''],
             'large'  : [20, 0.7, ''],
             'color'  : [20, 0.3, 'COLZ']}

class Plot2D:
    """A wrapper class for multiple TH2Fs"""

    #-----------------------------------------------------------
    def __init__(self, Name, xNBins, xlow, xhigh, yNBins, ylow, yhigh, xaxis = '', yaxis = ''):
        """Constructor"""
        self._n = 0
        self._TH2Fs = []
        self._Name = Name
        self._Title = Name
        self._xNBins = xNBins
        self._xlow  = xlow
        self._xhigh = xhigh
        self._yNBins = yNBins
        self._ylow  = ylow
        self._yhigh = yhigh
        self._xaxis = xaxis
        self._yaxis = yaxis

        self._styleList = [] #0: color, 1: fillstyle, 2: linestyle, 3: line width,
        self._nameList  = []
        self._labels = []
        self._labelPositions = []



    #-----------------------------------------------------------
    def add(self, color, style, name):
        """Add a new TH1F with color and style"""
        
        #Create and add the histogram
        newTH2F = TH2F(self._Name + str(self._n), name, self._xNBins, self._xlow, self._xhigh, self._yNBins, self._ylow, self._yhigh)
        self._TH2Fs.append(newTH2F)

        #Figure out the color
        newColor = TColor.GetColor(color)

        #Figure out and add the style
        self._styleList.append([newColor] + styleDict[style])

        #Add the name to namelist
        self._nameList.append(name)

        #Increment the count of histograms
        self._n += 1

        return 0



    #-----------------------------------------------------------
    def addLabel(self, text, x = -1, y = -1):
        """Add a new label under the legend"""
        self._labels.append(text)
        self._labelPositions.append([x, y])

        return 0



    #-----------------------------------------------------------
    def fill(self, xValue, yValue, index, weight = 1.0):
        """Fill one of the histograms"""

        try:
            self._TH2Fs[index].Fill(xValue, yValue, weight)
        except IndexError:
            print 'Wrong histogram index'

        return 0



    #-----------------------------------------------------------
    def applySettings(self):
        """Apply histograms settings"""

        # Load ATLAS style

        for i in range(0, len(self._TH2Fs)):
            self._TH2Fs[i].SetMarkerColor(self._styleList[i][0])
            self._TH2Fs[i].SetMarkerStyle(self._styleList[i][1])
            self._TH2Fs[i].SetMarkerSize(self._styleList[i][2])

            x = self._TH2Fs[i].GetXaxis()
            x.SetTitle(self._xaxis)
            x.SetTitleOffset(1.0)
            x.SetTitleSize(0.06)

            y = self._TH2Fs[i].GetYaxis()
            y.SetTitle(self._yaxis)
            y.SetTitleOffset(0.8)
            y.SetTitleSize(0.06)

        return 0




    #-----------------------------------------------------------
    def draw(self, option=''):
        """Print the histogram to file"""

        self.applySettings()
            
        #Deploy the canvas
        canvas = TCanvas('c1', 'c1', 10, 10, 800, 600)

        #Draw the histograms
        self._TH2Fs[0].Draw(self._styleList[0][3])
        
        for i in range(1, len(self._TH2Fs)):
            self._TH2Fs[i].Draw('SAME' + self._styleList[i][3])

        #Get the legend
        #if option.find('L') > -1:
        #    legend = self.makeLegend(option)
        #    if option.find('E') > -1:
        #        legend.AddEntry(error, '', 'F')
        #    legend.Draw()

        TL = TLatex()

        for i in range(0, len(self._labels)):
            if self._labels[i] == 'ATLASPreliminary':
                ATLASLabel(self._labelPositions[i][0], self._labelPositions[i][1], 'Preliminary', TColor.GetColor(black))
            else:
                TL.SetNDC()
                TL.SetTextSize(0.045)
                TL.DrawLatex(self._labelPositions[i][0], self._labelPositions[i][1], self._labels[i])

        #Print the canvas
        canvas.Print(self._Name + '.png')


    #-----------------------------------------------------------
    def getTH1F(self, index):
        return self._TH2Fs[index]
