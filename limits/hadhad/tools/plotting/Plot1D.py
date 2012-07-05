#Plot1D class

from ROOT import *
import math
from palette import *
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()

class Plot1D:
    """A wrapper class for multiple TH1Fs"""

    #-----------------------------------------------------------
    def __init__(self, Name, nBins, xlow, xhigh, xaxis = '', yaxis = ''):
        """Constructor"""
        self._n = 0
        self._TH1Fs = []
        self._Name = Name
        self._Title = Name
        self._nBins = nBins
        self._xlow  = xlow
        self._xhigh = xhigh
        self._xaxis = xaxis
        self._yaxis = yaxis

        self._styleList = [] #0: color, 1: fillstyle, 2: linestyle, 3: line width,
        self._nameList  = []
        self._stack = []
        self._labels = []
        self._labelPositions = []
    

    #-----------------------------------------------------------
    def add(self, color, style, name, stack = True):
        """Add a new TH1F with color and style"""
        
        #Create and add the histogram
        newTH1F = TH1F(self._Name + str(self._n), name, self._nBins, self._xlow, self._xhigh)
        self._TH1Fs.append(newTH1F)

        #Stack and don't stack
        self._stack.append(stack)

        #Figure out the color
        newColor = TColor.GetColor(color)

        #Figure out and add the style
        if style == 'dashLeft':
            self._styleList.append([newColor, 3004, 1, 2, '', newColor])

        elif style == 'dashRight':
            self._styleList.append([newColor, 3005, 1, 2, '', newColor])

        elif style == 'fill':
            self._styleList.append([newColor, 1001, 1, 1, '', newColor])

        elif style == 'points':
            self._styleList.append([newColor, 0, 0, 0, 'P0E', newColor])

        else:
            self._styleList.append([newColor, 3004, 1, 2, '', newColor])

        

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
    def fill(self, value, index, weight = 1.0):
        """Fill one of the histograms"""

        try:
            self._TH1Fs[index].Fill(value, weight)
        except IndexError:
            print 'Wrong histogram index'

        return 0


    #-----------------------------------------------------------
    def applySettings(self):
        """Apply histograms settings"""

        # Load ATLAS style

        for i in range(0, len(self._TH1Fs)):
            self._TH1Fs[i].SetLineColor(self._styleList[i][5])
            self._TH1Fs[i].SetMarkerColor(self._styleList[i][0])
            self._TH1Fs[i].SetMarkerStyle(20)
            self._TH1Fs[i].SetMarkerSize(1.0);
            self._TH1Fs[i].SetLineStyle(self._styleList[i][2])
            self._TH1Fs[i].SetLineWidth(self._styleList[i][3])
            self._TH1Fs[i].SetFillStyle(self._styleList[i][1])
            self._TH1Fs[i].SetFillColor(self._styleList[i][0])

            x = self._TH1Fs[i].GetXaxis()
            x.SetTitle(self._xaxis)
            x.SetTitleOffset(1.0)
            x.SetTitleSize(0.06)

            y = self._TH1Fs[i].GetYaxis()
            y.SetTitle(self._yaxis)
            y.SetTitleOffset(0.8)
            y.SetTitleSize(0.06)

        return 0


    #-----------------------------------------------------------
    def makeLegend(self, option=''):
        """Figures out the location and size of the legend"""

        #Figure out the position of the legend corners
        binMax = [0]*self._nBins
        allMax = 0
        
        for h in self._TH1Fs:
            hMax = h.GetMaximum()
            if hMax > allMax:
                allMax = hMax
            for i in range(1, self._nBins):
                binContent = h.GetBinContent(i)
                if binContent > binMax[i-1]:
                    binMax[i-1] = binContent

        position = 0.1
        interval = int(0.3*self._nBins)

        minInterval = 1000000000000
        for i in range(1, len(binMax) - interval):
            sumBins = 0
            for j in range(i, i + interval):
                sumBins += binMax[j]
            if sumBins < minInterval:
                minInterval = sumBins
                position = float(i)/float(self._nBins)

        position = 0.18 + 0.73*position

        nLabels = 0

        if option.find('E') > -1:
            nLabels += 1

        #Instantiate the legend
        leg = TLegend(position, 0.84 - (self._n + nLabels)*0.05, position + 0.20, 0.84)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)

        for i in range(0, len(self._TH1Fs)):
            if self._styleList[i][4] != '':
                leg.AddEntry(self._TH1Fs[i], '', 'P')
            else:
                leg.AddEntry(self._TH1Fs[i], '', 'F')

        #Position labels
        for i in range(0, len(self._labels)):
            labelx = self._labelPositions[i][0]
            labely = self._labelPositions[i][1]

            if labelx == -1 and labely == -1:
                self._labelPositions[i][0] = position
                self._labelPositions[i][1] = 0.84 - (self._n + nLabels + 0.4)*0.06

                nLabels += 1
                
        return leg
        

    #-----------------------------------------------------------
    def draw(self, option=''):
        """Print the histogram to file"""

        self.applySettings()

        #Figure out the yaxis range:
        maxBin = 0
        for h in self._TH1Fs:
            currentMax = h.GetMaximum()
            if currentMax > maxBin:
                maxBin = currentMax

        for h in self._TH1Fs:
            h.SetMaximum(maxBin*1.2)
            
        #Deploy the canvas
        canvas = TCanvas('c1', 'c1', 10, 10, 800, 600)

        error = TH1F()

        #Draw the histograms
        if option.find('S') > -1:
            stack = THStack('DaStak', 'Titre of Da Stak')
            if option.find('E') > -1:
                error = TH1F(self._Name + 'Error', 'Stat. Unc.', self._nBins, self._xlow, self._xhigh)

            for i in range(0, len(self._TH1Fs)):
                if self._stack[i]:
                    stack.Add(self._TH1Fs[i])
                    if option.find('E') > -1:
                        error.Add(self._TH1Fs[i])

            stack.Draw()

            x = stack.GetXaxis()
            x.SetTitle(self._xaxis)
            x.SetTitleOffset(1.0)
            x.SetTitleSize(0.06)

            y = stack.GetYaxis()
            y.SetTitle(self._yaxis)
            y.SetTitleOffset(0.8)
            y.SetTitleSize(0.06)

            if option.find('L') > -1:
                stack.SetMaximum(stack.GetMaximum()*(1.2 + 0.08*self._n ))

            stack.Draw()
            
            if option.find('E') > -1:
                error.SetFillStyle(3154)
                error.SetFillColor(TColor.GetColor(darkred))
                error.SetLineColor(0)
                error.SetMarkerStyle(0)
                error.SetMarkerSize(0)

                error.Draw('SAMEE2')

            for i in range(0, len(self._TH1Fs)):
                if not self._stack[i]:
                     self._TH1Fs[i].Draw('SAME' + self._styleList[i][4])

        else:
            self._TH1Fs[0].Draw(self._styleList[0][4])

            for i in range(1, len(self._TH1Fs)):
                self._TH1Fs[i].Draw('SAME' + self._styleList[i][4])

        #Get the legend
        if option.find('L') > -1:
            legend = self.makeLegend(option)
            if option.find('E') > -1:
                legend.AddEntry(error, '', 'F')
            legend.Draw()

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
        return self._TH1Fs[index]
