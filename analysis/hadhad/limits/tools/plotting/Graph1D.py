#Graph1D class

from ROOT import *
import math
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()

styleDict = {}

styleDict['point0']   = [21, 1, 1, 'P']
styleDict['point1']   = [20, 1, 1, 'P']
styleDict['style2']   = [22, 1, 1, 'P']
styleDict['style3']   = [23, 1, 1, 'P']
styleDict['band']     = [23, 1, 1, '3']
styleDict['dotline']  = [20, 1, 2, 'LP']
styleDict['dashline'] = [20, 2, 2, 'C']

class Graph1D:
    """A wrapper class for multiple TGraphs"""

    #-----------------------------------------------------------
    def __init__(self, Name, xaxis = '', yaxis = ''):
        """Constructor"""
        self.__n = 0
        self.__TGraphs  = []
        self.__TLines   = []

        self.__Name   = Name
        self.__Title  = Name

        self.__xlow   = []
        self.__xhigh  = []
        self.__ylow   = []
        self.__yhigh  = []

        self.__xmin   = 0.0
        self.__xmax   = 1.0
        self.__ymin   = 0.0
        self.__ymax   = 1.0

        self.__xaxis  = xaxis
        self.__yaxis  = yaxis

        self.__pointList = []

        self.__canvas = TCanvas('c1', 'c1', 10, 10, 800, 600)

        self.__styleList = [] #0: color, 1: markerstyle, 2: linestyle, 3: size, 4: draw option
        self.__nameList  = []
        self.__labels = []
        self.__labelPositions = []


    #-----------------------------------------------------------
    def add(self, color, style, name):
        """Add a new TGraphAsymmErrors with color and style"""

        #Create and add the histogram
        newTGraph = TGraphAsymmErrors()
        self.__TGraphs.append(newTGraph)
        self.__pointList.append([])

        #Figure out the color
        newColor = TColor.GetColor(color)

        #Figure out and add the style
        self.__styleList.append([newColor] + styleDict[style])

        #Add the name to namelist
        self.__nameList.append(name)

        #Increment the count of histograms
        self.__n += 1

        return


    #----------------------------------------------------------
    def addLine(self, x1, y1, x2, y2):
        """ Add a line on the graph """

        newTLine = TLine(x1, y1, x2, y2)
        self.__TLines.append(newTLine)

        return

    #-----------------------------------------------------------
    def addLabel(self, text, x = -1, y = -1):
        """Add a new label under the legend"""
        self.__labels.append(text)
        self.__labelPositions.append([x, y])

        return 0


    #-----------------------------------------------------------
    def addBayesDivide(self, passed, total, color, style, name):
        """Makes a TGraph from two TH1Fs"""
        newTGraph = TGraphAsymmErrors()
        newTGraph.BayesDivide(passed, total)
        self.__TGraphs.append(newTGraph)
        self.__pointList.append([])

        #Figure out the color
        newColor = TColor.GetColor(color)

        #Figure out and add the style
        self.__styleList.append([newColor] + styleDict[style])

        #Add the name to namelist
        self.__nameList.append(name)

        #Increment the count of histograms
        self.__n += 1

        return

    #-----------------------------------------------------------
    def fill(self, index, xvalue, yvalue, xerrorlow = 0, xerrorhigh = 0, yerrorlow = 0, yerrorhigh = 0):
        """Fill one of the histograms"""

        try:
            self.__pointList[index].append([xvalue, yvalue, xerrorlow, xerrorhigh, yerrorlow, yerrorhigh])
        except IndexError:
            print 'Wrong histogram index'

        return


    #-----------------------------------------------------------
    def transferPoints(self):
        """Tranfers the points from the internal pointList to the TGraphs"""

        for i in range(0, self.__n):
            if len(self.__pointList[i]) > 0:
                self.__pointList[i].sort(key=lambda point: point[0])
                for j in range(0, len(self.__pointList[i])):
                    x = self.__pointList[i][j][0]
                    y = self.__pointList[i][j][1]
                    exlow  = self.__pointList[i][j][2]
                    exhigh = self.__pointList[i][j][3]
                    eylow  = self.__pointList[i][j][4]
                    eyhigh = self.__pointList[i][j][5]

                    self.__TGraphs[i].SetPoint(j, x, y)
                    self.__TGraphs[i].SetPointEXlow(j, exlow)
                    self.__TGraphs[i].SetPointEXhigh(j, exhigh)
                    self.__TGraphs[i].SetPointEYlow(j, eylow)
                    self.__TGraphs[i].SetPointEYhigh(j, eyhigh)

                self.__xlow.append(self.__pointList[i][0][0] - self.__pointList[i][0][2])
                self.__xhigh.append(self.__pointList[i][-1][0] + self.__pointList[i][-1][3])

                self.__pointList[i].sort(key=lambda point: point[1])

                self.__ylow.append(self.__pointList[i][0][1] - self.__pointList[i][0][4])
                self.__yhigh.append(self.__pointList[i][-1][1] + self.__pointList[i][-1][5])

            else:
                nPoints = self.__TGraphs[i].GetN()
                for j in range(0, nPoints):

                    x = Double(0)
                    y = Double(0)

                    self.__TGraphs[i].GetPoint(j, x, y)
                    exlow  = self.__TGraphs[i].GetErrorXlow(j)
                    exhigh = self.__TGraphs[i].GetErrorXhigh(j)
                    eylow  = self.__TGraphs[i].GetErrorYlow(j)
                    eyhigh = self.__TGraphs[i].GetErrorYhigh(j)


                    self.__pointList[i].append([x, y, exlow, exhigh, eylow, eyhigh])

                self.__pointList[i].sort(key=lambda point: point[0])

                self.__xlow.append(self.__pointList[i][0][0] - self.__pointList[i][0][2])
                self.__xhigh.append(self.__pointList[i][-1][0] + self.__pointList[i][-1][3])

                self.__pointList[i].sort(key=lambda point: point[1])

                self.__ylow.append(self.__pointList[i][0][1] - self.__pointList[i][0][4])
                self.__yhigh.append(self.__pointList[i][-1][1] + self.__pointList[i][-1][5])

        return


    #-----------------------------------------------------------
    def getLegendCoord(self, xlow = -1, xhigh = -1, ylow = -1, yhigh = -1):
        """Set the ranges of the TGraphs"""

        if xlow > 0 and xhigh > 0 and ylow > 0 and yhigh > 0:
            return xlow, xhigh, ylow, yhigh

        legXlow  = 0.18
        legXhigh = 0.95
        legYlow  = 0.10
        legYhigh = 0.84

        xran = self.__xmax - self.__xmin
        yran = self.__ymax - self.__ymin

        cx = 0
        cy = 0

        #Figure out the corner which is the less busy for the legend
        corners = [2, 2, 2, 2] # 0: UpperLeft, 1: UpperRight, 2: LowerRight, 3: LowerLeft
        pointX  = [0, 1, 1, 0]
        pointY  = [1 ,1 ,0 ,0]

        for i in range(0, self.__n):

            for j in range(0, len(self.__pointList[i])):
                x = (self.__pointList[i][j][0] - self.__xmin)/xran
                y = (self.__pointList[i][j][1] - self.__ymin)/yran

                exlow  = self.__pointList[i][j][2]/xran
                exhigh = self.__pointList[i][j][3]/xran
                eylow  = self.__pointList[i][j][4]/yran
                eyhigh = self.__pointList[i][j][5]/yran

                dUpperLeft  = self.distance(x - exlow, y + eyhigh, 0, 1)
                dUpperRight = self.distance(x + exhigh, y + eyhigh, 1, 1)
                dLowerRight = self.distance(x + exhigh, y - eylow, 1, 0)
                dLowerLeft  = self.distance(x - exlow, y - eylow, 0, 0)

                if dUpperLeft < corners[0]: corners[0] = dUpperLeft
                if dUpperRight < corners[1]: corners[1] = dUpperRight
                if dLowerRight < corners[2]: corners[2] = dLowerRight
                if dLowerLeft < corners[3]: corners[3] = dLowerLeft


        cornerIndex = corners.index(max(corners))

        #Figure out the other corner of the legend
        if cornerIndex == 0:
            cx = 0
            cy = 1

        if cornerIndex == 1:
            cx = 1
            cy = 1

        if cornerIndex == 2:
            cx = 1
            cy = 0

        maxd = 0

        bestX = 0
        bestXlow = 0
        bestXhigh = 0
        bestY = 0
        bestYlow = 0
        bestYhigh = 0

        for i in range(0, self.__n):
            for j in range(0, len(self.__pointList[i])):
                x = (self.__pointList[i][j][0] - self.__xmin)/xran
                y = (self.__pointList[i][j][1] - self.__ymin)/yran

                exlow  = self.__pointList[i][j][2]/xran
                exhigh = self.__pointList[i][j][3]/xran
                eylow  = self.__pointList[i][j][4]/yran
                eyhigh = self.__pointList[i][j][5]/yran

                d = self.area(x,y,cx,cy)

                if d > maxd:
                    pointIsIn = False
                    for point in self.__pointList[i]:
                        px = point[0]/xran
                        py = point[1]/yran

                        if (px < max([x,cx])) and (px > min([x,cx])) and (py < max([y,cy])) and (py > min([y,cy])):
                            pointIsIn = True
                            break

                    if not pointIsIn:
                        maxd = d
                        bestX = x
                        bextXlow = exlow
                        bestXhigh = exhigh
                        bestY = y
                        bestYlow = eylow
                        bestYhigh = eyhigh

        if cornerIndex == 0:
                legXhigh = legXlow + 0.23
                legYlow  = legYhigh - 0.08*self.__n

        if cornerIndex == 1:
                legXlow  = legXhigh - 0.23
                legYlow  = legYhigh - 0.08*self.__n

        if cornerIndex == 2:
                legXlow  = legXhigh - 0.23
                legYhigh = legYlow + 0.08*self.__n

        if cornerIndex == 3:
                legXhigh = legXlow + 0.23
                legYhigh = legYlow + 0.08*self.__n

        return legXlow, legXhigh, legYlow, legYhigh



    #-----------------------------------------------------------
    def applySettings(self):
        """Apply histograms settings"""

        # Load ATLAS style

        for i in range(0, len(self.__TGraphs)):
            self.__TGraphs[i].SetMarkerColor(self.__styleList[i][0])
            self.__TGraphs[i].SetLineColor(self.__styleList[i][0])
            self.__TGraphs[i].SetFillColor(self.__styleList[i][0])
            self.__TGraphs[i].SetMarkerStyle(self.__styleList[i][1])
            self.__TGraphs[i].SetMarkerSize(self.__styleList[i][3])
            self.__TGraphs[i].SetLineWidth(self.__styleList[i][3])
            self.__TGraphs[i].SetLineStyle(self.__styleList[i][2])

            x = self.__TGraphs[i].GetXaxis()
            x.SetTitle(self.__xaxis)
            x.SetTitleOffset(1.0)
            x.SetTitleSize(0.06)

            y = self.__TGraphs[i].GetYaxis()
            y.SetTitle(self.__yaxis)
            y.SetTitleOffset(1.2)
            y.SetTitleSize(0.06)

        return


    #-----------------------------------------------------------
    def makeLegend(self, option=''):
        """Figures out the location and size of the legend"""

        Xlow, Xhigh, Ylow, Yhigh = self.getLegendCoord()
        Xlow = .2
        Xhigh = .5
        Ylow = .6
        Yhigh = .8
        #Instantiate the legend
        leg = TLegend(Xlow, Ylow, Xhigh, Yhigh)
        leg.SetFillColor(TColor.GetColor('#ffffff'))
        leg.SetBorderSize(0)

        for i in range(0, len(self.__TGraphs)):
            legStyle =  self.__styleList[i][4]
            if self.__styleList[i][4] == '3':
                legStyle = 'F'
            if self.__styleList[i][4] == 'C':
                legStyle = 'L'

            leg.AddEntry(self.__TGraphs[i], self.__nameList[i] , legStyle)

        return leg


    #-----------------------------------------------------------
    def draw(self, option=''):
        """Print the histogram to file"""

        self.transferPoints()
        self.getMinMax()
        self.applySettings()

        #Figure out y range
        for i in range(0, self.__n):
            self.__TGraphs[i].SetMinimum(self.__ymin)
            self.__TGraphs[i].SetMaximum(self.__ymax)
            self.__TGraphs[i].GetXaxis().SetRangeUser(self.__xmin, self.__xmax)

        #Draw the histograms
        self.__TGraphs[0].Draw('A' + self.__styleList[0][4])

        for i in range(1, self.__n):
            self.__TGraphs[i].Draw('SAME' + self.__styleList[i][4])

        #Draw the lines
        for i in range(0, len(self.__TLines)):
            self.__TLines[i].Draw('SAME')


        if 'L' in option:
            #Get the legend
            legend = self.makeLegend()

            #Draw the legend
            legend.Draw('SAME')

        #Draw the labels
        TL = TLatex()

        for i in range(0, len(self.__labels)):
            if self.__labels[i] == 'ATLASPreliminary':
                ATLASLabel(self.__labelPositions[i][0], self.__labelPositions[i][1], 'Preliminary', kBlack)
            else:
                TL.SetNDC()
                TL.SetTextSize(0.045)
                TL.DrawLatex(self.__labelPositions[i][0], self.__labelPositions[i][1], self.__labels[i])

        #Print the canvas
        self.__canvas.Print(self.__Name + '.png')

        return


    #-----------------------------------------------------------
    def getTGraph(self, index):
        return self.__TGraphs[index]


    #-----------------------------------------------------------
    def distance(self, x1, y1, x2, y2):
        """Calculates the cartesian distance between two points"""
        x = x2-x1
        y = y2-y1

        return sqrt(x**2 + y**2)


    #-----------------------------------------------------------
    def area(self, x1, y1, x2, y2):
        """Calculates the cartesian distance between two points"""
        x = abs(x2-x1)
        y = abs(y2-y1)

        if abs(x-y) != 0:
            return (x*y)/abs(x-y)
        else:
            return 1


    #-----------------------------------------------------------
    def getMinMax(self):
        """Get Maxmimum and Minimum coordinates for the graphs"""

        xmax = max(self.__xhigh)
        xmin = min(self.__xlow)
        ymax = max(self.__yhigh)
        ymin = min(self.__ylow)

        yran = ymax - ymin
        xran = xmax - xmin

        self.__xmax = xmax + 0.05*xran
        self.__xmin = xmin - 0.05*xran
        self.__ymax = ymax + 0.05*yran
        self.__ymin = ymin - 0.05*yran

        return


    #-----------------------------------------------------------
    def convertToRelCanvasCoord(self, x, y):
        """Convert absolute points on the plots to canvas relative coordinates"""

        canvasX1 = self.__canvas.GetX1()
        canvasX2 = self.__canvas.GetX2()
        canvasY1 = self.__canvas.GetY1()
        canvasY2 = self.__canvas.GetY2()

        print canvasX1, canvasY1, canvasX2, canvasY2

        canvasXrange = canvasX2 - canvasX1
        canvasYrange = canvasY2 - canvasY1

        plotX1 = self.__canvas.GetFrame().GetX1()
        plotX2 = self.__canvas.GetFrame().GetX2()
        plotY1 = self.__canvas.GetFrame().GetY1()
        plotY2 = self.__canvas.GetFrame().GetY2()

        print plotX1, plotY1, plotX2, plotY2

        plotXrange = plotX2 - plotX1
        plotYrange = plotY2 - plotY1

        xAbs = x*plotXrange + plotX1
        yAbs = y*plotYrange + plotY1

        xRel = (xAbs - canvasX1)/canvasXrange
        yRel = (yAbs - canvasY1)/canvasYrange

        return xRel, yRel

