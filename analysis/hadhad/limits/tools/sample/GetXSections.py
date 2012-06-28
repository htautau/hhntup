from subprocess import Popen, PIPE
from Samples import *
import sys

useAMIXsections = False

output = open('Xsections.txt', 'w')

#Load the Real cross-sections and k-factors

twikiFile = open('XsectionsTwiki.txt', 'r')

tlines = twikiFile.readlines()

kFactor = 1.0
Data = {}

BR = {}

BR[100] = 8.36E-02
BR[105] = 8.25E-02
BR[110] = 8.02E-02
BR[115] = 7.65E-02
BR[120] = 7.10E-02
BR[125] = 6.37E-02
BR[130] = 5.48E-02
BR[135] = 4.52E-02
BR[140] = 3.54E-02
BR[145] = 2.61E-02
BR[150] = 1.78E-02

for line in tlines:
    params = line.split(' ')

    if line.find('kfac=') != -1:
        kFactor = float(params[1].split('=')[1])

    if line.find('mc11_7TeV') != -1:
        sampleName = params[2].split('.')[2]
        Xsect      = float(params[1].split('*')[0])
        filterEff  = float(params[1].split('*')[1])

        if (line.find('PowHegPythia_ggH') > -1) or (line.find('PowHegPythia_VBFH') > -1):
            mass = int(params[2].split('_')[2].lstrip('ggVBFH'))

            Xsect *= BR[mass]*4.562e-01
            

        Data[sampleName] = [Xsect, filterEff, kFactor]

Total = len(sampleList)
print 'Samples to process : ', Total
count = 0

for sample in sampleList:
    print count, sample.name
    command = 'ami dataset info ' + sample.D3PDName

    child = Popen(args = command, stdout=PIPE, shell = True)
    out = child.communicate()[0]

    lines = out.split('\n')

    nEvents = 0
    Xsection = 0.0

    for line in lines:
        if line.find('totalEvents') > -1:
            nEvents = int(line.split(' ')[-1])
        if line.find('approx_crossSection') > -1:
            Xsection = float(line.split(' ')[-1])
                

    if not useAMIXsections:
        Xsection  = Data[sample.name][0]
    filterEff = Data[sample.name][1]
    kFactor   = Data[sample.name][2]

    output.write(sample.name + ' ' + str(nEvents) + ' ' + str(Xsection) + ' ' + str(filterEff) + ' ' + str(kFactor) + '\n')
    count += 1

output.close

    
