#!/usr/bin/python

import os
import glob
import sys
sys.path.append( os.path.join( os.getcwd(), '..' ) )
from progressbar import ProgressBar, Percentage, AnimatedMarker
from Configuration import *

#Parse channel argument
DoE  = False
DoMu = False

try:
    if sys.argv[1] == 'e' : DoE  = True
    if sys.argv[1] == 'mu': DoMu = True
except IndexError:
    print 'Please provide a channel (\'e\' or \'mu\')'
    sys.exit()

#Load cross sections
LoadXsections = True

Xsections = {}

if LoadXsections:
    f = open('Xsections.txt')
    lines=  f.readlines()
    for line in lines:
        arguments = line.split(' ')
        Xsections[arguments[0]] = [arguments[1].rstrip(' '), arguments[2].rstrip('\n'), arguments[3].rstrip('\n'), arguments[4].rstrip('\n')]

#Specify where the files are, and what they are called
path = None
dataPath = None
pattern = None
lepton = None

if DoE:
    path     = ePathToProcessorFiles
    dataPath = ePathToMCSkims
    pattern  = eFileNamePattern
    lepton   = eLeptonTag
 
if DoMu:
    path     = muPathToProcessorFiles
    dataPath = muPathToMCSkims
    pattern  = muFileNamePattern
    lepton   = muLeptonTag
    

#Specigy other parameters
tab = '    '
mcType = 'mc11c'

#Retrieve the files
cwd = os.getcwd()
os.chdir(path)
ProcessorFiles = glob.glob(pattern)
ProcessorFiles.sort()
os.chdir(cwd)

#Specify the output file
output = open(cwd + '/' + lepton + 'Samples.py', 'w')
sampleList = []

#Write file header
output.write('from Sample import Sample\n\n')

os.chdir(dataPath)

#Progress Bar
pbar = ProgressBar(widgets=['Progress : ', Percentage(), ' ', AnimatedMarker(markers='_,.-~*\'`\'*~-.,')], maxval=len(ProcessorFiles)).start()
counter = 0

#Add samples
for f in ProcessorFiles:
    tag = f.partition('Processor.')[2].rstrip('.root')
    sampleName = tag.rstrip('.data')
    sampleTreeNameTest = tag + '_' + lepton + 'lh_test'
    sampleTreeNameTrain = tag + '_' + lepton + 'lh_train'
    samplePath = path + '/' + f
    sampleType = 'DATA'
    if f.find(mcType) != -1:
        sampleType = 'MC'
        sampleName = tag.rstrip('.' + mcType)

    if sampleType == 'DATA': continue

    #Get search string for pyAMI
    inputSkimCandidates = glob.glob('*' + sampleName + '*')

    # Get highest dataset version
    D3PDName = ''
    SkimNames = ''
    version = 0
    
    for c in inputSkimCandidates:
        D3PDNameItems = c.partition('mtm.LHSkim.')[2].split('.v')
        SkimVersion = int(D3PDNameItems[1][0])
        if SkimVersion > version:
            D3PDName = D3PDNameItems[0]
            version = SkimVersion
            SkimName = 'user.mtm.LHSkim.' + D3PDName + '.v' + str(version) + '/'

    sampleList.append(sampleName)

    output.write(sampleName + ' = Sample(\n')
    output.write(tab + 'path          = \'' + samplePath + '\',\n')
    output.write(tab + 'treeNameTest  = \'' + sampleTreeNameTest + '\',\n')
    output.write(tab + 'treeNameTrain = \'' + sampleTreeNameTrain + '\',\n')
    output.write(tab + 'sampleType    = \'' + sampleType + '\',\n')
    output.write(tab + 'D3PDName      = \'' + D3PDName + '\',\n')
    output.write(tab + 'SkimName      = \'' + SkimName + '\',\n')
    output.write(tab + 'name          = \'' + sampleName + '\',\n')
    if LoadXsections:
        output.write(tab + 'D3PDnEvents    = ' + Xsections[sampleName][0] + ',\n')
        output.write(tab + 'Xsection       = ' + Xsections[sampleName][1] + ',\n')
        output.write(tab + 'FilterEff      = ' + Xsections[sampleName][2] + ',\n')
        output.write(tab + 'kFactor        = ' + Xsections[sampleName][3] + ',\n')
        
    output.write(tab + ')\n\n')

    pbar.update(counter + 1)
    counter += 1

output.write('Data = Sample(\n')
output.write(tab + 'path          = \'' + path + '/' + lepton + 'LHProcessor.data.root\',\n')
output.write(tab + 'treeNameTest  = \'data_' + lepton + 'lh_test\',\n')
output.write(tab + 'treeNameTrain = \'data_' + lepton + 'lh_train\',\n')
output.write(tab + 'sampleType    = \'DATA\',\n')
output.write(tab + 'name          = \'data\',\n')
if LoadXsections:
    output.write(tab + 'Lumi        = 1,\n')
    output.write(tab + 'D3PDnEvents = 1,\n')
    output.write(tab + 'Xsection    = 1,\n')
output.write(tab + ')\n\n')

output.write('sampleList = [\n')

# Print the sample list
for sample in sampleList:
    output.write(tab + sample + ',\n')

output.write(tab + ']')
    

output.close()
        
