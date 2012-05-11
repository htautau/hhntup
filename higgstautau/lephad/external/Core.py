#Core of functions used in External
#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

import os
filedir = os.path.dirname(os.path.abspath(__file__))

def compileC(filename, target="libs"):
   from ROOT import gSystem
   if not gSystem.CompileMacro(filename, "k-", "", filedir+"/"+target):
       print "@@@@@ ERROR: Failed to compile", filename
       sys.exit()
