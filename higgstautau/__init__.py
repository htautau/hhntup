import rootpy
import logging
import os
import ROOT

if not os.getenv('MVA_NO_BATCH', False):
    ROOT.gROOT.SetBatch(True)

rootpy.log.basic_config_colorized()

log = logging.getLogger('higgstautau')
if not os.environ.get("DEBUG", False):
    log.setLevel(logging.INFO)

if hasattr(logging, 'captureWarnings'):
    logging.captureWarnings(True)
