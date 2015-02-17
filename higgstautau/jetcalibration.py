# --
# February 16, 2015: Converted to XAOD !
# --
"""
See instructions here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ApplyJetCalibration2014
"""
from rootpy.tree.filtering import EventFilter

from ROOT import JetCalibrationTool

from . import log; log = log[__name__]
from . import datasets
from . import store_helper


class JetCalibration(EventFilter):
    """
    This class applies jet calibrations
    """
    def __init__(self, datatype, 
                 passthrough=False, **kwargs):
        super(JetCalibration_xAOD, self).__init__(
            passthrough=passthrough,
            **kwargs)

        if not passthrough:
            self.isdata = datatype in (datasets.DATA, datasets.EMBED)
            from ROOT import JetCalibrationTool
            self.jet_calib_tool = JetCalibrationTool('jet_calib_tool')
            self.jet_calib_tool.setProperty('std::string')('JetCollection', 'AntiKt4LCTopo')
            self.jet_calib_tool.setProperty('std::string')('ConfigFile', 'JES_Full2012dataset_May2014.config')
            self.jet_calib_tool.setProperty('std::string')('CalibSequence', 'JetArea_Residual_Origin_EtaJES_GSC')
            self.jet_calib_tool.setProperty('bool')('IsData', self.isdata)
            self.jet_calib_tool.initialize()

    def passes(self, event):
        jets_copy = store_helper.shallowCopyJetContainer(event.jets.collection)
        for jet in jets_copy:
            self.jet_calib_tool.applyCalibration(jet)
        event.jets.collection = jets_copy
        return True

