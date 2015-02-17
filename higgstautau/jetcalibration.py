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
        super(JetCalibration, self).__init__(
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


class JetResolution(EventFilter):
    """
    This class applies the jet smearing in MonteCarlo samples.
    It needs to be called after the JetCalibration since a shallow copy
    is needed to apply the ASG tool.
    """
    def __init__(self, passthrough=False, **kwargs):
        super(JetResolution, self).__init__(
            passthrough=passthrough,
            **kwargs)

        if not passthrough:
            from ROOT import JERTool, JERSmearingTool
            self.jer_tool = JERTool('JERTool')
            self.jer_tool.setProperty('std::string')("PlotFileName", "JetResolution/JERProviderPlots_2012.root")
            self.jer_tool.setProperty('std::string')("CollectionName", "AntiKt4LCTopoJets")
            self.jer_tool.setProperty('std::string')("BeamEnergy", "8TeV")
            self.jer_tool.setProperty('std::string')("SimulationType", "FullSim")
            self.jer_tool.initialize()

            self.jer_smearing_tool = JERSmearingTool('JERSmearingTool')
            self.jer_smearing_tool.setProperty('std::string')('JERToolName', 'JERTool')
            self.jer_smearing_tool.setJERTool(self.jer_tool)
            self.jer_smearing_tool.setNominalSmearing(True)
            self.jer_smearing_tool.initialize()
            
    def passes(self, event):

        for jet in event.jets:
            self.jer_smearing_tool.applyCorrection(jet)
        return True
