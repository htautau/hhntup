"""
See instructions here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Jets_NEW
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetCalibrationToolsWinter2011
"""
from atlastools import datasets
from rootpy.tree.filtering import EventFilter

# ATLAS tools imports
from externaltools import ApplyJetCalibration
from ROOT import JetCalibrationTool


class JetCalibration(EventFilter):
    """
    This class applies jet calibrations and allows all events to pass through
    """
    def __init__(self, datatype, year, verbose=False, **kwargs):

        super(JetCalibration, self).__init__(**kwargs)
        self.year = year
        self.datatype = datatype
        self.isdata = datatype in (datasets.DATA, datasets.EMBED)
        self.verbose = verbose
        self.algo = "AntiKt4LCTopo"
        if year == 2011:
            self.config_file = ApplyJetCalibration.get_resource(
                    'CalibrationConfigs/InsituJES_2011_Preliminary.config')
        elif year == 2012:
            self.config_file = ApplyJetCalibration.get_resource(
                    'CalibrationConfigs/Rel17_JES.config')
        else:
            raise ValueError("No JES calibration defined for year %d" % year)
        self.jes_tool = JetCalibrationTool(
                self.algo,
                self.config_file,
                self.isdata)

    def passes(self, event):

        # For the pile-up correction, we need mu and NPV(2+ tracks)
        mu = event.averageIntPerXing
        NPV = 0 # count the number of vertices with 2 or more tracks
        for vertex in event.vertices:
            if vertex.nTracks >= 2:
                NPV += 1

        if self.verbose:
            print "JETS BEFORE RECALIBRATION"
            self.print_jets(event)

        for jet in event.jets:
            Eraw    = jet.constscale_E
            eta_det = jet.constscale_eta
            eta     = jet.EtaOrigin
            phi     = jet.PhiOrigin
            m       = jet.MOrigin
            # Pile-up, origin, EtaJES correction applied, i.e. to OFFSET_ORIGIN_ETAJES scale
            calib_jet = self.jes_tool.ApplyOffsetEtaJES(
                Eraw, eta_det, eta, phi, m, mu, NPV)
            jet.E = calib_jet.E()
            jet.m = calib_jet.M()
            jet.pt = calib_jet.Pt()
            jet.eta = calib_jet.Eta()
            jet.phi = calib_jet.Phi()

        if self.verbose:
            print "JETS AFTER RECALIBRATION"
            self.print_jets(event)

        return True

    def print_jets(self, event):

        for i in xrange(event.jet_n):
            # use raw access to D3PD branches as an extra check
            print "E: %.3f pT: %.3f eta: %.5f phi: %.5f M: %.3f" % (
                    event.jet_E[i],
                    event.jet_pt[i],
                    event.jet_eta[i],
                    event.jet_phi[i],
                    event.jet_m[i])
