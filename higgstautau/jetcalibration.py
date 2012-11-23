"""
See instructions here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Jets_NEW
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetCalibrationToolsWinter2011
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetCalibrationToolsForPhysicsAnalyses
"""
from atlastools import datasets
from rootpy.tree.filtering import EventFilter

# ATLAS tools imports
from externaltools import ApplyJetCalibration
from ROOT import JetCalibrationTool


class JetCalibration(EventFilter):
    """
    This class applies jet calibrations
    """
    def __init__(self, datatype, year,
            verbose=False,
            passthrough=False,
            **kwargs):

        super(JetCalibration, self).__init__(
                passthrough=passthrough,
                **kwargs)

        if not passthrough:
            self.year = year
            self.datatype = datatype
            self.isdata = datatype in (datasets.DATA, datasets.EMBED)
            self.verbose = verbose
            self.algo = 'AntiKt4LCTopo'
            if year == 2011:
                config = 'InsituJES_2011_Preliminary.config'
                # use Rel17_JES_AFII.config for fastsim MC
            elif year == 2012:
                config = 'JES_August2012.config'
                # use JES_August2012_AFII.config for fastsim MC
            else:
                raise ValueError(
                        "No JES calibration defined for year %d" % year)
            print "Using JES config %s" % config
            self.config_file = ApplyJetCalibration.get_resource(
                    'CalibrationConfigs/%s' % config)
            self.jes_tool = JetCalibrationTool(
                    self.algo,
                    self.config_file,
                    self.isdata)

    def passes(self, event):

        # For the pile-up correction, we need mu and NPV(2+ tracks)
        mu = event.averageIntPerXing

        # count the number of vertices with 2 or more tracks
        NPV = 0
        for vertex in event.vertices:
            if vertex.nTracks >= 2:
                NPV += 1

        if self.verbose:
            print "JETS BEFORE RECALIBRATION"
            self.print_jets(event)

        if self.year == 2011:
            for jet in event.jets:
                Eraw    = jet.constscale_E
                eta_det = jet.constscale_eta
                eta     = jet.EtaOrigin
                phi     = jet.PhiOrigin
                m       = jet.MOrigin
                ## Record old values
                jet.EtaOriginEM = jet.pt
                jet.PhiOriginEM = jet.E
                #calib_jet = self.jes_tool.ApplyOffsetEtaJES(
                #    Eraw, eta_det, eta, phi, m, mu, NPV)
                calib_jet = self.jes_tool.ApplyEtaJES(
                    Eraw, eta_det, eta, phi, m)
                jet.E = calib_jet.E()
                jet.m = calib_jet.M()
                jet.pt = calib_jet.Pt()
                jet.eta = calib_jet.Eta()
                jet.phi = calib_jet.Phi()

        elif self.year == 2012:
            for jet in event.jets:
                Eraw    = jet.constscale_E
                eta_det = jet.constscale_eta
                eta     = jet.constscale_eta
                phi     = jet.constscale_phi
                m       = jet.constscale_m
                ## Record old values
                jet.EtaOriginEM = jet.pt
                jet.PhiOriginEM = jet.E
                calib_jet = self.jes_tool.ApplyOffsetEtaJES(
                    Eraw, eta_det, eta, phi, m, mu, NPV)
                jet.E = calib_jet.E()
                jet.m = calib_jet.M()
                jet.pt = calib_jet.Pt()
                jet.eta = calib_jet.Eta()
                jet.phi = calib_jet.Phi()
        else:
            raise ValueError('Invalid year in jet calibration: %d' % self.year)

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
