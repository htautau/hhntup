"""
See instructions here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetCalibrationToolsForPhysicsAnalyses
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ApplyJetCalibration2011
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ApplyJetCalibration2012
"""
from rootpy.tree.filtering import EventFilter
from rootpy import ROOTError

from externaltools import ApplyJetCalibration
from ROOT import JetCalibrationTool

from . import log; log = log[__name__]
from . import datasets


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
                config = 'JES_Full2012dataset_Preliminary_Jan13.config'
                # use JES_Full2012dataset_Preliminary_AFII_Jan13.config for fastsim MC
            else:
                raise ValueError(
                    "No JES calibration defined for year %d" % year)
            log.info("Using JES config %s" % config)
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
                calib_jet = self.jes_tool.ApplyOffsetEtaJES(
                    Eraw, eta_det, eta, phi, m, mu, NPV)
                jet.E = calib_jet.E()
                jet.m = calib_jet.M()
                jet.pt = calib_jet.Pt()
                jet.eta = calib_jet.Eta()
                jet.phi = calib_jet.Phi()

        elif self.year == 2012:
            rho = event.Eventshape_rhoKt4LC
            for jet in event.jets:
                Eraw    = jet.constscale_E
                eta     = jet.constscale_eta
                phi     = jet.constscale_phi
                m       = jet.constscale_m
                Ax      = jet.ActiveAreaPx
                Ay      = jet.ActiveAreaPy
                Az      = jet.ActiveAreaPz
                Ae      = jet.ActiveAreaE
                try:
                    calib_jet = self.jes_tool.ApplyJetAreaOffsetEtaJES(
                        Eraw, eta, phi, m, Ax, Ay, Az, Ae, rho, mu, NPV)
                    jet.E = calib_jet.E()
                    jet.m = calib_jet.M()
                    jet.pt = calib_jet.Pt()
                    jet.eta = calib_jet.Eta()
                    jet.phi = calib_jet.Phi()
                except ROOTError as e:
                    print "JET ERROR"
                    print "Run: ", event.RunNumber
                    print "Event: ", event.EventNumber
                    print "Eraw", Eraw
                    print "eta", eta
                    print "phi", phi
                    print "m", m
                    print "Ax", Ax
                    print "Ay", Ay
                    print "Az", Az
                    print "Ae", Ae
                    raise e
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
