# branches to ignore in the input tree and not include in the output tree

REMOVE = [
    "vxp_trk_weight",

    "cl_*",
    "ph_*",

    # these large mc branches are useless since
    # they contain barcodes and not indices
    # use mc_parent_index and mc_child_index
    "mc_children",
    "mc_parents",
    "mc_barcode",
    # don't need vertices
    "mc_vx_*",

    #"jet_AntiKt4TopoEM_*", <== jet cleaning recommendation is with TopoEM jets
    #"jet_AntiKt4LCTopo_*",  <== need these for MET systematics
    "jet_AntiKt6*",
    #"jet_flavor_*",  <== need these for systematics...
    "jet_*Assoc*",

    "tau_otherTrk_*",
    "tau_cell_*",
    "tau_cluster_*",

    # don't need the following trigger info
    "EF_2e*",
    "EF_2mu*",
    "EF_2j*",
    "EF_xe*",
    "EF_xs*",
    "EF_e*",
    "EF_mu*",
    "EF_MU*",
    "EF_g*",
    "EF_j*",
    "EF_g*",
    "L1_*",
    "L2_*",
    "trig_L2_*",
    "trig_EF_trigmuonef*",
    "trig_EF_trigmugirl*",
    "EF_b*",
    "EF_2b*",
    "*_TileMu_*",
    "trig_L1_*",
    "trig_EF_tau_seedCalo*",
    "trig_EF_tau_otherTrk*",
    "trig_EF_tau_track*",
    "trig_RoI_L2_*",
    "trig_roidescriptor_*",
    "trig_EF_topocl_*",
    "trig_EF_jet*",
    "trig_EF_feb_*",
    "trig_EF_bjet_*",
    "trig_RoI_EF_j_*",
    "trig_RoI_EF_b_*",
    "trig_RoI_EF_mu_*",

    #"tau_jet_*",

    "efo_*",

    "isCalibration",
    "isSimulation",
    "isTestBeam",

    "egtruth_*",

    # don't need pantau info
    "tau_pantau_*",

    # extra tau into we don't need
    "tau_seedCalo_track_*",
    "tau_seedCalo_wideTrk_*",

    "muonTruth*",
    # needed for JVF syst ONLY IN THE SIGNAL SAMPLES
    #"jet_antikt4truth_*",
    "collcand_*",

    "el_*",
    "mu_*",
    "MET_*Reg*",
    # corrupt branch in some d3pds
    "MET_CorrTopo_*",
    "trk_z0_wrtBL",
    "trk_err_z0_wrtBS",
    "trig_EF_met_*",
    "MET_LocHadTopo_*",
    "trk_theta_qoverp_err_wrtBL",
    "trk_err_phi_wrtBS",
    "tau_jet_e_TileExt0",

    #"mcevt_pdf*",
]

# override globs above
KEEP = [
    "el_n",
    "el_cl_E",
    "el_tracketa",
    "el_trackphi",
    "el_author",
    "el_charge",
    "el_loosePP",
    "el_mediumPP",
    "el_tightPP",
    "el_OQ",
    # required for electron ID recalc
    "el_cl_eta",
    "el_cl_phi",
    "el_m",
    "el_deltaeta1",
    "el_deltaphi2",
    "el_Emax2",
    "el_emaxs1",
    "el_etas2",
    "el_Ethad",
    "el_Ethad1",
    "el_expectHitInBLayer",
    "el_f1",
    "el_f3",
    "el_isEM",
    "el_nBLayerOutliers",
    "el_nBLHits",
    "el_nPixelOutliers",
    "el_nPixHits",
    "el_nSCTOutliers",
    "el_nSiHits",
    "el_nTRTHits",
    "el_nTRTOutliers",
    "el_reta",
    "el_trackd0_physics",
    "el_trackqoverp",
    "el_TRTHighTOutliersRatio",
    "el_weta2",
    "el_wstot",

    "mu_staco_n",
    "mu_staco_E",
    "mu_staco_pt",
    "mu_staco_eta",
    "mu_staco_phi",
    "mu_staco_loose",
    "mu_staco_medium",
    "mu_staco_tight",
    "mu_staco_isCombinedMuon",
    "mu_staco_isSegmentTaggedMuon",
    "mu_staco_expectBLayerHit",
    "mu_staco_nBLHits",
    "mu_staco_nPixHits",
    "mu_staco_nPixelDeadSensors",
    "mu_staco_nSCTHits",
    "mu_staco_nSCTDeadSensors",
    "mu_staco_nPixHoles",
    "mu_staco_nSCTHoles",
    "mu_staco_nTRTHits",
    "mu_staco_nTRTOutliers",

    "tau_jet_jvtxf",
]

# additional branches to remove only in the output tree

REMOVE_OUTPUT = [
    'trk_*',
    'tau_jet_*',
    'mc_pt',
    'mc_phi',
    'mc_eta',
    'mc_m',
    'mc_child_index',
    'mc_parent_index',
    'mc_pdgId',
    'mc_charge',
    'mc_status',
    'mc_n',
    'jet_antikt4truth_*',
    'mc_event_weight',
    'mcevt_*',
]

# override REMOVE_OUTPUT above
KEEP_OUTPUT = []
