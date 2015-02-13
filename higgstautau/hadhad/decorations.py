
SF_DEFAULT = 1.

def decorate_tau(tau):

    # truth matching
    tau.matched = False
    tau.matched_dR = 9999.

    tau.id = 0
    
    tau.centrality = 0.
    tau.centrality_boosted = 0.

    # vertex association
    tau.vertex_prob = 0.

    # overlap checking
    tau.min_dr_jet = 9999.

    # efficiency scale factor if matches truth
    tau.id_sf = SF_DEFAULT
    tau.id_sf_high = SF_DEFAULT
    tau.id_sf_low = SF_DEFAULT
    tau.id_sf_stat_high = SF_DEFAULT
    tau.id_sf_stat_low = SF_DEFAULT
    tau.id_sf_sys_high = SF_DEFAULT
    tau.id_sf_sys_low = SF_DEFAULT

    # trigger efficiency
    tau.trigger_sf = SF_DEFAULT
    tau.trigger_sf_high = SF_DEFAULT
    tau.trigger_sf_low = SF_DEFAULT
    tau.trigger_sf_mc_stat_high = SF_DEFAULT
    tau.trigger_sf_mc_stat_low = SF_DEFAULT
    tau.trigger_sf_data_stat_high = SF_DEFAULT
    tau.trigger_sf_data_stat_low = SF_DEFAULT
    tau.trigger_sf_stat_high = SF_DEFAULT
    tau.trigger_sf_stat_low = SF_DEFAULT
    tau.trigger_sf_sys_high = SF_DEFAULT
    tau.trigger_sf_sys_low = SF_DEFAULT

    tau.trigger_sf_stat_scale_high = SF_DEFAULT
    tau.trigger_sf_stat_scale_low = SF_DEFAULT

    tau.trigger_eff = SF_DEFAULT
    tau.trigger_eff_high = SF_DEFAULT
    tau.trigger_eff_low = SF_DEFAULT
    tau.trigger_eff_stat_high = SF_DEFAULT
    tau.trigger_eff_stat_low = SF_DEFAULT
    tau.trigger_eff_sys_high = SF_DEFAULT
    tau.trigger_eff_sys_low = SF_DEFAULT

    tau.trigger_eff_stat_scale_high = SF_DEFAULT
    tau.trigger_eff_stat_scale_low = SF_DEFAULT

    # fake rate scale factor for taus that do not match truth
    tau.fakerate_sf = SF_DEFAULT
    tau.fakerate_sf_high = SF_DEFAULT
    tau.fakerate_sf_low = SF_DEFAULT

    # fake rate reco scale factor for taus that do not match truth
    tau.fakerate_sf_reco = SF_DEFAULT
    tau.fakerate_sf_reco_high = SF_DEFAULT
    tau.fakerate_sf_reco_low = SF_DEFAULT

    # colliniear mass approx
    tau.collinear_momentum_fraction = -9999.

    # track recounting
    tau.numTrack_recounted = -1

    # BCH cleaning
    tau.BCHMedium = False
    tau.BCHTight = False
