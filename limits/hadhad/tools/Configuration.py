from VariableContainer import VariableContainer, PlotInfoContainer
import math

"""
===========================================================================
Paths and patterns
===========================================================================
"""

#Muon Channel
muPathToProcessorFiles = '/cluster/data05/michel/ProcessorFiles/17-06-2012'
muPathToMCSkims = '/global/mtm/data/SKIMS/MC'
muFileNamePattern = 'muLHProcessor*root'
muLeptonTag = 'mu'

#Electron Channel
ePathToProcessorFiles = '/cluster/data03/andres/higgs/higgspy/ProcessorFiles/2012_06_02/'
ePathToMCSkims = '/global/mtm/data/SKIMS/MC'
eFileNamePattern = 'eLHProcessor*root'
eLeptonTag = 'e'



"""
===========================================================================
BDT training parameters
===========================================================================
"""
VBFTrainingParameters = '!H:!V:UseYesNoLeaf=False:NTrees=40:MaxDepth=100000:BoostType=AdaBoost:AdaBoostBeta=0.2:SeparationType=GiniIndex:nCuts=50:PruneMethod=CostComplexity:PruneStrength=-1:PruneBeforeBoost=True:nEventsMin=150'

ggFTrainingParameters = '!H:!V:UseYesNoLeaf=False:NTrees=40:MaxDepth=100000:BoostType=AdaBoost:AdaBoostBeta=0.2:SeparationType=GiniIndex:nCuts=50:PruneMethod=CostComplexity:PruneStrength=-1:PruneBeforeBoost=True:nEventsMin=150'



"""
===========================================================================
BDT input variables
===========================================================================
"""
# Specify Variable name as it appears in the processor ntuple, plus 'F' or 'I' for Float or Integer
        

####################################################################
# VBF BDT variables

#Muon Channel
muVBFVars = VariableContainer()
muVBFVars.Add('dr_tau_muon',              'F')
muVBFVars.Add('higgs_pt',                 'F')
muVBFVars.Add('mass_transverse_met_muon', 'F')
muVBFVars.Add('met_phi_centrality',       'F')
muVBFVars.Add('sphericity',               'F')
muVBFVars.Add('mass_j1_j2',               'F')
muVBFVars.Add('eta_product_j1_j2',        'F')
muVBFVars.Add('muon_centrality_j1_j2',    'F')
muVBFVars.Add('eta_delta_j1_j2',          'F')
muVBFVars.Add('mass_mmc_tau_muon',        'F')

#Electron Channel
eVBFVars = VariableContainer()
eVBFVars.Add('dr_tau_electron',          'F')
eVBFVars.Add('higgs_pt',                 'F')
eVBFVars.Add('mass_transverse_met_electron','F')
eVBFVars.Add('met_phi_centrality',       'F')
eVBFVars.Add('sphericity',               'F')
eVBFVars.Add('mass_j1_j2',               'F')
eVBFVars.Add('eta_product_j1_j2',        'F')
eVBFVars.Add('electron_centrality_j1_j2','F')
eVBFVars.Add('eta_delta_j1_j2',          'F')



####################################################################
# ggF BDT variables

#Muon Channel
muggFVars = VariableContainer()
muggFVars.Add('dr_tau_muon',              'F')
muggFVars.Add('higgs_pt',                 'F')
muggFVars.Add('mass_transverse_met_muon', 'F')
muggFVars.Add('met_phi_centrality',       'F')
muggFVars.Add('HT',                       'F')
muggFVars.Add('mass_mmc_tau_muon',        'F')

#Electron Channel
eggFVars = VariableContainer()
eggFVars.Add('dr_tau_electron',              'F')
eggFVars.Add('higgs_pt',                     'F')
eggFVars.Add('mass_transverse_met_electron', 'F')
eggFVars.Add('met_phi_centrality',           'F')
eggFVars.Add('HT',                           'F')



"""
===========================================================================
VBF Category definitions
===========================================================================
"""

def muVBFCategory(tree):
    """
    Defines the VBF category, the complement of which is the ggF category
    """
    #Jet multiplicity, VBF cut
    if not tree.numJets >= 2: return False

    #Eta product
    #if not tree.eta_product_j1_j2 < 0: return False

    #Leptons between the jets
    if not (tree.muon_centrality_j1_j2 > 1/math.e): return False

    #Delta Eta of the jets
    #if not tree.eta_delta_j1_j2 > 2.0: return False

    #Jet mass
    if not tree.mass_j1_j2 > 100000: return False

    return True


def eVBFCategory(tree):
    """
    Defines the VBF category, the complement of which is the ggF category
    """
    #Jet multiplicity, VBF cut
    if not tree.numJets >= 2: return False

    #Eta product
    #if not tree.eta_product_j1_j2 < 0: return False

    #Leptons between the jets
    if not (tree.electron_centrality_j1_j2 > 1/math.e): return False

    #Delta Eta of the jets
    #if not tree.eta_delta_j1_j2 > 2.0: return False

    #Jet mass
    if not tree.mass_j1_j2 > 100000: return False

    return True



"""
===========================================================================
Variables to Plot
===========================================================================
"""

####################################################################
# VBF variables

#Muon Channel
muVBFPlots = PlotInfoContainer()
muVBFPlots.Add('HT',                       25, 30,    550,  'H_{T} [GeV]', 0.001)
muVBFPlots.Add('mass_j1_j2',               25, 0,     1000, 'm_{j1j2} [GeV]', 0.001)
muVBFPlots.Add('mass_all_jets',            25, 0,     1000, 'm_{all j} [GeV]', 0.001)
muVBFPlots.Add('eta_product_j1_j2',        25, -15,   15,   '#eta_{j1}#times#eta_{j2}')
muVBFPlots.Add('eta_delta_j1_j2',          25, 0,     10,   '#Delta#eta_{j1j2}')
muVBFPlots.Add('tau_centrality_j1_j2',     20, 0,     1,    '#tau #eta centrality j1j2')
muVBFPlots.Add('muon_centrality_j1_j2',    20, 0,     1,    '#mu #eta centrality j1j2')
muVBFPlots.Add('met_phi_centrality',       25, -1.42, 1.42, 'MET #phi centrality')
muVBFPlots.Add('numJets',                  10, 0,     10,   'N. Jets')
muVBFPlots.Add('leadJetPt',                25, 25,    250, 'lead Jet p_{T} [GeV]', 0.001)
muVBFPlots.Add('mass_transverse_met_muon', 25, 0,     200,  'm_{T MET#mu} [GeV]', 0.001)
muVBFPlots.Add('mass_transverse_met_tau',  25, 0,     200,  'm_{T MET#tau} [GeV]', 0.001)
muVBFPlots.Add('MET',                      25, 0,     170,  'missing E_{T} [GeV]', 0.001)
muVBFPlots.Add('cos_theta_tau_muon',       50, -2.1,  1.1,  'cos#theta_{#tau#mu}')
muVBFPlots.Add('tau_BDTJetScore',          25, 0.40,  1,    '#tau jet BDT score')
muVBFPlots.Add('sphericity',               25, 0,     0.8,  'sphericity')
muVBFPlots.Add('aplanarity',               25, 0,     0.2,  'aplanarity')
muVBFPlots.Add('tau_chpiemeovercaloeme',   25, -2,    2,    '(E_{track sys}^{#tau} - E_{HAD}^{#tau})/E_{EM}^{#tau}')
muVBFPlots.Add('tau_j1_j2_phi_centrality', 25, -1.42, 1.42, '#tau #phi centrality')
muVBFPlots.Add('mass_collinear_tau_muon',  25, 0,     400, 'm_{coll}^{#tau#mu} [GeV]', 0.001)
muVBFPlots.Add('mass2_vis_tau_muon',       25, 0,     400, 'm_{vis}^{#tau#mu} [GeV]', 0.001)
muVBFPlots.Add('mass_mmc_tau_muon',        25, 0,     400, 'm_{mmc}^{#tau#mu} [GeV]')
muVBFPlots.Add('pt_mmc_tau_muon',          25, 0,     250, 'Higgs p_{T, mmc} [GeV]')
muVBFPlots.Add('met_mmc_tau_muon',         25, 0,     170,  'missing E_{T, mmc} [GeV]')
muVBFPlots.Add('nvtx',                     25, 0,     25,  'number of vertices')
muVBFPlots.Add('tau_numTrack',             10, 0,     10,  'N_{tracks}^{#tau}')
muVBFPlots.Add('muon_charge',              10, -2,     8,  'muon charge')
muVBFPlots.Add('ddr_tau_muon',             25, 0,     2.0, '#Delta#DeltaR_{#tau#mu}')
muVBFPlots.Add('dr_tau_muon',              25, 0,     6.0, '#DeltaR_{#tau#mu}')
muVBFPlots.Add('higgs_pt',                 25, 0,     250, 'Higgs p_{T} [GeV]')

#Electron Channel
eVBFPlots = PlotInfoContainer()
eVBFPlots.Add('HT',                       25, 30,    550,  'H_{T} [GeV]', 0.001)
eVBFPlots.Add('mass_j1_j2',               25, 0,     1000, 'm_{j1j2} [GeV]', 0.001)
eVBFPlots.Add('mass_all_jets',            25, 0,     1000, 'm_{all j} [GeV]', 0.001)
eVBFPlots.Add('eta_product_j1_j2',        25, -15,   15,   '#eta_{j1}#times#eta_{j2}')
eVBFPlots.Add('eta_delta_j1_j2',          25, 0,     10,   '#Delta#eta_{j1j2}')
eVBFPlots.Add('tau_centrality_j1_j2',     25, 0,     1,    '#tau #eta centrality j1j2')
eVBFPlots.Add('electron_centrality_j1_j2',20, 0,     1,    'e #eta centrality j1j2')
eVBFPlots.Add('met_phi_centrality',       25, -1.42, 1.42, 'MET #phi centrality')
eVBFPlots.Add('numJets',                  10, 0,     10,   'N. Jets')
eVBFPlots.Add('leadJetPt',                25, 25,    250, 'lead Jet p_{T} [GeV]', 0.001)
eVBFPlots.Add('mass_transverse_met_electron', 25, 0,     200,  'm_{T METe} [GeV]', 0.001)
eVBFPlots.Add('mass_transverse_met_tau',  25, 0,     200,  'm_{T MET#tau} [GeV]', 0.001)
eVBFPlots.Add('MET',                      25, 0,     170,  'missing E_{T} [GeV]', 0.001)
eVBFPlots.Add('cos_theta_tau_electron',   50, -2.1,  1.1,  'cos#theta_{#taue}')
eVBFPlots.Add('tau_BDTJetScore',          25, 0.40,  1,    '#tau jet BDT score')
eVBFPlots.Add('sphericity',               25, 0,     0.8,  'sphericity')
eVBFPlots.Add('aplanarity',               25, 0,     0.2,  'aplanarity')
eVBFPlots.Add('tau_chpiemeovercaloeme',   25, -2,    2,    '(E_{track sys}^{#tau} - E_{HAD}^{#tau})/E_{EM}^{#tau}')
eVBFPlots.Add('tau_j1_j2_phi_centrality', 25, -1.42, 1.42, '#tau #phi centrality')
eVBFPlots.Add('mass_collinear_tau_electron',  25, 0,     400, 'm_{coll}^{#taue} [GeV]', 0.001)
eVBFPlots.Add('mass2_vis_tau_electron',       25, 0,     400, 'm_{vis}^{#taue} [GeV]', 0.001)
eVBFPlots.Add('mass_mmc_tau_electron',        25, 0,     400, 'm_{mmc}^{#taue} [GeV]')
eVBFPlots.Add('pt_mmc_tau_muon',          25, 0,     250, 'Higgs p_{T, mmc} [GeV]')
eVBFPlots.Add('met_mmc_tau_muon',         25, 0,     170,  'missing E_{T, mmc} [GeV]')
eVBFPlots.Add('nvtx',                     25, 0,     25,  'number of vertices')
eVBFPlots.Add('tau_numTrack',             10, 0,     10,  'N_{tracks}^{#tau}')
eVBFPlots.Add('muon_charge',              10, -2,     8,  'muon charge')
eVBFPlots.Add('ddr_tau_electron',             25, 0,     2.0, '#Delta#DeltaR_{#taue}')
eVBFPlots.Add('dr_tau_electron',              25, 0,     6.0, '#DeltaR_{#taue}')
eVBFPlots.Add('higgs_pt',                 25, 0,     250, 'Higgs p_T [GeV]')


####################################################################
# ggF variables

#Muon Channel
muggFPlots = PlotInfoContainer()
muggFPlots.Add('HT',                       50, 30,    550,  'H_{T} [GeV]', 0.001)
muggFPlots.Add('mass_j1_j2',               50, 0,     1000, 'm_{j1j2} [GeV]', 0.001)
muggFPlots.Add('mass_all_jets',            50, 0,     1000, 'm_{all j} [GeV]', 0.001)
muggFPlots.Add('eta_product_j1_j2',        50, -15,   15,   '#eta_{j1}#times#eta_{j2}')
muggFPlots.Add('eta_delta_j1_j2',          50, 0,     10,   '#Delta#eta_{j1j2}')
muggFPlots.Add('tau_centrality_j1_j2',     20, 0,     1,    '#tau #eta centrality j1j2')
muggFPlots.Add('muon_centrality_j1_j2',    20, 0,     1,    '#mu #eta centrality j1j2')
muggFPlots.Add('met_phi_centrality',       50, -1.42, 1.42, 'MET #phi centrality')
muggFPlots.Add('numJets',                  10, 0,     10,   'N. Jets')
muggFPlots.Add('leadJetPt',                50, 25,    250, 'lead Jet p_{T} [GeV]', 0.001)
muggFPlots.Add('mass_transverse_met_muon', 50, 0,     200,  'm_{T MET#mu} [GeV]', 0.001)
muggFPlots.Add('mass_transverse_met_tau',  50, 0,     200,  'm_{T MET#tau} [GeV]', 0.001)
muggFPlots.Add('MET',                      50, 0,     170,  'missing E_{T} [GeV]', 0.001)
muggFPlots.Add('cos_theta_tau_muon',       75, -2.1,  1.1,  'cos#theta_{#tau#mu}')
muggFPlots.Add('tau_BDTJetScore',          50, 0.40,  1,    '#tau jet BDT score')
muggFPlots.Add('sphericity',               50, 0,     0.8,  'sphericity')
muggFPlots.Add('aplanarity',               50, 0,     0.2,  'aplanarity')
muggFPlots.Add('tau_chpiemeovercaloeme',   50, -2,    2,    '(E_{track sys}^{#tau} - E_{HAD}^{#tau})/E_{EM}^{#tau}')
muggFPlots.Add('tau_j1_j2_phi_centrality', 50, -1.42, 1.42, '#tau #phi centrality')
muggFPlots.Add('mass_collinear_tau_muon',  50, 0,     400, 'm_{coll}^{#tau#mu} [GeV]', 0.001)
muggFPlots.Add('mass2_vis_tau_muon',       50, 0,     400, 'm_{vis}^{#tau#mu} [GeV]', 0.001)
muggFPlots.Add('mass_mmc_tau_muon',        50, 0,     400, 'm_{mmc}^{#tau#mu} [GeV]')
muggFPlots.Add('pt_mmc_tau_muon',          50, 0,     250, 'Higgs p_{T, mmc} [GeV]')
muggFPlots.Add('met_mmc_tau_muon',         50, 0,     170,  'missing E_{T, mmc} [GeV]')
muggFPlots.Add('nvtx',                     25, 0,     25,  'number of vertices')
muggFPlots.Add('tau_numTrack',             10, 0,     10,  'N_{tracks}^{#tau}')
muggFPlots.Add('muon_charge',              10, -2,     8,  'muon charge')
muggFPlots.Add('ddr_tau_muon',             50, 0,     2.0, '#Delta#DeltaR_{#tau#mu}')
muggFPlots.Add('dr_tau_muon',              50, 0,     6.0, '#DeltaR_{#tau#mu}')
muggFPlots.Add('higgs_pt',                 50, 0,     250, 'Higgs p_T [GeV]')

#Electron Channel
eggFPlots = PlotInfoContainer()
eggFPlots.Add('HT',                       50, 30,    550,  'H_{T} [GeV]', 0.001)
eggFPlots.Add('mass_j1_j2',               50, 0,     1000, 'm_{j1j2} [GeV]', 0.001)
eggFPlots.Add('mass_all_jets',            50, 0,     1000, 'm_{all j} [GeV]', 0.001)
eggFPlots.Add('eta_product_j1_j2',        50, -15,   15,   '#eta_{j1}#times#eta_{j2}')
eggFPlots.Add('eta_delta_j1_j2',          50, 0,     10,   '#Delta#eta_{j1j2}')
eggFPlots.Add('tau_centrality_j1_j2',     20, 0,     1,    '#tau #eta centrality j1j2')
eggFPlots.Add('electron_centrality_j1_j2',    20, 0,     1,    '#mu #eta centrality j1j2')
eggFPlots.Add('met_phi_centrality',       50, -1.42, 1.42, 'MET #phi centrality')
eggFPlots.Add('numJets',                  10, 0,     10,   'N. Jets')
eggFPlots.Add('leadJetPt',                50, 25,    250, 'lead Jet p_{T} [GeV]', 0.001)
eggFPlots.Add('mass_transverse_met_electron', 50, 0,     200,  'm_{T METe} [GeV]', 0.001)
eggFPlots.Add('mass_transverse_met_tau',  50, 0,     200,  'm_{T MET#tau} [GeV]', 0.001)
eggFPlots.Add('MET',                      50, 0,     170,  'missing E_{T} [GeV]', 0.001)
eggFPlots.Add('cos_theta_tau_electron',       75, -2.1,  1.1,  'cos#theta_{#taue}')
eggFPlots.Add('tau_BDTJetScore',          50, 0.40,  1,    '#tau jet BDT score')
eggFPlots.Add('sphericity',               50, 0,     0.8,  'sphericity')
eggFPlots.Add('aplanarity',               50, 0,     0.2,  'aplanarity')
eggFPlots.Add('tau_chpiemeovercaloeme',   50, -2,    2,    '(E_{track sys}^{#tau} - E_{HAD}^{#tau})/E_{EM}^{#tau}')
eggFPlots.Add('tau_j1_j2_phi_centrality', 50, -1.42, 1.42, '#tau #phi centrality')
eggFPlots.Add('mass_collinear_tau_electron',  50, 0,     400, 'm_{coll}^{#taue} [GeV]', 0.001)
eggFPlots.Add('mass2_vis_tau_electron',       50, 0,     400, 'm_{vis}^{#taue} [GeV]', 0.001)
eggFPlots.Add('mass_mmc_tau_electron',        50, 0,     400, 'm_{mmc}^{#taue} [GeV]')
eggFPlots.Add('pt_mmc_tau_muon',          50, 0,     250, 'Higgs p_{T, mmc} [GeV]')
eggFPlots.Add('met_mmc_tau_muon',         50, 0,     170,  'missing E_{T, mmc} [GeV]')
eggFPlots.Add('nvtx',                     25, 0,     25,  'number of vertices')
eggFPlots.Add('tau_numTrack',             10, 0,     10,  'N_{tracks}^{#tau}')
eggFPlots.Add('muon_charge',              10, -2,     8,  'muon charge')
eggFPlots.Add('ddr_tau_electron',             50, 0,     2.0, '#Delta#DeltaR_{#taue}')
eggFPlots.Add('dr_tau_electron',              50, 0,     6.0, '#DeltaR_{#taue}')
eggFPlots.Add('higgs_pt',                 50, 0,     250, 'Higgs p_T [GeV]')
