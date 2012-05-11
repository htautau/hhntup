#include <iostream>
#include <cmath>

#include "CaloIsoCorrection.h"


//-----------------------------------------------------------------------
//
// Internal Methods
//
namespace AnalysisFramework
{
namespace External
{

using namespace CaloIsoCorrection;

const float GeV = 1000.;

namespace CaloIsoCorrection{
  float GetPtCorrection       (float energy, float etaS2, float etaPointing, float etaCluster, 
                               float radius, bool is_mc=true, 
                               bool isConversion=false, ParticleType parttype = ELECTRON, Version ver = REL17);
  float GetPtCorrectionTopo   (float energy, float etaS2, float etaPointing, float etaCluster, 
                               float radius, bool is_mc=true, 
                               bool isConversion=false, ParticleType parttype = ELECTRON, Version ver = REL17);
  float GetPtCorrectionError  (float energy, float etaS2, float etaPointing, float etaCluster, 
                               float radius, bool is_mc=true, ParticleType parttype = ELECTRON, Version ver = REL17);

  int GetRadius(float radius);
  int GetEtaBinFine(float eta);
  int GetEtaBinCoarse(float eta);
  
  float GetPtCorrectionFactor(float eta,  
                            const float* mc_leakage_corrections_ptr = 0, 
                            const float* data_leakage_corrections_ptr = 0);
  float GetPtCorrectionValue(float energy, float etaPointing, float etaCluster, float scale_factor);

  float GetNPVCorrection (float eta, float radius, unsigned int nPV, bool is_mc, ParticleType parttype = ELECTRON, Version ver = REL17);
  float GetNPVCorrectionError (float eta, float radius, unsigned int nPV, bool is_mc, ParticleType parttype = ELECTRON, Version ver = REL17);

  float GetEDCorrection (float Etcone40, float Etcone40_ED_corrected, float radius);

  // array sizes
  const unsigned int nBinsEtaFine = 10;
  const unsigned int nBinsEtaCoarse = 5;
  
  // eta bins
  const float eta_bins_fine[nBinsEtaFine+1] = {0.0, 0.10, 0.60, 0.80, 1.15, 1.37, 1.52, 1.81, 2.01, 2.37, 2.47};
  const float eta_bins_coarse[nBinsEtaCoarse+1] = {0.0, 0.60, 1.37, 1.52, 1.81, 2.47};

// ----------------------------------------------------------------------------
// ------------- 2010 leakage corrections -------------------------------------
  // OLD isolation corrections: fine grained in eta, derived from MC08
//   const float mc_leakage_corrections_photon_20[] = {0.014582 , 0.013887 , 0.013786, 0.017952 , 0.020799 , 0.022892 , 0.022947 , 0.027274, 0.034607 , 0.04111  , 0.043765 };
//   const float mc_leakage_corrections_photon_30[] = {0.0149955, 0.0143672, 0.01439 , 0.0189413, 0.0221763, 0.0251387, 0.0244524, 0.028928, 0.0373025, 0.044983 , 0.0450644};
//   const float mc_leakage_corrections_photon_40[] = {0.0150641, 0.0144416, 0.014464, 0.0191335, 0.0224708, 0.0256225, 0.0248054, 0.029223, 0.0378147, 0.0457925, 0.0457752};
// 
//   const float mc_leakage_corrections_electron_20[] = {0.0146779, 0.0140476, 0.0142736, 0.0187919, 0.0242244, 0.0327271, 0.0236067, 0.0273382, 0.0347824, 0.0405783, 0.0488965};
//   const float mc_leakage_corrections_electron_30[] = {0.0151632, 0.0145416, 0.0149607, 0.0198783, 0.0258172, 0.0404185, 0.0257027, 0.0289142, 0.0375322, 0.0451606, 0.0555181};
//   const float mc_leakage_corrections_electron_40[] = {0.0152929, 0.0146146, 0.0150541, 0.0200918, 0.0261464, 0.043179 , 0.0262671, 0.0292846, 0.0380784, 0.0462397, 0.056945 };

  // OLD additional isolation corrections: coarse grained in eta, measured with 2010 data
//  const float data_leakage_corrections_photon_40  [] = {.0030, .0044, 0., .0012, .0041};
//   const float data_leakage_corrections_electron_20[] = {.0040, .0074, 0., .0075, .0064};
//   const float data_leakage_corrections_electron_30[] = {.0041, .0072, 0., .0094, .0023};
//   const float data_leakage_corrections_electron_40[] = {.0042, .0072, 0., .0111, .0038};

  // OLD additional isolation correction errors: coarse grained in eta, measured with 2010 data
//   const float data_leakage_corrections_error_photon_40  [] = {.0013, .0012, 0., .0019, .0014};
//   const float data_leakage_corrections_error_electron_20[] = {.0010, .0011, 0., .0047, .0019};
//   const float data_leakage_corrections_error_electron_30[] = {.0018, .0016, 0., .0069, .0026};
//   const float data_leakage_corrections_error_electron_40[] = {.0026, .0033, 0., .0062, .0038};
// ----------------------------------------------------------------------------
  
// -------------------------------------------------------------------------------------------
// ------------- early 2011 (rel 16) leakage corrections -------------------------------------
  // OLD isolation corrections: fine grained in eta, derived from MC10
  const float mc_rel16_leakage_correction_slopes_electron_15[] = {0.0134623, 0.0125705, 0.0127194, 0.016388 , 0.0192139, 0.0191176, 0.019365 , 0.0224323, 0.0276738, 0.0306204, 0.0395363};
  const float mc_rel16_leakage_correction_slopes_electron_20[] = {0.0146765, 0.0137598, 0.0141393, 0.0185445, 0.0222504, 0.0303701, 0.022529 , 0.0255273, 0.0322996, 0.0367724, 0.0488965};
  const float mc_rel16_leakage_correction_slopes_electron_25[] = {0.0150446, 0.0141208, 0.0146059, 0.0192928, 0.0233765, 0.0355231, 0.0239549, 0.0266049, 0.0340472, 0.0395156, 0.0533984};
  const float mc_rel16_leakage_correction_slopes_electron_30[] = {0.0151549, 0.0142345, 0.0147978, 0.0195921, 0.0238402, 0.0374921, 0.0245635, 0.0270403, 0.0347653, 0.0407656, 0.0555181};
  const float mc_rel16_leakage_correction_slopes_electron_35[] = {0.0151991, 0.0142886, 0.0148807, 0.0197355, 0.0240916, 0.0389585, 0.024845 , 0.0272151, 0.0350951, 0.0414101, 0.0565319};
  const float mc_rel16_leakage_correction_slopes_electron_40[] = {0.0152502, 0.014324 , 0.014922 , 0.019808 , 0.0242069, 0.0398394, 0.0250683, 0.0273238, 0.0352379, 0.0417484, 0.056945 };

  const float mc_rel16_leakage_correction_slopes_photon_converted_15[] = {0.0124144, 0.0117418, 0.0122051, 0.0159856, 0.0189189, 0.0190803, 0.0195599, 0.0230894, 0.0287442, 0.0335009, 0.040739 };
  const float mc_rel16_leakage_correction_slopes_photon_converted_20[] = {0.0134282, 0.0127686, 0.0134831, 0.0180318, 0.0218801, 0.030747 , 0.0227755, 0.0262832, 0.0335691, 0.0401866, 0.043765 };
  const float mc_rel16_leakage_correction_slopes_photon_converted_25[] = {0.0136567, 0.0130497, 0.0138831, 0.0187274, 0.0229529, 0.0360431, 0.0242219, 0.0274141, 0.0354121, 0.0432104, 0.0450998};
  const float mc_rel16_leakage_correction_slopes_photon_converted_30[] = {0.0137181, 0.0131058, 0.0140011, 0.0189988, 0.023384 , 0.038064 , 0.0248426, 0.0278485, 0.0361593, 0.044579 , 0.0450644};
  const float mc_rel16_leakage_correction_slopes_photon_converted_35[] = {0.0137829, 0.0131081, 0.01406  , 0.019133 , 0.0235932, 0.0396043, 0.0251306, 0.0280451, 0.0365044, 0.0452737, 0.0452159};
  const float mc_rel16_leakage_correction_slopes_photon_converted_40[] = {0.0137756, 0.0131087, 0.0140897, 0.0192113, 0.0237036, 0.0405692, 0.0253519, 0.0281581, 0.0366884, 0.0456197, 0.0457752};

  const float mc_rel16_leakage_correction_slopes_photon_unconverted_15[] = {0.0136592, 0.0127361, 0.0122583, 0.0152527, 0.0170659, 0.0157067, 0.0195983, 0.0242494, 0.0287865, 0.033302 , 0.040739 };
  const float mc_rel16_leakage_correction_slopes_photon_unconverted_20[] = {0.0148505, 0.0139408, 0.0135913, 0.0171323, 0.0196672, 0.0254563, 0.021963 , 0.0274632, 0.0335511, 0.039913 , 0.043765 };
  const float mc_rel16_leakage_correction_slopes_photon_unconverted_25[] = {0.0152095, 0.0143039, 0.0140069, 0.0177627, 0.0206302, 0.029655 , 0.0229128, 0.0285693, 0.0353589, 0.0428317, 0.0450998};
  const float mc_rel16_leakage_correction_slopes_photon_unconverted_30[] = {0.0153482, 0.0144291, 0.0141718, 0.018018 , 0.0210249, 0.0312614, 0.0232916, 0.028999 , 0.0361068, 0.0442073, 0.0450644};
  const float mc_rel16_leakage_correction_slopes_photon_unconverted_35[] = {0.0153847, 0.0144704, 0.0142212, 0.0181336, 0.0212057, 0.032474 , 0.0234635, 0.0291766, 0.0364465, 0.0448465, 0.0452159};
  const float mc_rel16_leakage_correction_slopes_photon_unconverted_40[] = {0.0154202, 0.0144914, 0.0142582, 0.0181816, 0.0213156, 0.0331942, 0.0235738, 0.0292826, 0.0365929, 0.0451769, 0.0457752};

  // OLD additional isolation corrections: coarse grained in eta, measured with 2011 periods B-E
  const float data_rel16_leakage_correction_slopes_photon_20  [] = {.0032, .0073, 0., .0046, .0082};
  const float data_rel16_leakage_correction_slopes_photon_30  [] = {.0036, .0086, 0., .0044, .0086};
  const float data_rel16_leakage_correction_slopes_photon_40  [] = {.0040, .0101, 0., .0048, .0079};
   
  const float data_rel16_leakage_correction_slopes_electron_20[] = {.0042, .0087, 0., .0128, .0114};
  const float data_rel16_leakage_correction_slopes_electron_30[] = {.0058, .0119, 0., .0142, .0130};
  const float data_rel16_leakage_correction_slopes_electron_40[] = {.0062, .0168, 0., .0170, .0147};

  // OLD additional isolation correction errors: coarse grained in eta, measured with 2011 periods B-E
  const float data_rel16_leakage_correction_slopes_error_photon_20  [] = {.0005, .0009, 0., .0012, .0013};
  const float data_rel16_leakage_correction_slopes_error_photon_30  [] = {.0006, .0012, 0., .0016, .0017};
  const float data_rel16_leakage_correction_slopes_error_photon_40  [] = {.0010, .0019, 0., .0020, .0019};
  
  const float data_rel16_leakage_correction_slopes_error_electron_20[] = {.0008, .0012, 0., .0024, .0016};
  const float data_rel16_leakage_correction_slopes_error_electron_30[] = {.0029, .0018, 0., .0032, .0025};
  const float data_rel16_leakage_correction_slopes_error_electron_40[] = {.0040, .0026, 0., .0049, .0027};  
// -------------------------------------------------------------------------------------------
  
// -------------------------------------------------------------------------------------------
// ------------- full 2011 (rel 17) leakage corrections --------------------------------------
  // CURRENT isolation corrections: fine grained in eta, derived from MC11
  const float mc_rel17_leakage_correction_slopes_electron_15[] = { 0.01466, 0.01421, 0.01427, 0.01872, 0.02008, 0.02181, 0.02141, 0.02636, 0.03285, 0.03564 };
  const float mc_rel17_leakage_correction_slopes_electron_20[] = { 0.01616, 0.01572, 0.01596, 0.02139, 0.02339, 0.03350, 0.02499, 0.03038, 0.03870, 0.04337 };
  const float mc_rel17_leakage_correction_slopes_electron_25[] = { 0.01663, 0.01616, 0.01667, 0.02223, 0.02445, 0.03852, 0.02611, 0.03174, 0.04093, 0.04692 };
  const float mc_rel17_leakage_correction_slopes_electron_30[] = { 0.01689, 0.01635, 0.01684, 0.02256, 0.02485, 0.04223, 0.02660, 0.03232, 0.04182, 0.04846 };
  const float mc_rel17_leakage_correction_slopes_electron_35[] = { 0.01695, 0.01642, 0.01693, 0.02268, 0.02501, 0.04403, 0.02685, 0.03254, 0.04223, 0.04928 };
  const float mc_rel17_leakage_correction_slopes_electron_40[] = { 0.01701, 0.01646, 0.01702, 0.02272, 0.02517, 0.04550, 0.02698, 0.03267, 0.04242, 0.04964 };
  const float mc_rel17_leakage_correction_slopes_electron_tp40[] = { 0.01470, 0.01420, 0.01450, 0.02080, 0.02400, 0.0, 0.02630, 0.03020, 0.03980, 0.04850 };
  
  const float mc_rel17_leakage_correction_offsets_electron_15[] = { 21.71, 36.00, 132.56, 191.64, 263.46, 619.58, 288.75, 121.92, 102.35, 175.02 };
  const float mc_rel17_leakage_correction_offsets_electron_20[] = { 76.67, 85.35, 184.90, 276.72, 384.97, 595.60, 657.99, 231.88, 170.48, 312.30 };
  const float mc_rel17_leakage_correction_offsets_electron_25[] = { 90.44, 105.02, 267.95, 420.38, 555.09, 1014.50, 765.83, 283.50, 224.18, 357.30 };
  const float mc_rel17_leakage_correction_offsets_electron_30[] = { 119.55, 127.09, 279.48, 430.96, 571.81, 846.86, 968.01, 354.46, 263.11, 455.21 };
  const float mc_rel17_leakage_correction_offsets_electron_35[] = { 138.79, 161.87, 371.47, 572.08, 754.17, 1249.38, 1000.44, 389.22, 295.72, 464.28 };
  const float mc_rel17_leakage_correction_offsets_electron_40[] = { 180.08, 187.89, 363.31, 553.46, 707.60, 1006.20, 1105.73, 434.21, 312.78, 535.90 };
  const float mc_rel17_leakage_correction_offsets_electron_tp40[] = { -371.80, -338.90, -284.00, -322.10, -238.80, 0.0, -437.70, -477.50, -538.10, -491.90 };
  
  // photons: last eta bin isn't used
  const float mc_rel17_leakage_correction_slopes_photon_converted_15[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_converted_20[] = { 0.01600, 0.01560, 0.01580, 0.02130, 0.02320, 0.0, 0.02450, 0.03000, 0.03840, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_converted_25[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_converted_30[] = { 0.01630, 0.01600, 0.01630, 0.02240, 0.02460, 0.0, 0.02610, 0.03190, 0.04150, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_converted_35[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_converted_40[] = { 0.01610, 0.01590, 0.01620, 0.02250, 0.02480, 0.0, 0.02650, 0.03220, 0.04210, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_converted_tp40[] = { 0.01490, 0.01430, 0.01450, 0.02090, 0.02370, 0.0, 0.02590, 0.02970, 0.03930, 0.0 };
  
  const float mc_rel17_leakage_correction_offsets_photon_converted_15[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_converted_20[] = { 93.30, 101.40, 270.60, 369.10, 514.70, 0.0, 586.10, 160.80, 193.80, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_converted_25[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_converted_30[] = { 286.60, 241.60, 501.60, 570.70, 739.50, 0.0, 860.00, 264.50, 270.40, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_converted_35[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_converted_40[] = { 478.70, 383.80, 679.20, 725.70, 938.70, 0.0, 999.30, 322.80, 316.20, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_converted_tp40[] = { -432.10, -368.20, -284.10, -294.40, -141.20, 0.0, -455.70, -470.10, -474.10, 0.0 };
  
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_15[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_20[] = { 0.01630, 0.01560, 0.01560, 0.02060, 0.02240, 0.0, 0.02480, 0.03060, 0.03830, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_25[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_30[] = { 0.01690, 0.01630, 0.01640, 0.02170, 0.02360, 0.0, 0.02630, 0.03260, 0.04140, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_35[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_40[] = { 0.01710, 0.01640, 0.01650, 0.02190, 0.02390, 0.0, 0.02660, 0.03290, 0.04200, 0.0 };
  const float mc_rel17_leakage_correction_slopes_photon_unconverted_tp40[] = { 0.01470, 0.01410, 0.01410, 0.01940, 0.02230, 0.0, 0.02550, 0.02990, 0.03860, 0.0 };
  
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_15[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_20[] = { -17.70, 12.60, 97.80, 126.50, 186.20, 0.0, 200.80, 24.00, 62.30, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_25[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_30[] = { 5.40, 33.80, 141.60, 199.50, 295.40, 0.0, 336.50, 64.80, 90.40, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_35[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_40[] = { 13.30, 62.00, 177.00, 267.10, 406.20, 0.0, 419.80, 89.40, 105.90, 0.0 };
  const float mc_rel17_leakage_correction_offsets_photon_unconverted_tp40[] = { -373.70, -347.40, -308.10, -372.90, -394.30, 0.0, -461.00, -469.00, -493.30, 0.0 };
// ------------------------------------------------------------------------------------------
  
// -------------------------------------------------------------------------------------------
// ------------- early 2011 (rel 16) pileup corrections --------------------------------------

  // nPV pileup linear corrections, derived on data (2011, periods D,E)
  // units are MeV
  // >= 2 tracks/vertex
  const float rel16_npv_correction_slopes_20[] = { 53.,  59., 0.,  70.,  76.};
  const float rel16_npv_correction_slopes_30[] = {145., 163., 0., 175., 189.};
  const float rel16_npv_correction_slopes_40[] = {280., 306., 0., 331., 351.};

  // errors (relative) for nPV pileup corrections
  // - dominant source of errors are: run-to-run and bunch-to-bunch variation in 2011 data
  const float rel16_npv_correction_slopes_relative_error_20[] = {0.38, 0.38, 0.0, 0.38, 0.38};
  const float rel16_npv_correction_slopes_relative_error_30[] = {0.28, 0.38, 0.0, 0.38, 0.38};
  const float rel16_npv_correction_slopes_relative_error_40[] = {0.30, 0.30, 0.0, 0.30, 0.30};
// -------------------------------------------------------------------------------------------
  
  
// -------------------------------------------------------------------------------------------
// ------------- full 2011 (rel 17) pileup corrections --------------------------------------
  // all corrections use >= 2 tracks/vertex
  
  // linear corrections from Zee rel 17 full 2011 data
  // - not used in favor of quadratic corrections below
//   const float rel17_npv_correction_slopes_20[] = {066., 075., 0., 065., 075.};
//   const float rel17_npv_correction_slopes_30[] = {163., 191., 0., 173., 196.};
//   const float rel17_npv_correction_slopes_40[] = {315., 357., 0., 333., 366.};
 
  // quadratic corrections from rel 17 full 2011 data  
  const float rel17_npv_correction_slopes_20[] = {  42.,  52., 0.0,  58.,  45. };
  const float rel17_npv_correction_slopes_30[] = { 140., 137., 0.0, 152., 119. };
  const float rel17_npv_correction_slopes_40[] = { 253., 248., 0.0, 287., 227. };

  const float rel17_npv_correction_quads_20[] = { 1.65, 1.66, 0.0, 0.467, 2.19 };
  const float rel17_npv_correction_quads_30[] = { 1.74, 3.63, 0.0, 1.61, 5.39 };
  const float rel17_npv_correction_quads_40[] = { 4.26, 7.05, 0.0, 3.50, 9.55 };

// -------------------------------------------------------------------------------------------
  
}
//-----------------------------------------------------------------------
  

//-----------------------------------------------------------------------
// User function
// Returns the nPV corrected isolation
//
float CaloIsoCorrection::GetNPVCorrectedIsolation(unsigned int nPV,
                                                  float etaS2,
                                                  float radius,
                                                  bool is_mc,
                                                  float Etcone_value,
                                                  ParticleType parttype,
                                                  Version ver) {
  
  float npv_correction = GetNPVCorrection(etaS2, radius, nPV, is_mc, parttype, ver);
  return Etcone_value - npv_correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the ED corrected isolation
//
float CaloIsoCorrection::GetEDCorrectedIsolation(float Etcone40,
                                                 float Etcone40_ED_corrected,
                                                 float radius,
                                                 float Etcone_value,
                                                 ParticleType /*parttype*/,
                                                 Version /*ver*/) {
  
  float ED_correction = GetEDCorrection(Etcone40, Etcone40_ED_corrected, radius);
  return Etcone_value - ED_correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the pt leakage corrected isolation
//
float CaloIsoCorrection::GetPtCorrectedIsolation(float energy, 
                                                 float etaS2, 
                                                 float etaPointing, 
                                                 float etaCluster, 
                                                 float radius, 
                                                 bool is_mc, 
                                                 float Etcone_value,
                                                 bool isConversion,
                                                 ParticleType parttype,
                                                 Version ver){
  
  float pt_correction = GetPtCorrection(energy, etaS2, etaPointing, etaCluster, radius, is_mc, isConversion, parttype, ver);
  return Etcone_value - pt_correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the pt leakage corrected isolation
//
float CaloIsoCorrection::GetPtCorrectedTopoIsolation(float energy, 
                                                     float etaS2, 
                                                     float etaPointing, 
                                                     float etaCluster, 
                                                     float radius, 
                                                     bool is_mc, 
                                                     float Etcone_value,
                                                     bool isConversion,
                                                     ParticleType parttype,
                                                     Version ver){
  
  float pt_correction = GetPtCorrectionTopo(energy, etaS2, etaPointing, etaCluster, radius, is_mc, isConversion, parttype, ver);
  return Etcone_value - pt_correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the pt leakage + nPV corrected isolation
//
float CaloIsoCorrection::GetPtNPVCorrectedIsolation(unsigned int nPV,
                                                    float energy, 
                                                    float etaS2, 
                                                    float etaPointing, 
                                                    float etaCluster, 
                                                    float radius, 
                                                    bool is_mc, 
                                                    float Etcone_value,
                                                    bool isConversion,
                                                    ParticleType parttype,
                                                    Version ver){
  
  float pt_correction = GetPtCorrection(energy, etaS2, etaPointing, etaCluster, radius, is_mc, isConversion, parttype, ver);
  float npv_correction = GetNPVCorrection(etaS2, radius, nPV, is_mc, parttype);
  return Etcone_value - pt_correction - npv_correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the pt leakage + energy density corrected isolation
//
float CaloIsoCorrection::GetPtEDCorrectedIsolation(float Etcone40,
                                                   float Etcone40_ED_corrected,
                                                   float energy, 
                                                   float etaS2, 
                                                   float etaPointing, 
                                                   float etaCluster, 
                                                   float radius, 
                                                   bool is_mc, 
                                                   float Etcone_value,
                                                   bool isConversion,
                                                   ParticleType parttype,
                                                   Version ver){
  
  float pt_correction = GetPtCorrection(energy, etaS2, etaPointing, etaCluster, radius, is_mc, isConversion, parttype, ver);
  float ED_correction = GetEDCorrection(Etcone40, Etcone40_ED_corrected, radius);
  return Etcone_value - pt_correction - ED_correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the error on the nPV pileup correction
//
float CaloIsoCorrection::GetNPVCorrectedIsolationError(unsigned int nPV,
                                                       float etaS2,
                                                       float radius,
                                                       bool is_mc,
                                                       ParticleType parttype,
                                                       Version ver){
  
  return GetNPVCorrectionError(etaS2, radius, nPV, is_mc, parttype, ver);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// User function
// Returns the error on the data-derived pT leakage correction
//
float CaloIsoCorrection::GetPtCorrectedIsolationError(float energy, 
                                                      float etaS2, 
                                                      float etaPointing, 
                                                      float etaCluster, 
                                                      float radius, 
                                                      bool is_mc,
                                                      ParticleType parttype,
                                                      Version ver){
  
  return GetPtCorrectionError(energy, etaS2, etaPointing, etaCluster, radius, is_mc, parttype, ver);
}
//-----------------------------------------------------------------------




//***************************************************************************************************************************
//***************************************************************************************************************************
// Internal Functions Below.
//***************************************************************************************************************************
//***************************************************************************************************************************

//-----------------------------------------------------------------------
// Internal function
// Gets the pt correction factor
//
float CaloIsoCorrection::GetPtCorrection(float energy, 
                                             float etaS2, float etaPointing, float etaCluster, 
                                             float radius, 
                                             bool is_mc, 
                                             bool isConversion, ParticleType parttype, Version ver){
                                             
  int newrad = GetRadius(radius);
  const float* mc_correction_slopes_ptr(0);
  const float* mc_correction_offsets_ptr(0);
  const float* data_correction_slopes_ptr(0);
  switch(newrad){
    case 15: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_15; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_15;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_converted_15; 
          }
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_15; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_15;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_unconverted_15; 
          }
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_15;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_15;
        } else if (ver == REL16) {
          mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_electron_15;
        }
        
      }
      break;
    case 20: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_20; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_20;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_converted_20; 
          }
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_20; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_20;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_unconverted_20; 
          }
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_20;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_20;
        } else if (ver == REL16) {
          mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_electron_20;
        }
        
      }
      if (!is_mc && (ver == REL16)) {
        if (parttype == PHOTON) {
          // photons
          data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_photon_20; 
        } else {
          // electrons
          data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_electron_20; 
        }
      }
      break;
    case 25: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_25; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_25;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_converted_25; 
          }
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_25; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_25;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_unconverted_25; 
          }
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_25;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_25;
        } else if (ver == REL16) {
          mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_electron_25;
        }
        
      }
      break;
    case 30: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_30; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_30;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_converted_30; 
          }
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_30; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_30;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_unconverted_30; 
          }
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_30;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_30;
        } else if (ver == REL16) {
          mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_electron_30;
        }
        
      }
      if (!is_mc && (ver == REL16)) {
        if (parttype == PHOTON) {
          // photons
          data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_photon_30; 
        } else {
          // electrons
          data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_electron_30; 
        }
      }
      break;
    case 35: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_35; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_35;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_converted_35; 
          }
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_35; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_35;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_unconverted_35; 
          }
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_35;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_35;
        } else if (ver == REL16) {
          mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_electron_35;
        }
        
      }
      break;
    case 40: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_40; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_40;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_converted_40; 
          }
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_40; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_40;
          } else if (ver == REL16) {
            mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_photon_unconverted_40; 
          }
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_40;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_40;
        } else if (ver == REL16) {
          mc_correction_slopes_ptr = (float*)mc_rel16_leakage_correction_slopes_electron_40;
        }
        
      }
      if (!is_mc && (ver == REL16)) {
        if (parttype == PHOTON) {
          // photons
          data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_photon_40; 
        } else {
          // electrons
          data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_electron_40; 
        }
      }
      break;
    default:
      std::cerr << "Unable to retrieve leakage correction for cone with radius = " << radius << "." << std::endl
                << "--- Radii must be one of {.15, .20, .25, .30, .35, .40} OR {15, 20, 25, 30, 35, 40}." << std::endl;
      return 0.;
  }

  float scale_factor = GetPtCorrectionFactor(etaS2, mc_correction_slopes_ptr, data_correction_slopes_ptr);
  float offset = 0.;
  if (mc_correction_offsets_ptr) offset = GetPtCorrectionFactor(etaS2, mc_correction_offsets_ptr);
  
  return offset + GetPtCorrectionValue(energy, etaPointing, etaCluster, scale_factor);

}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Internal function
// Gets the pt correction factor for topoIso
//
float CaloIsoCorrection::GetPtCorrectionTopo(float energy, 
                                             float etaS2, float etaPointing, float etaCluster, 
                                             float radius, 
                                             bool, 
                                             bool isConversion, ParticleType parttype, Version ver){
  if (ver != REL17) return 0.;
  int newrad = GetRadius(radius);
  const float* mc_correction_slopes_ptr(0);
  const float* mc_correction_offsets_ptr(0);
  switch(newrad){
    case 40: 
      if (parttype == PHOTON) {
        // photons
        if(isConversion) {
          // converted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_converted_tp40; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_converted_tp40;
          } 
        } else {
          // unconverted
          if (ver == REL17) {
            mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_photon_unconverted_tp40; 
            mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_photon_unconverted_tp40;
          } 
        }
      } else {
        // electrons
        if (ver == REL17) {
          mc_correction_slopes_ptr = (float*)mc_rel17_leakage_correction_slopes_electron_tp40;
          mc_correction_offsets_ptr = (float*)mc_rel17_leakage_correction_offsets_electron_tp40;
        } 
        
      }
      break;
    default:
      std::cerr << "Unable to retrieve leakage correction for topoIso cone with radius = " << radius << "." << std::endl
                << "--- Only size 40 currently supported." << std::endl;
      return 0.;
  }

  float scale_factor = GetPtCorrectionFactor(etaS2, mc_correction_slopes_ptr);
  float offset = 0.;
  if (mc_correction_offsets_ptr) offset = GetPtCorrectionFactor(etaS2, mc_correction_offsets_ptr);
  
  return offset + GetPtCorrectionValue(energy, etaPointing, etaCluster, scale_factor);

}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Internal function
// Gets the photon correction
//
float CaloIsoCorrection::GetPtCorrectionError(float energy, 
                                                  float etaS2, float etaPointing, float etaCluster, 
                                                  float radius, 
                                                  bool is_mc,
                                                  ParticleType parttype, Version ver){
  if (is_mc || (ver != REL16)) return 0.;
  int newrad = GetRadius(radius);
  const float* data_correction_slopes_ptr(0);
  switch(newrad){
    case 20: 
      if (parttype == PHOTON) data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_error_photon_20; 
      else data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_error_electron_20; 
      break;
    case 30: 
      if (parttype == PHOTON) data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_error_photon_30; 
      else data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_error_electron_30; 
      break;
    case 40: 
      if (parttype == PHOTON) data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_error_photon_40; 
      else data_correction_slopes_ptr = (float*)data_rel16_leakage_correction_slopes_error_electron_40; 
      break;
    default:
      std::cerr << "Unable to retrieve leakage correction errors for cone with radius = " << radius << "." << std::endl
                << "--- Radii must be one of {.20, .30, .40} OR {20, 30, 40}." << std::endl;
      return 0.;
  }

  float scale_factor = GetPtCorrectionFactor(etaS2, 0, data_correction_slopes_ptr);
  
  return GetPtCorrectionValue(energy, etaPointing, etaCluster, scale_factor);

}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Internal function
// Used to retrieve the correct radius
int CaloIsoCorrection::GetRadius(float radius){
  int newrad = 0;
  // avoid roundoff errors by adding 0.1
  if(radius < 1)  newrad = (int)(radius * 100 + 0.1); 
  else            newrad = (int)radius;
  return newrad;
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Internal function
// Used to retrieve the correct (fine) eta bin number
//
int CaloIsoCorrection::GetEtaBinFine(float eta){
  int eta_bin=-1;
  float fabs_eta = fabs(eta);
  for (unsigned int i=0; i < nBinsEtaFine; ++i) {
    if ((fabs_eta >= eta_bins_fine[i]) && (fabs_eta < eta_bins_fine[i+1])) {
      eta_bin = i;
      break;
    }
  }
  return eta_bin;
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Internal function
// Used to retrieve the correct (coarse) eta bin number
//
int CaloIsoCorrection::GetEtaBinCoarse(float eta){
  int eta_bin=-1;
  float fabs_eta = fabs(eta);
  for (unsigned int i=0; i < nBinsEtaCoarse; ++i) {
    if ((fabs_eta >= eta_bins_coarse[i]) && (fabs_eta < eta_bins_coarse[i+1])) {
      eta_bin = i;
      break;
    }
  }
  return eta_bin;
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Internal function
// Returns the appropriate corrections value
//
float CaloIsoCorrection::GetPtCorrectionFactor(float eta,  
                                             const float* mc_leakage_corrections_ptr, 
                                             const float* data_leakage_corrections_ptr){
  if(!mc_leakage_corrections_ptr && !data_leakage_corrections_ptr) 
    return 0.;

  int eta_bin_fine = GetEtaBinFine(eta);
  int eta_bin_coarse = GetEtaBinCoarse(eta);
  
  float correction = 0.;
  if (mc_leakage_corrections_ptr && (eta_bin_fine >= 0)) correction += mc_leakage_corrections_ptr[eta_bin_fine];
  if (data_leakage_corrections_ptr && (eta_bin_coarse >= 0)) correction += data_leakage_corrections_ptr[eta_bin_coarse];

  return correction;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Internal function
// Does the final pt scaling
float CaloIsoCorrection::GetPtCorrectionValue(float energy, float etaPointing, float etaCluster, float scale_factor){
  // apply the correction to et
  double etaForPt = ((fabs(etaPointing - etaCluster) < 0.15) ? etaPointing : etaCluster);
  double et = (fabs(etaForPt)<99.) ? energy/cosh(etaForPt) : 0.;
  
  return scale_factor * et;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Internal function
// Returns an absolute correction based on number of vertices
float CaloIsoCorrection::GetNPVCorrection(float eta, float radius, unsigned int nPV, bool /*is_mc*/, ParticleType /*parttype*/, Version ver){
  
  int newrad = GetRadius(radius);
  const float* npv_correction_slopes_ptr(0);
  const float* npv_correction_quads_ptr(0);
  switch(newrad){
    case 20: 
      if (ver == REL17) {
        npv_correction_slopes_ptr = (float*)rel17_npv_correction_slopes_20; 
        npv_correction_quads_ptr = (float*)rel17_npv_correction_quads_20; 
      } else if (ver == REL16) {
        npv_correction_slopes_ptr = (float*)rel16_npv_correction_slopes_20;
      }
      break;
    case 30: 
      if (ver == REL17) {
        npv_correction_slopes_ptr = (float*)rel17_npv_correction_slopes_30; 
        npv_correction_quads_ptr = (float*)rel17_npv_correction_quads_30; 
      } else if (ver == REL16) {
        npv_correction_slopes_ptr = (float*)rel16_npv_correction_slopes_30;
      }
      break;
    case 40: 
      if (ver == REL17) {
        npv_correction_slopes_ptr = (float*)rel17_npv_correction_slopes_40; 
        npv_correction_quads_ptr = (float*)rel17_npv_correction_quads_40; 
      } else if (ver == REL16) {
        npv_correction_slopes_ptr = (float*)rel16_npv_correction_slopes_40;
      }
      break;
    default:
      std::cerr << "Unable to retrieve npv corrections for electron with radius = " << radius << "." << std::endl
                << "--- Radii must be one of {.20, .30, .40} OR {20, 30, 40}." << std::endl;
      return 0.;
  }
  
  int eta_bin_coarse = GetEtaBinCoarse(eta);
  
  float correction = 0.;
  if (npv_correction_slopes_ptr && (eta_bin_coarse >= 0)) {
    correction = npv_correction_slopes_ptr[eta_bin_coarse]*float(nPV);
    if (npv_correction_quads_ptr) {
      correction += npv_correction_quads_ptr[eta_bin_coarse]*float(nPV*nPV);
//    std::cout << "Applying npv correction for cone size " << newrad << ", eta: " << eta << ": " << npv_correction_slopes_ptr[eta_bin_coarse] << " * npv + " 
//              << npv_correction_quads_ptr[eta_bin_coarse] << " * npv^2 for npv == " << nPV << " gives " << correction << " MeV" << std::endl;
    }
  }
  
  return correction;
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Internal function
// Returns error on correction based on number of vertices
float CaloIsoCorrection::GetNPVCorrectionError(float eta, float radius, unsigned int nPV, bool is_mc, ParticleType parttype, Version ver){
  
  if (ver != REL16) return 0.;
    
  int newrad = GetRadius(radius);
  const float* npv_correction_slopes_err_ptr(0);
  switch(newrad){
    case 20: 
      npv_correction_slopes_err_ptr = (float*)rel16_npv_correction_slopes_relative_error_20; 
      break;
    case 30: 
      npv_correction_slopes_err_ptr = (float*)rel16_npv_correction_slopes_relative_error_30; 
      break;
    case 40: 
      npv_correction_slopes_err_ptr = (float*)rel16_npv_correction_slopes_relative_error_40; 
      break;
    default:
      std::cerr << "Unable to retrieve npv correction errors for radius = " << radius << "." << std::endl
                << "--- Radii must be one of {.20, .30, .40} OR {20, 30, 40}." << std::endl;
      return 0.;
  }
  
  int eta_bin_coarse = GetEtaBinCoarse(eta);
  
  float error = 0.;
  if (npv_correction_slopes_err_ptr && (eta_bin_coarse >= 0)) error = npv_correction_slopes_err_ptr[eta_bin_coarse]*GetNPVCorrection(eta,radius,nPV,is_mc,parttype);
  
  return error;
}
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
// Internal function
// Returns an absolute correction based on Energy Density for cone size 40
float CaloIsoCorrection::GetEDCorrection (float Etcone40, float Etcone40_ED_corrected, float radius){
  int newrad = GetRadius(radius);
  float floatrad = float(newrad)/100.;
  
  // energy density in a cone of 0.4 (minus the core!)
  float ED40 = Etcone40 - Etcone40_ED_corrected;
  if (newrad == 40) return ED40;
  const int neta = 5, nphi = 7;
  const float deta = 0.0125*neta*2;
  const float dphi = M_PI/256.*nphi*2;
  float area40 = M_PI * pow(0.4,2) - (deta * dphi);
  float areanew = M_PI * pow(floatrad,2) - (deta * dphi);
  // use the ratio of areas to calculate the energy density in new cone size
  float EDnew = ED40 * (areanew/area40);  
  return EDnew;
}
//-----------------------------------------------------------------------

}
}

