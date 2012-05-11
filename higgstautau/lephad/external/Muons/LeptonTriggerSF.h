#ifndef MYLEPTONTRIGGERSF_H__
#define MYLEPTONTRIGGERSF_H__

#include <vector>
#include <map>
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "../Electrons/egammaSFclass.h"

const double commonSystMTSG = 0.01;

namespace AnalysisFramework
{
namespace External
{

namespace TrigMuonEff{
  enum SFDataPeriod {
    perUnDefined=-1,
    per2011B_I=0,
    per2011J_MwoL3_L4=1,
    per2011L3_L4=2,

    //for HSG3 specific use
    per2011J=3,
    per2011K=4,
    per2011J_K=5,
    per2011J_M=6,
    per2011L_M=7
  };
}

enum muon_quality{
  loose=0,
  combined=1
};

enum electron_quality{
  loosepp=0,
  mediumpp=1,
  tightpp=2
};

class LeptonTriggerSF {

 public:

  ~LeptonTriggerSF();

  LeptonTriggerSF(std::string path="");

  void closefile();

  void initialize();

  std::map<TString,TH2*> _EfficiencyMap;

  std::pair<double, double> GetTriggerSF(int runnumber, bool useGeV, std::vector<TLorentzVector> muons, muon_quality q, std::vector<TLorentzVector> electrons, electron_quality p, int var = 0);

  std::pair<double, double> MuEff(TrigMuonEff::SFDataPeriod period, bool isData, TLorentzVector muon, int mu_quality) const;
  std::pair<double, double> ElEff_MC (TLorentzVector electron, int set_mc) const;
  std::pair<double, double> ElEff_Data (TLorentzVector electron, int set_mc, int set_data) const;

  void setThresholds(bool useGeV, int runnumber);
  TrigMuonEff::SFDataPeriod getDataPeriod(Int_t runNumber);
  double check_Phi_range(double phi) const;
  double SFerror(double a, double b, double c, double d, double e, double f);

  TFile* file_Muon_Trig_Eff;


 private:

  std::string path_to_root_files;

  bool inGeV;
  double muThreshold;
  double elThreshold;
  double m_phi_boundary_barrel;
  double m_phi_boundary_endcap;

  TString decide_mu_quality(int mu_q) const;
  int decide_el_quality(int runnumber, int el_quality, bool isSF);

  egammaSFclass *obj_egammaSF;

};

}
}
#endif
