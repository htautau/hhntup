/////////////////////////////////////////////////////////////////////////////////////
//  Trigger Scale Factor (SF) Tool for multilepton (muons and electrons) analyses  //
/////////////////////////////////////////////////////////////////////////////////////
//  Marilyn Marx (marx@cern.ch)                                                    //
//  Junjie Zhu (junjie@umich.edu)                                                  //
/////////////////////////////////////////////////////////////////////////////////////
// How to use:
//
// pair<double, double> GetTriggerSF(int runnumber, bool useGeV, vector<TLorentzVector> muons, muon_quality q, vector<TLorentzVector> electrons, electron_quality p, int var){
//
// where runnumber is an int associating to a real data runnumber in mc,
// useGeV is a bool you set depending on whether quantities are in MeV (useGeV=false) or GeV (useGeV=true) in your code,
// vector<TLorentzVector> should be filled (separately for muons and electrons) with all good leptons passing your object selection cuts,
// (N.B. make sure CLUSTER eta is used in the electron TLorentzVector)
// muon_quality is an enum, for q enter either loose or combined,
// electron_quality is an enum, for p enter either loosepp, mediumpp or tightpp.
// var is an int for systematic studies, it can be set to 1 or 2 for +1 or -1 sigma variation on the SF
//
// GetTriggerSF().first returns the event SF
// GetTriggerSF().second returns the uncertainty (stat + sys) on the event SF
/////////////////////////////////////////////////////////////////////////////////////

#include "LeptonTriggerSF.h"

#include "TFile.h"
#include "TROOT.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace AnalysisFramework::External::TrigMuonEff;
using AnalysisFramework::External::TrigMuonEff::SFDataPeriod;

namespace AnalysisFramework
{
namespace External
{

void LeptonTriggerSF::closefile(){

  file_Muon_Trig_Eff->Close();
  delete file_Muon_Trig_Eff;

}

void LeptonTriggerSF::initialize() {

  cout << "LeptonTriggerSF : Initializing" << endl;

  file_Muon_Trig_Eff = new TFile((path_to_root_files+"muon_trigger_sf.root").c_str());
  if(file_Muon_Trig_Eff) cout << "LeptonTriggerSF : File " << file_Muon_Trig_Eff->GetName() << " was found." << endl;
  else cout << "LeptonTriggerSF : ERROR Rootfile failed to open." << endl;

  //unsigned int nperiods = 3;
  unsigned int nperiods = 8;
  const TString quality[2] = {"loose","combined"};
  const TString bins[2] = {"fine","coarse"};
  const TString type[2] = {"data","mc"};
  const TString region[2] = {"barrel","endcap"};
  //const TString trigger[3] = {"mu18_MG","mu18_MG_medium","mu18_MG_medium"};
  const TString trigger[8] = {"mu18_MG","mu18_MG_medium","mu18_MG_medium","mu18_MG_medium","mu18_MG_medium","mu18_MG_medium","mu18_MG_medium","mu18_MG_medium"};
  //const TString periods[3] = {"BtoI","JtoMwoL3toL4","L3toL4"};
  const TString periods[8] = {"BtoI","JtoMwoL3toL4","L3toL4","J","K","JtoK","JtoM","LtoM"};
  for(unsigned int iqu = 0 ; iqu < 2; iqu++){
    for(unsigned int ibins = 0; ibins < 2; ibins++){
      for(unsigned int iperiod = 0; iperiod < nperiods; iperiod++){
	for(unsigned int iregion = 0; iregion < 2; iregion++){
	  for(unsigned int itype = 0; itype < 2; itype++){
	    TString histname = ("_MuonTrigEff_" + periods[iperiod] + "_" + quality[iqu] + "_EtaPhi_" + bins[ibins] + "_" + region[iregion] + "_" + type[itype]);
	    _EfficiencyMap[histname] = (TH2F *)(file_Muon_Trig_Eff->Get((quality[iqu] + "/" + bins[ibins] + "/" + trigger[iperiod]+ "_" + periods[iperiod]+"/etaphi_" + region[iregion] + "_eff_" + type[itype] + "_period" + periods[iperiod] + "_EF_" + trigger[iperiod]).Data()));
	  }}}}}

  bool failedinitialise = false;
  m_phi_boundary_barrel = TMath::Pi();
  m_phi_boundary_endcap = TMath::Pi();
  for(map<TString,TH2*>::iterator mit = _EfficiencyMap.begin(); mit !=_EfficiencyMap.end(); mit++){
    if(!mit->second) failedinitialise = true;
    // setting lower phi to account for unusual binning
    if(mit->first.Contains("_MuonTrigEff_JtoMwoL3toL4_loose_EtaPhi_fine_barrel_data")) m_phi_boundary_barrel = mit->second->GetYaxis()->GetXmin();
    if(mit->first.Contains("_MuonTrigEff_JtoMwoL3toL4_loose_EtaPhi_fine_endcap_data")) m_phi_boundary_endcap = mit->second->GetYaxis()->GetXmin();
    if(m_phi_boundary_barrel < -TMath::Pi()) m_phi_boundary_barrel += 2.*TMath::Pi();
    if(m_phi_boundary_endcap < -TMath::Pi()) m_phi_boundary_endcap += 2.*TMath::Pi();
  }

  if(failedinitialise){
    cout << "Not all histograms could be initialized, exiting." << endl;
    exit (1);
  }
  else cout << "LeptonTriggerSF : Initialization of LeptonTriggerSF successful " << endl;

  obj_egammaSF = new egammaSFclass();

  return;
}

/// Constructor
LeptonTriggerSF::LeptonTriggerSF(const string path) {
  path_to_root_files = path;
  if (path_to_root_files.empty()) { // default to InstallArea/share for the files if running in Athena
    char *tmparea = getenv("TestArea");
    if (tmparea != NULL) {
      path_to_root_files =  string(tmparea) + "/InstallArea/share/";
      cout << "LeptonTriggerSF : Using default dir: "<< path_to_root_files << endl;
    }
    else cout << "You are not running in Athena but also did not set a path to the root file." << endl;
  }
  else cout << "LeptonTriggerSF : Using user defined path: " << path_to_root_files << endl;
  inGeV = false;
  initialize();
}

/// Destructor
LeptonTriggerSF::~LeptonTriggerSF() {
  closefile();

}


std::pair<double, double> LeptonTriggerSF::GetTriggerSF(int runnumber, bool useGeV, std::vector<TLorentzVector> muons, muon_quality q, std::vector<TLorentzVector> electrons, electron_quality p, int var){

  SFDataPeriod period = getDataPeriod(runnumber);
  if(period == perUnDefined) {
    cout << "RunNumber is not in 2011 dataset. No scale factors calculated. Please use Runnumber between 178044-191933" << endl;
    return (make_pair(0.,0.));
  }

  setThresholds(useGeV, runnumber);

  int set_data = decide_el_quality(runnumber, p, true);
  int set_mc = decide_el_quality(runnumber, p, false);

  double rate_not_fired_data = 1.;
  double rate_not_fired_mc = 1.;

  // needed for uncertainty calculation
  double sq_err_eff_data = 0.;
  double sq_err_eff_mc = 0.;
  double inv_sq_eff_data = 0.;
  double inv_sq_eff_mc = 0.;

  for(unsigned int ielec = 0; ielec < electrons.size(); ielec ++) {
    double eff_data = 0., eff_mc = 0.;
    double err_data = 0., err_mc = 0.;

    if(electrons[ielec].Pt() < elThreshold){ // add 2.47 eta cut for electrons here plus information message for user?
      eff_data=0;
      eff_mc=0;
      err_data=0;
      err_mc=0;
    }

    else{
      // get efficiency from data
      eff_data= ElEff_Data(electrons[ielec],set_mc, set_data).first;
      err_data= ElEff_Data(electrons[ielec],set_mc, set_data).second;
      // get efficiency from MC
      eff_mc=  ElEff_MC(electrons[ielec],set_mc).first;
      err_mc=  ElEff_MC(electrons[ielec],set_mc).second;
    }

    rate_not_fired_data *= (1-eff_data);
    rate_not_fired_mc *= (1-eff_mc);

    // needed for uncertainty calculation
    sq_err_eff_data += pow(err_data,2);
    sq_err_eff_mc += pow(err_mc,2);
    inv_sq_eff_data += pow((1-eff_data),-2);
    inv_sq_eff_mc += pow((1-eff_mc),-2);
  }

  for(unsigned int imuon = 0; imuon < muons.size(); imuon ++) {
    double eff_data = 0., eff_mc = 0.;
    double err_data = 0., err_mc = 0.;

    if(muons[imuon].Pt() < muThreshold){
      eff_data=0;
      eff_mc=0;
      err_data=0;
      err_mc=0;
    }

    else{
      // get efficiency from data
      eff_data= MuEff(period,true, muons[imuon], q).first;
      err_data= MuEff(period,true, muons[imuon], q).second;
      // get efficiency from MC
      eff_mc=  MuEff(period,false, muons[imuon], q).first;
      err_mc=  MuEff(period,false, muons[imuon], q).second;
    }

    rate_not_fired_data *= (1-eff_data);
    rate_not_fired_mc *= (1-eff_mc);

    // needed for uncertainty calculation
    sq_err_eff_data += err_data*err_data;
    sq_err_eff_mc += err_mc*err_mc;
    inv_sq_eff_data += pow((1-eff_data),-2);
    inv_sq_eff_mc += pow((1-eff_mc),-2);
  }

  double event_SF = 1.;
  double event_SF_err = 0.;

  // prevent events with no triggered electrons or muons
  if ((electrons.size() != 0 || muons.size() != 0) && ((1-rate_not_fired_mc) != 0)){
    event_SF = (1-rate_not_fired_data)/(1-rate_not_fired_mc);
    event_SF_err = SFerror(rate_not_fired_data, rate_not_fired_mc, sq_err_eff_data, sq_err_eff_mc, inv_sq_eff_data, inv_sq_eff_mc);
  }

  if(var == 1) return std::pair<double, double> (event_SF + event_SF_err, event_SF_err);
  else if(var == 2) return std::pair<double, double> (event_SF - event_SF_err, event_SF_err);
  else return std::pair<double, double> (event_SF, event_SF_err);
}

std::pair<double, double> LeptonTriggerSF::MuEff(SFDataPeriod period, bool isData, TLorentzVector muon, int mu_quality) const{

  const double mu_eta = muon.Eta();
  double mu_phi = check_Phi_range(muon.Phi());
  //double mu_pt=muon.Pt();
  //if(!inGeV) mu_pt=mu_pt/1000.;

  // fix phi range for unusual binning
  if(fabs(mu_eta) < 1.05){
    if(m_phi_boundary_barrel < 0) {
      if(mu_phi < m_phi_boundary_barrel) mu_phi += 2*TMath::Pi();
    } else {
      if(mu_phi >= m_phi_boundary_barrel) mu_phi -= 2*TMath::Pi();
    }
  }
  else{
    if(m_phi_boundary_endcap < 0) {
      if(mu_phi < m_phi_boundary_endcap) mu_phi += 2*TMath::Pi();
    } else {
      if(mu_phi >= m_phi_boundary_endcap) mu_phi -= 2*TMath::Pi();
    }
  }

  TString type, region, dataperiod, bins, quality;
  type = (isData) ? "_data" : type = "_mc";
  region = (fabs(mu_eta) < 1.05) ? "_barrel" :  region = "_endcap" ;
  if(period == per2011B_I){
    dataperiod = "BtoI_";
    bins = "fine";
  }
  if(period == per2011J_MwoL3_L4){
    dataperiod = "JtoMwoL3toL4_";
    bins = "fine";
  }
  if(period == per2011L3_L4){
    dataperiod = "L3toL4_";
    bins = "coarse";
  }
  if(period == per2011J){
    dataperiod = "J_";
    bins = "coarse";
  }
  if(period == per2011K){
    dataperiod = "K_";
    bins = "coarse";
  }
  if(period == per2011J_K){
    dataperiod = "JtoK_";
    bins = "coarse";
  }
  if(period == per2011J_M){
    dataperiod = "JtoM_";
    bins = "coarse";
  }
  if(period == per2011L_M){
    dataperiod = "LtoM_";
    bins = "coarse";
  }
  quality = decide_mu_quality(mu_quality);

  Int_t bin = -1;
  double eff = 0;
  double stat = 0;

  TString hist = "_MuonTrigEff_" + dataperiod + quality + "_EtaPhi_" + bins + region + type;
  map<TString,TH2*>::const_iterator mapit = _EfficiencyMap.find(hist);
  if(mapit == _EfficiencyMap.end()){
    cout << "Could not find what you are looking for in the efficiency map. This is a bug." << endl;
    exit(1);
  }
  bin = mapit->second->FindBin(mu_eta, mu_phi);
  eff = mapit->second->GetBinContent(bin);
  stat = mapit->second->GetBinError(bin);

  double syst = commonSystMTSG*eff;
  double err  = sqrt(stat*stat + syst*syst);

  return std::pair<double, double> ( eff, err );

}

std::pair<double, double> LeptonTriggerSF::ElEff_MC(TLorentzVector electron, int set_mc) const{

  const double el_eta = electron.Eta();
  double el_ET = electron.Pt();
  if(inGeV) el_ET *= 1000.0;

  std::pair<float,float> mc_eff_error = obj_egammaSF->scaleFactor(el_eta, el_ET, set_mc, 0, 6, 1);

  double eff = mc_eff_error.first;
  double err = mc_eff_error.second;

  return std::pair<double, double> ( eff, err );
}

std::pair<double, double> LeptonTriggerSF::ElEff_Data(TLorentzVector electron, int set_mc, int set_data) const{

  const double el_eta = electron.Eta();
  double el_ET = electron.Pt();
  if(inGeV) el_ET *= 1000.0;

  std::pair<float,float> mc_eff_error = obj_egammaSF->scaleFactor(el_eta, el_ET, set_mc, 0, 6, 1);
  std::pair<float,float> sf_error = obj_egammaSF->scaleFactor(el_eta, el_ET, set_data, 0, 6, 1);

  double eff = mc_eff_error.first * sf_error.first;
  double err = eff*( sqrt(pow(mc_eff_error.second/mc_eff_error.first,2) + pow(sf_error.second/sf_error.first,2)) );

  return std::pair<double, double> ( eff, err );
}

TString LeptonTriggerSF::decide_mu_quality(int mu_quality) const{

  TString mu_q;
  if(mu_quality == 0) mu_q = "loose";
  else if(mu_quality == 1) mu_q = "combined";
  else{
    cout << "mu_quality has to be 0 for loose or 1 for combined muons." << endl;
    exit(1);
  }

  return mu_q;
}

int LeptonTriggerSF::decide_el_quality(int runnumber, int el_quality, bool isSF){

  int set = 0;

  if(runnumber <= 186755){ // e20_medium
    // loose++
    if(el_quality == 0){
      if(isSF) set = 8;
      else set = 20;
    }
    // medium++
    if(el_quality == 1){
      if(isSF) set = 8;
      else set = 9;
    }
    // tight++
    if(el_quality == 2){
      if(isSF) set = 10;
      else set = 11;
    }
  }
  else if(runnumber <= 187815){ // e22_medium
    if(el_quality == 0){
      if(isSF) set = 12;
      else set = 21;
    }
    if(el_quality == 1){
      if(isSF) set = 12;
      else set = 13;
    }
    if(el_quality == 2){
      if(isSF) set = 14;
      else set = 15;
    }
  }
  else if(runnumber <= 191933){ // e22vh_medium1
    if(el_quality == 0){
      if(isSF) set = 16;
      else set = 22;
    }
    if(el_quality == 1){
      if(isSF) set = 16;
      else set = 17;
    }
    if(el_quality == 2){
      if(isSF) set = 18;
      else set = 19;
    }
  }
  else cout << "Sorry, not a good runnumber!" << endl;

  return set;

}

void LeptonTriggerSF::setThresholds(bool useGeV, int runnumber){

  muThreshold = 20000.;
  if(runnumber <= 186755) elThreshold = 21000.;
  else if(runnumber <= 191933) elThreshold = 23000.;
  else cout << "This is not a good runnumber for setting electron plateau thresholds." << endl;

  if(useGeV){
    inGeV = true;
    muThreshold = 20.;
    if(runnumber <= 186755) elThreshold = 21.;
    else if(runnumber <= 191933) elThreshold = 23.;
    else cout << "This is not a good runnumber for setting electron plateau thresholds." << endl;
  }
  return;
}

SFDataPeriod LeptonTriggerSF::getDataPeriod(Int_t runNumber){

  if(runNumber < 178044) return perUnDefined;
  if(runNumber <= 186493) return per2011B_I;
  if(runNumber <= 189090 || (runNumber >= 189639 && runNumber <= 191933)) return per2011J_MwoL3_L4;
  if(runNumber >= 189184 && runNumber <= 189610) return per2011L3_L4;
  else cout << "Cannot associate this run number to a period." << endl;

  return perUnDefined;
}

// Helper function to deal with possible phi ambiguity
double LeptonTriggerSF::check_Phi_range(double phi) const  {

  double newphi = phi;
  const double pi = TMath::Pi();

  if (newphi > pi) {
    printf("<muon_SF>: WARNING: Muon phi %4.2f > pi! Using (phi-2*pi) \n", phi);
    newphi -= 2*pi;
  }
  if (newphi < -pi) {
    printf("<muon_SF>: WARNING: Muon phi %4.2f < -pi! Using (phi+2*pi) \n", phi);
    newphi += 2*pi;
  }

  return newphi;
}

double LeptonTriggerSF::SFerror(double a, double b, double c, double d, double e, double f){

  double err_SF = 0;
  err_SF = sqrt(((a*a)/(pow(1-b,2))) * c * e + ((b*b*pow(1-a,2))/(pow(1-b,4))) * d * f);

  return err_SF;
}

}
}
