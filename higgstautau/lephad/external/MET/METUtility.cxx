//see header for instruction
#include<cmath>
#include "METUtility.h"

using std::fabs;


#ifdef METUTIL_STANDALONE
ClassImp(METObject);
ClassImp(PhysicsObject);
ClassImp(METUtility);
#endif //METUTIL_STANDALONE

namespace AnalysisFramework
{
namespace External
{

///////////////////////////////////////////////////////////////////////////////////////////////
///functions for PhysicsObjects helper class///////////////////////////////////////////////////   
///////////////////////////////////////////////////////////////////////////////////////////////


void PhysicsObject::setMomenergy(float _pT, float _eta, float _phi, float _E, string caseFlag) {

  // TLorentzVector::SetPtEtaPhiE() applies a fabs() to the pT
  // This screws up the angle in the case of -ve pT clusters/objects
  double pi = atan(1)*4;
  if(_pT < 0) _phi += pi;
 
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets")
    m_secondaryMomenergy.SetPtEtaPhiE(_pT, _eta, _phi, _E);
  
  else m_momenergy.SetPtEtaPhiE(_pT, _eta, _phi, _E);

}//end of setMomenergy


void PhysicsObject::setEnergyUncertainty(float energyUp, float energyDown, string /*caseFlag*/) {
  //  if(caseFlag=="spectro" || caseFlag=="Spectro") {
  //  -- right now the muon group doesn't provide different uncertainties for each half of the muon
  //  m_relEnergyUncertaintySpectro.first = energyUp; 
  //  m_relEnergyUncertaintySpectro.second = energyDown;
  //  }
  m_relEnergyUncertainty.first = energyUp; 
  m_relEnergyUncertainty.second = energyDown;


}//end of setEnergyUncertainty

void PhysicsObject::setResShift(float resUp, float resDown, string caseFlag) {
  
   if(caseFlag=="spectro" || caseFlag=="Spectro") {
    m_relativeResolutionSpectro.first = resUp; 
    m_relativeResolutionSpectro.second = resDown;    
  }
  else if(caseFlag=="comboMS" || caseFlag=="comboMs" || caseFlag=="ComboMS" || caseFlag=="comboms") {
    m_relativeResolutionComboMS.first = resUp; 
    m_relativeResolutionComboMS.second = resDown;
  }

  else{
    m_relativeResolutionShift.first = resUp; 
    m_relativeResolutionShift.second = resDown;
  }
  
}//end of setResUncertainty


float PhysicsObject::E(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.E();
  }
 
  return m_momenergy.E();
}//end of pt

float PhysicsObject::Et(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.Et();
  }
  
  return m_momenergy.Et();
}//end of pt


float PhysicsObject::Pt(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.Pt();
  }
  return m_momenergy.Pt();
}//end of pt

float PhysicsObject::Px(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.Px();
  }

  return m_momenergy.Px();
}//end of px

float PhysicsObject::Py(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.Py();
  }
 
  return m_momenergy.Py();
}//end of py

float PhysicsObject::phi(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.Phi();
  }

  return m_momenergy.Phi();
}//end of phi


float PhysicsObject::eta(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro" || caseFlag == "LCJets" || caseFlag == "lcjets") {
    return m_secondaryMomenergy.Eta();
  }
 
  return m_momenergy.Eta();
}//end of phi



pair<float, float> PhysicsObject::resShift(string caseFlag) {
  if(caseFlag == "secondary" || caseFlag == "spectro" || caseFlag=="Spectro") {
    return m_relativeResolutionSpectro;
  }
   else if(caseFlag=="comboms" || caseFlag=="comboMS" || caseFlag=="ComboMS" || caseFlag=="Comboms" || caseFlag=="ComboMs") {
    return m_relativeResolutionComboMS;
  }   
  else if(caseFlag=="clusterStandard") { 
   float  _resUp = 0.0;
   float  _resDown = 0.0;
   pair<float, float> _resolution;
  _resolution.first = _resUp;
  _resolution.second = _resDown;
  return _resolution;
  }

  return m_relativeResolutionShift;
}//end of resShift


pair<float, float> PhysicsObject::energyShift(string caseFlag) {
   if(caseFlag=="clusterStandard") { //use if MET ever gets a default for cluster resolution uncertainty
    float a = (fabs(m_momenergy.Eta()) < 3.2 ) ? 0.03 : 0.1;
    float shift =  a * (1200.0/fabs(Pt()*wet()));
    pair<float, float> uncertainty;
    uncertainty.first = shift;
    uncertainty.second = -1.0*shift;
    return uncertainty;
  }
  return m_relEnergyUncertainty;
}//end of energyShift


float PhysicsObject::resolution() {
  

  return m_relativeResolution; 
}//end of resolution


void PhysicsObject::setResolution(float res) {
 
  m_relativeResolution = res;
}//end of set resolution


///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////METUtility classes////////////////////
///////////////////////////////////////////////////////////////////////////////





METUtility::METUtility(bool verbose) {


  m_verbose = verbose;
  m_softJetCut = 20000.0;
  m_pileUpUncertainty = 0.066;
  m_isMuid = false;

  m_doForwardEtaCut = false;
  m_doCellFix = false;

  //m_doCellFixOverride = true;

  m_doRefEle = true;
  m_doRefGamma = true;
  m_doRefTau = true;
  m_doRefJet = true;
  m_doRefMuon = true;
  m_doMuonTotal = true;
  m_doSoftJets = true;
  m_doCellOut = false;
  m_doCellOutEflow = true;

  m_useStandardClusterRes = false;
  m_useStandardClusterEnergySigma = false;


  m_doSignificance = false;

  reset();
} //end of

METUtility::METUtility(bool doRefEle, bool doRefGamma, bool doRefTau, bool doRefJet, bool doSoftJets, bool doRefMuon, bool doMuonTotal, bool doCellOut, bool doCellOutEflow, bool isMuid, double softJetCut, bool verbose) {

  m_verbose = verbose;
  m_pileUpUncertainty = 0.066;
  
  m_doForwardEtaCut = false;
  m_doCellFix = false;
  //m_doCellFixOverride = true;


  m_doRefEle = doRefEle;
  m_doRefGamma = doRefGamma;
  m_doRefTau = doRefTau;
  m_doRefJet = doRefJet;
  m_doSoftJets = doSoftJets;
  m_doRefMuon = doRefMuon;
  m_doMuonTotal = doMuonTotal;
  m_doCellOut = doCellOut;
  m_doCellOutEflow = doCellOutEflow;
  
  m_useStandardClusterRes = false;
  m_useStandardClusterEnergySigma = false;
  
  m_softJetCut = softJetCut;
  
  if(m_doSoftJets == false) m_softJetCut = 0.0; // so smeared jets don't get lose from the RefJet term
  
  m_isMuid = isMuid;
  
  m_doSignificance = false;

  reset();
}//end of creation function with arguments

void METUtility::reset() {
  
   m_jets.clear();
   m_electrons.clear();
   m_muons.clear();
   m_taus.clear();
   m_clusters.clear();
   m_tracks.clear();
   m_photons.clear();

   m_cellOut.setBase(0.,0.,0.);
   m_cellOutEflow.setBase(0.,0.,0.);
   m_softJets.setBase(0.,0.,0.);
   m_refMuon.setBase(0.,0.,0.);
   m_refJet.setBase(0.,0.,0.);
   m_refEle.setBase(0.,0.,0.);
   m_refTau.setBase(0.,0.,0.);
   m_refGamma.setBase(0.,0.,0.);
   m_muonTotal.setBase(0.,0.,0.);

   m_jetsMomentaSet = false;
   m_photonsMomentaSet = false;
   m_electronsMomentaSet = false;
   m_muonsMomentaSet = false;
   m_tausMomentaSet = false;
   m_clustersMomentaSet = false;
   m_trackMomentaSet = false;


   m_jetsUncertaintiesSet = false;
   m_photonsUncertaintiesSet = false;
   m_electronsUncertaintiesSet = false;
   m_muonsUncertaintiesSet = false;
   m_tausUncertaintiesSet = false;
   m_clustersUncertaintiesSet = false;
   m_trackUncertaintiesSet = false;
   m_muonsComboIDUncertaintiesSet = false;
   m_muonsComboMSUncertaintiesSet = false;
 
   m_jetsResShiftsSet = false;
   m_photonsResShiftsSet = false;
   m_electronsResShiftsSet = false;
   m_muonsResShiftsSet = false;
   m_tausResShiftsSet = false;
   m_clustersResShiftsSet = false;
   m_tracksResShiftsSet = false;

   m_jetsResolutionsSet = false;
   m_photonsResolutionsSet = false;
   m_electronsResolutionsSet = false;
   m_muonsResolutionsSet = false;
   m_tausResolutionsSet = false;
   m_clustersResolutionsSet = false;
   m_tracksResolutionsSet = false;
}//end of 

METObject METUtility::deltaMissingET(string term, string systematics) {

  //use commas as delimiters

  if(systematics=="All"  || systematics=="ALL"  || systematics=="all") {
    systematics = "";
    if(m_doRefEle) systematics.append("EES,EER,");
    if(m_doRefGamma) systematics.append("PES,PER,");    
    if(m_doRefTau) systematics.append("TES,TER,");
    if(m_doRefJet) systematics.append("JES,JER");
    if(m_doMuonTotal) systematics.append("MES,MERID,MERMS");
    if(m_doSoftJets) systematics.append("AllClusters");
    if(m_doCellOut) systematics.append("AllClusters");
    if(m_doCellOutEflow) systematics.append("AllClusters");
  }//end of all check

  for(unsigned int i=0; i < systematics.length(); ++i)
    if(systematics[i] == ' ') systematics.erase(i,1);

  vector<string> sources;
  string::size_type lastPos = systematics.find_first_not_of(",", 0);
  string::size_type pos = systematics.find_first_of(",", lastPos);
  while (string::npos != pos || string::npos != lastPos) {
    sources.push_back(systematics.substr(lastPos, pos - lastPos));
    lastPos = systematics.find_first_not_of(",", pos);
    pos = systematics.find_first_of(",",lastPos);
  }//end of loop

  Float_t summedSquares_et = 0.0;
  Float_t summedSquares_etx = 0.0;
  Float_t summedSquares_ety = 0.0;
  Float_t summedSquares_sumet = 0.0;
  METObject metTerm;
  string systo = "";
  for(unsigned int j = 0; j < sources.size(); ++j) {
    systo = sources.at(j);
    systo.append("Diff");
    metTerm = getMissingET(term, systo);
    summedSquares_et += metTerm.et()*metTerm.et();
    summedSquares_etx += metTerm.etx()*metTerm.etx();
    summedSquares_ety += metTerm.ety()*metTerm.ety();
    summedSquares_sumet += metTerm.sumet()*metTerm.sumet();
  }
  summedSquares_et = sqrt(summedSquares_et);
  summedSquares_etx = sqrt(summedSquares_etx);
  summedSquares_ety = sqrt(summedSquares_ety);
  summedSquares_sumet = sqrt(summedSquares_sumet);
  
  METObject deltaMET(summedSquares_etx, summedSquares_ety, summedSquares_sumet);
  deltaMET.setEt(summedSquares_et); //substitutes our et diff for what would be summed in quadrature from etx/ety
  return deltaMET;
  }//end of deltaMissingET

METObject METUtility::absDeltaMissingET(string term, string systematics) {
  METObject relMET = deltaMissingET(term, systematics);
  METObject theMET = getMissingET(term);

  Float_t _et = relMET.et()*theMET.et();
  Float_t _etx = relMET.etx()*theMET.etx();
  Float_t _ety = relMET.ety()*theMET.ety();
  Float_t _sumet = relMET.sumet()*theMET.sumet();
  
  METObject absDeltaMET(_etx, _ety, _sumet);
  absDeltaMET.setEt(_et);
  return absDeltaMET;

}//end of absDeltaMissingET


METObject METUtility::getMissingET(string term, string systematic) {
  

 
  //  m_doCellFix = false;
  // if(m_doCellFix && (systematic.find("Up") != string::npos || systematic.find("up") != string::npos || systematic.find("down") != string::npos || systematic.find("Down") != string::npos)) m_doCellFix = true;


  size_t pos = systematic.find("diff");
  if(pos == string::npos) pos = systematic.find("Diff");
  
  if(pos != string::npos) {
    string systo = systematic.substr(0, pos);
    

    float _etx = MissingETHelper(term).etx();
    if(_etx != 0) _etx = (fabs(_etx - MissingETHelper(term, systo+"Up").etx()) + fabs(_etx - MissingETHelper(term, systo+"Down").etx()))/(2*_etx);
    float _ety = MissingETHelper(term).ety();
    if(_ety != 0) _ety = (fabs(_ety - MissingETHelper(term, systo+"Up").ety()) + fabs(_ety - MissingETHelper(term, systo+"Down").ety()))/(2*_ety);
    float _sumet = MissingETHelper(term).sumet();
    if(_sumet != 0) _sumet = (fabs(_sumet - MissingETHelper(term, systo+"Up").sumet()) + fabs(_sumet - MissingETHelper(term, systo+"Down").sumet()))/(2*_sumet);
    float _et = MissingETHelper(term).et();
    if(_et != 0) _et = (fabs(_et - MissingETHelper(term, systo+"Up").et()) + fabs(_et - MissingETHelper(term, systo+"Down").et()))/(2*_et);
    

    METObject metDiff(_etx, _ety, _sumet);
    metDiff.setEt(_et); //substitutes our et diff for what would be summed in quadrature from etx/ety
    return metDiff;
    
  }//end of diff logic
  else{
    return MissingETHelper(term, systematic);
  }//end of else

  
}//end of MissingET





//most of the excitment happens hre
METObject METUtility::MissingETHelper(string term, string systematic) {

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  if((term == "RefFinal" || term == "RefEle") && m_doRefEle) {
    m_scaled_refEle = RefEle(systematic);
    if(term == "RefEle") return m_scaled_refEle;
    _etx += m_scaled_refEle.etx();
    _ety += m_scaled_refEle.ety();
    _sumet += m_scaled_refEle.sumet();
  }
  if((term == "RefFinal" || term == "RefGamma") && m_doRefGamma) {
    m_scaled_refGamma = RefGamma(systematic);
    if(term == "RefGamma") return m_scaled_refGamma;
    _etx += m_scaled_refGamma.etx();
    _ety += m_scaled_refGamma.ety();
    _sumet += m_scaled_refGamma.sumet();
  }
  if((term == "RefFinal" || term == "RefTau") && m_doRefTau) {
    m_scaled_refTau = RefTau(systematic);
    if(term == "RefTau") return m_scaled_refTau;
    _etx += m_scaled_refTau.etx();
    _ety += m_scaled_refTau.ety();
    _sumet += m_scaled_refTau.sumet();
  }
  if((term == "RefFinal" || term == "RefJet") && m_doRefJet) {
    m_scaled_refJet = RefJet(systematic);
     if(term == "RefJet") return m_scaled_refJet;    
    _etx += m_scaled_refJet.etx();
    _ety += m_scaled_refJet.ety();
    _sumet += m_scaled_refJet.sumet();
  }
  if((term == "RefFinal" || term == "RefMuon") && m_doRefMuon) {
    m_scaled_refMuon = RefMuon(systematic);
    if(term == "RefMuon") return m_scaled_refMuon;
    _etx += m_scaled_refMuon.etx(); 
    _ety += m_scaled_refMuon.ety();
    _sumet += m_scaled_refMuon.sumet();
  }
  if((term == "RefFinal" || term == "MuonTotal") && m_doMuonTotal) {
    m_scaled_muonTotal = MuonTotal(systematic);
    if(term == "MuonTotal") return m_scaled_muonTotal;
    _etx += m_scaled_muonTotal.etx();
    _ety += m_scaled_muonTotal.ety();
    //muons do not contribute to sumet except for refmuon, it is a calorimeter term
  }
  if((term == "RefFinal" || term == "CellOut") && m_doCellOut) {
    m_scaled_cellOut = CellOut(systematic);
    if(term == "CellOut") return m_scaled_cellOut;
    _etx += m_scaled_cellOut.etx();
    _ety += m_scaled_cellOut.ety();
    _sumet += m_scaled_cellOut.sumet();
  }
  if((term == "RefFinal" || term == "CellOutEflow") && m_doCellOutEflow) {
    m_scaled_cellOutEflow = CellOutEflow(systematic);
    if(term == "CellOutEflow") return m_scaled_cellOutEflow;
    _etx += m_scaled_cellOutEflow.etx(); 
    _ety += m_scaled_cellOutEflow.ety();
    _sumet += m_scaled_cellOutEflow.sumet();
  } 
  if((term == "RefFinal" || term == "SoftJets") && m_doSoftJets) {
    m_scaled_softJets = SoftJets(systematic);
    if(term == "SoftJets") return m_scaled_softJets;
    _etx += m_scaled_softJets.etx();
    _ety += m_scaled_softJets.ety();
    _sumet += m_scaled_softJets.sumet();
  }
  
  METObject refFinal(_etx, _ety, _sumet);
  //refFinal.setEt(-1.0);
  float siggy = 0.0;
  if(m_doSignificance) {
    siggy = METSignificance(systematic);
  }//end of do sig
  
  refFinal.setSignificance(siggy);
 
 return refFinal;

}//end of MissingET main function


void METUtility::setObjectsHelper(string type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, vector<vector<float> > wet, vector<vector<float> > wex, vector<vector<float> > wey, vector<vector<unsigned int> > statusWord) {
 vector<AnalysisFramework::External::PhysicsObject> _objectVector;

  //why do it like this? instead of giving an object multiple weights, I'll just do a different PhysicsObject and have multiple muons, 
  //such that I have combined muons, spectro muons, track muons, and refmuon muons and all have different weights, and the statusWord 
  //tells how to use them. It's sort of how it works when METPerformance uses the composition map at AOD level... 
  //Of course the filler for uncertainties and resolutions has to be adjusted to match...
  for(unsigned int i = 0; i < wet.size(); ++i) {
    for(unsigned int j = 0; j < wet.at(i).size(); ++j) {
      PhysicsObject _object;// = new PhysicsObject;
      _object.setMomenergy(pT.at(i), eta.at(i), phi.at(i), E.at(i));
      float _wex = wex.at(i).at(j);
      float _wey = wey.at(i).at(j);
      float _wet = wet.at(i).at(j);
      unsigned int _statusWord = statusWord.at(i).at(j);
      _object.setWeights(_wex, _wey, _wet);
      _object.setStatusCode(_statusWord);
      _object.setIndex(i);
      _objectVector.push_back(_object);
    }//end of inner weight loop
  }//end of outer weight loop
  
  if(type == "Electrons" || type == "electrons") {
    m_electrons.clear();
    m_electrons = _objectVector;
    m_electronsMomentaSet = true;
  }
  else if(type == "Photons" || type == "photons") {
    m_photons.clear();
    m_photons = _objectVector;
    m_photonsMomentaSet = true;
  }
  else if(type == "Taus" || type == "taus") {
    m_taus.clear();
    m_taus = _objectVector;
    m_tausMomentaSet = true;
  }
  else if(type == "Jets" || type == "jets") {
    m_jets.clear();
    m_jets = _objectVector;
    m_jetsMomentaSet = true;
  }
  else if(type == "Muons" || type == "muons") {
    m_muons.clear();
    m_muons = _objectVector;
    m_muonsMomentaSet = true;
  }
  else if(type == "Clusters" || type == "clusters") {
    m_clusters.clear();
    m_clusters = _objectVector;
    m_clustersMomentaSet = true;
  }
  else if(type == "Tracks" || type == "tracks") {
    m_tracks.clear();
    m_tracks = _objectVector;
    m_trackMomentaSet = true;
  } 

}//end of setobjecthelper



void METUtility::setObjectsHelper(string type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, vector<float> wet, vector<float> wex, vector<float> wey, vector<unsigned short> statusWord) {


  vector<vector<float> > newWex;// = new vector<vector<float> >;
  vector<vector<float> > newWey;// = new vector<vector<float> >;
  vector<vector<float> > newWet;// = new vector<vector<float> >;
  vector<vector<unsigned int> > newStatusWord;// = new vector<vector<unsigned int> >;


  for(unsigned int i = 0; i < statusWord.size(); ++i) {

    vector<float> wexHold;
    vector<float> weyHold;
    vector<float> wetHold;
    vector<unsigned int> statusWordHold;

    wexHold.push_back(wex.at(i));
    weyHold.push_back(wey.at(i));
    wetHold.push_back(wet.at(i));
    statusWordHold.push_back(static_cast<unsigned int>(statusWord.at(i)));

    newWex.push_back(wexHold);
    newWey.push_back(weyHold);
    newWet.push_back(wetHold);
    newStatusWord.push_back(statusWordHold);

   
  }//end of loop


  setObjectsHelper(type, (pT), (eta), (phi), (E), newWet, newWex, newWey, newStatusWord);


}//end of setObjectshelper 




void METUtility::setObjects(string type, const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {

  setObjectsHelper(type, (*pT), (*eta), (*phi), (*E), (*wet), (*wex), (*wey), (*statusWord));
}//adapter to use vector<float> instead of vector<vector<float> >




void METUtility::setObjects(string type, const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  
  setObjectsHelper(type, (*pT), (*eta), (*phi), (*E), (*wet), (*wex), (*wey), (*statusWord));

}//end of setObjects



void METUtility::setObjectMomenta(string type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E) {

  if(type == "Electrons" || type == "electrons") { 
  setObjectMomentaHelper(m_electrons, pT, eta, phi, E);
  }
  else if(type == "Photons" || type == "photon") { 
    setObjectMomentaHelper(m_photons, pT, eta, phi, E);
  }
  else if(type == "Taus" || type == "taus") { 
    setObjectMomentaHelper(m_taus, pT, eta, phi, E);
  }
  else if(type == "Jets" || type == "jets") { 
    setObjectMomentaHelper(m_jets, pT, eta, phi, E);
  }
  else if(type == "Muons" || type == "muons") { 
    setObjectMomentaHelper(m_muons, pT, eta, phi, E);
  }
  else if(type == "Spectromuons" || type == "spectroMuons" || type == "SpectroMuons" || type == "spectromuons") { 
    setObjectMomentaHelper(m_muons, pT, eta, phi, E, "secondary");
  }
  else if(type == "LCJets" || type == "lcjets") { ///this won't work
    setObjectMomentaHelper(m_jets, pT, eta, phi, E, "secondary");
  }
  else if(type == "Clusters" || type == "clusters") { 
    setObjectMomentaHelper(m_clusters, pT, eta, phi, E);
  }
  else if(type == "Tracks" || type == "tracks") { 
    setObjectMomentaHelper(m_tracks, pT, eta, phi, E);
  }
  
}//end of setMomenta

void METUtility::setObjectMomentaHelper(vector<AnalysisFramework::External::PhysicsObject> &object, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, string caseFlag) {

  for(unsigned int i = 0; i < object.size(); ++i) {
     object.at(i).setMomenergy(pT.at(object.at(i).index()), eta.at(object.at(i).index()), phi.at(object.at(i).index()), E.at(object.at(i).index()), caseFlag);
  }//end of object loops


}//end of setMomentaHelper




void METUtility::setObjectEnergyUncertainties(string type, const vector<float> energyUp, const vector<float> energyDown) {

  if(type == "Electrons" || type == "electrons") { 
    setObjectEnergyUncertaintiesHelper(m_electrons,energyUp, energyDown);
    m_electronsUncertaintiesSet = true;
  }
  else if(type == "Photons" || type == "photon") { 
    setObjectEnergyUncertaintiesHelper(m_photons, energyUp, energyDown);
    m_photonsUncertaintiesSet = true;
  }
  else if(type == "Taus" || type == "taus") { 
    setObjectEnergyUncertaintiesHelper(m_taus, energyUp, energyDown);
    m_tausUncertaintiesSet = true;
  }
  else if(type == "Jets" || type == "jets") { 
    setObjectEnergyUncertaintiesHelper(m_jets, energyUp, energyDown);
    m_jetsUncertaintiesSet = true;
  }
  else if(type == "Muons" || type == "muons") { 
    setObjectEnergyUncertaintiesHelper(m_muons, energyUp, energyDown);
    m_muonsUncertaintiesSet = true;
  }
  else if(type == "Spectromuons" || type == "spectroMuons" || type == "SpectroMuons" || type == "spectromuons") {
    setObjectEnergyUncertaintiesHelper(m_muons, energyUp, energyDown, "spectro");
    m_muonsSpectroUncertaintiesSet = true;
  }
  else if(type=="comboms" || type=="comboMS" || type=="ComboMS" || type=="Comboms" || type=="ComboMs") {
   setObjectEnergyUncertaintiesHelper(m_muons, energyUp, energyDown, "comboMS");
   m_muonsComboMSUncertaintiesSet = true;
 }
 else if(type=="comboid" || type=="comboID" || type=="ComboID" || type=="Comboid" || type=="ComboId") {
   setObjectEnergyUncertaintiesHelper(m_muons, energyUp, energyDown, "comboID");
   m_muonsComboIDUncertaintiesSet = true;
 }
  else if(type == "Clusters" || type == "clusters") { 
    setObjectEnergyUncertaintiesHelper(m_clusters, energyUp, energyDown);
    m_clustersUncertaintiesSet = true;
  }
  else if(type == "Tracks" || type == "tracks") { 
    setObjectEnergyUncertaintiesHelper(m_tracks, energyUp, energyDown);
    m_trackUncertaintiesSet = true;
  }
  
}//end of setObjectUncertainties


void METUtility::setObjectEnergyUncertaintiesHelper(vector<AnalysisFramework::External::PhysicsObject> &object, vector<float> energyUp, vector<float> energyDown, string caseFlag) {
  if(object.size() == 0 && m_verbose) cout << "Set object four-momentum and weights first!" << endl;
   //don't assume number of objects matches number of uncertainties (muons and clusters won't), use index to match
  for(unsigned int i = 0; i < object.size(); ++i) {
    object.at(i).setEnergyUncertainty(energyUp.at(object.at(i).index()), energyDown.at(object.at(i).index()), caseFlag);
  }//end of object loops

}//end of setObjectUncertaintiesHelper
 


void METUtility::setObjectResolutionShift(string type, const vector<float> resUp, const vector<float> resDown) {

  if(type == "Electrons" || type == "electrons") { 
    setObjectResolutionShiftHelper(m_electrons,resUp, resDown);
    m_electronsResShiftsSet = true;
  }
  else if(type == "Photons" || type == "photon") { 
    setObjectResolutionShiftHelper(m_photons, resUp, resDown);
    m_photonsResShiftsSet = true;
  }
  else if(type == "Taus" || type == "taus") { 
    setObjectResolutionShiftHelper(m_taus, resUp, resDown);
    m_tausResShiftsSet = true;
  }
  else if(type == "Jets" || type == "jets") { 
    setObjectResolutionShiftHelper(m_jets, resUp, resDown);
    m_jetsResShiftsSet = true;
  }
  else if(type == "Muons" || type == "muons") { 
    setObjectResolutionShiftHelper(m_muons, resUp, resDown);
    m_muonsResShiftsSet = true;
  }
  else if(type == "Spectromuons" || type == "spectroMuons" || type == "SpectroMuons" || type == "spectromuons") {
    setObjectResolutionShiftHelper(m_muons, resUp, resDown, "spectro");
    m_muonsSpectroResShiftsSet = true;
  }
  else if(type=="comboms" || type=="comboMS" || type=="ComboMS" || type=="Comboms" || type=="ComboMs") {
   setObjectResolutionShiftHelper(m_muons, resUp, resDown, "comboMS");
   m_muonsComboMSResShiftsSet = true;
  }
  else if(type=="comboid" || type=="comboID" || type=="ComboID" || type=="Comboid" || type=="ComboId") {
    setObjectResolutionShiftHelper(m_muons, resUp, resDown, "comboID");
    m_muonsComboIDResShiftsSet = true;
  }
  else if(type == "Clusters" || type == "clusters") { 
    setObjectResolutionShiftHelper(m_clusters, resUp, resDown);
    m_clustersResShiftsSet = true;
  }
  else if(type == "Tracks" || type == "tracks") { 
    setObjectResolutionShiftHelper(m_tracks, resUp, resDown);
    m_tracksResShiftsSet = true;
  }
  
}//end of setObjectResolutions


void METUtility::setObjectResolutionShiftHelper(vector<AnalysisFramework::External::PhysicsObject> &object, vector<float> resUp, vector<float> resDown, string caseFlag) {
  if(object.size() == 0 && m_verbose) cout << "Set object four-momentum and weights first!" << endl;
  
  //don't assume number of objects matches number of resolutions (muons and clusters won't), use index to match
  for(unsigned int i = 0; i < object.size(); ++i) {
    object.at(i).setResShift(resUp.at(object.at(i).index()), resDown.at(object.at(i).index()), caseFlag);
  }//end of object loops

}//end of setObjectUncertaintiesHelper
 




void METUtility::setObjectResolutions(string type, const vector<float> resolutions) {
  
  if(type == "Electrons" || type == "electrons") { 
    setObjectResolutionsHelper(m_electrons, resolutions);
    m_electronsResolutionsSet = true;
  }
  else if(type == "Photons" || type == "photon") { 
    setObjectResolutionsHelper(m_photons,resolutions);
    m_photonsResolutionsSet = true;
  }
  else if(type == "Taus" || type == "taus") { 
    setObjectResolutionsHelper(m_taus, resolutions);
    m_tausResolutionsSet = true;
  }
  else if(type == "Jets" || type == "jets") { 
    setObjectResolutionsHelper(m_jets, resolutions);
    m_jetsResolutionsSet = true;
  }
  else if(type == "Muons" || type == "muons") { 
    setObjectResolutionsHelper(m_muons, resolutions);
    m_muonsResolutionsSet = true;
  }
  else if(type == "Spectromuons" || type == "spectroMuons" || type == "SpectroMuons" || type == "spectromuons") { 
    setObjectResolutionsHelper(m_muons, resolutions, "spectro");
    m_muonsSpectroResolutionsSet = true;
  }
   //  else if(type == "ComboMS" || type == "Comboms" || type == "comboMS" || type == "comboms") { 
//     setObjectResolutionsHelper(m_muons, resolutions, "comboMS");
//     m_muonsSpectroResolutionsSet = true;
//   }
//   else if(type == "ComboID" || type == "Comboid" || type == "comboID" || type == "comboid") { 
//     setObjectResolutionsHelper(m_muons, resolutions, "comboID");
//     m_muonsComboIDResolutionsSet = true;
//   }
  else if(type == "Clusters" || type == "clusters") { 
    setObjectResolutionsHelper(m_clusters, resolutions);
    m_muonsComboMSResolutionsSet = true;
  }
  else if(type == "Tracks" || type == "tracks") { 
    setObjectResolutionsHelper(m_tracks, resolutions);
    m_tracksResolutionsSet = true;
  }
  
}//end of setObjectResolutions


void METUtility::setObjectResolutionsHelper(vector<AnalysisFramework::External::PhysicsObject> &objects, vector<float> resolutions, string caseFlag) {
  caseFlag = "whoopty";
  if(objects.size() == 0 && m_verbose) cout << "Set object four-momentum and weights first!" << endl;
  for(unsigned int i = 0; i < objects.size(); ++i) {
    objects.at(i).setResolution(resolutions.at(objects.at(i).index()));
  }//end of object loops
  
}//end of setObjectResolutionsHelper


void METUtility::setMETTerm(string term, float _etx, float _ety, float _sumet) {

  if(term == "CellOut") {
    m_cellOut.setEtx(_etx);
    m_cellOut.setEty(_ety);
    m_cellOut.setSumet(_sumet);
    m_cellOut.setIsValid(true);
  }  
  else if(term == "CellOutEflow") {
    m_cellOutEflow.setEtx(_etx);
    m_cellOutEflow.setEty(_ety);
    m_cellOutEflow.setSumet(_sumet);
    m_cellOutEflow.setIsValid(true);
  }  
  else if(term == "SoftJets") {
    m_softJets.setEtx(_etx);
    m_softJets.setEty(_ety);
    m_softJets.setSumet(_sumet);
    m_softJets.setIsValid(true);
  }  
  else if(term == "RefMuon") {
    m_refMuon.setEtx(_etx);
    m_refMuon.setEty(_ety);
    m_refMuon.setSumet(_sumet);
    m_refMuon.setIsValid(true);
  }
  else if(term == "RefJet") {
    m_refJet.setEtx(_etx);
    m_refJet.setEty(_ety);
    m_refJet.setSumet(_sumet);
    m_refJet.setIsValid(true);
  }
  else if(term == "RefEle") {
    m_refEle.setEtx(_etx);
    m_refEle.setEty(_ety);
    m_refEle.setSumet(_sumet);
    m_refEle.setIsValid(true);
  }
  else if(term == "RefGamma") {
    m_refGamma.setEtx(_etx);
    m_refGamma.setEty(_ety);
    m_refGamma.setSumet(_sumet);
    m_refGamma.setIsValid(true);
  }
  else if(term == "RefTau") {
    m_refTau.setEtx(_etx);
    m_refTau.setEty(_ety);
    m_refTau.setSumet(_sumet);
    m_refTau.setIsValid(true);
  }
  else if(term == "MuonTotal") {
    m_muonTotal.setEtx(_etx);
    m_muonTotal.setEty(_ety);
    m_muonTotal.setSumet(_sumet);
    m_muonTotal.setIsValid(true);
  }
}


void METUtility::defineMissingET(bool doRefEle, bool doRefGamma, bool doRefTau, bool doRefJet, bool doSoftJets, bool doRefMuon, bool doMuonTotal, bool doCellOut, bool doCellOutEflow) {

  m_doRefEle = doRefEle;
  m_doRefGamma = doRefGamma;
  m_doRefTau = doRefTau;
  m_doRefJet = doRefJet;
  m_doSoftJets = doSoftJets;
  m_doRefMuon = doRefMuon;
  m_doMuonTotal = doMuonTotal;
  m_doCellOut = doCellOut;
  m_doCellOutEflow = doCellOutEflow;

  if(m_doSoftJets == false) m_softJetCut = 0.0; // so smeared jets don't get lost from the RefJet term

}


METObject METUtility::RefEle(string systematic) {

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_electrons.size(); ++i) {
    scale = 1.0;
   
    if(m_doCellFix && m_electrons.at(i).wet() <= .5) continue;
    if(m_doForwardEtaCut && fabs(m_electrons.at(i).eta()) > 4.5) continue;
 
    if(systematic == "None" || systematic == "none") scale = 1.0;
    else if(systematic == "EERUp" && m_electronsResShiftsSet) scale = 1.0 + m_electrons.at(i).resShift().first;
    else if(systematic == "EERDown" && m_electronsResShiftsSet) scale = 1.0 + m_electrons.at(i).resShift().second;
    else if(systematic == "EESUp" && m_electronsUncertaintiesSet) scale = 1.0 + m_electrons.at(i).energyShift().first;
    else if(systematic == "EESDown" && m_electronsUncertaintiesSet) scale = 1.0 + m_electrons.at(i).energyShift().second;
    
  
    _etx -= m_electrons.at(i).Px()*m_electrons.at(i).wex()*scale;
    _ety -= m_electrons.at(i).Py()*m_electrons.at(i).wey()*scale;
    _sumet += m_electrons.at(i).Pt()*m_electrons.at(i).wet()*scale;
   
  }//end of loop
 
  if(m_electronsMomentaSet == false && m_refEle.isValid()) {
    _etx += m_refEle.etx();
    _ety += m_refEle.ety();
    _sumet += m_refEle.sumet();
  }

  METObject _refEle(_etx, _ety, _sumet);
  return _refEle;
}//end of RefEle term remaker



METObject METUtility::RefGamma(string systematic) {

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_photons.size(); ++i) {
    scale = 1.0;
  
    if(m_doCellFix && m_photons.at(i).wet() <= .5) continue;
    if(m_doForwardEtaCut && fabs(m_photons.at(i).eta()) > 4.5) continue;
 
    if(systematic == "None" || systematic == "none") scale = 1.0;
    else if(systematic == "PERUp" && m_photonsResShiftsSet) scale = 1.0 + m_photons.at(i).resShift().first;
    else if(systematic == "PERDown" && m_photonsResShiftsSet) scale = 1.0 + m_photons.at(i).resShift().second;
    else if(systematic == "PESUp" && m_photonsUncertaintiesSet) scale = 1.0 + m_photons.at(i).energyShift().first;
    else if(systematic == "PESDown" && m_photonsUncertaintiesSet) scale = 1.0 + m_photons.at(i).energyShift().second;
    
 
    _etx -= m_photons.at(i).Px()*m_photons.at(i).wex()*scale;
    _ety -= m_photons.at(i).Py()*m_photons.at(i).wey()*scale;
    _sumet += m_photons.at(i).Pt()*m_photons.at(i).wet()*scale;

  }//end of loop

  if(m_photonsMomentaSet == false  && m_refGamma.isValid()) {
    _etx += m_refGamma.etx();
    _ety += m_refGamma.ety();
    _sumet += m_refGamma.sumet();
  }

  METObject _refGamma(_etx, _ety, _sumet);
  return _refGamma;
}//end of RefGamma

METObject METUtility::RefTau(string systematic) {

  
  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_taus.size(); ++i) {
    scale = 1.0;

    if(m_doCellFix && m_taus.at(i).wet() <= 0.5) continue;
    if(m_doForwardEtaCut && fabs(m_taus.at(i).eta()) > 4.5) continue;
 
    if(systematic == "None" || systematic == "none") scale = 1.0;
    else if(systematic == "TERUp" && m_tausResShiftsSet) scale = 1.0 + m_taus.at(i).resShift().first;
    else if(systematic == "TERDown" && m_tausResShiftsSet) scale = 1.0 + m_taus.at(i).resShift().second;
    else if(systematic == "TESUp" && m_tausUncertaintiesSet) scale = 1.0 + m_taus.at(i).energyShift().first;
    else if(systematic == "TESDown" && m_tausUncertaintiesSet) scale = 1.0 + m_taus.at(i).energyShift().second;

    _etx -= m_taus.at(i).Px()*m_taus.at(i).wex()*scale;
    _ety -= m_taus.at(i).Py()*m_taus.at(i).wey()*scale;
    _sumet += m_taus.at(i).Pt()*m_taus.at(i).wet()*scale;

  }//end of loop

  if(m_tausMomentaSet == false  && m_refTau.isValid()) {
    _etx += m_refTau.etx();
    _ety += m_refTau.ety();
    _sumet += m_refTau.sumet();
  }

   const int nBins = 34;
   double refTauFracUncert[nBins] = {0.136409, 0.153857, 0.153689, 0.149694, 0.138453, 0.139799, 0.133924, 0.139531, 0.138844, 0.134789, 0.138834, 0.132662, 0.142114, 0.13908, 0.13661, 0.128625, 0.133791, 0.125408, 0.122293, 0.133035, 0.132411, 0.131334, 0.128113, 0.13843, 0.12271, 0.129, 0.125528, 0.125005, 0.130233, 0.125374, 0.120071, 0.125658, 0.121132, 0.121722}; //CellOut Sum Et Bins, adapted to reftau at same scaling
   double refTauSumEtBounds[nBins] = {-5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160};//Cell Out Fraction uncertainty for bin ranges. adapted to reftau

  double slope = -1.19302e-4;//error 1.9e-05
  double intercept = 0.141673; //error 0.0022
  
  
  //strictly speaking this is for CellOutEflow, but we added it so analyses without taus have something
  int bin = -1;
  
  if(systematic == "RefTauUp" || systematic == "RefTauDown" ) {
    for(int i = 0; i < nBins-1; ++i) {
      if(fabs(_sumet)/1000.0 >= refTauSumEtBounds[i] && fabs(_sumet)/1000.0 < refTauSumEtBounds[i+1]) bin = i;
    }//end of loop
    float shift = refTauFracUncert[bin];
    if(fabs(_sumet)/1000.0 >= refTauSumEtBounds[nBins-1]) shift = slope*fabs(_sumet) + intercept;
    
    if(systematic == "RefTauUp" ) scale = 1.0 + shift;
    else if(systematic == "RefTauDown") scale = 1.0 - shift;
    _etx *= scale;
    _ety *= scale;
    _sumet *= scale;

  }//end of RefTauSystematic


//   if(systematic == "PileUpUp" || systematic == "PileUpDown") {
//     if(systematic == "PileUpUp") scale = 1.0 + m_pileUpUncertainty;
//     else if(systematic == "PileUpDown") scale = 1.0 - m_pileUpUncertainty;
//     _etx *= scale;
//     _ety *= scale;
//     _sumet *= scale;
//   }//end of pileup systematic
  

  METObject _refTau(_etx, _ety, _sumet);
return _refTau;
}//end of RefTau

METObject METUtility::RefJet(string systematic) {

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_jets.size(); ++i) {
    if(m_jets.at(i).Pt() > m_softJetCut) {
      scale = 1.0;

      if(m_doCellFix && m_jets.at(i).wet() <= .5) {continue;}
      if(m_doForwardEtaCut && fabs(m_jets.at(i).eta()) > 4.5) continue;
 
      if(systematic == "None" || systematic == "none") scale = 1.0;
      else if(systematic == "JERUp" && m_jetsResShiftsSet) scale = 1.0 + m_jets.at(i).resShift().first;
      else if(systematic == "JERDown" && m_jetsResShiftsSet) scale = 1.0 + m_jets.at(i).resShift().second;
      else if(systematic == "JESUp" && m_jetsUncertaintiesSet) scale = 1.0 + m_jets.at(i).energyShift().first;
      else if(systematic == "JESDown" && m_jetsUncertaintiesSet) scale = 1.0 + m_jets.at(i).energyShift().second;
       _etx -= m_jets.at(i).Px()*m_jets.at(i).wex()*scale;
       _ety -= m_jets.at(i).Py()*m_jets.at(i).wey()*scale;
       _sumet += m_jets.at(i).Pt()*m_jets.at(i).wet()*scale;
    }//end of softjet cut
  }//end of loop


  if(m_jetsMomentaSet == false  && m_refJet.isValid()) {
    _etx += m_refJet.etx();
    _ety += m_refJet.ety();
    _sumet += m_refJet.sumet();
  }

//   if(systematic == "PileUpUp" || systematic == "PileUpDown") {
//   if(systematic == "PileUpUp") scale = 1.0 + m_pileUpUncertainty;
//   else if(systematic == "PileUpDown") scale = 1.0 - m_pileUpUncertainty;
   
//   _etx *= scale;
//   _ety *= scale;
//   _sumet *= scale;
//   }//end of pileup systematic

   METObject _refJet(_etx, _ety, _sumet);
  return _refJet;
}//end of jet term remaker



METObject METUtility::RefMuon(string systematic) {
  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;
  
  if(m_refMuon.isValid()) {
    _etx += m_refMuon.etx();
    _ety += m_refMuon.ety();
    _sumet += m_refMuon.sumet();
  } else {
    float scale = 1.0;
    for(unsigned int i = 0; i < m_muons.size(); ++i) {
      scale = 1.0;

      if(systematic == "None" || systematic == "none") scale = 1.0;
      if(m_doForwardEtaCut && fabs(m_muons.at(i).eta()) > 4.5) continue;
 
      if(MissingETTags::usesREFMUON(m_muons.at(i).statusWord())) {
	_etx -= m_muons.at(i).Px()*m_muons.at(i).wex()*scale;
	_ety -= m_muons.at(i).Py()*m_muons.at(i).wey()*scale;
	_sumet += m_muons.at(i).Pt()*m_muons.at(i).wet()*scale;
      }//end of statusWord check

    }//end of loop
  }

  METObject _refMuon(_etx, _ety, _sumet);
  return _refMuon;
}//end of refmuon


METObject METUtility::MuonTotal(string systematic) {

 
  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  float scale = 1.0;
  
  for(unsigned int i = 0; i < m_muons.size(); ++i) {
    scale = 1.0;

    if(MissingETTags::usesREFMUON(m_muons.at(i).statusWord())) continue;
    if(!m_isMuid && MissingETTags::usesMUID(m_muons.at(i).statusWord())) continue;
    if(m_isMuid && !(MissingETTags::usesMUID(m_muons.at(i).statusWord()))) continue;

    if(MissingETTags::usesDEFAULT(m_muons.at(i).statusWord()) ||
       MissingETTags::usesTRACK(m_muons.at(i).statusWord())) {

      if(systematic == "None" || systematic == "none")
	scale = 1.0;
      else if(systematic == "MERIDUp" && m_muonsComboIDResShiftsSet)
	scale = 1.0 + m_muons.at(i).resShift("comboID").first;
      else if(systematic == "MERIDDown" && m_muonsComboIDResShiftsSet)
	scale = 1.0 + m_muons.at(i).resShift("comboID").second;
      else if(systematic == "MERMSUp" && m_muonsComboMSResShiftsSet)
	scale = 1.0 + m_muons.at(i).resShift("comboMS").first;
      else if(systematic == "MERMSDown" && m_muonsComboMSResShiftsSet)
	scale = 1.0 + m_muons.at(i).resShift("comboMS").second;
      else if(systematic == "MESUp" && m_muonsUncertaintiesSet)
	scale = 1.0 + m_muons.at(i).energyShift("default").first;
      else if(systematic == "MESDown" && m_muonsUncertaintiesSet)
	scale = 1.0 + m_muons.at(i).energyShift("default").second;

       _etx -= m_muons.at(i).Px()*m_muons.at(i).wex()*scale;
       _ety -= m_muons.at(i).Py()*m_muons.at(i).wey()*scale;
       _sumet += m_muons.at(i).Pt()*m_muons.at(i).wet()*scale;
       
    }//end of default muons
    else if(MissingETTags::usesSPECTRO(m_muons.at(i).statusWord())) {
     
      if(systematic == "None" || systematic == "none")
	scale = 1.0;
      else if(systematic == "MERMSUp" && m_muonsComboMSResShiftsSet)
	scale = 1.0 + m_muons.at(i).resShift("spectro").first;
      else if(systematic == "MERMSDown" && m_muonsComboMSResShiftsSet)
	scale = 1.0 + m_muons.at(i).resShift("spectro").second;
      else if(systematic == "MESUp" && m_muonsUncertaintiesSet)
	scale = 1.0 + m_muons.at(i).energyShift().first;
      else if(systematic == "MESDown" && m_muonsUncertaintiesSet)
	scale = 1.0 + m_muons.at(i).energyShift().second;
      
      _etx -= m_muons.at(i).Px("spectro")*m_muons.at(i).wex()*scale;
      _ety -= m_muons.at(i).Py("spectro")*m_muons.at(i).wey()*scale;
      _sumet += m_muons.at(i).Pt("spectro")*m_muons.at(i).wet()*scale;
    }//end of spectro
    // else if(MissingETTags::usesTRACK(m_muons.at(i).statusWord())) {
//        if(systematic == "None" || systematic == "none") scale = 1.0;
//        else if(systematic == "MERIDUp" && m_muonsResShiftsSet) scale = 1.0 + m_muons.at(i).resShift("comboID").first;
//        else if(systematic == "MERIDDown" && m_muonsResShiftsSet) scale = 1.0 + m_muons.at(i).resShift("comboID").second;
//        else if(systematic == "MERMSUp" && m_muonsResShiftsSet) scale = 1.0 + m_muons.at(i).resShift("comboMS").first;
//        else if(systematic == "MERMSDown" && m_muonsResShiftsSet) scale = 1.0 + m_muons.at(i).resShift("comboMS").second;
//        else if(systematic == "MESUp" && m_muonsUncertaintiesSet) scale = 1.0 + m_muons.at(i).energyShift().first;
//        else if(systematic == "MESDown" && m_muonsUncertaintiesSet) scale = 1.0 + m_muons.at(i).energyShift().second;

       
//        //this isn't any really different from default, but I keep them separate
//       _etx -= m_muons.at(i).Px()*m_muons.at(i).wex()*scale;
//       _ety -= m_muons.at(i).Py()*m_muons.at(i).wey()*scale;
//       _sumet += m_muons.at(i).Pt()*m_muons.at(i).wet()*scale;
//    }//end of track
    
    
  }//end of loop


  if(m_muonsMomentaSet == false && m_muonTotal.isValid()) {
    _etx += m_muonTotal.etx();
    _ety += m_muonTotal.ety();
    _sumet += m_muonTotal.sumet();
  }

  
  METObject _MuonTotal(_etx, _ety, _sumet);
  return _MuonTotal;
}//end of Muons


METObject METUtility::CellOut(string systematic) {

  float scale = 1.0;

  const int nBins = 34;
  double cellOutSumEtBounds[nBins] = {-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0}; //CellOut Sum Et Bins
  double cellOutFracUncert[nBins] = {0.136409, 0.153857, 0.153689, 0.149694, 0.138453, 0.139799, 0.133924, 0.139531, 0.138844, 0.134789, 0.138834, 0.132662, 0.142114, 0.13908, 0.13661, 0.128625, 0.133791, 0.125408, 0.122293, 0.133035, 0.132411, 0.131334, 0.128113, 0.13843, 0.12271, 0.129, 0.125528, 0.125005, 0.130233, 0.125374, 0.120071, 0.125658, 0.121132, 0.121722};//Cell Out Fraction uncertainty for bin ranges. 
  
  double slope = -1.19302e-4;//error 1.9e-05
  double intercept = 0.141673; //error 0.0022
  
//   if(!m_clustersMomentaSet) {
//     if(m_verbose) cout << "Clusters momenta and weights not loaded" << endl;
//     return m_cellOut;
//   }

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  string clusterEnergyCase = "default";
  string clusterResCase = "default";
  if(m_useStandardClusterRes) clusterResCase = "clusterStandard";
  if(m_useStandardClusterEnergySigma) clusterEnergyCase = "clusterStandard";

  for(unsigned int i = 0; i < m_clusters.size(); ++i) {

    if(MissingETTags::usesEFLOW_CLUSTER(m_clusters.at(i).statusWord())) {continue;}
    if(m_doForwardEtaCut && fabs(m_clusters.at(i).eta()) > 4.5) continue;
 
    scale = 1.0;
    if(systematic == "None" || systematic == "none") scale = 1.0;
    else if(systematic == "CERUp" && m_clustersResShiftsSet) scale = 1.0 + m_clusters.at(i).resShift(clusterResCase).first;
    else if(systematic == "CERDown" && m_clustersResShiftsSet) scale = 1.0 + m_clusters.at(i).resShift(clusterResCase).second;
    else if(systematic == "CESUp" && m_clustersUncertaintiesSet) scale = 1.0 + m_clusters.at(i).energyShift(clusterEnergyCase).first;
    else if(systematic == "CESDown" && m_clustersUncertaintiesSet) scale = 1.0 + m_clusters.at(i).energyShift(clusterEnergyCase).second;
    
    
    _etx -= m_clusters.at(i).Px()*m_clusters.at(i).wex()*scale;
    _ety -= m_clusters.at(i).Py()*m_clusters.at(i).wey()*scale;
    _sumet += m_clusters.at(i).Pt()*m_clusters.at(i).wet()*scale;
    
  }//end of loop
  
  if(m_clustersMomentaSet == false && m_cellOut.isValid()) {
    _etx += m_cellOut.etx();
    _ety += m_cellOut.ety();
    _sumet += m_cellOut.sumet();
  }


  //////////////////////////Cell corrections for object overlap --- ie, if the weight is <= .5, it's due to overlap between two objects. The cluster got assigned to one object, but cells remained attached to the other. They should have gone to CellOut, but didn't.
  if(m_doCellFix) {

    if(m_doRefEle) {

     for(unsigned int i = 0; i < m_electrons.size(); ++i) {
      if(fabs(m_electrons.at(i).wet()) > 0 && m_electrons.at(i).wet() <= 0.5) {
	if(m_doForwardEtaCut && fabs(m_electrons.at(i).eta()) > 4.5) continue;
 
	_etx -= m_electrons.at(i).Px()*m_electrons.at(i).wex();
	_ety -= m_electrons.at(i).Py()*m_electrons.at(i).wey();
	_sumet += m_electrons.at(i).Pt()*m_electrons.at(i).wet();

      }//weight <=.5
    }//end of loop
    }//end of m_doRefEle

    if(m_doRefGamma) {
    for(unsigned int i = 0; i < m_photons.size(); ++i) {
      if(fabs(m_photons.at(i).wet()) > 0 && m_photons.at(i).wet() <= 0.5) {
	if(m_doForwardEtaCut && fabs(m_photons.at(i).eta()) > 4.5) continue;
 
	_etx -= m_photons.at(i).Px()*m_photons.at(i).wex();
	_ety -= m_photons.at(i).Py()*m_photons.at(i).wey();
	_sumet += m_photons.at(i).Pt()*m_photons.at(i).wet();
      }//weight <=.5
    }//end of loop
    }//end of do RefGamma

    if(m_doRefTau) {
    for(unsigned int i = 0; i < m_taus.size(); ++i) {
      if(fabs(m_taus.at(i).wet()) > 0 && m_taus.at(i).wet() <= 0.5) {
	if(m_doForwardEtaCut && fabs(m_taus.at(i).eta()) > 4.5) continue;
 
	_etx -= m_taus.at(i).Px()*m_taus.at(i).wex();
	_ety -= m_taus.at(i).Py()*m_taus.at(i).wey();
	_sumet += m_taus.at(i).Pt()*m_taus.at(i).wet();
      }//weight <=.5
    }//end of loop
    }//end of if do tau

    if(m_doRefJet || m_doSoftJets) {
    for(unsigned int i = 0; i < m_jets.size(); ++i) {
      if(fabs(m_jets.at(i).wet()) > 0 && m_jets.at(i).wet() <= 0.5) {
	if(m_doForwardEtaCut && fabs(m_jets.at(i).eta()) > 4.5) continue;
 
	if(!m_doSoftJets && m_jets.at(i).Pt() <= m_softJetCut) continue;

	if(!m_doRefJet && m_jets.at(i).Pt() > m_softJetCut) continue;

	_etx -= m_jets.at(i).Px("secondary")*m_jets.at(i).wex();
	_ety -= m_jets.at(i).Py("secondary")*m_jets.at(i).wey();
	_sumet += m_jets.at(i).Pt("secondary")*m_jets.at(i).wet();
      }//weight <= .5
    }//end of loop
  }//end of if-jet

  }//end of cell fix 



  //strictly speaking this is for CellOutEflow, but I added it to have something
  int bin = -1;
  
  if(systematic == "CellOutUp" || systematic == "CellOutDown" || systematic == "AllClustersUp" || systematic == "AllClustersDown") {
    for(int i = 0; i < nBins-1; ++i) {
      if(fabs(_sumet)/1000.0 >= cellOutSumEtBounds[i] && fabs(_sumet)/1000.0 < cellOutSumEtBounds[i+1]) bin = i;
    }//end of loop
    float shift = cellOutFracUncert[bin];
    if(fabs(_sumet)/1000.0 >= cellOutSumEtBounds[nBins-1]) shift = slope*fabs(_sumet)/1000.0 + intercept;
    if(systematic == "CellOutUp"|| systematic == "AllClustersUp") scale = 1.0 + shift;
    else if(systematic == "CellOutDown" || systematic == "AllClustersDown") scale = 1.0 - shift;
    _etx *= scale;
    _ety *= scale;
    _sumet *= scale;
    
    
  }//end of MET PLHC CellOut Systematic



  if(systematic == "PileUpUp" || systematic == "PileUpDown") {
    if(systematic == "PileUpUp") scale = 1.0 + m_pileUpUncertainty;
    else if(systematic == "PileUpDown") scale = 1.0 - m_pileUpUncertainty;
    _etx *= scale;
    _ety *= scale;
    _sumet *= scale;
  }//end of pileup systematic
  
  

  METObject _cellOut(_etx, _ety, _sumet);
  return _cellOut;
    
  

}//end of CellOut



METObject METUtility::CellOutEflow(string systematic) {
  float scale = 1.0;

   const int nBins = 34;
  double cellOutSumEtBounds[nBins] = {-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0}; //CellOut Sum Et Bins
  double cellOutFracUncert[nBins] = {0.136409, 0.153857, 0.153689, 0.149694, 0.138453, 0.139799, 0.133924, 0.139531, 0.138844, 0.134789, 0.138834, 0.132662, 0.142114, 0.13908, 0.13661, 0.128625, 0.133791, 0.125408, 0.122293, 0.133035, 0.132411, 0.131334, 0.128113, 0.13843, 0.12271, 0.129, 0.125528, 0.125005, 0.130233, 0.125374, 0.120071, 0.125658, 0.121132, 0.121722};//Cell Out Fraction uncertainty for bin ranges. 
  
  double slope = -1.19302e-4;//error 1.9e-05
  double intercept = 0.141673; //error 0.0022
 
  //if(!m_clustersMomentaSet || !m_trackMomentaSet) {
  // if(m_verbose) cout << "Clusters momenta and weights not loaded" << endl;
  //return m_cellOutEflow;
    //}

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;

  string clusterEnergyCase = "default";
  string clusterResCase = "default";
  if(m_useStandardClusterRes) clusterResCase = "clusterStandard";
  if(m_useStandardClusterEnergySigma) clusterEnergyCase = "clusterStandard";

  for(unsigned int i = 0; i < m_clusters.size(); ++i) {

    if(!(MissingETTags::usesEFLOW_CLUSTER(m_clusters.at(i).statusWord()))) {continue;}

    if(m_doForwardEtaCut && fabs(m_clusters.at(i).eta()) > 4.5) continue;

    scale = 1.0;
    if(systematic == "None" || systematic == "none") scale = 1.0;
    else if(systematic == "CERUp" && m_clustersResShiftsSet) scale = 1.0 + m_clusters.at(i).resShift(clusterResCase).first;
    else if(systematic == "CERDown" && m_clustersResShiftsSet) scale = 1.0 + m_clusters.at(i).resShift(clusterResCase).second;
    else if(systematic == "CESUp" && m_clustersUncertaintiesSet) scale = 1.0 + m_clusters.at(i).energyShift(clusterEnergyCase).first;
    else if(systematic == "CESDown" && m_clustersUncertaintiesSet) scale = 1.0 + m_clusters.at(i).energyShift(clusterEnergyCase).second;
    
    _etx -= m_clusters.at(i).Px()*m_clusters.at(i).wex()*scale;
    _ety -= m_clusters.at(i).Py()*m_clusters.at(i).wey()*scale;
    _sumet += m_clusters.at(i).Pt()*m_clusters.at(i).wet()*scale;
    
  }//end of loop

  for(unsigned int i = 0; i < m_tracks.size(); ++i) {

    scale = 1.0;
    if(systematic == "None" || systematic == "none") scale = 1.0;
    else if(systematic == "TrkERUp" && m_tracksResShiftsSet) scale = 1.0 + m_tracks.at(i).resShift().first;
    else if(systematic == "TrkERDown" && m_tracksResShiftsSet) scale = 1.0 + m_tracks.at(i).resShift().second;
    else if(systematic == "TrkESUp" && m_trackUncertaintiesSet) scale = 1.0 + m_tracks.at(i).energyShift().first;
    else if(systematic == "TrkESDown" && m_trackUncertaintiesSet) scale = 1.0 + m_tracks.at(i).energyShift().second;
    _etx -= m_tracks.at(i).Px()*m_tracks.at(i).wex()*scale;
    _ety -= m_tracks.at(i).Py()*m_tracks.at(i).wey()*scale;
    _sumet += m_tracks.at(i).Pt()*m_tracks.at(i).wet()*scale;

  }//end of track loop

  if((m_trackMomentaSet == false || m_clustersMomentaSet == false) && m_cellOutEflow.isValid()) {
    if(m_verbose) cout << "Clusters momenta and weights not loaded" << endl;
    _etx = m_cellOutEflow.etx();
    _ety = m_cellOutEflow.ety();
    _sumet = m_cellOutEflow.sumet();
  }


  //////////////////////////Cell corrections for object overlap --- ie, if the weight is <= .5, it's due to overlap between two objects. The cluster got assigned to one object, but cells remained attached to the other. They should have gone to CellOut, but didn't.
  if(m_doCellFix) {
    if(m_doRefEle) {
      for(unsigned int i = 0; i < m_electrons.size(); ++i) {
	if(fabs(m_electrons.at(i).wet()) > 0 && m_electrons.at(i).wet() <= 0.5) {
	  if(m_doForwardEtaCut && fabs(m_electrons.at(i).eta()) > 4.5) continue;
	  
	  _etx -= m_electrons.at(i).Px()*m_electrons.at(i).wex();
	  _ety -= m_electrons.at(i).Py()*m_electrons.at(i).wey();
	  _sumet += m_electrons.at(i).Pt()*m_electrons.at(i).wet();
	}//weight <=.5
      }//end of loop
    }//end of m_doRefEle
    
    if(m_doRefGamma) {
      for(unsigned int i = 0; i < m_photons.size(); ++i) {
	if(fabs(m_photons.at(i).wet()) > 0 && m_photons.at(i).wet() <= 0.5) {
	  if(m_doForwardEtaCut && fabs(m_photons.at(i).eta()) > 4.5) continue;
	  
	  _etx -= m_photons.at(i).Px()*m_photons.at(i).wex();
	  _ety -= m_photons.at(i).Py()*m_photons.at(i).wey();
	  _sumet += m_photons.at(i).Pt()*m_photons.at(i).wet();
	}//weight <=.5
      }//end of loop
    }//end of do RefGamma
    
    if(m_doRefTau) {
      for(unsigned int i = 0; i < m_taus.size(); ++i) {
      if(fabs(m_taus.at(i).wet()) > 0 && m_taus.at(i).wet() <= 0.5) {
	if(m_doForwardEtaCut && fabs(m_taus.at(i).eta()) > 4.5) continue;
 
	_etx -= m_taus.at(i).Px()*m_taus.at(i).wex();
	_ety -= m_taus.at(i).Py()*m_taus.at(i).wey();
	_sumet += m_taus.at(i).Pt()*m_taus.at(i).wet();
      }//weight <=.5
    }//end of loop
    }//end of if do tau

    if(m_doRefJet || m_doSoftJets) {
    for(unsigned int i = 0; i < m_jets.size(); ++i) {
      if(fabs(m_jets.at(i).wet()) > 0 && m_doCellFix && m_jets.at(i).wet() <= 0.5) {
	if(m_doForwardEtaCut && fabs(m_jets.at(i).eta()) > 4.5) continue;
 
	if(!m_doSoftJets && m_jets.at(i).Pt() <= m_softJetCut) continue;

	if(!m_doRefJet && m_jets.at(i).Pt() > m_softJetCut) continue;
	//	cout << m_jets.at(i).Px("secondary") << " vs " << m_jets.at(i).Px() << " for the secondary" << endl;
	_etx -= m_jets.at(i).Px("secondary")*m_jets.at(i).wex();
	_ety -= m_jets.at(i).Py("secondary")*m_jets.at(i).wey();
	_sumet += m_jets.at(i).Pt("secondary")*m_jets.at(i).wet();
      }//weight <= .5
    }//end of loop
  }//end of if-jet

  }//end of cell fix 

 int bin = -1;

  if(systematic == "CellOutEflowUp" || systematic == "CellOutEflowDown" || systematic == "AllClustersUp" || systematic == "AllClustersDown") {
    for(int i = 0; i < nBins-1; ++i) {
      if(fabs(_sumet)/1000.0 >= cellOutSumEtBounds[i] && fabs(_sumet)/1000.0 < cellOutSumEtBounds[i+1]) bin = i;
    }//end of loop
    float shift = cellOutFracUncert[bin];
    if(fabs(_sumet)/1000.0 >= cellOutSumEtBounds[nBins-1]) shift = slope*fabs(_sumet)/1000.0 + intercept;
     
    if(systematic == "CellOutEflowUp" || systematic == "AllClustersUp") scale = 1.0 + shift;
    else if(systematic == "CellOutEflowDown" || systematic == "AllClustersDown") scale = 1.0 - shift;
    _etx *= scale;
    _ety *= scale;
    _sumet *= scale;
        
    
  }//end of MET PLHC CellOut Systematic




  if(systematic == "PileUpUp" || systematic == "PileUpDown") {
    if(systematic == "PileUpUp") scale = 1.0 + m_pileUpUncertainty;
    else if(systematic == "PileUpDown") scale = 1.0 - m_pileUpUncertainty;
    _etx *= scale;
    _ety *= scale;
    _sumet *= scale;
  }//end of pileup systematic
  
  
  METObject _cellOutEflow(_etx, _ety, _sumet);
  return _cellOutEflow;
  
}//end of CellOutEflow





METObject METUtility::SoftJets(string systematic) {
  float scale = 1.0;

   const int nBins = 62;
   double softJetsSumEtBounds[nBins] = {-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0}; //SoftJets Sum Et Bins

   double softJetsFracUncert[nBins] = {0, 0.00871549, 0.0897552, 0.0892624, 0.0841947, 0.0872207, 0.0894101, 0.090524, 0.0902687, 0.0909246, 0.0926338, 0.0917935, 0.0932977, 0.0892082, 0.093984, 0.0974575, 0.100168, 0.100101, 0.0979492, 0.100863, 0.102229, 0.102693, 0.108248, 0.10524, 0.100301, 0.104939, 0.110576, 0.103886, 0.0982754, 0.102644, 0.105031, 0.0989276, 0.101526, 0.0975329, 0.102695, 0.112649, 0.108279, 0.0992751, 0.111916, 0.107718, 0.115603, 0.115638, 0.107733, 0.111464, 0.103725, 0.102329, 0.117397, 0.112757, 0.103333, 0.0925, 0.085404, 0.114804, 0.150765, 0.117442, 0.104841, 0.112848, 0.109667, 0.116017, 0.113594, 0.107885, 0.111078, 0.112407}; //Cell Out Fraction uncertainty for bin ranges. 

   double slope = 3.88910e-05; //error = 1.76722e-05
   double intercept = 0.0993247;//error 3.03451e-03

  float _etx = 0;
  float _ety = 0;
  float _sumet = 0;
  
  for(unsigned int i = 0; i < m_jets.size(); ++i) {
    if(m_jets.at(i).Pt() <= m_softJetCut) {
      scale = 1.0;
      
      if(m_doCellFix && m_jets.at(i).wet() <= .5) continue;
      if(m_doForwardEtaCut && fabs(m_jets.at(i).eta()) > 4.5) continue;
 
      if(systematic == "None" || systematic == "none") scale = 1.0;
      //not applying these since they usually don't go down to 7 TeV. Use SoftJetsUp SoftJetsDown
      //else if(systematic == "JERUp" && m_jetsResShiftsSet) scale = 1.0 + m_jets.at(i).resShift().first;
      //else if(systematic == "JERDown" && m_jetsResShiftsSet) scale = 1.0 + m_jets.at(i).resShift().second;
      //else if(systematic == "JESUp" && m_jetsUncertaintiesSet) scale = 1.0 + m_jets.at(i).energyShift().first;
      //else if(systematic == "JESDown" && m_jetsUncertaintiesSet) scale = 1.0 + m_jets.at(i).energyShift().second;
      _etx -= m_jets.at(i).Px()*m_jets.at(i).wex()*scale;
      _ety -= m_jets.at(i).Py()*m_jets.at(i).wey()*scale;
      _sumet += m_jets.at(i).Pt()*m_jets.at(i).wet()*scale;
    }//end of softjet cut
  }//end of loop

  if(m_jetsMomentaSet == false && m_softJets.isValid()) {
    _etx = m_softJets.etx();
    _ety = m_softJets.ety();
    _sumet = m_softJets.sumet();
  }

 int bin = -1;

  if(systematic == "SoftJetsUp" || systematic == "SoftJetsDown" || systematic == "AllClustersUp" || systematic == "AllClustersDown") {
  for(int j = 0; j < nBins-1; ++j) {
    if(fabs(_sumet)/1000.0 >= softJetsSumEtBounds[j] && fabs(_sumet)/1000.0 < softJetsSumEtBounds[j+1]) bin = j;
  }//end of loop
  float shift = softJetsFracUncert[bin];
  if(fabs(_sumet)/1000.0 >= softJetsSumEtBounds[nBins-1]) shift = slope*fabs(_sumet)/1000.0 + intercept;
  if(systematic == "SoftJetsUp" || systematic == "AllClustersUp") scale = 1.0 + shift;
  else if(systematic == "SoftJetsDown" || systematic == "AllClustersDown") scale = 1.0 - shift;
  _etx *= scale;
  _ety *= scale;
  _sumet *= scale;

  }//end of MET PLHC SoftJets systemati



  if(systematic == "PileUpUp" || systematic == "PileUpDown") {
    if(systematic == "PileUpUp") scale = 1.0 + m_pileUpUncertainty;
    else if(systematic == "PileUpDown") scale = 1.0 - m_pileUpUncertainty;
    _etx *= scale;
    _ety *= scale;
    _sumet *= scale;
  }//end of pileup systematic
  
  

  METObject _softJets(_etx, _ety, _sumet);
  return _softJets;


}//end of SoftJets


float METUtility::METSignificance(string systematic) {

  float denominator = 0;
  float shift = 0.0;
  float deltaPhi = 0;
  if(m_doSignificance) {
    
    if(m_doRefEle) {
      if(m_electronsResolutionsSet) {
	for(unsigned int i = 0; i < m_electrons.size(); ++i) {
	  shift = 0.0;
	  if(systematic == "None" || systematic == "none") shift = 0.0;
	  else if(systematic == "EERUp" && m_electronsResShiftsSet) shift = m_electrons.at(i).resShift().first;
	  else if(systematic == "EERDown" && m_electronsResShiftsSet) shift = m_electrons.at(i).resShift().second;
	  
	  deltaPhi = m_refFinal.phi() - m_electrons.at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_electrons.at(i).resolution()*(m_electrons.at(i).resolution() + shift);
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_refEle.sumet();
      }
    }//end of if doRefEle


    if(m_doRefGamma) {
      if(m_photonsResolutionsSet) {
	for(unsigned int i = 0; i < m_photons.size(); ++i) {
	  shift = 0.0;
	  if(systematic == "None" || systematic == "none") shift = 0.0;
	 
	  deltaPhi = m_refFinal.phi() - m_photons.at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_photons.at(i).resolution()*(m_photons.at(i).resolution() + shift);
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_refGamma.sumet();
      }
    }//end of if doRefGamma
    

    if(m_doRefTau) {
      if(m_tausResolutionsSet) {
	for(unsigned int i = 0; i < m_taus.size(); ++i) {
	  shift = 0.0;
	  if(systematic == "None" || systematic == "none") shift = 0.0;
	  
	  deltaPhi = m_refFinal.phi() - m_taus.at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	  denominator += deltaPhi*deltaPhi*m_taus.at(i).resolution()*(m_taus.at(i).resolution() + shift);
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_refTau.sumet();
      }
    }//end of if doRefTau
    

    if(m_doRefJet) {
      if(m_jetsResolutionsSet) {
	for(unsigned int i = 0; i < m_jets.size(); ++i) {
	  if(m_jets.at(i).Pt() > m_softJetCut) {
	    shift = 0.0;
	    if(systematic == "None" || systematic == "none") shift = 0.0;
	    
	    deltaPhi = m_refFinal.phi() - m_jets.at(i).phi();
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	    deltaPhi = TMath::Cos(deltaPhi);
	    
	    denominator += deltaPhi*deltaPhi*m_jets.at(i).resolution()*(m_jets.at(i).resolution() + shift);
	  }//soft jet check	
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_refJet.sumet();
      }
    }//end of if doRefJet
    

    if(m_doMuonTotal) {
      if(m_muonsResolutionsSet) {
	for(unsigned int i = 0; i < m_muons.size(); ++i) {
	    shift = 0.0;

	    if(MissingETTags::usesDEFAULT(m_muons.at(i).statusWord())) {
	      if(systematic == "None" || systematic == "none") shift = 0.0;
	      
	    }//end of default muons
	    else if(MissingETTags::usesSPECTRO(m_muons.at(i).statusWord())) {
	      if(systematic == "None" || systematic == "none") shift = 0.0;
	     
	    }//end of spectro
	    else if(MissingETTags::usesTRACK(m_muons.at(i).statusWord())) {
	      if(systematic == "None" || systematic == "none") shift = 0.0;
	      
	    }//end of spectro
	    
	    deltaPhi = m_refFinal.phi() - m_muons.at(i).phi();
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	    deltaPhi = TMath::Cos(deltaPhi);
	    
	    denominator += deltaPhi*deltaPhi*m_muons.at(i).resolution()*(m_muons.at(i).resolution() + shift);
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_muonTotal.sumet();
      }
    }//end of if doMuonTotal
    
    if(m_doRefMuon) denominator += 0.5*0.5*m_scaled_refMuon.sumet();
    

    if(m_doSoftJets) {
      if(m_jetsResolutionsSet) {
	for(unsigned int i = 0; i < m_jets.size(); ++i) {
	  if(m_jets.at(i).Pt() <= m_softJetCut) {
	    shift = 0.0;
	    if(systematic == "None" || systematic == "none") shift = 0.0;
	    
	    deltaPhi = m_refFinal.phi() - m_jets.at(i).phi();
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	    deltaPhi = TMath::Cos(deltaPhi);
	    
	    denominator += deltaPhi*deltaPhi*m_jets.at(i).resolution()*(m_jets.at(i).resolution() + shift);
	  }//soft jet check	
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_refJet.sumet();
      }
    }//end of if doSoftJets
    

    if(m_doCellOut) {
      if(m_clustersResolutionsSet) {
	for(unsigned int i = 0; i < m_clusters.size(); ++i) {
	  if(MissingETTags::usesEFLOW_CLUSTER(m_clusters.at(i).statusWord())) {continue;}
	  
	  string clusterResCase = "default";
	  if(m_useStandardClusterRes) clusterResCase = "clusterStandard";
	  
	  shift = 0.0;
	  if(systematic == "None" || systematic == "none") shift = 0.0;
	  
	  deltaPhi = m_refFinal.phi() - m_clusters.at(i).phi();
	  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
	  deltaPhi = TMath::Cos(deltaPhi);
	  
	    denominator += deltaPhi*deltaPhi*m_clusters.at(i).resolution()*(m_clusters.at(i).resolution() + shift);
	}//end of loop
      }//end of resolutions filled requirement
      else{
	denominator += 0.5*0.5*m_scaled_cellOut.sumet();
      }
    }//end of if doCellOut
    



  }//end of if doSignificance

  return m_refFinal.et()/TMath::Sqrt(denominator);


}//end of METSignificance



vector<AnalysisFramework::External::PhysicsObject> METUtility::GetObjects(std::string objectName) {
  if(objectName == "Electrons" || objectName == "electrons") {
    return m_electrons;
   
  }
  else if(objectName == "Photons" || objectName == "photon") {
    return m_photons;
   
  }
  else if(objectName == "Taus" || objectName == "taus") {
    return m_taus;
    
  }
  else if(objectName == "Jets" || objectName == "jets") {
    return m_jets;
    
  }
  else if(objectName == "Muons" || objectName == "muons") {
    return m_muons;
    
  }
  else if(objectName == "Clusters" || objectName == "clusters") {
    return m_clusters;
    
  }
  else if(objectName == "Tracks" || objectName == "tracks") {
    return m_tracks;
  } 

 return m_jets;//by default, and so the function won't complain during compiling
}//end of GetObjects

vector<AnalysisFramework::External::PhysicsObject> METUtility::getPhysicsObject(std::string type) {
  if(type == "jets" || type == "Jets") {
    return m_jets;
  }
  if(type == "electrons" || type == "Electrons") {
    return m_electrons;
  }
  if(type == "photons" || type == "Photons") {
    return m_photons;
  }
  if(type == "taus" || type == "Taus") {
    return m_taus;
  }
  if(type == "muons" || type == "Muons") {
    return m_muons;
  }
  if(type == "clusters" || type == "Clusters") {
    return m_clusters;
  }
  if(type == "Tracks" || type == "tracks") {
    return m_tracks;
  }
  
  return m_jets;//default return

}//end of getPhysicsObject


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///Further helper function to deal with the fact that not every Performance group dumps the same variables///////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void METUtility::setJetParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  setObjects("jets", pT, eta, phi, E, wet, wex, wey, statusWord);
}//end of setJetParameters

void METUtility::setElectronParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  vector<float> E;// = new const vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.51099891);
    E.push_back(_hlv.E());
  }
  setObjectsHelper("electrons", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setElectronParameters


void METUtility::setPhotonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.0);
    E.push_back(_hlv.E());
  }
  setObjectsHelper("photons", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setPhotonParameters


void METUtility::setTauParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {

  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 1776.8);
    E.push_back(_hlv.E());
  }
  
  setObjectsHelper("taus", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setTauParameters


void METUtility::setClusterParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  setObjects("clusters", pT, eta, phi, E, wet, wex, wey, statusWord);
}//end of setClusterParameters


void METUtility::setTrackParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  setObjects("tracks", pT, eta, phi, E, wet, wex, wey, statusWord);
}//end of setTrackParameters


void METUtility::setMuonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord) {
  
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 105.658367);
    E.push_back(_hlv.E());
  }
  
  setObjectsHelper("muons", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setJetParameters

void METUtility::setExtraMuonParameters(const vector<float> *qOverPSpectro, const vector<float> *thetaSpectro, const vector<float> *phiSpectro, const vector<float> *charge) {

  vector<float> mu_staco_ms_eta;// = new vector<float>;
  vector<float> mu_staco_ms_pt;// = new vector<float>;
  vector<float> mu_staco_ms_E;// = new vector<float>;
  
  for(unsigned int iMu = 0; iMu < phiSpectro->size(); ++iMu) {
    float p = 0.0;
    if(qOverPSpectro->at(iMu) != 0.0) p = charge->at(iMu)/qOverPSpectro->at(iMu);
    
    float E = TMath::Sqrt(105.658367*105.658367 + p*p);
    mu_staco_ms_E.push_back(E);
    
    float p_T = p*TMath::Sin(thetaSpectro->at(iMu));
    mu_staco_ms_pt.push_back(p_T);
    
    float eta = -1.0*TMath::Log(TMath::Tan(thetaSpectro->at(iMu))/2.0);
    mu_staco_ms_eta.push_back(eta);
  }//end of spectro muon loop

  setObjectMomenta("spectroMuons", mu_staco_ms_pt, mu_staco_ms_eta, (*phiSpectro), mu_staco_ms_E);
  
}//end of setMuonExtraParameters




void METUtility::setExtraMuonParameters(const vector<float> *mu_staco_ms_pt, const vector<float> *thetaSpectro, const vector<float> *phiSpectro) {

  vector<float> mu_staco_ms_eta;// = new vector<float>;
  vector<float> mu_staco_ms_E;// = new vector<float>;
    
  float E = 0.0;
  float eta = 0.0;
   
   for(unsigned int iMu = 0; iMu < phiSpectro->size(); ++iMu) {
    float p = 0.0;
    if(TMath::Sin(thetaSpectro->at(iMu)) != 0) p = mu_staco_ms_pt->at(iMu)/TMath::Sin(thetaSpectro->at(iMu));  
    
    E = TMath::Sqrt(105.658367*105.658367 + p*p);
    mu_staco_ms_E.push_back(E);
    
    //float eta = -1.0*TMath::Log(TMath::Tan(thetaSpectro->at(iMu))/2.0);
    eta = -1.0*TMath::Log(fabs(TMath::Tan(thetaSpectro->at(iMu))/2.0))*TMath::Tan(thetaSpectro->at(iMu))/fabs(TMath::Tan(thetaSpectro->at(iMu)));
    mu_staco_ms_eta.push_back(eta);
    
  }//end of spectro muon loop

   setObjectMomenta("spectroMuons", (*mu_staco_ms_pt), mu_staco_ms_eta, (*phiSpectro), mu_staco_ms_E);
 
}//end of setMuonExtraParameters




void METUtility::setExtraJetParameters(const vector<float> *moment, const vector<float> */*mass*/, const vector<float> *eta, const vector<float> *phi) {

  vector<float> lc_pt;
  vector<float> lc_E;
   
   float scale = 1.0;
   float oldE = 0.0;
   float oldPt = 0.0;
   for(unsigned int i = 0; i < moment->size(); ++i) {
     scale = 1.0;
     if(moment->at(i) > 0) scale = 1.0/moment->at(i);
     //scale = 1.0;
     for(unsigned int j = 0; j < m_jets.size(); ++j) {
       if(i == m_jets.at(j).index()) {
	 oldE = m_jets.at(j).E();
	 oldPt = m_jets.at(j).Pt();
	 break;
       }//index matching if
     }//m_jets loop

     float newE = oldE*scale;
//      float newMass = mass->at(i)*scale;
//     float theta = 2.0*tanh(exp(-1.0*eta->at(i)));
     float newPt = oldPt*scale;//sqrt(newE*newE - newMass*newMass)*sin(theta);
          
     lc_pt.push_back(newPt);
     lc_E.push_back(newE);
     
   }//end of loop
  
   setObjectMomenta("lcjets", lc_pt, (*eta), (*phi),lc_E);
   
 }//end of function
 


void METUtility::setJetParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {
  setObjects("jets", pT, eta, phi, E, wet, wex, wey, statusWord);
}//end of setJetParameters

void METUtility::setElectronParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.51099891);
    E.push_back(_hlv.E());
  }
  setObjectsHelper("electrons", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setElectronParameters


void METUtility::setPhotonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 0.0);
    E.push_back(_hlv.E());
  }
  setObjectsHelper("photons", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setPhotonParameters


void METUtility::setTauParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {

  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 1776.8);
    E.push_back(_hlv.E());
  }
  
  setObjectsHelper("taus", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setTauParameters


void METUtility::setClusterParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {
  setObjects("clusters", pT, eta, phi, E, wet, wex, wey, statusWord);
}//end of setClusterParameters


void METUtility::setTrackParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {
  setObjects("tracks", pT, eta, phi, E, wet, wex, wey, statusWord);
}//end of setTrackParameters


void METUtility::setMuonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord) {
  
  vector<float> E;// = new vector<float>;
  for(unsigned int i = 0; i < pT->size(); ++i) {
    TLorentzVector _hlv;
    _hlv.SetPtEtaPhiM(pT->at(i), eta->at(i), phi->at(i), 105.658367);
    E.push_back(_hlv.E());
  }
setObjectsHelper("muons", (*pT), (*eta), (*phi), E, (*wet), (*wex), (*wey), (*statusWord));
}//end of setJetParameters

}
}
