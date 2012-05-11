//////////////////////////////////////////
//METUtility.h
//Authors: Jet Goodson
//First version: 5/5/2011
//
//The METUtility is intended to provide the user at D3PD level 
//with the correct recipe and information for rebuilding MET to include
//scaling and smearing, and for calculating the systematic 
//uncertainty on MET. 
//
//While it's designed for D3PD/Ntuples, we hope to keep it usable by AOD users
////////////////////////////////////////////
//Instructions:
//Prior to root macro
// .L METUtility.h+
// 
////////////////////////////////////////////

#ifndef _MYMETUTILITY_
#define _MYMETUTILITY_

#include "TNamed.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <string>
#include <utility>

using namespace std;

namespace AnalysisFramework
{
namespace External
{

struct MissingETTags {//copied from MissingETEvent/MissingETComposition -- changed tag to Tag to avoid potential conflicts

  enum Tags {
    UNKNOWN           = 0x0000,
    DEFAULT           = 0x0001,
    SPECTRO           = 0x0002,
    TRACK             = 0x0004,
    REFMUON           = 0x0008,
    MUID              = 0x0010,
    EFLOW_CLUSTER     = 0x0020
    //NEXT_COMES      = 0x0040,
    //THEN            = 0x0080,
    //AND_THEN        = 0x0100,
  };

//isDEFAULT means use default pt and no other bits are set
static bool isDEFAULT(unsigned short tag)
{ return tag == DEFAULT; }

//usesDEFAULT means use default pt but other bits could also be set
static bool usesDEFAULT(unsigned short tag)
{ return tag & DEFAULT; }

static bool usesSPECTRO(unsigned short tag)
{ return tag & SPECTRO; }

static bool usesTRACK(unsigned short tag)
{ return tag & TRACK; }

static bool usesREFMUON(unsigned short tag)
{ return tag & REFMUON; }

static bool usesMUID(unsigned short tag)
{ return tag & MUID; }

static bool usesEFLOW_CLUSTER(unsigned short tag)
{ return tag & EFLOW_CLUSTER; }

static void setUsesDEFAULT(unsigned short &tag)
{ tag = tag | DEFAULT;}

static void setUsesSPECTRO(unsigned short &tag)
{ tag = tag | SPECTRO;}

static void setUsesTRACK(unsigned short &tag)
{ tag = tag | TRACK;}

static void setUsesREFMUON(unsigned short &tag)
{ tag = tag | REFMUON;}

static void setUsesMUID(unsigned short &tag)
{ tag = tag | MUID;}

static void setUsesEFLOW_CLUSTER(unsigned short &tag)
{ tag = tag | EFLOW_CLUSTER;}
};

class METObject : public TNamed {
 public:
  METObject(){m_etx = 0.0; m_ety = 0.0; m_sumet = 0.0; m_et = -1.0; m_significance = 0; m_isValid = false;}//will return false until changed, preferably after etx/ety/sumet are set
  METObject(float _etx, float _ety, float _sumet){m_significance = 0; m_sumet = _sumet; m_etx = _etx; m_ety = _ety; m_et = -1.0; m_isValid = true;}
  ~METObject(){}

  float etx() {return m_etx;}
  float ety() {return m_ety;}
  float sumet() {return m_sumet;}
  float phi() {return atan2(m_ety, m_etx);}
  float et() {return ((m_et < 0) ? sqrt(m_etx*m_etx + m_ety*m_ety) : m_et);}
  float sig() {return et()/(.5*TMath::Sqrt(m_sumet));}
  float significance(){return m_significance;} //for the more complex MET significance
  bool isValid(){return m_isValid;}

  void setEtx(float _etx){m_etx = _etx;}
  void setEty(float _ety){m_ety = _ety;}
  void setSumet(float _sumet){m_sumet = _sumet;}
  void setEt(float _et){m_et = _et;} //don't do this unless you know what you're doing, it's mostly for systematic diffs. if you're doing actually MET let etx/ety compute et
  void setBase(float _etx, float _ety, float _sumet){m_sumet = _sumet; m_etx = _etx; m_ety = _ety;}
  void setSignificance(float _sig){m_significance = _sig;}
  void setIsValid(bool status){m_isValid = status;}
 private:
  float m_sumet;
  float m_etx;
  float m_ety;
  float m_et; //this is mostly for return systematic diffs, i.e., (up - down)/none
  float m_significance;
  bool m_isValid;

  friend class METUtility;


//#ifdef METUTIL_STANDALONE
  ClassDef(METObject,1);
//#endif //METUTIL_STANDALONE

};

///////////////////////////////////////////////////////////////////////////////////

class PhysicsObject : public TNamed {
 public:
  PhysicsObject(){m_relativeResolution = -1.0;}
  ~PhysicsObject(){}

  void setMomenergy(float _pT, float _eta, float _phi, float _E, string caseFlag="default");//{m_momenergy->SetPtEtaPhiE(pT, eta, phi, E);}
  void setWeights(float _wex, float _wey, float _wet){m_weights[0] = _wex; m_weights[1] = _wey; m_weights[2] = _wet;}
  void setStatusCode(unsigned int status){m_statusWord = status;}
  void setEnergyUncertainty(float energyUp, float energyDown, string caseFlag="default");//{m_relEnergyUncertainty.first = energyUp; m_relEnergyUncertainty.second = energyDown;}
  void setResolution(float res);
  void setResShift(float resUp, float resDown, string caseFlag="default");
  // void setIsMuon(bool status){m_isMuon = status;}
  void setIndex(int _index){m_index = _index;}
 
  float E(string caseFlag="default");
  float Et(string caseFlag="default");
  float Pt(string caseFlag="default");
  float Px(string caseFlag="default");
  float Py(string caseFlag="default");
  float phi(string caseFlag="default");
  float eta(string caseFlag="default");
  float wex(){return m_weights[0];}
  float wey(){return m_weights[1];}
  float wet(){return m_weights[2];}
  unsigned int statusWord(){return m_statusWord;}
  pair<float, float> resShift(string caseFlag="default");
  pair<float, float> energyShift(string caseFlag="default");
  float resolution();
  unsigned int index(){return m_index;}
 
 private:
  TLorentzVector m_momenergy;
  TLorentzVector m_secondaryMomenergy;
  float m_weights[3];
  unsigned int m_statusWord;
  float m_relativeResolution;
  unsigned int m_index;   //used for keeping track of uncertainties and weights for in vector<vector<> >


  pair<float, float> m_relEnergyUncertainty; 
  pair<float, float> m_relativeResolutionShift;
 

  //bool m_isMuon; //muons have different pts, combined, track, spectro
  //let's just ComboID  in m_relativeResolutionShift
  pair<float, float> m_relativeResolutionComboMS;
  pair<float, float> m_relativeResolutionSpectro;
  
  
  friend class METUtility;

//#ifdef METUTIL_STANDALONE
  ClassDef(PhysicsObject,1);
//#endif //METUTIL_STANDALONE

};

class METUtility : public TNamed {

 public:
  METUtility(bool verbose = false);
  METUtility(bool doRefEle, bool doRefGamma, bool doRefTau, bool doRefJet, bool doSoftJets, bool doRefMuon, bool doMuonTotal, bool doCellOut, bool doCellOutEflow, bool isMuid, double softJetCut, bool verbose=false);
  ~METUtility(){}

//reset all data
  void reset();

//objects to read things out
  METObject getMissingET(string term, string systematic="none");
  METObject deltaMissingET(string term, string systematics = "All");
  METObject absDeltaMissingET(string term, string systematics = "All");
  METObject MissingETHelper(string term, string systematic="none");
  

  
  //objects to set things up
   void setObjects(string type, const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);
   void setObjects(string type, const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);

   void setObjectsHelper(string type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, vector<vector<float> > wet, vector<vector<float> > wex, vector<vector<float> > wey, vector<vector<unsigned int> > statusWord);
   void setObjectsHelper(string type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E,vector<float> wet, vector<float> wex, vector<float> wey, vector<unsigned short> statusWord);

  void setObjectEnergyUncertainties(string type, const vector<float> energyUp, const vector<float> energyDown); 
  void setObjectEnergyUncertaintiesHelper(vector<PhysicsObject> &objects, vector<float> energyUp, vector<float> energyDown, string caseFlag="default");
 
  void setObjectResolutionShift(string type, const vector<float> resUp, const vector<float> resDown); 
  void setObjectResolutionShiftHelper(vector<PhysicsObject> &objects, vector<float> resUp, vector<float> resDown, string caseFlag="default");
  
  void setObjectResolutions(string type, const vector<float> resolutions);
  void setObjectResolutionsHelper(vector<PhysicsObject> &objects, vector<float> resolutions, string caseFlag="default");
 
  void setObjectMomenta(string type, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E);
  void setObjectMomentaHelper(vector<PhysicsObject> &objects, vector<float> pT, vector<float> eta, vector<float> phi, vector<float> E, string caseFlag ="default");

  void setMETTerm(string term, float _etx, float _ety, float _sumet);
  void setSoftJetCut(float cut){m_softJetCut = cut;}
  void setIsMuid(bool status){m_isMuid = status;}
  void setVerbosity(bool status){m_verbose = status;}
  void setCellFix(bool status){m_doCellFix = status;}
  void doSignificance(bool status){m_doSignificance = status;}
  void doForwardEtaCut(bool status){m_doForwardEtaCut = status;}
  void setPileUpUncertainty(double uncert){m_pileUpUncertainty = uncert;}


  void setClusterDefaults(bool energyStat, bool resStat=false){m_useStandardClusterRes = resStat; m_useStandardClusterEnergySigma = energyStat;}

  void defineMissingET(bool doRefEle, bool doRefGamma, bool doRefTau, bool doRefJet, bool doSoftJets, bool doRefMuon, bool doMuonTotal, bool doCellOut, bool doCellOutEflow);

  float METSignificance(string systematic);
 
  vector<PhysicsObject> GetObjects(std::string objectName);

  
  METObject RefEle(string systematic);
  METObject RefGamma(string systematic);
  METObject RefTau(string systematic);
  METObject RefJet(string systematic);
  METObject RefMuon(string systematic);
  METObject MuonTotal(string systematic);
  METObject CellOut(string systematic);
  METObject CellOutEflow(string systematic);
  METObject SoftJets(string systematic);

   

  ///helper functions to deal with odd variables
  void setJetParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);
  void setElectronParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);
  void setPhotonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);
  void setTauParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);

  void setClusterParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);
  void setTrackParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);

  void setMuonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<vector<float> > *wet, const vector<vector<float> > *wex, const vector<vector<float> > *wey, const vector<vector<unsigned int> > *statusWord);
  void setExtraMuonParameters(const vector<float> *qOverPSpectro, const vector<float> *thetaSpectro, const vector<float> *phiSpectro, const vector<float> *charge);
  void setExtraMuonParameters(const vector<float> *mu_staco_ms_pt, const vector<float> *thetaSpectro, const vector<float> *phiSpectro);



  void setExtraJetParameters(const vector<float> *moment, const vector<float> *mass, const vector<float> *eta, const vector<float> *phi);
 


  void setJetParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);
  void setElectronParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);
  void setPhotonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);
  void setTauParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);

  void setClusterParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);
  void setTrackParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *E, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);

  void setMuonParameters(const vector<float> *pT, const vector<float> *eta, const vector<float> *phi, const vector<float> *wet, const vector<float> *wex, const vector<float> *wey, const vector<unsigned short> *statusWord);

  vector<PhysicsObject> getPhysicsObject(std::string type="jets");

 private:
  vector<PhysicsObject> m_jets;
  vector<PhysicsObject> m_electrons;
  vector<PhysicsObject> m_muons;
  vector<PhysicsObject> m_taus;
  vector<PhysicsObject> m_clusters;
  vector<PhysicsObject> m_tracks;
  vector<PhysicsObject> m_photons;
  


  METObject m_cellOut; //these hold the original inputs, for instance cell out we don't want to feed a systematically altered term to the next systematic
  METObject m_cellOutEflow;
  METObject m_softJets;
  METObject m_refMuon;
  METObject m_refJet;
  METObject m_refEle;
  METObject m_refTau;
  METObject m_refGamma;
  METObject m_muonTotal;
  

  METObject m_scaled_cellOut; //these hold the current scaled MET term
  METObject m_scaled_cellOutEflow;
  METObject m_scaled_softJets;
  METObject m_scaled_refMuon;
  METObject m_scaled_refJet;
  METObject m_scaled_refEle;
  METObject m_scaled_refTau;
  METObject m_scaled_refGamma;
  METObject m_scaled_muonTotal;
  
  METObject m_refFinal;
 
  //control variables

  //check whether momenta are set
  bool m_jetsMomentaSet;
  bool m_photonsMomentaSet;
  bool m_electronsMomentaSet;
  bool m_muonsMomentaSet;
  bool m_tausMomentaSet;
  bool m_clustersMomentaSet;
  bool m_trackMomentaSet;


  //check whether uncertainties are set
  bool m_jetsUncertaintiesSet;
  bool m_photonsUncertaintiesSet;
  bool m_electronsUncertaintiesSet;
  bool m_muonsUncertaintiesSet;
  bool m_muonsTrackUncertaintiesSet;
  bool m_muonsSpectroUncertaintiesSet;
  bool m_tausUncertaintiesSet;
  bool m_clustersUncertaintiesSet;
  bool m_trackUncertaintiesSet;
  bool m_muonsComboIDUncertaintiesSet;
  bool m_muonsComboMSUncertaintiesSet;

  //check whether resShifts are set
  bool m_jetsResShiftsSet;
  bool m_photonsResShiftsSet;
  bool m_electronsResShiftsSet;
  bool m_muonsResShiftsSet;
  bool m_muonsTrackResShiftsSet;
  bool m_muonsSpectroResShiftsSet;
  bool m_muonsComboIDResShiftsSet;
  bool m_muonsComboMSResShiftsSet;
  bool m_tausResShiftsSet;
  bool m_clustersResShiftsSet;
  bool m_tracksResShiftsSet;

  //check whether resolutions are set
  bool m_jetsResolutionsSet;
  bool m_photonsResolutionsSet;
  bool m_electronsResolutionsSet;
  bool m_muonsResolutionsSet;
  bool m_muonsTrackResolutionsSet;
  bool m_muonsSpectroResolutionsSet;
  bool m_muonsComboIDResolutionsSet;
  bool m_muonsComboMSResolutionsSet;
  bool m_tausResolutionsSet;
  bool m_clustersResolutionsSet;
  bool m_tracksResolutionsSet;

  //which terms to do
  bool m_doRefEle;
  bool m_doRefGamma;
  bool m_doRefTau;
  bool m_doRefJet;
  bool m_doRefMuon;
  bool m_doMuonTotal;
  bool m_doSoftJets;
  bool m_doCellOut;
  bool m_doCellOutEflow;
 
  
  double m_pileUpUncertainty;
  bool m_doForwardEtaCut;
  bool m_doCellFix;
  bool m_doCellFixOverride;//this is the manual shutoff, if true cellfix is allowed, if false, cellfix is not allowed at all this way we can ovreride it for release 17 even if a new tag of this tool isn't immediatly available

  float m_softJetCut; //the cut used to put jets in RefJet or SoftJets
  bool m_isMuid; //are muons muid

  bool m_useStandardClusterRes;        //switches to do the default cluster resolution and energy
  bool m_useStandardClusterEnergySigma;
  
  bool m_doSignificance;

  bool m_verbose;

  

//#ifdef METUTIL_STANDALONE
  ClassDef(METUtility,1);
//#endif //METUTIL_STANDALONE
 
};//end of METUtility class

}
}

#endif // _METUTILITY_
