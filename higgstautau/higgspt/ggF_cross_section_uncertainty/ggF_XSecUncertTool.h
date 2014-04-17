#ifndef _GGF_XSECUNCERTTOOL_
#define _GGF_XSECUNCERTTOOL_

/**********************************************
 *  ggF_XSecUncertTool
 *
 *  Tool to evaluate the Njet dependent ggF cross section uncertainty
 *
 *  Author: Dag Gillberg
 *
 **********************************************
 */

#include <iostream>
#include <vector>

#include <TString.h>
#include <TH1F.h>
#include <TEnv.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>

using namespace std;

typedef TString Str;
typedef std::vector<Str> StrV;
typedef std::vector<double> VecD;
typedef std::vector<VecD> VVecD;

class ggF_XSecUncertTool : public TNamed {

 public:
  /*  Constructor
   *  Parameters:
   *    1. method JVE or BLPTW
   *       (Jet-Veto Efficinecy or the combined resummed calculation by BLPTW)
   *    2. configFile - name of TEnv file with input values
   */  
  ggF_XSecUncertTool( Str method, Str configFile="data/ggF_uncert_ATLAS_jets.config", 
		      bool verbose=true );
  
  // Returns the relative uncertainty amplitude for each of the
  // uncertainty components
  double relUncert(int iComp, int NtruthJets);
  
  double yieldUncert(int Ntruthjets) { return relUncert(0,Ntruthjets); }
  double migration01Uncert(int Ntruthjets) { return relUncert(1,Ntruthjets); }
  double migration12Uncert(int Ntruthjets) { return relUncert(2,Ntruthjets); }
  // this method is will abort with an error message with the current (Jan 2014) inputs
  // as we don't yet have any 3-jet input values
  double migration23Uncert(int Ntruthjets) { return relUncert(3,Ntruthjets); }


  // Sum in quadrature of the above, DON'T USE IN ANALYSIS
  double totalXsecUncert(int Ntruthjets);
  VecD totalXsecUncert();
  
  // Event weights for uncertainty propagation
  // Simply obtained by 1+uncert for "up" variation
  // and 1-uncert for "down" variation, where
  // uncert is the relative uncert. amplitude of the three above methods
  double uncertShift(int iComp, int Ntruthjets, bool isUp) {
    return isUp ? 1.0+relUncert(iComp,Ntruthjets) : 1.0-relUncert(iComp,Ntruthjets);
  }
  double uncertShiftUp(int iComp, int Ntruthjets) { return uncertShift(iComp,Ntruthjets,true); }
  double uncertShiftDown(int iComp, int Ntruthjets) { return uncertShift(iComp,Ntruthjets,false); }

  double yieldShiftUp(int Ntruthjets)   { return uncertShift(0,Ntruthjets,true); }
  double yieldShiftDown(int Ntruthjets) { return uncertShift(0,Ntruthjets,false); }
  double migration01ShiftUp(int Ntruthjets)   { return uncertShift(1,Ntruthjets,true); }
  double migration01ShiftDown(int Ntruthjets) { return uncertShift(1,Ntruthjets,false); }
  double migration12ShiftUp(int Ntruthjets)   { return uncertShift(2,Ntruthjets,true); }
  double migration12ShiftDown(int Ntruthjets) { return uncertShift(2,Ntruthjets,false); }

  // dump the NP uncertainty amplitudes to the screen
  void printNPs();

  // print correlation details to the screen
  void printCorrelationDetails();

  // A few methods for fun
  // double ggFsigmaTot()
 private:
  
  void init(Str configFile);
  
  void fatal(Str msg) { printf("\n  ggF_XSecUncertTool\n  FATAL:\n\n    %s\n\n",msg.Data()); abort(); }

  // Get value from the TEnv file
  double getValue(Str key);
  VecD getValues(Str key);

  // helpers
  void printNP(Str uName, const VecD &relU);
  //double addInQuad(double a, double b) { return sqrt(a*a+b*b); }
  //double addInQuad(double a, double b, double c) { return sqrt(a*a+b*b+c*c); }
  void add(VecD &vec, double a) { vec.push_back(a); }
  //void add(VecD &vec, double a, double b) { add(vec,a); add(vec,b); }
  //void add(VecD &vec, double a, double b, double c) { add(vec,a,b); add(vec,c); }

  // uncertaintiy covariance and correlation between jet bins
  double cov(int ijet, int jjet);
  double corr(int ijet, int jjet);

  /*
   *  Private variables
   */

 private:
  
  // pointer to - and name of - input file
  TEnv *m_input;
  Str m_inputFN;

  // Which method to use: JVE, ST or RIST ?
  Str m_method;

  // How many jet bins?
  int m_Njetbins;

  // jet definition
  double m_pTcut, m_jetR;

  // verbose or not?
  bool m_verb;

  // vector of vectors of relative uncertainty amplitudes
  // (Njets x Njets) elements
  VVecD m_relUncert;

  // cross section central values 
  // associated with uncertainties
  VecD m_sigExcl, m_sigIncl;
};

#endif
