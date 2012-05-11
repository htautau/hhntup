
///***************************************************
///
/// Class: JERProviderMMCMoriond
/// Author: Gaston Romeo <glromeo@cern.ch>
/// Plots by: Gaston Romeo <glromeo@cern.ch>
///
/// Provide the JER and its uncertainty.
/// Created: Jan/19
/// Version (see cmt/version.cmt)
///
/// (The uncertainty corresponds to the systematic error and it is 100% correlated point by point.
/// The statistical error is found to be negligible)
///
/// Based on ATL-COM-PHYS-2011-240
///
/// 
/// ----- Usage: https://twiki.cern.ch/twiki/bin/view/Main/JetEnergyResolutionProvider
///
/// 1) Link the final library, located in: JERProviderMMCMoriond/StandAlone/libJERProviderMMCMoriond.so
/// 2) Create an instance, ie:
///    JERProviderMMCMoriond myJER("AntiKt6TopoJES","Truth","JERProviderMMCMoriondPlots.root");
/// 3) Initialize the instance: myJER.init();
/// 3) Call myJER.getSigma(pt,y) to get the resolution (pt in GeV)
/// 4) Call myJER.getUncert(pt,y) to get its uncertainty (pt in GeV)   
///     
///   
///   
///***************************************************
#ifndef _JERProviderMMCMoriondMMCMORIOND_
#define _JERProviderMMCMoriondMMCMORIOND_

#include "TNamed.h"
#include "TFile.h"
#include "TF1.h"
#include <TRandom.h>
#include <iostream>
//#include <stream.h>
#include <cmath>
#include <map>

using std::cout;
using std::endl;
using std::map;

namespace AnalysisFramework
{
namespace External
{

class JERProviderMMCMoriond : public TNamed
{

 public: 

  // Constructor, destructor
  JERProviderMMCMoriond() { }
  JERProviderMMCMoriond(std::string CollectionName, std::string MethodName = "Truth", std::string FileName="JERProviderPlots.root");
  virtual ~JERProviderMMCMoriond();
  
  //initialization
  void init();

  // Read the parametrization
  TF1* getParam(double y);

  // Read sigma
  float getSigma(double pt, double y);

  // Read the uncertainty
  float getUncert(double pt, double y);

  // inputCollection 
  void setInputCollection(std::string CollectionName, std::string MethodName);
  
/*#ifdef JER_STANDALONE*/
	ClassDef(JERProviderMMCMoriond,1);
/*#endif*/

 private:
  
  static const int m_nParam = 6;
  //static const double m_GeV = 1000.;
  static const int m_nY = 6;

  // strings
  std::string m_collectionName;
  std::string m_methodName;
  std::string m_fileName;
  
  //initialization flag
  bool m_isInit;

  // Input file
  TFile* m_inputFile; //!

  // Pointer to the correct graph
  map<int,TF1*> m_jerFunc;

  // Parameters and their errors
  map<int,float> m_param;

  // Uncertainty from Data/MC comparison
  float m_uncert;

  // Offset from Data/MC comparison
  float m_offset;

};

}
}
#endif
