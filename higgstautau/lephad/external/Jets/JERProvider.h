
///***************************************************
///
/// Class: JERProvider
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
/// 1) Link the final library, located in: JERProvider/StandAlone/libJERProvider.so
/// 2) Create an instance, ie:
///    JERProvider myJER("AntiKt6TopoJES","Truth","JERProviderPlots.root");
/// 3) Initialize the instance: myJER.init();
/// 3) Call myJER.getSigma(pt,y) to get the resolution (pt in GeV)
/// 4) Call myJER.getUncert(pt,y) to get its uncertainty (pt in GeV)
///
///
///
///***************************************************
#ifndef _MYJERPROVIDER_
#define _MYJERPROVIDER_

#include "TNamed.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
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

class JERProvider : public TNamed
{

 public:

  // Constructor, destructor
  JERProvider() { }
  JERProvider(std::string CollectionName, std::string MethodName = "Truth", std::string FileName="JERProviderPlots.root");
  virtual ~JERProvider();

  // filename -- only allowed before initialisation
  void setFileName(std::string fileName);

  //initialization
  void init();

  // Read the parametrization
  TF1* getParam(double y);

  // Read the offset from data/MC
  TGraphErrors* getOffset(double y);



  // Read the error from data/MC
  TGraphErrors* getSyst(double y);

  // inputCollection
  void setInputCollection(std::string CollectionName, std::string MethodName);

  // Read the offset from data/MC for the pT bin considered
  float getOffset(double pt, double y);

  // Read the resolution from the MC parameterisation
  float getRelResolutionMC(double pt, double y);

  // Obtain the resolution for data (sigma_MC + sigma_MC*offset)
  float getRelResolutionData(double pt, double y);

  // Obtain the uncertainty on the resolution for data
  float getResolutionUncert(double pt, double y);

  // The smearing extra factor for shifting the central value of MC
  float getSmearingFactorMC(double pt, double y);

  // The smearing extra factor for the systematic studies
  float getSmearingFactorSyst(double pt, double y);

/*#ifdef JER_STANDALONE*/
	ClassDef(JERProvider,1);
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

  // Pointer to the correct graph
  map<int,TGraphErrors*> m_jerOffset;

  // Pointer to the correct graph
  map<int,TGraphErrors*> m_jerSyst;

  // Uncertainty from Data/MC comparison
  float m_uncert;

  // Offset from Data/MC comparison
  float m_offset;

  // Systematics
  float m_syst;


};

}
}

#endif

