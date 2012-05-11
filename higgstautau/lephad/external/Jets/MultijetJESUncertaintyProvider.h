#ifndef _MYMULTIJETJESUNCERTAINTYPROVIDER_
#define _MYMULTIJETJESUNCERTAINTYPROVIDER_

#include "JESUncertaintyProvider.h"

namespace AnalysisFramework
{
namespace External
{

class MultijetJESUncertaintyProvider : public JESUncertaintyProvider{

 public:

  // Constructors, destructor
  MultijetJESUncertaintyProvider(std::string CollectionName="AntiKt6TopoJetsEM", std::string AnalysisFileName="MJESUncertainty.root", std::string MJESFileName="MJESUncertainty.root", std::string FileName="JESUncertainty.root");
   ~MultijetJESUncertaintyProvider();
 

  // Initialize the provider
  virtual void init();
  void includeFlavorComposition(bool include = true);

  // Open input file
  //virtual void openInputFile(std::string FileName);




  // Relative positive uncertainty - by default no pile-up term, no close-by jets systematic (dR = 3.0)
  double getRelPosUncert(double pT, double Eta, double dRmin=3.0, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);
  // Absolute positive uncertainty - by default no pile-up term, no close-by jets systematic (dR = 3.0)
  double getAbsPosUncert(double pT, double Eta, double dRmin=3.0, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);
  // Relative negative uncertainty - by default no pile-up term, no close-by jets systematic (dR = 3.0)
  double getRelNegUncert(double pT, double Eta, double dRmin=3.0, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);
  // Absolute negative uncertainty - by default no pile-up term, no close-by jets systematic (dR = 3.0)
  double getAbsNegUncert(double pT, double Eta, double dRmin=3.0, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);

  // Relative uncertainty - by default no pile-up term, no close-by jets systematic (dR = 3.0)
  // isUp = true if relative positive uncertainty is requested
  // isUp = false if relative negative uncertainty is requested
  double getRelUncert(double pT, double Eta, double dRmin, bool isUp, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);

  // Get a copy of the 2D Graph containing the Uncertainties
  // isUp = true if relative positive uncertainty is requested
  // isUp = false if relative negative uncertainty is requested
  // to be implemented
//  TH2D* getUncGraphCopy(bool isUp, Components UncertComps, unsigned int nVtx);

	//#ifdef JES_STANDALONE
	ClassDef(MultijetJESUncertaintyProvider,1);
	//#endif //JES_STANDALONE

  private:
  
  //Pointers to the flavor composition histograms
  TH2D*  m_flavorCompGluGraph; 
  TH2D*  m_flavorCompLightGraph; 

  //Pointer to the deltaR histogram
  TH2D*  m_deltaRGraph; 

  //Pointers to the sample related histograms
  //Pointer to the gluon fraction histogram
  TH2D*  m_gluonFraction; 
  //Pointer to the gluon fraction error histogram
  TH2D*  m_gluonFractionError; 
  //Pointer to the response histogram
  TH2D*  m_responseSample; 

  // Input File with sample related histograms(gluon fraction, gluon fraction error, sample response)
  TFile* m_analysisInputFile; //! transient, don't save
  // Input File with the MJES related plots
  TFile* m_MJESInputFile; //! transient, don't save  


  // Name of analysis file
  std::string m_analysisFileName;
  // Name of MJES file
  std::string m_MJESFileName;

  
  
  // Decide whether to include the flavor compositio term
  bool m_includeFlavorComp;
  

  // Count how often error message is printed
  static int m_counter;

  // Helper function to increment m_counter
  static void incrementCounter(){
  m_counter++;
 }

  // Change the inputCollection
  virtual bool setInputCollection(std::string CollectionName);

  // Helper function to check requested components (global uncertainty from JESUncertaintyProvider) and adding sample specific flavor composition uncertainty on top of it
  double getComponents(double pT, double Eta, double dRmin, bool isUp, Components UncertComps, unsigned int nVtx);

  // Helper function to calculate prefactor alphaC
  double getAlphaC(int currentBin, bool isUp);

  // Helper function to calculate prefactor alphaR
  double getAlphaR(int currentBin, bool isUp, double alphaC, double avgResponse);

  // Helper funtion to return average response - protect against r_avg>1+Delta+ and r_avg < 1-Delta-
  double getResponseSample(int currentBin);

};

}
}
#endif
