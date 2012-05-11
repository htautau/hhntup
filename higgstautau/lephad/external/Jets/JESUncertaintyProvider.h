///***********************************************************************
///
/// JESUncertaintyProvider
/// Authors: C. Doglioni, P.O. DeViveiros
/// Plots by: C. Doglioni, M. Duehrssen, D. Gillberg, D. Miller
/// PROOF-compatibility: B. Butler
///
/// Provide a jet-by-jet uncertainty on the Jet Energy Scale (assumed 100% correlated bin by bin).
/// Each unique component of the uncertainty can be accessed separately.
/// 
/// See the Components enum below for a list of components of the JES 
/// which can be accessed separately.
/// 
///#### Basic usage:
///#### More detailed instructions on this twiki:
///#### https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JESUncertaintyProvider
///
/// 1) Include as a header in a C++ based analysis
/// 2) Create an instance, ie:
///       JESUncertaintyProvider myJES;
/// 3) Initialize it:
///       myJES.init()  
/// 3) Call one of the functions
///       myJES.getRelUncert(myPt, myEta);
/// 4) By default, one gets uncertainties without 
///    the pile-up term. That can be called separately, 
///    and the number of vertices needs to be specified as an additional argument
///
///#### Using separate terms
///
/// One can use the Components enum to access different components of the
/// JES or different combinations.
///
/// 1) Create the list of components
///     JESUncertaintyProvider::Components myComps =
///        JESUncertaintyProvider::Components(
///        JESUncertaintyProvider::CALORIMETER +
///        JESUncertaintyProvider::NOISETHRESHOLDS )
/// 2) Pass the list as a last parameter:
///     myJES.getRelUncert(myPt, myEta, myComps);
/// 3) One can also used pre-defined combinations
///     myJES.getRelUncert(myPt, myEta, 
///        JESUncertaintyProvider::NOPILEUP)
///
///#### Obtaining the 2D plots
///
/// Calling getUncGraphCopy will return a copy of a plot with all
/// uncertainties in pT, Eta bins. One can pass the components requested
/// as an argument.
///
///***********************************************************************
#ifndef _MYJESUNCERTAINTYPROVIDER_
#define _MYJESUNCERTAINTYPROVIDER_

#include "TNamed.h"
#include "TFile.h"
#include "TH2D.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
using std::cout;
using std::map;

namespace AnalysisFramework
{
namespace External
{

class JESUncertaintyProvider : public TNamed
{

public:

  enum Components { 
    // Basic components
    CALORIMETER = 1,                  // Calorimeter uncertainty from single particle propagation
    NOISETHRESHOLDS = 2,              // Use topocluster noise thresholds from data
    PERUGIATUNE = 4,                  // Perugia 2010 Pythia tune 
    ALPGENHERWIGJIMMY = 8,            // Alpgen+Herwig+Jimmy
    ETAINTERCALIBRATION = 16,         // Uncertainty due to intercalibration (endcap wrt central region)
    CLOSURE = 32,                     // Non-closure of numerical inversion constants
    PILEUP = 64,                      // Uncertainty due to in-time pile-up
    
    // Suggested combinations 
    NOCLOSURE = 95,                   // Allow users to include their own non-closure
    NOPILEUP = 63,                    // For data with no pile-up
    NOPILEUPNOCLOSURE = 31,           // Data with no pile-up, no non-closure
    ALL = 127                         // Everything
  }; 
  

  // Constructors, destructor
  JESUncertaintyProvider(std::string CollectionName="AntiKt6TopoJetsEM", std::string FileName="JESUncertainty.root");
   ~JESUncertaintyProvider();
  
  // Read the uncertainties off of the graph - By default, no Pile-Up term!
  double getRelUncert(double pT, double Eta, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);
  double getAbsUncert(double pT, double Eta, Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);
  
  // Pass a pointer to the correct graph
  //TH2D* getUncGraphCopy(Components UncertComps = JESUncertaintyProvider::NOPILEUP, unsigned int nVtx=1);

  // Initialize the provider
  virtual void init();

  // Open input file
  virtual TFile* openInputFile(std::string FileName = "JESUncertainty.root");

	//#ifdef JES_STANDALONE
	ClassDef(JESUncertaintyProvider,1);
	//#endif //JES_STANDALONE

 protected:

  // Change the inputCollection
  bool setInputCollection(std::string CollectionName);
    
  // Maximum number of components
  static const unsigned int m_nUncertainties_EMJES = 7;
  static const unsigned int m_nUncertainties_inSitu = 3;
  
  // Maximum number of vertices available
  static const unsigned int m_nVertices = 7;
  
  // For standalone compatibility
  const double m_GeV;
  
  // Pointers to the uncertainty histograms 
  map<int,TH2*> m_uncGraph;

  // Pointers to the pileup histograms
  map<int,TH2*> m_pileupUncGraph;

  // Input File
  TFile* m_inputFile; //! transient, do not save

  // Helper function to check requested components
  double getComponents(int currentBin, Components UncertComps, unsigned int nVtx = 1);

  // Name of jet collection
  std::string m_collectionName;

  // Name of file
  std::string m_fileName; 
  
  // Initialization flag
  bool m_isInit;
  
  // Deal with in-situ uncertainties (temporary)
  bool m_doInSitu;

  
};

}
}

#endif
