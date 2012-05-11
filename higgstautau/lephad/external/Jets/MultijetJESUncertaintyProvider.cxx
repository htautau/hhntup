#include "MultijetJESUncertaintyProvider.h"

namespace AnalysisFramework
{
namespace External
{


int MultijetJESUncertaintyProvider::m_counter=0;

// Constructor
MultijetJESUncertaintyProvider::MultijetJESUncertaintyProvider(std::string CollectionName, std::string AnalysisFileName, std::string MJESFileName, std::string FileName):
JESUncertaintyProvider(CollectionName, FileName), m_analysisFileName(AnalysisFileName), m_MJESFileName(MJESFileName)
{
  // Include Flavor Composition Uncertainty by default
  m_includeFlavorComp = true;
  
  m_flavorCompGluGraph = NULL; 
  m_flavorCompLightGraph = NULL; 
  m_deltaRGraph = NULL; 
  m_gluonFraction = NULL; 
  m_gluonFractionError = NULL; 
  m_responseSample = NULL; 
}

// Destructor
MultijetJESUncertaintyProvider::~MultijetJESUncertaintyProvider()
{

	//delete histograms, we own them
	if(m_flavorCompGluGraph) delete m_flavorCompGluGraph;
	if(m_flavorCompLightGraph) delete m_flavorCompLightGraph;  
	if(m_deltaRGraph) delete m_deltaRGraph;
	if(m_gluonFraction) delete m_gluonFraction;
	if(m_gluonFractionError) delete m_gluonFractionError;  
	if(m_responseSample) delete m_responseSample;  

}



void MultijetJESUncertaintyProvider::init()
{

  // Prevent multiple initializations 
  if (m_isInit == false)
    {

    m_inputFile = 0;
    m_MJESInputFile = 0;  
    m_analysisInputFile = 0;

    // Uncertainty Provider input file
    //JESUncertaintyProvider::openInputFile(m_fileName);
    m_inputFile = openInputFile(m_fileName);
    // MJES input file
    m_MJESInputFile = openInputFile(m_MJESFileName);  
    // Analysis input file
    m_analysisInputFile = openInputFile(m_analysisFileName);

    if(!m_inputFile)      
      {
	Error( "MultijetJESUncertaintyProvider::init()", "ERROR: Input File %s could not be found in the current directory", m_fileName.c_str());
      }
     if(!m_MJESInputFile)      
      {
	Error( "MultijetJESUncertaintyProvider::init()", "ERROR: Input File %s could not be found in the current directory", m_MJESFileName.c_str());
      }     
    if(!m_analysisInputFile)
      {
	Error( "MultijetJESUncertaintyProvider::init()", "ERROR: Input File %s could not be found in the current directory", m_analysisFileName.c_str());
      }
    else 
      {
      // The flag will be set as initialized if everything goes right
      m_isInit = setInputCollection(m_collectionName);
      }

	// Close the files
	m_analysisInputFile->Close();
	delete m_analysisInputFile; 
	
	m_MJESInputFile->Close();
	delete m_MJESInputFile; 	
	
	m_inputFile->Close();
	delete m_inputFile;
	
    Info( "MultijetJESUncertaintyProvider::init()", "===================================" );
    Info( "MultijetJESUncertaintyProvider::init()", "Initializing the MultijetJESUncertaintyProvider tool");
    Info( "MultijetJESUncertaintyProvider::init()", ("Uncertainty for "+m_collectionName).c_str() );
    if (!m_doInSitu) {
      Info( "MultijetJESUncertaintyProvider::init()", "Using uncertainty in ATLAS-CONF-2011-032 from JESUncertaintyProvider tool" );
      Info( "MultijetJESUncertaintyProvider::init()", "Updated pile-up uncertainty from 2011 data" );
    }
    else {
      Info( "MultijetJESUncertaintyProvider::init()", "Using combination of in-situ techniques from JESUncertaintyProvider tool" );
      Info( "MultijetJESUncertaintyProvider::init()", "Separate JES uncertainty components not available" );
      Info( "MultijetJESUncertaintyProvider::init()", "Pile-up uncertainty not available" );
    }
    Info( "MultijetJESUncertaintyProvider::init()", "For data produced with rel 16.6 (full 2011 dataset and MC)");
    Info( "MultijetJESUncertaintyProvider::init()", "(baseline uncertainty breakdown not available for 2011)");
    Info( "MultijetJESUncertaintyProvider::init()", "===================================");
    }
  else {
    Warning( "MultijetJESUncertaintyProvider::init()", "WARNING: MultijetJESUncertaintyProvider already initialized, skipping re-initialization");
  }
  
}



// Open the analysis file 
//void MultijetJESUncertaintyProvider::openInputFile(std::string FileName, TFile analysisInputFilePtr)
//{
  // Open Input File
  // The ROOT file containing the uncertainty plots must be placed in the current directory for the standalone version
//  analysisInputFilePtr = new TFile(FileName.c_str());
//}


// Read the plots from the chosen Jet Collection
bool MultijetJESUncertaintyProvider::setInputCollection(std::string CollectionName)
{


  // Jet Collection used
  std::string suffixName;

  if(CollectionName == "AntiKt6EMJESTopoJets" || CollectionName == "AntiKt6TopoJetsEM" )
    suffixName = "_AntiKt6Topo_EMJES";
  else if(CollectionName == "AntiKt4EMJESTopoJets" || CollectionName == "AntiKt4TopoJetsEM") 
    suffixName = "_AntiKt4Topo_EMJES"; 
  else if(CollectionName == "AntiKt4LCTopoJets") 
    suffixName = "_AntiKt4Topo_LCJES";
  else if(CollectionName == "AntiKt6LCTopoJets") 
    suffixName = "_AntiKt6Topo_LCJES";
  else if(CollectionName == "AntiKt4GCWTopoJets" || CollectionName == "AntiKt4TopoGCWJets") 
    suffixName = "_AntiKt4Topo_GCWJES";
  else if(CollectionName == "AntiKt6GCWTopoJets" || CollectionName == "AntiKt6TopoGCWJets") 
    suffixName = "_AntiKt6Topo_GCWJES";

  else
  {
    Warning("JESUncertaintyProvider::setInputCollection()", "ERROR: Name not recognized, using default AntiKt6EMJESTopoJets");
    m_collectionName = "AntiKt6TopoJetsEM";
    suffixName = "_AntiKt6Topo_EMJES";
  }

  //set JESUncertaintyProvider input collection...and continue with multijet plots if everything goes right
  if(!(JESUncertaintyProvider::setInputCollection(CollectionName))) return false;

  // Pull the correct flavor composition graph

  // Gluon flavor composition graph
  m_MJESInputFile->GetObject(("flavorCompGlu"+suffixName).c_str(),m_flavorCompGluGraph);

  if(!m_flavorCompGluGraph) 
    {
	  Error("MultijetJESUncertaintyProvider::SetInputCollection()", "ERROR: Problem finding Required Input Graph: flavorCompGlu%s",suffixName.c_str());
          return false;
    }
    
  m_flavorCompGluGraph->SetName(TString("the").Append("flavorCompGlu").Append(suffixName));
  m_flavorCompGluGraph->SetDirectory(0);

  // Light quark flavor composition graph
  m_MJESInputFile->GetObject(("flavorCompLight"+suffixName).c_str(),m_flavorCompLightGraph);
  
  if(!m_flavorCompLightGraph)
    {
	  Error("MultijetJESUncertaintyProvider::SetInputCollection()","ERROR: Problem finding Required Input Graph: flavorCompLight%s",suffixName.c_str());
          return false;
    }

  m_flavorCompLightGraph->SetName(TString("the").Append("flavorCompLight").Append(suffixName));
  m_flavorCompLightGraph->SetDirectory(0);

  // Pull the correct deltaR graph
  m_MJESInputFile->GetObject(("deltaR_new"+suffixName).c_str(),m_deltaRGraph);

  if(!m_deltaRGraph) 
    {
	  Error("MultijetJESUncertaintyProvider::SetInputCollection()","ERROR: Problem finding Required Input Graph: deltaR_new%s",suffixName.c_str());
          return false;
    }
    
  m_deltaRGraph->SetName(TString("the").Append("deltaR_new").Append(suffixName));
  m_deltaRGraph->SetDirectory(0);


  // Analysis sample graphs

  // Pull the correct gluon fraction graph
  m_analysisInputFile->GetObject(("gluonFraction"+suffixName).c_str(),m_gluonFraction);

  if(!m_gluonFraction && m_includeFlavorComp)
    {
	  Warning("MultijetJESUncertaintyProvider::SetInputCollection()","Problem finding Required Input Graph gluonFraction%s in input file %s", suffixName.c_str(), m_analysisFileName.c_str());

    }


  // Pull the correct gluon fraction error graph
  m_analysisInputFile->GetObject(("gluonFractionError"+suffixName).c_str(),m_gluonFractionError);
  
  if(!m_gluonFractionError && m_includeFlavorComp)
    {
	  Warning("MultijetJESUncertaintyProvider::SetInputCollection()","Problem finding Required Input Graph gluonFractionError%s in input file %s", suffixName.c_str(), m_analysisFileName.c_str());

    }


  // If gluon fraction or gluon fraction error graph/s is/are missing -> use default gluonFraction and gluonFractionError graphs provided in input file "JESUncertainty.root" instead -> this will set gluon fraction = gluon fraction error = 0.5 resulting in alphaC=alphaR=1 when calling getAlphaC() and getAlphaR() later in the code
  if(!m_gluonFraction || !m_gluonFractionError)
    {

	  if(m_includeFlavorComp) Warning("MultijetJESUncertaintyProvider::SetInputCollection()","Not enough information to calculate prefactor alphaC: Setting alphaC = 1!");

  		m_MJESInputFile->GetObject(("gluonFractionDefault"+suffixName).c_str(),m_gluonFraction); 
  		m_MJESInputFile->GetObject(("gluonFractionErrorDefault"+suffixName).c_str(),m_gluonFractionError);

    }


  m_gluonFraction->SetName(TString("the").Append("gluonFraction").Append(suffixName));
  m_gluonFraction->SetDirectory(0);
  m_gluonFractionError->SetName(TString("the").Append("gluonFractionError").Append(suffixName));
  m_gluonFractionError->SetDirectory(0);

  // Pull the correct response graph
  m_analysisInputFile->GetObject(("rSample_rIncl"+suffixName).c_str(),m_responseSample);


  // If the response graph is not available pull the default graph provided in input file "JESUncertainty.root" -> this will set r_sample/r_incl = 1
  if(!m_responseSample) 
    {
      	if(m_includeFlavorComp) 
	  Warning("MultijetJESUncertaintyProvider::SetInputCollection()","Problem finding Required Input Graph rSample_rIncl%s in input file %s! Setting response_sample/r_incl = 1!", suffixName.c_str(), m_analysisFileName.c_str());

    // Use default response graph
  	m_MJESInputFile->GetObject(("rSample_rInclDefault"+suffixName).c_str(),m_responseSample);

    }

  m_responseSample->SetName(TString("the").Append("rSample_rIncl").Append(suffixName));
  m_responseSample->SetDirectory(0);

  return true;
}


// Include/exclude flavor composition term
void MultijetJESUncertaintyProvider::includeFlavorComposition(bool include)
{
  m_includeFlavorComp = include;
}


// Relative positive uncertainty
double MultijetJESUncertaintyProvider::getRelPosUncert(double pT, double Eta, double dRmin, Components UncertComps, unsigned int nVtx)
{
  return getRelUncert(pT, Eta, dRmin, true, UncertComps, nVtx);
}

// Absolute positive uncertainty 
double MultijetJESUncertaintyProvider::getAbsPosUncert(double pT, double Eta, double dRmin, Components UncertComps, unsigned int nVtx)
{
  return getRelPosUncert(pT, Eta, dRmin, UncertComps, nVtx)*pT;
}

// Relative negative uncertainty
double MultijetJESUncertaintyProvider::getRelNegUncert(double pT, double Eta, double dRmin, Components UncertComps, unsigned int nVtx)
{
  return getRelUncert(pT, Eta, dRmin, false, UncertComps, nVtx);
}

// Absolute negative uncertainty
double MultijetJESUncertaintyProvider::getAbsNegUncert(double pT, double Eta, double dRmin, Components UncertComps, unsigned int nVtx)
{
  return getRelNegUncert(pT, Eta, dRmin, UncertComps, nVtx)*pT;
}



// Relative Uncertainty (returns relative positive/negative uncertainty according to updown=true/false)
double MultijetJESUncertaintyProvider::getRelUncert(double pT, double Eta, double dRmin, bool isUp, Components UncertComps, unsigned int nVtx)
{


  // Convert units
  pT /= m_GeV;

  // Protect against jets in the wrong range

  if(pT <= 15 || pT >= 7000)
    {
    Warning("MultijetJESUncertaintyProvider::getRelUncert()", "pT outside of covered range (15-7000): Returning -1");
    return -1;
    }

  // Protect against jets in the wrong region
  if(fabs(Eta) >= 4.5)
    {
    Warning("MultijetJESUncertaintyProvider::getRelUncert()", "Eta outside of covered range (0.0<|eta|<4.5): Returning -1");
    return -1;
    }


  // Use the last filled value in the histogram
  if(pT > 2500)
    pT = 2400.;


  // Add uncertainties
  return getComponents(pT, Eta, dRmin, isUp, UncertComps, nVtx);


}


/*
TH2D* MultijetJESUncertaintyProvider::getUncGraphCopy(bool isUp, Components UncertComps, unsigned int nVtx)
{
  TH2D* myPlot = (TH2D*)m_uncGraph[0]->Clone();
  myPlot->Reset("ICE");
  
  
  // Set to default values to avoid compiler warnings
  isUp = true; 
  UncertComps=JESUncertaintyProvider::NOPILEUP; 
  nVtx=1;
  

    Warning("MultijetJESUncertaintyProvider::getUncGraphCopy()", "MultijetJESUncertaintyProvider::getUncGraphCopy() not implemented... Returning empty graph.");  

  
  return myPlot;
}
*/

// Check requested components (global uncertainty) and add up flavor composition --> moved close-by jets contribution to separate function
double MultijetJESUncertaintyProvider::getComponents(double pT, double Eta, double dRmin, bool isUp, Components UncertComps, unsigned int nVtx)
{

  // Uncertainty components to be added
  double globalUnc(0);
  double flavComp(0);
  double closeBy(0);

  // This works since all MJES graphs (EM+JES, LC, GCW) follow the same binning conventions as the EM+JES plots from the JESUncertaintyProvider
  int currentBin = m_flavorCompLightGraph->FindBin(pT, fabs(Eta));

  // Pick up uncertainty components from the JESUncertaintyProvider tool
  // EM+JES (using getComponents() as before)
  if (!m_doInSitu){
  globalUnc = JESUncertaintyProvider::getComponents(currentBin, UncertComps, nVtx);
  }
  // Pick up the LC uncertainty - make sure to re-convert units to MeV before calling JESUncertaintyProvider::getRelUncert()
  else{
  globalUnc = JESUncertaintyProvider::getRelUncert(pT*1000,Eta);
  }


  // Flavor composition uncertainty

  if(m_includeFlavorComp){
  
    // ratio r_sample/r_incl (response in dijet sample)
    double responseSample =  getResponseSample(currentBin);

    // Prefactor alphaC (->flavor composition uncertainty)
    double alphaC = getAlphaC(currentBin, isUp);
 
    // Choose the correct flavor composition graph (light quark flavor composition graph if positive uncertainty (isUp) is requested, else choose gluon composition graph)
    TH2D* flavorCompGraph = isUp ? m_flavorCompLightGraph  : m_flavorCompGluGraph;

    flavComp = alphaC*(fabs(flavorCompGraph->GetBinContent(currentBin)-responseSample))/responseSample;


  }


  // Additional contribution for non-isolated jets

  // Make sure we're not running out of the dR-range covered in the close-by histograms
  if(dRmin<1.5)
  {
  
    closeBy = m_deltaRGraph->Interpolate(pT, dRmin);
  }


  // Add the uncertainties in quadrature
  return sqrt(pow(globalUnc,2)+pow(flavComp,2)+pow(closeBy,2));


}




// Prefactor alphaC
double MultijetJESUncertaintyProvider::getAlphaC(int currentBin, bool isUp){

  // Get the correct light quark (isUp) or gluon (!isUp) fraction
  double fraction = isUp ? (1-m_gluonFraction->GetBinContent(currentBin)) : m_gluonFraction->GetBinContent(currentBin);
  // Get the correct fraction error
  double fracUnc = m_gluonFractionError->GetBinContent(currentBin);


  // Make sure (1-fraction) != 0 -> avoid division by zero by setting alphaC=1
  if((1-fraction)==0){
     Warning("MultijetJESUncertaintyProvider::getAlphaC()", "gluon fraction = 0/1: setting alphaC = 1 to avoid division by zero!");  
    return 1;
  }

  // Make sure fracUnc <= (1-fraction) (otherwise alphaC > 1!)
  if(fracUnc > (1-fraction))
    {
    fracUnc = (1-fraction);
    }

    return fracUnc/(1-fraction);


}


// Prefactor alphaR (->flavor response uncertainty)  -- currently, we don't need alphaR since the flavor response terms were dropped
double MultijetJESUncertaintyProvider::getAlphaR(int currentBin, bool isUp, double alphaC, double responseSample){


  double alphaR = alphaC;


  // Modify alphaR if average response >= 1 and positive error is requested (isUp)
  if(responseSample >=1 && isUp)
    {
    alphaR=alphaR+(1-alphaR)*(responseSample-1)/(m_flavorCompLightGraph->GetBinContent(currentBin));
    }

  // Modify alphaR if average response <= 1 and negative error is requested (!isUp)
  else if(responseSample <=1 && !isUp)
    {
    alphaR=alphaR+(1-alphaR)*(1-responseSample)/(m_flavorCompGluGraph->GetBinContent(currentBin));
    }	

  return alphaR;

}

// Pick up r_sample/r_incl from input histos and make sure the values are within a reasonable range
double MultijetJESUncertaintyProvider::getResponseSample(int currentBin){

  double responseSample =  m_responseSample->GetBinContent(currentBin);

  // Protect against response_sample/r_incl > response_lightQuark/r_incl (= Delta_pos)
  if(responseSample > m_flavorCompLightGraph->GetBinContent(currentBin))
    {
    if (m_counter < 10)
      {
      Warning("MultijetJESUncertaintyProvider::getResponseSample()", "analysis sample response larger than light quark jet response: setting response(sample)/r_incl = response(light quarks)/r_incl --> Check input graph responseSample in file %s !", m_analysisFileName.c_str()) ;
      }
    if (m_counter == 10)
      {
      Warning("MultijetJESUncertaintyProvider::getResponseSample()", ("analysis sample response larger than light quark jet response: setting response(sample)/r_incl = response(light quarks)/r_incl --> Check input graph responseSample in file " +  m_analysisFileName + "! (last message)").c_str()) ;
      }
    incrementCounter();
    responseSample =  m_flavorCompLightGraph->GetBinContent(currentBin);
    }

  // Protect against response_sample/r_incl < response_gluon/r_incl (= Delta_neg)
  else if(responseSample < m_flavorCompGluGraph->GetBinContent(currentBin))
    {
    if (m_counter < 10)
      {
      Warning("MultijetJESUncertaintyProvider::getResponseSample()", "analysis sample response larger than gluon jet response: setting response(sample)/r_incl = response(gluons)/r_incl --> Check input graph responseSample in file %s !", m_analysisFileName.c_str()) ;
      }
    if (m_counter == 10)
      {
      Warning("MultijetJESUncertaintyProvider::getResponseSample()", "analysis sample response larger than gluon jet response: setting response(sample)/r_incl = response(gluons)/r_incl --> Check input graph responseSample in file %s !", m_analysisFileName.c_str()) ;
      }
    incrementCounter();
    responseSample =  m_flavorCompGluGraph->GetBinContent(currentBin);
    }

  return responseSample;

}

}
}
