#include "JESUncertaintyProvider.h"


namespace AnalysisFramework
{
namespace External
{

// Constructor
JESUncertaintyProvider::JESUncertaintyProvider(std::string CollectionName, std::string FileName):
  m_GeV(1000.0), m_collectionName(CollectionName), m_fileName(FileName), m_isInit(false), m_doInSitu(false)
{

  m_inputFile = 0;
  m_uncGraph.clear();
  m_pileupUncGraph.clear();

}


// Destructor
JESUncertaintyProvider::~JESUncertaintyProvider()
{

  // delete the histograms (we own them now)
  map<int, TH2* >::iterator coll = m_uncGraph.begin();
  for(;coll!=m_uncGraph.end();coll++) { 
    //harmless: no pile-up for in-situ
    //delete also the pileup histo
    //if(coll->first==6) continue; //last unc pointer is a to pileup histo
    delete coll->second;
  }
  //only non in-situ has pileup (???? still needed ?????)
  if (!m_doInSitu) {
    coll = m_pileupUncGraph.begin(); 
    for(;coll!=m_pileupUncGraph.end();coll++) {
      delete coll->second;
    }
  }

}

void JESUncertaintyProvider::init() 
{
  
  //prevent multiple initializations
  if (m_isInit == false) {

    m_inputFile = openInputFile(m_fileName);

    if(!m_inputFile)
    {
      Error( "JESUncertaintyProvider::init()", "ERROR: Input File %s could not be found", m_fileName.c_str());
    } else {
      //the flag will be set as initialized if everything goes right
      m_isInit = setInputCollection(m_collectionName);

      m_inputFile->Close();
      //is someone trying to delete this afterwards?
      //delete m_inputFile;
    }
    
    Info( "JESUncertaintyProvider::init()", "===================================" );
    Info( "JESUncertaintyProvider::init()", "Initializing the JESUncertaintyProvider tool");
    Info( "JESUncertaintyProvider::init()", "Uncertainty for %s jets", m_collectionName.c_str() );
    if (!m_doInSitu) 
      Info( "JESUncertaintyProvider::init()", "Using uncertainty in ATLAS-CONF-2011-032" );
    else {
      Info( "JESUncertaintyProvider::init()", "Using combination of in-situ techniques" );
      Info( "JESUncertaintyProvider::init()", "Separate JES uncertainty components not available" );
      Info( "JESUncertaintyProvider::init()", "--> total uncertainty will always be returned even if components requested" );
    }
    Info( "JESUncertaintyProvider::init()", "For data produced with rel 16 (full 2010 dataset and MC)");
    Info( "JESUncertaintyProvider::init()", "===================================");

  }
  
  else {
    Warning( "JESUncertaintyProvider::init()", "WARNING: JESUncertaintyProvider already initialized, skipping re-initialization");
  }

}

// Open the file 
TFile* JESUncertaintyProvider::openInputFile(std::string FileName) 
{
  // Open Input File
  // The ROOT file containing the uncertainty plots must be placed in the current directory for the standalone version
  return new TFile(FileName.c_str()); //this is to avoid the file pointer to be invalid when it gets out of the function, TFile::Open wouldn't work - root magic... 
}


// Read the plots from the chosen Jet Collection
bool JESUncertaintyProvider::setInputCollection(std::string CollectionName)
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
    suffixName = "_AntiKt6Topo_EMJES";
  }

  if (suffixName.find("EMJES") != std::string::npos) {
    
    std::string plotNames[m_nUncertainties_EMJES] = {"Calorimeter","NoiseThresholds","Perugia2010", "AlpgenJimmy","EtaIntercalibration", "NonClosure", "Pileup"};
    
    // Pull the correct uncertainty graphs
    for(unsigned int i=0; i<m_nUncertainties_EMJES; i++)
    {

      std::string currentPlot = "";

      // Pileup plot (=6) is handled separately
      if (i!=6)  
      {
	// Combine names to get the right plots - CINT compatibility issue
	currentPlot = plotNames[i]+suffixName;
	m_inputFile->GetObject(currentPlot.c_str(),m_uncGraph[i]);
	m_uncGraph[i]->SetName(TString("the").Append(currentPlot.c_str()));
	m_uncGraph[i]->SetDirectory(0);
      }

      else 
      {
	
	//2010: Default pile-up plot: no additional vertices (=empty)
	//currentPlot = "Pileup"+suffixName+"_NPV_1";
	
	//2011: Default pile-up plot
	currentPlot = "Pileup"+suffixName;
	m_inputFile->GetObject(currentPlot.c_str(),m_uncGraph[i]);
	m_uncGraph[i]->SetName(TString("the").Append(currentPlot.c_str()));
	m_uncGraph[i]->SetDirectory(0);
	
      }

      if(!m_uncGraph[i]) 
      {
	Error("JESUncertaintyProvider::SetInputCollection()", "ERROR: Problem finding Required Input Graph: %s", currentPlot.c_str());
	return false;
      }

    }
    /*
    // Pull the correct vertex-dependent pileup graph
    for (unsigned int i=0; i<m_nVertices; i++) 
    {

      // Turn the index into a string for the number of vertices 
      std::ostringstream osstream;
      osstream << i+1;
      std::string string_i = osstream.str();

      // Combine names to get the right plots - CINT compatibility issue
      std::string currentPlot = "Pileup"+suffixName+"_NPV_"+string_i;

      m_inputFile->GetObject(currentPlot.c_str(),m_pileupUncGraph[i]);
      m_pileupUncGraph[i]->SetName(TString("the").Append(currentPlot.c_str()));
      m_pileupUncGraph[i]->SetDirectory(0);

      if(!m_pileupUncGraph[i]) 
      {
	Error("JESUncertaintyProvider::SetInputCollection()", "ERROR: Problem finding Required Input Graph: %s", currentPlot.c_str());
	return false;
      }

    }
    */
  }//end if on EM+JES

  if (suffixName.find("LC") != std::string::npos || suffixName.find("GCW") != std::string::npos) {
    
    //set flag for later use
    m_doInSitu = true;
    
    std::string plotNames[m_nUncertainties_inSitu] = {"InSitu","EtaIntercalibration", "Pileup"};

    // Pull the correct uncertainty graphs
    for(unsigned int i=0; i<m_nUncertainties_inSitu; i++)
    {
      std::string currentPlot = "";
      // Pileup plot (=3) is handled separately
      if (i!=2) {
	// Combine names to get the right plots - CINT compatibility issue
	currentPlot = plotNames[i]+suffixName;
	
	m_inputFile->GetObject(currentPlot.c_str(),m_uncGraph[i]);
	
	m_uncGraph[i]->SetName(TString("the").Append(currentPlot.c_str()));
	m_uncGraph[i]->SetDirectory(0);
      }
      
      else {
        
	//2011: Default pile-up plot
	currentPlot = "Pileup"+suffixName;
	m_inputFile->GetObject(currentPlot.c_str(),m_uncGraph[i]);
	m_uncGraph[i]->SetName(TString("the").Append(currentPlot.c_str()));
	m_uncGraph[i]->SetDirectory(0);
	
      }
      
      if(!m_uncGraph[i]) 
      {
	Error("JESUncertaintyProvider::SetInputCollection()", "ERROR: Problem finding Required Input Graph: %s", currentPlot.c_str());
	return false;
      }

    }

  }//end if on LC

  return true;

}

// Absolute Uncertainty
double JESUncertaintyProvider::getAbsUncert(double pT, double Eta, Components UncertComps, unsigned int nVtx)
{
  return getRelUncert(pT, Eta, UncertComps, nVtx)*pT;
}

// Relative Uncertainty
double JESUncertaintyProvider::getRelUncert(double pT, double Eta, Components UncertComps, unsigned int nVtx)
{
  // Convert units
  pT /= m_GeV;

  // Protect against jets in the wrong range
  if(pT <= 15 || pT >= 7000)
  {
    Warning("JESUncertaintyProvider::getRelUncert()", "pT outside of covered range (15-7000): Returning -1");
    return -1;
  } 

  // Protect against jets in the wrong region
  if(fabs(Eta) >= 4.5)
  {
    Warning("JESUncertaintyProvider::getRelUncert()", "Eta outside of covered range (0.0<|eta|<4.5): Returning -1");
    return -1;
  }

  // Use the last filled value in the histogram
  if(pT > 2500)
    pT = 2400.;

  // Find the bin with the given pT, Eta value
  
  // Keep current code for EM+JES (temporary)
  if (!m_doInSitu) {
    int currentBin = m_uncGraph[0]->FindBin(pT, fabs(Eta));
    // add uncertainties in the proper way
    return getComponents(currentBin, UncertComps, nVtx);
  }
  
  //for the moment, the components/NVTX arguments are disregarded. Should think of a better idea.  
  
  else {
    
    // here hardwired: in-situ is the first component, intercalibration is the second component
    // (not that it would matter anyways...but just for naming purposes)
    // they have a different binning, so need to be taken into account separately
   
    int currentBinInSitu = m_uncGraph[0]->FindBin(pT, fabs(Eta));
    int currentBinIntercalib = m_uncGraph[1]->FindBin(pT, fabs(Eta));
    int currentBinPileup = m_uncGraph[2]->FindBin(pT, fabs(Eta));//the bins are the same but still...
    
    double inSituUnc= m_uncGraph[0]->GetBinContent(currentBinInSitu);
    double intercalibrationUnc= m_uncGraph[1]->GetBinContent(currentBinIntercalib);
    //add the 2011 as well - code duplication but just temporary dirty fix for a couple of days anyways...
    double pileupUnc = m_uncGraph[2]->GetBinContent(currentBinPileup);
        
    // add uncertainties in quadrature
    return sqrt(inSituUnc*inSituUnc + intercalibrationUnc*intercalibrationUnc + pileupUnc*pileupUnc);
    
  }
  
}

// Get a copy of the 2D Graph containing the Uncertainties
// 15/05: Function temporarily commmented out due to incompatibility with different binnings 
// This function will be substituted by the caching mechanism anyways

/*TH2D* JESUncertaintyProvider::getUncGraphCopy(Components UncertComps, unsigned int nVtx)
{
  TH2D* myPlot = (TH2D*)m_uncGraph[0]->Clone();
  myPlot->Reset("ICE");
  int nBinsX = myPlot->GetNbinsX();
  int nBinsY = myPlot->GetNbinsY();

  for(int g=0; g<nBinsX; g++)
  {
    for(int h=0; h<nBinsY; h++)
    {
      // Follow the math highlighted in the TH2D class
      int currentBin = (g+1)+(nBinsX+2)*(h+1);
      myPlot->SetBinContent(currentBin, getComponents(currentBin, UncertComps, nVtx));
    }
  }
  return myPlot;
}*/

// Construct the uncertainty from a set of components
double JESUncertaintyProvider::getComponents(int currentBin, Components UncertComps, unsigned int nVtx)
{
  // Terms which will be added in quadrature
  double quadrature(0);

  // Initialize the bitmask 
  // Use int multiplication to avoid potential machine 
  // precision mismatches when using the pow function
  int bitmask = 1;

  // Take care of picking up the right pileup contribution first:  

  // check we have enough vertices stored (return 7 if there are more vertices in the event)
  if (nVtx > m_nVertices) nVtx=7;
  // protection against zero-vertex events (thanks Dag!)
  if (nVtx == 0) nVtx=1;

  // 2010: substitute the pile-up plot on the fly
  //m_uncGraph[6] = m_pileupUncGraph[nVtx-1];///CHANGE HERE WHEN NUMBER OF COMPONENTS CHANGES
  // 2011: leave generic pile-up plot, disregard vertex argument (a bit of a hack...)

  //Loop on the uncertainties
  for(unsigned int i=0; i<m_nUncertainties_EMJES; i++)

  {
    double currentComponent = 0;

    // Check if current component is requested
    if(int(UncertComps) & bitmask)
      currentComponent = m_uncGraph[i]->GetBinContent(currentBin);

    // Now all uncertainties are added in quadrature
    quadrature += currentComponent*currentComponent;

    // Prepare the bitmask for the next iteration
    bitmask *=2;
  }

  return sqrt(quadrature); 

}

}
}
