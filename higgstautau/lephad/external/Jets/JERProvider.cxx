
///***************************************************
///
/// Class: JERProvider
/// Author: Gaston Romeo <glromeo@cern.ch>
/// Plots by: Gaston Romeo <glromeo@cern.ch>
///
/// Provide the JER and its uncertainty.
/// Created: Jan/19/2011
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
/// 3) Call myJER.getRelResolutionMC(pt,y) to get the resolution in Monte Carlo (pt in GeV)
/// 4) Call myJER.getRelResolutionData(pt,y) to get the resolution in Monte Carlo (pt in GeV)
/// 5) Call myJER.getResolutionUncert(pt,y) to get its uncertainty (pt in GeV)
///
///
///***************************************************

#include "JERProvider.h"

namespace AnalysisFramework
{
namespace External
{

JERProvider::JERProvider(std::string CollectionName, std::string MethodName, std::string FileName):
	m_collectionName(CollectionName), m_methodName(MethodName), m_fileName(FileName), m_isInit(false)
{
	//Nothing needed here
}

void JERProvider::setFileName(std::string fileName)
{

  if(!m_isInit) m_fileName = fileName;
  else cout << "ERROR: Input file cannot be changed once JERProvider is initialized!" << endl;
}

void JERProvider::init()
{


  // Open Input File
  if(!m_isInit) {
  	m_inputFile = new TFile(m_fileName.c_str());
  	if(!m_inputFile)
    	{
    		cout << "ERROR: Input File " << m_fileName << " could not be found." << endl;
    	} else {
    		setInputCollection(m_collectionName, m_methodName);
    		cout << "JERProvider initialized:  Collection name = " << m_collectionName.c_str() << endl;

    		m_isInit = true;
  	}
  	// Close the file
  	m_inputFile->Close();
  	delete m_inputFile;
  }
  else {
  	cout << "JERProvider already initialized!" << endl;
  }
}

JERProvider::~JERProvider()
{
  // delete the functions (we own them now)
  map<int, TF1* >::iterator func = m_jerFunc.begin();
  for(;func!=m_jerFunc.end();func++) {
    delete func->second;
  }
  // delete the functions (we own them now)
  map<int, TGraphErrors* >::iterator off = m_jerOffset.begin();
  for(;off!=m_jerOffset.end();off++) {
    delete off->second;
  }
  // delete the functions (we own them now)
  map<int, TGraphErrors* >::iterator syst = m_jerSyst.begin();
  for(;syst!=m_jerSyst.end();syst++) {
    delete syst->second;
  }
}

void JERProvider::setInputCollection(std::string CollectionName, std::string MethodName)
{

  std::string suffixName;

  if(CollectionName == "AntiKt6TopoJES")
    suffixName = "_AntiKt6TopoJES";
  else if(CollectionName == "AntiKt4TopoJES")
    suffixName = "_AntiKt4TopoJES";
  else if(CollectionName == "AntiKt6LCTopoJES")
    suffixName = "_AntiKt6LCTopoJES";
  else if(CollectionName == "AntiKt4LCTopoJES")
    suffixName = "_AntiKt4LCTopoJES";
  else
    {
      cout << "ERROR: " << CollectionName << " not implemented, using default AntiKt6TopoJES"
	   << endl;
      suffixName = "_AntiKt6TopoJES";
    }


  if(MethodName == "Truth")
    suffixName = MethodName+suffixName;
  /* else if (MethodName == "DijetBalance")  */
  /*   suffixName = MethodName+suffixName; */
  /* else if (MethodName == "Bisector")  */
  /*   suffixName = MethodName+suffixName; */
  else {
    cout << "ERROR: " << MethodName << " not implemented, using default Truth"
	 << endl;
    suffixName = "Truth_AntiKt6TopoJES";
  }

  std::string regions[m_nY] = {"_0","_1","_2","_3","_4","_5"};

   // Pull the correct graphs
  for(int i=0; i < m_nY; i++)
    {

      std::string currentPlot = suffixName+regions[i];

      m_inputFile->GetObject(currentPlot.c_str(),m_jerFunc[i]);
      if(!m_jerFunc[i])
	cout << " ERROR: Problem finding Required Input Graph: " <<  currentPlot.c_str() << endl;
      else {
	m_jerFunc[i]->SetName(TString("the").Append(currentPlot.c_str()));
      }
    }

  // Pull the correct offset from data/mc comparsion
  // (only up to 2.8 due to statistics) ---> m_nY-2
  for(int i=0; i < m_nY - 2 ; i++)
    {

      std::string currentPlot = "DataMCBisector_"+CollectionName+regions[i];

      m_inputFile->GetObject(currentPlot.c_str(),m_jerOffset[i]);
      if(!m_jerOffset[i])
	cout << " ERROR: Problem finding Required Input Graph: " <<  currentPlot.c_str() << endl;
      else {
	m_jerOffset[i]->SetName(TString("the").Append(currentPlot.c_str()));
      }
    }

  // Pull the correct uncertainty from data/mc comparsion
  // (only up to 2.8 due to statistics) ---> m_nY-2
  for(int i=0; i < m_nY - 2 ; i++)
    {

      std::string currentPlot = "DataMCBisectorUNCERT_"+CollectionName+regions[i];

      m_inputFile->GetObject(currentPlot.c_str(),m_jerSyst[i]);
      if(!m_jerSyst[i])
	cout << " ERROR: Problem finding Required Input Graph: " <<  currentPlot.c_str() << endl;
      else {
	m_jerSyst[i]->SetName(TString("the").Append(currentPlot.c_str()));
      }
    }



}

TF1* JERProvider::getParam(double y)
{

  if(!m_isInit) {
    cout << "JERProvider not initialized." << endl;
    return 0;
  }

  y = fabs(y);

  // Return the requested parametrization

  int index = -1;
  if (y <= 0.8) {
    index = 0;
    if (m_collectionName == "AntiKt6TopoJES" || m_collectionName == "AntiKt6LCTopoJES") m_syst = 0.08;
    if (m_collectionName == "AntiKt4TopoJES" || m_collectionName == "AntiKt4LCTopoJES") m_syst = 0.11;
  }
  else if (y >= 0.8 && y < 1.2) {
    index = 1;
    if (m_collectionName == "AntiKt6TopoJES" || m_collectionName == "AntiKt6LCTopoJES") m_syst = 0.10;
    if (m_collectionName == "AntiKt4TopoJES" || m_collectionName == "AntiKt4LCTopoJES") m_syst = 0.13;
  }
  else if (y >= 1.2 && y < 2.1) {
    index = 2;
    if (m_collectionName == "AntiKt6TopoJES" || m_collectionName == "AntiKt6LCTopoJES") m_syst = 0.09;
    if (m_collectionName == "AntiKt4TopoJES" || m_collectionName == "AntiKt4LCTopoJES") m_syst = 0.12;
  }
  else if (y >= 2.1 && y < 2.8) {
    index = 3;
    if (m_collectionName == "AntiKt6TopoJES" || m_collectionName == "AntiKt6LCTopoJES") m_syst = 0.10;
    if (m_collectionName == "AntiKt4TopoJES" || m_collectionName == "AntiKt4LCTopoJES") m_syst = 0.13;
  }
  else if (y >= 2.8 && y < 3.6) {
    index = 4;
    if (m_collectionName == "AntiKt6TopoJES" || m_collectionName == "AntiKt6LCTopoJES") { m_offset = 0.00; m_uncert = 0.20; m_syst = 0.11;}
    if (m_collectionName == "AntiKt4TopoJES" || m_collectionName == "AntiKt4LCTopoJES") { m_offset = 0.00; m_uncert = 0.20; m_syst = 0.14;}
  }
  else if (y >= 3.6 && y <= 4.5) {
    index = 5;
    if (m_collectionName == "AntiKt6TopoJES" || m_collectionName == "AntiKt6LCTopoJES") { m_offset = 0.00; m_uncert = 0.30; m_syst = 0.11;}
    if (m_collectionName == "AntiKt4TopoJES" || m_collectionName == "AntiKt4LCTopoJES") { m_offset = 0.00; m_uncert = 0.30; m_syst = 0.14;}
  }
  else {
    //Protect against jets in the wrong region
    cout << "WARNING: Y outside of covered range (0.0-4.5): Returning parametrization for 4.5" << endl;
    return m_jerFunc[5];
  }

  return m_jerFunc[index];

}

TGraphErrors* JERProvider::getOffset(double y)
{

  if(!m_isInit) {
    cout << "JERProvider not initialized." << endl;
    return 0;
  }

  y = fabs(y);

  // Return the requested region

  int index = -1;
  if (y <= 0.8) index = 0;
  else if (y >= 0.8 && y < 1.2) index = 1;
  else if (y >= 1.2 && y < 2.1) index = 2;
  else if (y >= 2.1 && y < 2.8) index = 3;
  else if (y >= 2.8 && y < 3.6) index = 3;  // Return [3] due to lacking of stats from data/mc beyond 2.8
  else if (y >= 3.6 && y <= 4.5) index = 3; // Return [3] due to lacking of stats from data/mc beyond 2.8
  else {
    //Protect against jets in the wrong region
    cout << "WARNING: Y outside of covered range (0.0-4.5): Returning parametrization for 4.5" << endl;
    return m_jerOffset[3];
  }

  return m_jerOffset[index];

}

TGraphErrors* JERProvider::getSyst(double y)
{

  if(!m_isInit) {
    cout << "JERProvider not initialized." << endl;
    return 0;
  }

  y = fabs(y);

  // Return the requested region

  int index = -1;
  if (y <= 0.8) index = 0;
  else if (y >= 0.8 && y < 1.2) index = 1;
  else if (y >= 1.2 && y < 2.1) index = 2;
  else if (y >= 2.1 && y < 2.8) index = 3;
  else if (y >= 2.8 && y < 3.6) index = 3;  //Return [3] due to lacking of stats from data/mc
  else if (y >= 3.6 && y <= 4.5) index = 3; //Return [3] due to lacking of stats from data/mc
  else {
    //Protect against jets in the wrong region
    cout << "WARNING: Y outside of covered range (0.0-4.5): Returning parametrization for 4.5" << endl;
    return m_jerSyst[3];
  }

  return m_jerSyst[index];

}


float JERProvider::getOffset(double pt, double y)
{

  TGraphErrors* jer_offset = getOffset(y);
  double offset;

  if ( fabs(y) < 2.8) offset = jer_offset->Eval(pt) / 100.;
  else offset = 0.0;

  return offset;

}

float JERProvider::getRelResolutionMC(double pt, double y)
{

  if(!m_isInit) {
    cout << "JERProvider not initialized." << endl;
    return 0;
  }

  y = fabs(y);

  // Protect against jets in the wrong region
  if(fabs(y) > 4.5)
    {
      //cout << "WARNING: Y outside of covered range (0.0-4.5): Y set to 4.5" << endl;
      y = 4.5;
    }

  // Protect against jets in the wrong range
  if(pt < 10.)
    {
      //cout << "WARNING: pt outside of covered range (10-2000): Pt set to 10 GeV" << endl;
      pt = 10.;
    }

  if(pt > 2000.)
    {
      //cout << "WARNING: pt outside of covered range (10-2000): Pt set to 5000 GeV" << endl;
      pt = 2000.;
    }


  // Get JER
  TF1* jer = getParam(y);

  //  Noise, Stochastic, Constant terms and their errors ---> [N,S,C,Ne,Se,Ce]
  m_param[0] = jer->GetParameter(0);   m_param[3] = jer->GetParError(0);
  m_param[1] = jer->GetParameter(1);   m_param[4] = jer->GetParError(1);
  m_param[2] = jer->GetParameter(2);   m_param[5] = jer->GetParError(2);

  // cout << "Parametrization: sqrt([0]*[0]/(x*x) + [1]*[1]/x + [2]*[2])" << endl;
  // for (int i = 0; i < m_nParam/2; i++) { cout << "Parameter " << i << ":\t" << m_param[i] <<"\t +/- \t" << m_param[i+3]<< endl;}

  double sigma_MC = sqrt(m_param[0]*m_param[0]/pt/pt + m_param[1]*m_param[1]/pt + m_param[2]*m_param[2]);

  return sigma_MC;

}

float JERProvider::getRelResolutionData(double pt, double y) {

  double sigma_MC = getRelResolutionMC(pt, y);
  double offset = getOffset(pt, y);

  double sigma_data = sigma_MC + offset*sigma_MC;

  return sigma_data;

}

float JERProvider::getResolutionUncert(double pt, double y)
{

	if(!m_isInit) {
		cout << "JERProvider not initialized." << endl;
		return 0;
	}

  y = fabs(y);

  // Protect against jets in the wrong region
  if(fabs(y) > 4.5)
    {
      //cout << "WARNING: Y outside of covered range (0.0-4.5): Y set to 4.5" << endl;
      y = 4.5;
    }

  // Protect against jets in the wrong range
  if(pt < 10.)
    {
      //cout << "WARNING: pt outside of covered range (10-2000): Pt set to 10 GeV" << endl;
      pt = 10.;
    }

  if(pt > 2000.)
    {
      //cout << "WARNING: pt outside of covered range (10-2000): Pt set to 2000 GeV" << endl;
      pt = 2000.;
    }


  TGraphErrors* jer_syst = getSyst(y);

  // Setting boundaries due to lack of statistics
  float ptmin = 10; float ptmax = 1000; float ptpiv = 1000;

  if ( y > 0.8 && y <= 1.2 ) { ptmin = 10; ptmax = 800; ptpiv = 800; }
  if ( y > 1.2 && y <= 2.1 ) { ptmin = 10; ptmax = 600; ptpiv = 600; }
  if ( y > 2.1 && y <= 2.8 ) { ptmin = 10; ptmax = 500; ptpiv = 500; }
  if ( y > 2.8 ) { ptmin = 10; ptmax = 500; ptpiv = 500;  }


  float factor = 0;

  if ( pt >= ptmin && pt <= ptmax) {
    if ( fabs(y) < 2.8 ) m_uncert = jer_syst->Eval(pt) / 100.;
    m_uncert = sqrt( m_uncert*m_uncert + m_syst*m_syst );
    m_uncert = m_uncert*getRelResolutionData(pt,y);
  }
  if ( pt > ptmax && pt <= 2*ptmax) {
    factor = (pt - ptmax)/ptpiv;
    if ( fabs(y) < 2.8 ) m_uncert = jer_syst->Eval(ptmax) / 100.;
    m_uncert = sqrt( m_uncert*m_uncert + m_syst*m_syst );
    m_uncert = (1+factor)*m_uncert*getRelResolutionData(pt,y); // Above validated range (no data available), we double the uncertainty at ptpiv to be conservative.
  }
  if ( pt > 2*ptmax) {
    if ( fabs(y) < 2.8 ) m_uncert = jer_syst->Eval(ptmax) / 100.;
    m_uncert = sqrt( m_uncert*m_uncert + m_syst*m_syst );
    m_uncert = 2.0*m_uncert*getRelResolutionData(pt,y);
  }

  return m_uncert;

}

float JERProvider::getSmearingFactorMC(double pt, double y) {

  //here we want to consider the intrinsic resolution of the MC and subtract it out
  //TMath::Sqrt( S*S - S_MC*S_MC );
  double S = getRelResolutionData(pt,y);
  double S_MC = getRelResolutionMC(pt,y);
  double smearingFactorMC = S > S_MC ? sqrt(S*S-S_MC*S_MC) : 0.0;
  return smearingFactorMC;

}

float JERProvider::getSmearingFactorSyst(double pt, double y) {

  //here we want to keep the MC resolution as a starting point, and oversmear it up to the worst possible value of the data resoln+uncertainty
  //TMath::Sqrt( (S+U)*(S+U) - (S*S) );
  double smearingFactorSyst = sqrt(pow(getRelResolutionData(pt,y)+getResolutionUncert(pt,y),2) - pow(getRelResolutionData(pt,y), 2));
  return smearingFactorSyst;

}

}
}
