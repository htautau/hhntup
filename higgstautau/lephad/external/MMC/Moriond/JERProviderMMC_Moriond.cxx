
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

#include "JERProviderMMC_Moriond.h"

namespace AnalysisFramework
{
namespace External
{

JERProviderMMCMoriond::JERProviderMMCMoriond(std::string CollectionName, std::string MethodName, std::string FileName):
	m_collectionName(CollectionName), m_methodName(MethodName), m_fileName(FileName), m_isInit(false)
{
	//Nothing needed here
}

void JERProviderMMCMoriond::init()
{
  // Open Input File
  if(!m_isInit) {
  	m_inputFile = new TFile(m_fileName.c_str());
  	if(!m_inputFile)
    	{
    		cout << "ERROR: Input File " << m_fileName << " could not be found." << endl;
    	} else {
    		setInputCollection(m_collectionName, m_methodName);
    		cout << "JERProviderMMCMoriond initialized:  Collection name = " << m_collectionName.c_str() << endl;
    
    		m_isInit = true;
  	}
  	// Close the file
  	m_inputFile->Close();
  	delete m_inputFile;
  }
  else {
  	cout << "JERProviderMMCMoriond already initialized!" << endl;
  }
}

JERProviderMMCMoriond::~JERProviderMMCMoriond()
{
  // delete the functions (we own them now)
	map<int, TF1* >::iterator func = m_jerFunc.begin();
	for(;func!=m_jerFunc.end();func++) { 
			delete func->second;
	}
}

void JERProviderMMCMoriond::setInputCollection(std::string CollectionName, std::string MethodName)
{
  
  std::string suffixName;
  
  if(CollectionName == "AntiKt6TopoJES") 
    suffixName = "_AntiKt6TopoJES";
  else if(CollectionName == "AntiKt4TopoJES") 
    suffixName = "_AntiKt4TopoJES";
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
			cout << " ERROR: Problem finding Required Input Graph: " <<  suffixName+regions[i] << endl;
		else {
    		m_jerFunc[i]->SetName(TString("the").Append(currentPlot.c_str()));
			//m_jerFunc[i]->SetDirectory(0);	
		}
	}
}

TF1* JERProviderMMCMoriond::getParam(double y)
{

	if(!m_isInit) {
		cout << "JERProviderMMCMoriond not initialized." << endl;
		return 0;
	}

  y = fabs(y); 
  
  // Pull uncertainty and offset from Data/MC comparison

  //cout << "getParam() called " << endl;
  //cout << " name = " << m_collectionName.c_str() << endl;
  //cout << " done " << endl;
  if (m_collectionName == "AntiKt6TopoJES") {
    if (y <= 0.8) { m_offset = 0.04; m_uncert = 0.05; }
    else if (y >= 0.8 && y < 1.2) { m_offset = 0.09; m_uncert = 0.1; }
    else if (y >= 1.2 && y < 2.1) { m_offset = 0.05; m_uncert = 0.1; }
    else if (y >= 2.1 && y < 2.8) { m_offset = 0.00; m_uncert = 0.14; }
    else if (y >= 2.8) { m_offset = 0.00; m_uncert = 0.14; }
    else { 
      cout << "ERROR: " << y << " not implemented, using default central region" << endl;
      m_offset = 0.04; m_uncert = 0.05;
    }
  }
  else if (m_collectionName == "AntiKt4TopoJES") {
    if (y <= 0.8) { m_offset = 0.03; m_uncert = 0.05; }
    else if (y >= 0.8 && y < 1.2) { m_offset = 0.09; m_uncert = 0.10; }
    else if (y >= 1.2 && y < 2.1) { m_offset = 0.07; m_uncert = 0.10; }
    else if (y >= 2.1 && y < 2.8) { m_offset = 0.07; m_uncert = 0.12; }
    else if (y >= 2.8) { m_offset = 0.07; m_uncert = 0.12; }
    else { 
      cout << "ERROR: " << y << " not implemented, using default central region" << endl;
      m_offset = 0.03; m_uncert = 0.05;
    }
  }
  else { cout << "ERROR: " << y << " not implemented " << endl;}

  // So far, agreeed that we should not smear the MC for the nominal result 
  // since the MC describes data within uncertainties.
  
  m_offset = 0.0; // offset set to zero!


  // Return the requested parametrization

  int index = -1;
  if (y <= 0.8) index = 0;
  else if (y >= 0.8 && y < 1.2) index = 1;
  else if (y >= 1.2 && y < 2.1) index = 2;
  else if (y >= 2.1 && y < 2.8) index = 3;
  else if (y >= 2.8 && y < 3.6) index = 4;
  else if (y >= 3.6 && y <= 4.5) index = 5;
  else {
    //Protect against jets in the wrong region
    cout << "WARNING: Y outside of covered range (0.0-4.5): Returning parametrization for 4.5" << endl;
    return m_jerFunc[0];
  }

  return m_jerFunc[index];

}


float JERProviderMMCMoriond::getUncert(double pt, double y)
{

	if(!m_isInit) {
		cout << "JERProviderMMCMoriond not initialized." << endl;
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
      //cout << "WARNING: pt outside of covered range (10-5000): Pt set to 10 GeV" << endl;
      pt = 10.;
    } 
  
  if(pt > 5000.)
    {
      //cout << "WARNING: pt outside of covered range (10-5000): Pt set to 5000 GeV" << endl;
      pt = 5000.;
    } 
  
  
  // Get JER
  TF1* jer = getParam(y); jer->GetParameter(0); 
  
  float ptmin = 30; float ptmax = 500; 
  float factor = 0;
  
  if ( pt < ptmin ) {
    factor = (ptmin-pt)/10.;
    m_uncert = (1+factor)*m_uncert*getSigma(pt,y); // factor agreed below 30 GeV 
  }
  if ( pt >= ptmin && pt <= ptmax) m_uncert = m_uncert*getSigma(pt,y); // Data / MC comparison
  if ( pt > ptmax && pt <= 1000.) {
    factor = (pt - ptmax)/500.;
    m_uncert = (1+factor)*m_uncert*getSigma(pt,y);  // Above validated range (no data available), we double the uncertainty at 1000 GeV to be conservative.
  }  
  if ( pt > 1000.) {
    m_uncert = 2.0*m_uncert*getSigma(pt,y);  
  }  
  
  return m_uncert;
}

float JERProviderMMCMoriond::getSigma(double pt, double y)
{

	if(!m_isInit) {
		cout << "JERProviderMMCMoriond not initialized." << endl;
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
      //cout << "WARNING: pt outside of covered range (10-5000): Pt set to 10 GeV" << endl;
      pt = 10.;
    } 
  
  if(pt > 5000.)
    {
      //cout << "WARNING: pt outside of covered range (10-5000): Pt set to 5000 GeV" << endl;
      pt = 5000.;
    } 
  
  
  // Get JER
  TF1* jer = getParam(y);

  //  Noise, Stochastic, Constant terms and their errors ---> [N,S,C,Ne,Se,Ce] 
  m_param[0] = jer->GetParameter(0);   m_param[3] = jer->GetParError(0); 
  m_param[1] = jer->GetParameter(1);   m_param[4] = jer->GetParError(1);
  m_param[2] = jer->GetParameter(2);   m_param[5] = jer->GetParError(2);
   
  // cout << "Parametrization: sqrt([0]*[0]/(x*x) + [1]*[1]/x + [2]*[2])" << endl;
  // for (int i = 0; i < m_nParam/2; i++) { cout << "Parameter " << i << ":\t" << m_param[i] <<"\t +/- \t" << m_param[i+3]<< endl;}

  double sigma = sqrt(m_param[0]*m_param[0]/pt/pt + m_param[1]*m_param[1]/pt + m_param[2]*m_param[2]);  

  sigma = sigma + sigma*m_offset; // offset == 0 for the time being as agreed.

  return sigma;

}

}
}
