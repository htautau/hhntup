#ifndef FAKERATESCALER_H_
#define FAKERATESCALER_H_

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TString.h>

namespace FakeRateScaleFactor{
  class FakeRateScaler{
  public:
    FakeRateScaler(TString filename="FakeRateScaleFactor.txt"){
      //read in the config file
      
      ifstream infile;
      infile.open(filename.Data(), ifstream::in);

      TString datamc, id, trigger;
      int ptbin;
      double ptcenter, nominal, stat, sysup, sysdown;

      while(!infile.eof()){
	infile >> datamc >> id >> trigger >> ptbin >> ptcenter >> nominal >> stat >> sysup >> sysdown;
	if(datamc!=""){
	  std::map<TString,std::vector<double> >* nominalmap=0;
	  std::map<TString,std::vector<double> >* statmap=0;
	  std::map<TString,std::vector<double> >* systupmap=0;
	  std::map<TString,std::vector<double> >* systdownmap=0;
	  if(datamc=="Data" && id=="Loose"){
	    nominalmap=&m_looseData;
	    statmap=&m_looseDatastatuncertainty;
	    systupmap=&m_looseDatasystuncertaintyup;
	    systdownmap=&m_looseDatasystuncertaintydown;
	  }
	  else if(datamc=="Data" && id=="Medium"){
	    nominalmap=&m_mediumData;
	    statmap=&m_mediumDatastatuncertainty;
	    systupmap=&m_mediumDatasystuncertaintyup;
	    systdownmap=&m_mediumDatasystuncertaintydown;
	  }
	  else if(datamc=="Data" && id=="Tight"){
	    nominalmap=&m_tightData;
	    statmap=&m_tightDatastatuncertainty;
	    systupmap=&m_tightDatasystuncertaintyup;
	    systdownmap=&m_tightDatasystuncertaintydown;
	  }
	  else if(datamc=="MC" && id=="Loose"){
	    nominalmap=&m_looseMC;
	    statmap=&m_looseMCstatuncertainty;
	  }
	  else if(datamc=="MC" && id=="Medium"){
	    nominalmap=&m_mediumMC;
	    statmap=&m_mediumMCstatuncertainty;
	  }
	  else if(datamc=="MC" && id=="Tight"){
	    nominalmap=&m_tightMC;
	    statmap=&m_tightMCstatuncertainty;
	  }
	  else{
	    std::cout << "ERROR: FakeRateScaler::FakeRateScaler(): unknown type (Data,MC)="<<datamc << " or id="<<id<<" (Loose,Medium,Tight)" << std::endl; 
	    return;
	  }
	  
	  std::map<TString, std::vector<double> >::iterator it=nominalmap->find(trigger);
	  if( it==nominalmap->end() ){
	    std::vector<double> vec(9, 1.);
	    nominalmap->insert(std::make_pair<TString,std::vector<double> >(trigger,vec));
	    it=nominalmap->find(trigger);
	  }

	  std::map<TString, std::vector<double> >::iterator itstat=statmap->find(trigger);
	  if( itstat==statmap->end() ){
	    std::vector<double> vec(9, 1.);
	    statmap->insert(std::make_pair<TString,std::vector<double> >(trigger,vec));
	    itstat=statmap->find(trigger);
	  }
	  
	  (it->second)[ptbin]=nominal;
	  (itstat->second)[ptbin]=stat;

	  //std::cout << datamc << " " << id << " " << trigger <<" "<< ptbin <<" "<< ptcenter << " " << nominal <<std::endl;

	  if(datamc=="Data"){
	    std::map<TString, std::vector<double> >::iterator itsystup=systupmap->find(trigger);
	    if( itsystup==systupmap->end() ){
	      std::vector<double> vec(9, 1.);
	      systupmap->insert(std::make_pair<TString,std::vector<double> >(trigger,vec));
	      itsystup=systupmap->find(trigger);
	    }
	    
	    std::map<TString, std::vector<double> >::iterator itsystdown=systdownmap->find(trigger);
	    if( itsystdown==systdownmap->end() ){
	      std::vector<double> vec(9, 1.);
	      systdownmap->insert(std::make_pair<TString,std::vector<double> >(trigger,vec));
	      itsystdown=systdownmap->find(trigger);
	    }

	    //std::cout << " aaa" << std::endl;
	    (itsystup->second)[ptbin]=sysup;
	    //std::cout << " bbb" << std::endl;
	    (itsystdown->second)[ptbin]=sysdown;
	    
	  }
	  
	  //std::cout << datamc << " " << id << " " << trigger <<" "<< ptbin <<" "<< ptcenter << " " << nominal << " " << stat << " " << sysup << " " << sysdown << std::endl;
	}
      }

      infile.close();

    };
    ~FakeRateScaler(){ };
    void printData(){
      std::vector<TString> ids;
      ids.push_back("Loose");
      ids.push_back("Medium");
      ids.push_back("Tight");

      std::vector<TString> triggers;
      triggers.push_back("EF_tau20_medium1");
      triggers.push_back("EF_tau29_medium1");
      triggers.push_back("EF_tau20T_medium1");
      triggers.push_back("EF_tau29T_medium1");

      std::vector<double> ptcenter;
      ptcenter.push_back(2.5);
      ptcenter.push_back(7.5);
      ptcenter.push_back(12.5);
      ptcenter.push_back(17.5);
      ptcenter.push_back(22.5);
      ptcenter.push_back(27.5);
      ptcenter.push_back(32.5);
      ptcenter.push_back(37.5);
      ptcenter.push_back(70.0);
      for(unsigned int iid=0; iid<ids.size(); iid++){
	for(unsigned int itrigger=0; itrigger<triggers.size(); itrigger++){
	  std::cout << "scale factors for " << ids[iid] << " ID and " << triggers[itrigger] <<std::endl;
	  for(unsigned int ipt=0; ipt<9; ipt++){
	    std::cout<< std::setw(5) << ptcenter[ipt] << "GeV" << std::setw(10) 
		     << getScaleFactor(ptcenter[ipt]*1000, ids[iid], triggers[itrigger])
		     << " + " 
		     << std::setw(10) << getScaleFactorUncertainty(ptcenter[ipt]*1000, ids[iid], triggers[itrigger], true)
		     << " - " 
		     << std::setw(10) << getScaleFactorUncertainty(ptcenter[ipt]*1000, ids[iid], triggers[itrigger], false)
		     << std::endl;
	  }
	  std::cout << std::endl;
	}
      }

    };

    //id can be Loose, Medium,Tight
    //trigger can be EF_tau20_medium1, EF_tau20T_medium1, EF_tau29_medium1, EF_tau29T_medium
    //tau pt in MeV
    double getScaleFactor(double taupt, TString id="Loose", TString trigger="EF_tau20_medium1"){
      std::map<TString,std::vector<double> >* currentmapData;
      std::map<TString,std::vector<double> >* currentmapMC;

      if(taupt<1000.){
	std::cout << "ERROR: FakeRateScale::getScaleFactor(pt="<<taupt<<"): pT should be given in MeV" << std::endl; 
	return -11111111111.;
      }

      if(id=="Loose"){
	currentmapData=&m_looseData;
	currentmapMC=&m_looseMC;
      }
      else if(id=="Medium"){
	currentmapData=&m_mediumData;
	currentmapMC=&m_mediumMC;
      }
      else if(id=="Tight"){
	currentmapData=&m_tightData;
	currentmapMC=&m_tightMC;
      }
      else{
	std::cout << "ERROR: FakeRateScale::getScaleFactor(id="<<id<<"): id unknown" << std::endl; 
	return -11111111111.;
      }
      
      std::map<TString, std::vector<double> >::iterator itData=currentmapData->find(trigger);
      std::map<TString, std::vector<double> >::iterator itMC=currentmapMC->find(trigger);
      if(itData!=currentmapData->end() && itMC!=currentmapMC->end()){
	int ptbin=getptbin(taupt);
	if((itMC->second)[ptbin]==0)
	  return 1.;
	return (itData->second)[ptbin]/(itMC->second)[ptbin];
      }
      else{
	std::cout << "ERROR: FakeRateScale::getScaleFactor(trigger="<<trigger<<"): trigger unknown" << std::endl; 
      }
      
      return -11111111111.;
    };

    double getScaleFactorUncertainty(double taupt, TString id="Loose", TString trigger="EF_tau20_medium1", bool up=true){
      std::map<TString,std::vector<double> >* currentData;
      std::map<TString,std::vector<double> >* currentMC;
      std::map<TString,std::vector<double> >* currentstatData;
      std::map<TString,std::vector<double> >* currentstatMC;
      std::map<TString,std::vector<double> >* currentsyst;

      if(taupt<1000.){
	std::cout << "ERROR: FakeRateScale::getScaleFactorUncertainty(pt="<<taupt<<"): pT should be given in MeV" << std::endl; 
	return -11111111111.;
      }

      if(id=="Loose"){
	currentData=&m_looseData;
	currentMC=&m_looseMC;
	currentstatData=&m_looseMCstatuncertainty;
	currentstatMC=&m_looseDatastatuncertainty;
	if(up)
	  currentsyst=&m_looseDatasystuncertaintyup;
	else
	  currentsyst=&m_looseDatasystuncertaintydown;
      }
      else if(id=="Medium"){
	currentData=&m_mediumData;
	currentMC=&m_mediumMC;
	currentstatData=&m_mediumDatastatuncertainty;
	currentstatMC=&m_mediumMCstatuncertainty;
	if(up)
	  currentsyst=&m_mediumDatasystuncertaintyup;
	else
	  currentsyst=&m_mediumDatasystuncertaintydown;
      }
      else if(id=="Tight"){
	currentData=&m_tightData;
	currentMC=&m_tightMC;
	currentstatData=&m_tightDatastatuncertainty;
	currentstatMC=&m_tightMCstatuncertainty;
	if(up)
	  currentsyst=&m_tightDatasystuncertaintyup;
	else
	  currentsyst=&m_tightDatasystuncertaintydown;
      }
      else{
	std::cout << "ERROR: FakeRateScale::getScaleFactorUncertainty(id="<<id<<"): id unknown" << std::endl; 
	return -11111111111.;
      }
      
      std::map<TString, std::vector<double> >::iterator itData=currentData->find(trigger);
      std::map<TString, std::vector<double> >::iterator itMC=currentMC->find(trigger);
      std::map<TString, std::vector<double> >::iterator itstatData=currentstatData->find(trigger);
      std::map<TString, std::vector<double> >::iterator itstatMC=currentstatMC->find(trigger);
      std::map<TString, std::vector<double> >::iterator itsyst=currentsyst->find(trigger);
      
      if(itData!=currentData->end() 
	 && itMC!=currentMC->end()
	 && itstatData!=currentstatData->end()
	 && itstatMC!=currentstatMC->end()
	 && itsyst!=currentsyst->end()
	 ){
	int ptbin=getptbin(taupt);
	if((itMC->second)[ptbin]==0)
	  return 0.;
	double uncertnum=sqrt(pow((itstatData->second)[ptbin],2.)+pow((itsyst->second)[ptbin],2.));
	double error2=pow(uncertnum/(itMC->second)[ptbin],2.)+pow((itstatMC->second)[ptbin]*(itData->second)[ptbin]/pow((itMC->second)[ptbin],2),2);
	
	//double error2=0.;
	
	//if((itData->second)[ptbin]!=0 && (itMC->second)[ptbin]!=0 )
	//std::cout << "relnum="<< (itstatData->second)[ptbin]/(itData->second)[ptbin] << ", reldenom=" << (itstatMC->second)[ptbin]/(itMC->second)[ptbin] <<std::endl;

	return sqrt(error2);
      }
      else{
	std::cout << "ERROR: FakeRateScale::getScaleFactorUncertainty(trigger="<<trigger<<"): trigger unknown" << std::endl; 
      }
      
      return -11111111111.;
    };
  private:
    int getptbin(double pt){
      if(pt>40000)
	return 8;
      else if(pt>35000)
	return 7;
      else if(pt>30000)
	return 6;
      else if(pt>25000)
	return 5;
      else if(pt>20000)
	return 4;
      else if(pt>15000)
	return 3;
      else if(pt>10000)
	return 2;
      else if(pt>5000)
	return 1;
      else
	return 0;
    };

    std::map<TString, std::vector<double> >  m_looseData;
    std::map<TString, std::vector<double> >  m_mediumData;
    std::map<TString, std::vector<double> >  m_tightData;

    std::map<TString, std::vector<double> >  m_looseMC;
    std::map<TString, std::vector<double> >  m_mediumMC;
    std::map<TString, std::vector<double> >  m_tightMC;

    std::map<TString, std::vector<double> >  m_looseDatastatuncertainty;
    std::map<TString, std::vector<double> >  m_mediumDatastatuncertainty;
    std::map<TString, std::vector<double> >  m_tightDatastatuncertainty;

    std::map<TString, std::vector<double> >  m_looseDatasystuncertaintyup;
    std::map<TString, std::vector<double> >  m_mediumDatasystuncertaintyup;
    std::map<TString, std::vector<double> >  m_tightDatasystuncertaintyup;

    std::map<TString, std::vector<double> >  m_looseDatasystuncertaintydown;
    std::map<TString, std::vector<double> >  m_mediumDatasystuncertaintydown;
    std::map<TString, std::vector<double> >  m_tightDatasystuncertaintydown;

    std::map<TString, std::vector<double> >  m_looseMCstatuncertainty;
    std::map<TString, std::vector<double> >  m_mediumMCstatuncertainty;
    std::map<TString, std::vector<double> >  m_tightMCstatuncertainty;

  };
}

#endif /*FAKERATESCALER_H_*/
