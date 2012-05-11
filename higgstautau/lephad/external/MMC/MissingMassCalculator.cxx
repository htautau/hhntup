/*
  Missing Mass Calculator 
*/

#include <iostream>
#include <fstream>
// #include "MissingMassCalculator/MissingMassCalculator.h" // this is for package
#include "MissingMassCalculator.h"
// #include "JetResolution/JERProviderAthena.h"
#include <TRandom.h>
#include <TObject.h>
//SpeedUp committed from revision 163876
#include <TStopwatch.h>

// #if !defined(__CINT__) | defined(__MAKECINT__)
// ClassImp(MissingMassCalculator);
// #endif

namespace AnalysisFramework
{
namespace External
{

// dTheta probability density function for hadronic taus
inline double myDelThetaHadFunc(double *x, double *par)
{
  double fitval=1.0E-10;
  if(x[0]>TMath::Pi() || x[0]<0.0) return fitval;
  const double arg=x[0];
  const double arg_L=arg;
  const double mean=par[1];
  const double sigmaG=par[2];
  const double mpv=par[3];
  const double sigmaL=par[4];
  double normL=par[5];
  if(normL<0.0) normL=0.0;
  //SpeedUp
  //const double f1=normL/(1.0+normL)*TMath::Gaus(arg,mean,sigmaG);
  //const double f2=TMath::Landau(arg_L,mpv,sigmaL)/(1.0+normL);
  //fitval=par[0]*(f1+f2);

   const double g1=normL*TMath::Gaus(arg,mean,sigmaG);
   const double g2=TMath::Landau(arg_L,mpv,sigmaL);
  fitval=par[0]*(g1+g2)/(1.0+normL);

  if(fitval<0.0) return 0.0;
  return fitval;
}

// dTheta probability density function for leptonic taus
inline double myDelThetaLepFunc(double *x, double *par)
{
  double fitval=1.0E-10;
  if(x[0]>TMath::Pi() || x[0]<0.0) return fitval;
  double arg=x[0]/par[1];

  double normL=par[5];
  if(normL<0.0) normL=0.0;

  if(arg<1) arg=sqrt(fabs(arg));
  else arg=arg*arg;
  const double arg_L=x[0];
  const double mean=1.0;
  const double sigmaG=par[2];
  const double mpv=par[3];
  const double sigmaL=par[4];
  // speedup save one division....
  // const double f1=normL/(1.0+normL)*TMath::Gaus(arg,mean,sigmaG);
  // const double f2=TMath::Landau(arg_L,mpv,sigmaL)/(1.0+normL);
  // fitval=par[0]*(f1+f2);

   const double g1=normL*TMath::Gaus(arg,mean,sigmaG);
   const double g2=TMath::Landau(arg_L,mpv,sigmaL);
   fitval=par[0]*(g1+g2)/(1.0+normL);


  if(fitval<0.0) return 0.0;
  return fitval;
}



//ClassImp(MissingMassCalculator)
//______________________________constructor________________________________
MissingMassCalculator::MissingMassCalculator(std::string JERProviderFile){
  fUseVerbose=0;
  fSpeedStudy=false;
  AlgorithmVersion=1; // use V9 by default
  beamEnergy=3500.0; // for now LHC default is sqrt(S)=7 TeV 
  Niter_fit1=20;
  Niter_fit2=30;
  Niter_fit3=10;
  dRmax_tau=0.2;
  SearchMode=0;
  Nsigma_METscan=3.0; // number of sigmas for MET scan
  dTheta3d_binMin=0.0025;
  dTheta3d_binMax=0.02;
  fJERsyst=0; // no JER systematics by default (+/-1: up/down 1 sigma)
  METresSyst=0; // no MET resolution systematics by default (+/-1: up/down 1 sigma) 
  fApplyMassScale=0; // don't apply mass scale correction by default
  InputInfo.dataType=1; // set to "data" by default
  fUseTailCleanup=1; // cleanup by default for lep-had Moriond 2012 analysis
  METScanScheme=1; // MET-scan scheme: 0- use JER; 1- use simple sumEt & missingHt for Njet=0 events in (lep-had winter 2012)
  ClearInputStuff(InputInfo);

  phi1cache=10.;
  phi2cache=10.;
  nuvecsol1.resize(4);
  nuvecsol2.resize(4);

  //--- define histograms for histogram method
  //--- upper limits need to be revisied in the future!!! It may be not enough for some analyses
  fMfit_all=new TH1F("h1","M",2500,0.0,beamEnergy/3.5); // all solutions 
  fPXfit_nu1=new TH1F("h2","Px1",10000,-beamEnergy/3.5,beamEnergy/3.5); // Px for nu1 
  fPYfit_nu1=new TH1F("h3","Py1",10000,-beamEnergy/3.5,beamEnergy/3.5); // Py for nu1 
  fPZfit_nu1=new TH1F("h4","Pz1",10000,-beamEnergy/3.5,beamEnergy/3.5); // Pz for nu1 
  fPXfit_nu2=new TH1F("h5","Px2",10000,-beamEnergy/3.5,beamEnergy/3.5); // Px for nu2 
  fPYfit_nu2=new TH1F("h6","Py2",10000,-beamEnergy/3.5,beamEnergy/3.5); // Py for nu2 
  fPZfit_nu2=new TH1F("h7","Pz2",10000,-beamEnergy/3.5,beamEnergy/3.5); // Pz for nu2 
  
  if(fUseVerbose==1)
    {
      gDirectory->pwd();
      gDirectory->ls();
      std::cout << "opening JER" << std::endl;
    }

  myJER=new JERProviderMMC("AntiKt4LCTopoJES","Truth",JERProviderFile,fUseVerbose);  
  myJER->init();

  if(fUseVerbose==1)
    {
      gDirectory->pwd();
      gDirectory->ls();
    }

  // [tau_type][parLG][par]
  // leptonic tau
  //-par[0]
  fit_param[0][0][0]=-9.82013E-1;
  fit_param[0][0][1]=9.09874E-1;
  fit_param[0][0][2]=0.0;
  fit_param[0][0][3]=0.0;
  //-par[1]
  fit_param[0][1][0]=9.96303E1;
  fit_param[0][1][1]=1.68873E1;
  fit_param[0][1][2]=3.23798E-2;
  fit_param[0][1][3]=0.0;
  //-par[2]
  fit_param[0][2][0]=0.0;
  fit_param[0][2][1]=0.0;
  fit_param[0][2][2]=0.0;
  fit_param[0][2][3]=0.3; // fit value is 2.8898E-1, I use 0.3
  //-par[3]
  fit_param[0][3][0]=9.70055E1; 
  fit_param[0][3][1]=6.46175E1; 
  fit_param[0][3][2]=3.20679E-2;
  fit_param[0][3][3]=0.0;
  //-par[4]
  fit_param[0][4][0]=1.56865;	  
  fit_param[0][4][1]=2.28336E-1;
  fit_param[0][4][2]=1.09396E-1;
  fit_param[0][4][3]=1.99975E-3;
  //-par[5]
  fit_param[0][5][0]=0.0;
  fit_param[0][5][1]=0.0;
  fit_param[0][5][2]=0.0;
  fit_param[0][5][3]=0.66;
  //-----------------------------------------------------------------
  // 1-prong hadronic tau
  //-par[0]
  fit_param[1][0][0]=-2.42674;  
  fit_param[1][0][1]=7.69124E-1;
  fit_param[1][0][2]=0.0;
  fit_param[1][0][3]=0.0;
  //-par[1]
  fit_param[1][1][0]=9.52747E1; 
  fit_param[1][1][1]=1.26319E1; 
  fit_param[1][1][2]=3.09643E-2;
  fit_param[1][1][3]=0.0;
  //-par[2]
  fit_param[1][2][0]=1.71302E1; 
  fit_param[1][2][1]=3.00455E1; 
  fit_param[1][2][2]=7.49445E-2;
  fit_param[1][2][3]=0.0;
  //-par[3]
  fit_param[1][3][0]=1.06137E2; 
  fit_param[1][3][1]=6.01548E1; 
  fit_param[1][3][2]=3.50867E-2;
  fit_param[1][3][3]=0.0;
  //-par[4]
  fit_param[1][4][0]=4.26079E-1;
  fit_param[1][4][1]=1.76978E-1;
  fit_param[1][4][2]=1.43419;   
  fit_param[1][4][3]=0.0;
  //-par[5]
  fit_param[1][5][0]=0.0;
  fit_param[1][5][1]=0.0;
  fit_param[1][5][2]=0.0;
  fit_param[1][5][3]=0.4;
  //-----------------------------------------------------------------
  // 3-prong hadronic tau
  //-par[0]
  fit_param[2][0][0]=-2.43533;  
  fit_param[2][0][1]=6.12947E-1;
  fit_param[2][0][2]=0.0;
  fit_param[2][0][3]=0.0;
  //-par[1]
  fit_param[2][1][0]=9.54202;	  
  fit_param[2][1][1]=2.80011E-1;
  fit_param[2][1][2]=2.49782E-1;
  fit_param[2][1][3]=0.0;
  //-par[2]
  fit_param[2][2][0]=1.61325E1; 
  fit_param[2][2][1]=1.74892E1; 
  fit_param[2][2][2]=7.05797E-2;
  fit_param[2][2][3]=0.0;
  //-par[3]
  fit_param[2][3][0]=1.17093E2; 
  fit_param[2][3][1]=4.70000E1; 
  fit_param[2][3][2]=3.87085E-2;
  fit_param[2][3][3]=0.0;
  //-par[4]
  fit_param[2][4][0]=4.16557E-1;
  fit_param[2][4][1]=1.58902E-1;
  fit_param[2][4][2]=1.53516;   
  fit_param[2][4][3]=0.0;
  //-par[5]
  fit_param[2][5][0]=0.0;
  fit_param[2][5][1]=0.0;
  fit_param[2][5][2]=0.0;
  fit_param[2][5][3]=0.95;


}
//______________________________destructor________________________________
MissingMassCalculator::~MissingMassCalculator(){

  if(fUseVerbose==1)
    {
      std::cout << " in MMC destructor " << std::endl;
      gDirectory->pwd();
      gDirectory->ls();
    }

  delete fMfit_all;
  delete fPXfit_nu1;  
  delete fPYfit_nu1;  
  delete fPZfit_nu1;  
  delete fPXfit_nu2;  
  delete fPYfit_nu2;  
  delete fPZfit_nu2;  
  delete myJER;
}

//_____________________________________________________________________________
// Main Method to run MissingMassCalculator
int MissingMassCalculator::RunMissingMassCalculator() {
  ClearOutputStuff(OutputInfo);
  FinalizeInputStuff(InputInfo);
  DoMetResolution(InputInfo);
  if(SearchMode==0) 
    {
      if(AlgorithmVersion==1) OutputInfo.FitStatus=DitauMassCalculatorV9(InputInfo.vistau1,InputInfo.type_visTau1, 
									 InputInfo.vistau2,InputInfo.type_visTau2, 
									 InputInfo.MetVec);
      if(AlgorithmVersion==2) OutputInfo.FitStatus=DitauMassCalculatorV10fast(InputInfo.vistau1,InputInfo.type_visTau1, 
									      InputInfo.vistau2,InputInfo.type_visTau2, 
									      InputInfo.MetVec);
      if(AlgorithmVersion==3) OutputInfo.FitStatus=DitauMassCalculatorV9fast(InputInfo.vistau1,InputInfo.type_visTau1, 
									 InputInfo.vistau2,InputInfo.type_visTau2, 
									 InputInfo.MetVec);

    }
  DoOutputInfo(OutputInfo);
  PrintResults(OutputInfo);
  ClearInputStuff(InputInfo);
  return 1;
}

// ---- input Met vector
void MissingMassCalculator::SetMetVec(TVector2 vec) {
  InputInfo.MetVec=vec;
  return;
}
// ----- input vis Tau vectors
void MissingMassCalculator::SetVisTauVec(int i, TLorentzVector vec) {
  if(i==0) InputInfo.vistau1=vec;
  if(i==1) InputInfo.vistau2=vec;
  return;
}
// ----- input vis Tau type
void MissingMassCalculator::SetVisTauType(int i, int tautype) {
  if(i==0) InputInfo.type_visTau1=tautype/10; 
  if(i==1) InputInfo.type_visTau2=tautype/10; 
  return;
}
// ----- input vis Tau N-prong
void MissingMassCalculator::SetNprong(int i, int nprong) {
  if(i==0) InputInfo.Nprong_tau1=nprong; 
  if(i==1) InputInfo.Nprong_tau2=nprong; 
  return;
}
// ----- input SumEt
void MissingMassCalculator::SetSumEt(double sumEt) {
  InputInfo.SumEt=sumEt;
  return;
}
// ----- input data type
void MissingMassCalculator::SetIsData(int val) {
  if(val==0 || val==1) InputInfo.dataType=val;
  return;
}
// ---- input number of jets with Et>25 GeV
void MissingMassCalculator::SetNjet25(int val) {
  if(val>-1) InputInfo.Njet25=val;
  return;
}

// ----- input Met Scan parameters
void MissingMassCalculator::SetMetScanParams(double phi, double sigmaP, double sigmaL) {
  InputInfo.phi_jet=phi;
  InputInfo.METsigmaP=sigmaP;
  InputInfo.METsigmaL=sigmaL;
  return;
}

// input is sumEt after electrons and taus have been removed
// data_code=0 for data and =1 for MC
void MissingMassCalculator::SetMetScanParamsUE(double sumEt, double phi_scan, int data_code) {
  InputInfo.phi_jet=phi_scan;
  if(sumEt>2.0*beamEnergy) sumEt=sumEt/1000.0; // it's likely that sumEt was entered in MeV; this fix won't work only for a very small fraction of events
  double sigma=1.0;
  double sigmaSyst=0.10; // 10% systematics for now (be conservative)
  double METresScale=0.7; // using inclusive number for winter 2012
  if(data_code==1) METresScale=0.7; // use the same for data & MC
  METresScale=METresScale*(1.0+METresSyst*sigmaSyst); 
  // MET resolution can't be perfect in presence of other objects (i.e., electrons, jets, taus), so assume minSumEt=5.0 for now
  sigma= sumEt>pow(3.0/METresScale,2) ? METresScale*sqrt(sumEt) : 3.0; // assume that MET resolution can't be better than 3 GeV 
  InputInfo.METsigmaP=sigma;
  InputInfo.METsigmaL=sigma;
  InputInfo.SumEt=sumEt;
  InputInfo.dataType=data_code; // Sasha added on 09/26/11
  return;
} 

// input jet vectors and for MET resolution due to jets
// only consider jets with Et>20 GeV
void MissingMassCalculator::SetMetScanParamsJets(std::vector<TLorentzVector> jets) {
  for(unsigned int i=0; i<jets.size(); i++)
    {
      if(jets[i].Pt()>20.0) InputInfo.jet4vecs.push_back(jets[i]);
    }
  // re-order jets
  if(InputInfo.jet4vecs.size()>1)
    {
      TLorentzVector jet1(0.0,0.0,0.0,0.0);
      for(unsigned int i=1; i<InputInfo.jet4vecs.size(); i++)
	{
	  if(InputInfo.jet4vecs[i].Pt()>InputInfo.jet4vecs[i-1].Pt())
	    {
	      jet1=InputInfo.jet4vecs[i-1];
	      InputInfo.jet4vecs[i-1]=InputInfo.jet4vecs[i];
	      InputInfo.jet4vecs[i]=jet1;
	    }
	}
    }
  return;
}

//-------- clearing ditau container
void MissingMassCalculator::ClearDitauStuff(DitauStuff &fStuff) {
  fStuff.Mditau_best=0.0; 
  fStuff.Sign_best=1.0E6; 
  fStuff.nutau1.SetPxPyPzE(0.0,0.0,0.0,0.0); 
  fStuff.nutau2.SetPxPyPzE(0.0,0.0,0.0,0.0); 
  fStuff.RMSoverMPV=0.0;

  return;
}

//------- clearing input stuff
void MissingMassCalculator::ClearInputStuff(InputInfoStuff &fStuff) {
  fStuff.MetVec.Set(0.0,0.0);
  fStuff.vistau1.SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.vistau2.SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.type_visTau1=-1;
  fStuff.type_visTau2=-1;
  fStuff.Nprong_tau1=-1;
  fStuff.Nprong_tau2=-1;
//   fStuff.dataType=-1;
  fStuff.phi_jet=0.0;
  fStuff.METsigmaP=0.0;
  fStuff.METsigmaL=0.0;
  fStuff.SumEt=0.0;
  fStuff.jet4vecs.clear();
  fStuff.Njet25=-1;
  fStuff.DelPhiTT=0.0;
  MHtSigma1=-1.0;
  MHtSigma2=-1.0;
  MHtGaussFr=-1.0;

  return;
}

// checks units of input variables, converts into [GeV] if needed
void MissingMassCalculator::FinalizeInputStuff(InputInfoStuff &fStuff) {

  if(fStuff.vistau1.P()>1.0 && fStuff.vistau2.P()>1.0) fStuff.DelPhiTT=fabs(TVector2::Phi_mpi_pi(fStuff.vistau1.Phi()-fStuff.vistau2.Phi()));
  if(fStuff.vistau1.P()>3000.0 || fStuff.vistau2.P()>3000.0) // if units are MeV
    {
      if(fUseVerbose==1) std::cout << "converting to GeV" << std::endl;
      double scale=1000.0;
      fStuff.MetVec=(1.0/scale)*fStuff.MetVec;
      fStuff.vistau1=(1.0/scale)*fStuff.vistau1;
      fStuff.vistau2=(1.0/scale)*fStuff.vistau2;
      fStuff.SumEt=fStuff.SumEt/scale;
      for(unsigned int i=0; i<fStuff.jet4vecs.size(); i++)
	{
	  fStuff.jet4vecs[i]=(1.0/scale)*fStuff.jet4vecs[i];
	  // correcting sumEt, give priority to SetMetScanParamsUE()
	  if(METScanScheme==0)
	    {
	      if((fStuff.METsigmaP<0.1 || fStuff.METsigmaL<0.1) 
		 && fStuff.SumEt>fStuff.jet4vecs[i].Pt() 
		 && fStuff.jet4vecs[i].Pt()>20.0) fStuff.SumEt=fStuff.SumEt-fStuff.jet4vecs[i].Pt();
	    }
	}
    }
  else // if units are GeV
    {
      // correcting sumEt, give priority to SetMetScanParamsUE()
      if(METScanScheme==0)
	{
	  if(fStuff.METsigmaP<0.1 || fStuff.METsigmaL<0.1) 
	    {
	      for(unsigned int i=0; i<fStuff.jet4vecs.size(); i++)
		{
		  if(fStuff.SumEt>fStuff.jet4vecs[i].Pt() && fStuff.jet4vecs[i].Pt()>20.0) fStuff.SumEt=fStuff.SumEt-fStuff.jet4vecs[i].Pt();
		}
	    }
	}
    }


  // correcting sumEt, give priority to SetMetScanParamsUE()
  if(fStuff.METsigmaP<0.1 || fStuff.METsigmaL<0.1)
    {
      // 	  // correct only for electrons and hadronic tau's
      // 	  if(fStuff.SumEt>fStuff.vistau1.Pt() && (fStuff.vistau1.M()<0.05 || fStuff.vistau1.M()>0.12)) fStuff.SumEt=fStuff.SumEt-fStuff.vistau1.Pt();
      // 	  if(fStuff.SumEt>fStuff.vistau2.Pt() && (fStuff.vistau2.M()<0.05 || fStuff.vistau2.M()>0.12)) fStuff.SumEt=fStuff.SumEt-fStuff.vistau2.Pt();
      // correct only for electrons
      if(fStuff.SumEt>fStuff.vistau1.Pt() && fStuff.vistau1.M()<0.05) fStuff.SumEt=fStuff.SumEt-fStuff.vistau1.Pt();
      if(fStuff.SumEt>fStuff.vistau2.Pt() && fStuff.vistau2.M()<0.05) fStuff.SumEt=fStuff.SumEt-fStuff.vistau2.Pt();
    }
  // give priority to SetVisTauType, only do this if type_visTau1 and type_visTau2 are not set
  if(fStuff.type_visTau1<0 && fStuff.type_visTau2<0 && fStuff.Nprong_tau1>-1 && fStuff.Nprong_tau2>-1)
    {
      if(fStuff.Nprong_tau1==0 || fStuff.Nprong_tau1==1 || fStuff.Nprong_tau1==3) fStuff.type_visTau1=fStuff.Nprong_tau1;
      if(fStuff.Nprong_tau2==0 || fStuff.Nprong_tau2==1 || fStuff.Nprong_tau2==3) fStuff.type_visTau2=fStuff.Nprong_tau2;
    }
  
  // checking input mass of hadronic tau-1
  if(fStuff.type_visTau1==1 && fStuff.vistau1.M()!=0.8) 
    {
      double _pt, _phi, _eta, _m;
      _pt=fStuff.vistau1.Pt();
      _phi=fStuff.vistau1.Phi();
      _eta=fStuff.vistau1.Eta();
      _m=0.8;
      fStuff.vistau1.SetPtEtaPhiM(_pt,_eta,_phi,_m);
    }
  if(fStuff.type_visTau1==3 && fStuff.vistau1.M()!=1.2) 
    {
      double _pt, _phi, _eta, _m;
      _pt=fStuff.vistau1.Pt();
      _phi=fStuff.vistau1.Phi();
      _eta=fStuff.vistau1.Eta();
      _m=1.2;
	  fStuff.vistau1.SetPtEtaPhiM(_pt,_eta,_phi,_m);
    }
  // checking input mass of hadronic tau-2
  if(fStuff.type_visTau2==1 && fStuff.vistau2.M()!=0.8) 
    {
      double _pt, _phi, _eta, _m;
      _pt=fStuff.vistau2.Pt();
      _phi=fStuff.vistau2.Phi();
      _eta=fStuff.vistau2.Eta();
      _m=0.8;
      fStuff.vistau2.SetPtEtaPhiM(_pt,_eta,_phi,_m);
    }
  if(fStuff.type_visTau2==3 && fStuff.vistau2.M()!=1.2) 
    {
      double _pt, _phi, _eta, _m;
      _pt=fStuff.vistau2.Pt();
      _phi=fStuff.vistau2.Phi();
      _eta=fStuff.vistau2.Eta();
      _m=1.2;
      fStuff.vistau2.SetPtEtaPhiM(_pt,_eta,_phi,_m);
    }
  
  // give priority to SetMetScanParamsUE()
  if(fStuff.METsigmaP<0.1 || fStuff.METsigmaL<0.1)
    {
      if(METScanScheme==1) // default for Winter 2012 
	{
	  if((fStuff.vistau1.M()<0.12 && fStuff.vistau2.M()>0.12) || (fStuff.vistau2.M()<0.12 && fStuff.vistau1.M()>0.12)) // lep-had case
	    {
	      if(fStuff.Njet25==0 && fStuff.MetVec.Mod()<20.0) // low MET lep-had category for Winter 2012
		{
		  // giving priority to external settings
		  if(MHtSigma1<0.0) MHtSigma1=5.89;
		  if(MHtSigma2<0.0) MHtSigma2=15.47;
		  if(MHtGaussFr<0.0) MHtGaussFr=0.48;
		  fStuff.METsigmaP=MHtSigma2; // sigma of 2nd Gaussian for missing Ht resolution
		  fStuff.METsigmaL=MHtSigma2;		      
		}
	      if(fStuff.Njet25==0 && fStuff.MetVec.Mod()>=20.0) // high MET lep-had category for Winter 2012
		{
		  // giving priority to external settings
		  if(MHtSigma1<0.0) MHtSigma1=6.47;
		  if(MHtSigma2<0.0) MHtSigma2=16.82;
		  if(MHtGaussFr<0.0) MHtGaussFr=0.4767;
		  fStuff.METsigmaP=MHtSigma2; // sigma of 2nd Gaussian for missing Ht resolution
		  fStuff.METsigmaL=MHtSigma2;		      
		}
	      if(fStuff.Njet25>0) // Inclusive 1-jet and VBF lep-had categories for Winter 2012
		{
		  double sigmaSyst=0.10; // 10% systematics for now (be conservative)
		  double METresScale=0.56*(1.0+METresSyst*sigmaSyst); // for events with jets & analysis cuts for winter 2012
		  double METoffset=3.73*(1.0+METresSyst*sigmaSyst); // for events with jets & analysis cuts for winter 2012
		  // MET resolution can't be perfect in presence of other objects (i.e., electrons, jets, taus), so assume minSumEt=5.0 for now
		  double sigma= fStuff.SumEt>0.0 ? METoffset+METresScale*sqrt(fStuff.SumEt) : METoffset;
		  fStuff.METsigmaP=sigma;
		  fStuff.METsigmaL=sigma;
		}		  
	    }
	  else // for now, use this for lep-lep & had-had events, Feb 6, 2012 
	    {
	      double sigmaSyst=0.10; // 10% systematics for now (be conservative)
	      double METresScale=0.56*(1.0+METresSyst*sigmaSyst); // for events with jets & analysis cuts for winter 2012
	      double METoffset=3.73*(1.0+METresSyst*sigmaSyst); // for events with jets & analysis cuts for winter 2012
	      // MET resolution can't be perfect in presence of other objects (i.e., electrons, jets, taus), so assume minSumEt=5.0 for now
	      double sigma= fStuff.SumEt>0.0 ? METoffset+METresScale*sqrt(fStuff.SumEt) : METoffset;
	      fStuff.METsigmaP=sigma;
	      fStuff.METsigmaL=sigma;		  
	    }
	}
      if(METScanScheme==0) // old scheme with JER
	{
	  if(fStuff.dataType==0 || fStuff.dataType==1) SetMetScanParamsUE(fStuff.SumEt,fStuff.phi_jet,fStuff.dataType);
	  else SetMetScanParamsUE(fStuff.SumEt,fStuff.phi_jet,0);
	}
    }
  
  return;
}

//------- clearing output stuff
void MissingMassCalculator::ClearOutputStuff(OutputInfoStuff &fStuff) {

  fStuff.FitStatus=0;
  fStuff.FitSignificance[0]=-1.0;
  fStuff.FitSignificance[1]=-1.0;
  fStuff.FittedMass[0]=0.0;
  fStuff.FittedMass[1]=0.0;
  fStuff.nuvec1[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.objvec1[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.nuvec1[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.objvec1[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.nuvec2[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.objvec2[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.nuvec2[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.objvec2[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.totalvec[0].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.FittedMetVec[0].Set(0.0,0.0);
  fStuff.totalvec[1].SetPxPyPzE(0.0,0.0,0.0,0.0);
  fStuff.FittedMetVec[1].Set(0.0,0.0);
  fStuff.RMS2MPV=0.0;

  return;
}

//---- finalizes MET resolution parameters
void MissingMassCalculator::DoMetResolution(InputInfoStuff &fStuff) {
  if(fStuff.jet4vecs.size()>0 && METScanScheme==0)
    {
      if(fStuff.jet4vecs[0].Pt()<20.0) return; // Sasha added on 09/26/11
      double pt=fStuff.jet4vecs[0].Pt(); // range [10,5000] in JERProvider-00-00-02
//       double eta=fabs(fStuff.jet4vecs[0].Eta()); // |eta|<4.5 in JERProvider-00-00-02
      double y=fabs(fStuff.jet4vecs[0].Rapidity()); // JERProvider-00-00-09 requires rapidity
//       double sigma1 = (fStuff.jet4vecs[0].Pt())*(myJER.getSigma(pt,eta));
      double sigma1=0.15; // dummy value for initialization
      if(fStuff.dataType==0) sigma1 = (fStuff.jet4vecs[0].Pt())*(myJER->getRelResolutionMC(pt,y));
      if(fStuff.dataType==1) sigma1 = (fStuff.jet4vecs[0].Pt())*(myJER->getRelResolutionData(pt,y));
      double dPhi=0.0;
      double uncert=0.0;
      if(fJERsyst!=0)
	{
// 	  uncert = myJER->getUncert(pt,eta);
	  uncert = myJER->getResolutionUncert(pt,y);
	  sigma1=sigma1+(fStuff.jet4vecs[0].Pt())*fJERsyst*uncert;
	}
      fStuff.phi_jet=fStuff.jet4vecs[0].Phi();
      fStuff.METsigmaL=sqrt(fStuff.METsigmaL*fStuff.METsigmaL+sigma1*sigma1);
      for(unsigned int i=1; i<fStuff.jet4vecs.size(); i++) // start from 2nd jet
	{
	  if(fStuff.jet4vecs[i].Pt()<20.0) continue; // Sasha added on 09/26/11
	  pt=fStuff.jet4vecs[i].Pt(); // range [10,5000] in JERProvider-00-00-02
// 	  eta=fabs(fStuff.jet4vecs[i].Eta()); // |eta|<4.5 in JERProvider-00-00-02
	  y=fabs(fStuff.jet4vecs[i].Rapidity()); // JERProvider-00-00-09 uses rapidity

// 	  sigma1 = (fStuff.jet4vecs[i].Pt())*(myJER->getSigma(pt,eta));
	  if(fStuff.dataType==0) sigma1 = (fStuff.jet4vecs[i].Pt())*(myJER->getRelResolutionMC(pt,y));
	  if(fStuff.dataType==1) sigma1 = (fStuff.jet4vecs[i].Pt())*(myJER->getRelResolutionData(pt,y));

	  if(fJERsyst!=0)
	    {
// 	      uncert = myJER->getUncert(pt,eta);
	      uncert = myJER->getResolutionUncert(pt,y);
	      sigma1=sigma1+(fStuff.jet4vecs[i].Pt())*fJERsyst*uncert;
	    }
	  dPhi=fStuff.jet4vecs[i].DeltaPhi(fStuff.jet4vecs[0]);
	  fStuff.METsigmaL=sqrt(fStuff.METsigmaL*fStuff.METsigmaL+sigma1*sigma1*cos(dPhi)*cos(dPhi));
	  fStuff.METsigmaP=sqrt(fStuff.METsigmaP*fStuff.METsigmaP+sigma1*sigma1*sin(dPhi)*sin(dPhi));	  
	}
    }
  return;
}


//---------------------------- Accessors to output parameters ------------------------
//
// return fit status
int MissingMassCalculator::GetFitStatus() {
  return OutputInfo.FitStatus;
}

// returns fit significance
double MissingMassCalculator::GetFitSignificance(int fitcode) {
  double signif=-1.0;
  if(fitcode>-1 && fitcode<2 && OutputInfo.FitStatus>0) signif=OutputInfo.FitSignificance[fitcode];
  return signif;
}

// returns RMS/MPV according to histogram method
double MissingMassCalculator::GetRms2Mpv() {
  return OutputInfo.RMS2MPV;
}

// returns fitted Mass
double MissingMassCalculator::GetFittedMass(int fitcode) {
  double mass=0.0;
  if(fitcode<0 || fitcode>2) return 0.0;
  if(fitcode>-1 && fitcode<2 && OutputInfo.FitStatus>0) mass=OutputInfo.FittedMass[fitcode];
  if(fitcode==2) mass=(OutputInfo.objvec1[1]+OutputInfo.objvec2[1]).M();
  return mass;
}

// returns neutrino 4-vec
TLorentzVector MissingMassCalculator::GetNeutrino4vec(int fitcode, int ind) {
  TLorentzVector vec(0.0,0.0,0.0,0.0);
  if(fitcode>-1 && fitcode<2 && OutputInfo.FitStatus>0) 
    {
      if(ind==0) vec=OutputInfo.nuvec1[fitcode];
      if(ind==1) vec=OutputInfo.nuvec2[fitcode];
    }
  return vec;
}

// returns full tau 4-vec
TLorentzVector MissingMassCalculator::GetTau4vec(int fitcode, int ind) {
  TLorentzVector vec(0.0,0.0,0.0,0.0);
  if(fitcode>-1 && fitcode<2 && OutputInfo.FitStatus>0) 
    {
      if(ind==0) vec=OutputInfo.objvec1[fitcode];
      if(ind==1) vec=OutputInfo.objvec2[fitcode];
    }
  return vec;
}

// returns 4-vec for resonance
TLorentzVector MissingMassCalculator::GetResonanceVec(int fitcode) {
  TLorentzVector vec(0.0,0.0,0.0,0.0);
  if(fitcode>-1 && fitcode<2 && OutputInfo.FitStatus>0) vec=OutputInfo.objvec1[fitcode]+OutputInfo.objvec2[fitcode];
  return vec;
}

// returns 2-vec for fitted MET
TVector2 MissingMassCalculator::GetFittedMetVec(int fitcode) {
  TVector2 vec(0.0,0.0);
  if(fitcode>-1 && fitcode<2) vec=OutputInfo.FittedMetVec[fitcode];
  return vec;
}

// filanizes output information
void MissingMassCalculator::DoOutputInfo(OutputInfoStuff &fStuff)
{
  if(SearchMode==0) // di-tau mode
    {
      if(fStuff.FitStatus>0)
	{
	  fStuff.FitSignificance[0]=fDitauStuffFit.Sign_best;
	  fStuff.FittedMass[0]=fDitauStuffFit.Mditau_best;
	  fStuff.nuvec1[0]=fDitauStuffFit.nutau1;
	  fStuff.objvec1[0]=fDitauStuffFit.nutau1+InputInfo.vistau1;
	  fStuff.nuvec2[0]=fDitauStuffFit.nutau2;
	  fStuff.objvec2[0]=fDitauStuffFit.nutau2+InputInfo.vistau2;
	  fStuff.totalvec[0]=fStuff.objvec1[0]+fStuff.objvec2[0];
	  fStuff.FittedMetVec[0]=InputInfo.MetVec; // temporary value=input value, to be replaced later
	  
	  double scale=MassScale(1,fDitauStuffHisto.Mditau_best,InputInfo.type_visTau1,InputInfo.type_visTau2); // only for histo method for now
	  
	  fStuff.FitSignificance[1]=fDitauStuffHisto.Sign_best;
	  fStuff.FittedMass[1]=scale*fDitauStuffHisto.Mditau_best; // changed on 09/22/11, multiply by scale
	  fStuff.nuvec1[1]=fDitauStuffHisto.nutau1;
	  fStuff.objvec1[1]=fDitauStuffHisto.nutau1+InputInfo.vistau1;
	  fStuff.nuvec2[1]=fDitauStuffHisto.nutau2;
	  fStuff.objvec2[1]=fDitauStuffHisto.nutau2+InputInfo.vistau2;
	  fStuff.totalvec[1]=fStuff.objvec1[1]+fStuff.objvec2[1];
	  fStuff.FittedMetVec[1]=InputInfo.MetVec; // temporary value=input value, to be replaced later
	  fStuff.RMS2MPV=fDitauStuffHisto.RMSoverMPV;
	}
      if(fStuff.FittedMass[1]<=(InputInfo.vistau1+InputInfo.vistau2).M()) ClearOutputStuff(OutputInfo); // to avoid cases when FitStatus=1 but mass is not reconstructed
    }
  return;
}

// Printout of final results
void MissingMassCalculator::PrintResults(OutputInfoStuff fStuff) {
  if(fUseVerbose==1)
    {
      std::cout<<"------------- Printing Input for MissingMassCalculator --------------"<<std::endl;
      std::cout<<"................................................................................."<<std::endl;
      std::cout<<"   SearchMode="<<SearchMode<<"  -- 0: ditau, 1: WW, 2: W->taunu"<<std::endl;
      std::cout<<" Beam energy ="<<beamEnergy<<"  sqrt(S) for collisions ="<<2.0*beamEnergy<<std::endl;
      std::cout<<" met_x="<<InputInfo.MetVec.Px()<<" met_y="<<InputInfo.MetVec.Py()<<" MET="<<InputInfo.MetVec.Mod()<<" met_phi="<<InputInfo.MetVec.Phi()<<std::endl;
      std::cout<<" sumEt="<<InputInfo.SumEt<<" METsigmaP="<<InputInfo.METsigmaP<<" METsigmaL="<<InputInfo.METsigmaL<<" phi_scan="<<InputInfo.phi_jet<<std::endl;
      std::cout<<" tau_type1="<<InputInfo.type_visTau1<<" tau_type2="<<InputInfo.type_visTau2<<std::endl;
      std::cout<<" 1st visible tau:  P="<<InputInfo.vistau1.P()<<" Pt="<<InputInfo.vistau1.Pt()<<" M="<<InputInfo.vistau1.M()
	       <<" Phi="<<InputInfo.vistau1.Phi()<<" Eta="<<InputInfo.vistau1.Eta()<<std::endl;
      std::cout<<" 2nd visible tau:  P="<<InputInfo.vistau2.P()<<" Pt="<<InputInfo.vistau2.Pt()<<" M="<<InputInfo.vistau2.M()
	       <<" Phi="<<InputInfo.vistau2.Phi()<<" Eta="<<InputInfo.vistau2.Eta()<<std::endl;
      if(InputInfo.jet4vecs.size()>0)
	{
	  for(unsigned int i=0; i<InputInfo.jet4vecs.size(); i++)
	    {
	      std::cout<<" Printing jets: jet-"<<i<<" E="<<InputInfo.jet4vecs[i].E()<<" Pt="
		       <<InputInfo.jet4vecs[i].Pt()<<" Phi="<<InputInfo.jet4vecs[i].Phi()<<" Eta="<<InputInfo.jet4vecs[i].Eta()<<std::endl;
	    }
	}

      std::cout<<"------------- Printing Final Results for MissingMassCalculator --------------"<<std::endl;
      std::cout<<"................................................................................."<<std::endl;
      std::cout<<"  Fit status="<<fStuff.FitStatus<<std::endl;
      std::cout<<"___  Results for Fit Method ___"<<std::endl;
      std::cout<<" sign="<<fStuff.FitSignificance[0]<<std::endl;
      std::cout<<" mass="<<fStuff.FittedMass[0]<<std::endl;
      if(fStuff.FitStatus>0)
	{
	  std::cout<<" Neutrino-1: P="<<fStuff.nuvec1[0].P()<<"  Pt="<<fStuff.nuvec1[0].Pt()<<"  Eta="<<fStuff.nuvec1[0].Eta()<<"  Phi="<<fStuff.nuvec1[0].Phi()<<std::endl;
	  std::cout<<" Neutrino-2: P="<<fStuff.nuvec2[0].P()<<"  Pt="<<fStuff.nuvec2[0].Pt()<<"  Eta="<<fStuff.nuvec2[0].Eta()<<"  Phi="<<fStuff.nuvec2[0].Phi()<<std::endl;
	  std::cout<<" Tau-1: P="<<fStuff.objvec1[0].P()<<"  Pt="<<fStuff.objvec1[0].Pt()<<"  Eta="<<fStuff.objvec1[0].Eta()<<"  Phi="<<fStuff.objvec1[0].Phi()<<std::endl;
	  std::cout<<" Tau-2: P="<<fStuff.objvec2[0].P()<<"  Pt="<<fStuff.objvec2[0].Pt()<<"  Eta="<<fStuff.objvec2[0].Eta()<<"  Phi="<<fStuff.objvec2[0].Phi()<<std::endl;
	  if(SearchMode==0) 
	    {
	      std::cout<<" dR(nu1-visTau1)="<<fStuff.nuvec1[0].DeltaR(InputInfo.vistau1)<<std::endl;
	      std::cout<<" dR(nu2-visTau2)="<<fStuff.nuvec2[0].DeltaR(InputInfo.vistau2)<<std::endl;
	      std::cout<<" Tau-1 mass="<<fStuff.objvec1[0].M()<<std::endl;
	      std::cout<<" Tau-2 mass="<<fStuff.objvec2[0].M()<<std::endl;
	    }
	
	  std::cout<<" Resonance: M="<<fStuff.totalvec[0].M()<<" P="<<fStuff.totalvec[0].P()
		   <<" Pt="<<fStuff.totalvec[0].Pt()<<" Phi="<<fStuff.totalvec[0].Phi()<<" Eta="<<fStuff.totalvec[0].Eta()<<std::endl;
	}
      std::cout<<"___  Mass solution for Histogram Method ___"<<std::endl;
      std::cout<<" sign="<<fStuff.FitSignificance[1]<<std::endl;
      std::cout<<" mass="<<fStuff.FittedMass[1]<<std::endl;
      std::cout<<" rms/mpv="<<fStuff.RMS2MPV<<std::endl;
      std::cout<<"___  All results for Histogram Method ___"<<std::endl;
      if(fStuff.FitStatus>0)
	{
	  std::cout<<" Neutrino-1: Pt="<<fStuff.nuvec1[1].Pt()<<"  Eta="<<fStuff.nuvec1[1].Eta()<<"  Phi="<<fStuff.nuvec1[1].Phi()<<std::endl;
	  std::cout<<" Neutrino-2: Pt="<<fStuff.nuvec2[1].Pt()<<"  Eta="<<fStuff.nuvec2[1].Eta()<<"  Phi="<<fStuff.nuvec2[1].Phi()<<std::endl;
	  std::cout<<" Tau-1: P="<<fStuff.objvec1[1].P()<<"  Pt="<<fStuff.objvec1[1].Pt()<<"  Eta="<<fStuff.objvec1[1].Eta()<<"  Phi="<<fStuff.objvec1[1].Phi()<<std::endl;
	  std::cout<<" Tau-2: P="<<fStuff.objvec2[1].P()<<"  Pt="<<fStuff.objvec2[1].Pt()<<"  Eta="<<fStuff.objvec2[1].Eta()<<"  Phi="<<fStuff.objvec2[1].Phi()<<std::endl;
	  if(SearchMode==0) 
	    {
	      std::cout<<" dR(nu1-visTau1)="<<fStuff.nuvec1[1].DeltaR(InputInfo.vistau1)<<std::endl;
	      std::cout<<" dR(nu2-visTau2)="<<fStuff.nuvec2[1].DeltaR(InputInfo.vistau2)<<std::endl;
	      std::cout<<" Tau-1 mass="<<fStuff.objvec1[1].M()<<std::endl;
	      std::cout<<" Tau-2 mass="<<fStuff.objvec2[1].M()<<std::endl;
	    }
	  std::cout<<" Resonance: M="<<fStuff.totalvec[1].M()<<" P="<<fStuff.totalvec[1].P()
		   <<" Pt="<<fStuff.totalvec[1].Pt()<<" Phi="<<fStuff.totalvec[1].Phi()<<" Eta="<<fStuff.totalvec[1].Eta()<<std::endl;
	}
      std::cout<<"__________________ End of Printout ___________________________________"<<std::endl;
    }
  return;
}

// returns P(nu1) & P(nu2)
int MissingMassCalculator::NuPsolution(TVector2 met_vec, double theta1, double phi1, 
					double theta2, double phi2, double &P1, double &P2) {
  int solution_code=0; // 0== no solution, 1==with solution
  P1=0.0;
  P2=0.0;
  double D=sin(theta1)*sin(theta2)*sin(phi2-phi1);
  if(fabs(D)>0.0) // matrix deteriminant is non-zero
    {
      P1=(met_vec.Px()*sin(phi2)-met_vec.Py()*cos(phi2))*sin(theta2)/D;
      P2=(met_vec.Py()*cos(phi1)-met_vec.Px()*sin(phi1))*sin(theta1)/D;
      if(P1>0.0 && P2>0.0) solution_code=1;
    }
  return solution_code;
}

// returns P1, P2, and theta1 & theta2 solutions
int MissingMassCalculator::NuPsolutionV2(const TVector2 & met_vec,const TLorentzVector & vec1,const  TLorentzVector & vec2, 
					 const double & mass1, const double & mass2, const double & phi1,const double & phi2, 
					 std::vector<TLorentzVector> &nu_vec1, std::vector<TLorentzVector> &nu_vec2) {

  int solution_code=0; // 0== no solution, 1==with solution
  //------ re-zeroing input
  nu_vec1.clear();
  nu_vec2.clear();

  std::vector<double> P1;
  std::vector<double> P2;
  std::vector<double> theta1;
  std::vector<double> theta2;
  TLorentzVector dummy(0.0,0.0,0.0,0.0);



  //------defining/calculating service parameters
  if(SearchMode==0 && (vec1.M()>1.5*mass1 || vec2.M()>1.5*mass2)) return 0; // no solutions if visM>1.5*M(tau)
  double Dsin=sin(phi2-phi1);
  if(fabs(Dsin)>0.0)
    {
      double e[2];

      e[0]=vec1.Px()*cos(phi1)+vec1.Py()*sin(phi1);
      e[1]=vec2.Px()*cos(phi2)+vec2.Py()*sin(phi2);
      double f[2];
      f[0]=vec1.Pz();
      f[1]=vec2.Pz();
      double Ev[2];
      Ev[0]=vec1.E();
      Ev[1]=vec2.E();
      double d[2];
      d[0]=fabs(mass1*mass1-vec1.M()*vec1.M()); // !!!! use abs value because mass of vis_tau can be slighty mis-measured 
      d[1]=fabs(mass2*mass2-vec2.M()*vec2.M()); // !!!! use abs value because mass of vis_tau can be slighty mis-measured 
      double R[2];
      R[0]=(met_vec.Px()*sin(phi2)-met_vec.Py()*cos(phi2))/Dsin;
      R[1]=(met_vec.Py()*cos(phi1)-met_vec.Px()*sin(phi1))/Dsin;

      double tA[2];
      double tB[2];
      double Determ[2];
      double Psolution1[2];
      double Psolution2[2];
      //--------------- obtaining P-solutions
      for(unsigned int i=0; i<2; i++)
	{
	  tA[i]=R[i]*e[i]+d[i]/2.0;
	  tB[i]=Ev[i]*Ev[i]-f[i]*f[i];
	  if(tB[i]<=0.0) return 0; // can't be zero
	  Determ[i]=Ev[i]*Ev[i]*tA[i]*tA[i]-(tA[i]*tA[i]+f[i]*f[i]*R[i]*R[i])*tB[i];
	  if(Determ[i]<0.0) return 0; // determinant is zero ==> no solutions
	  else Determ[i]=sqrt(Determ[i]);
	}
      Psolution1[0]=(Ev[0]*tA[0]-Determ[0])/tB[0];
      Psolution1[1]=(Ev[0]*tA[0]+Determ[0])/tB[0];
      Psolution2[0]=(Ev[1]*tA[1]-Determ[1])/tB[1];
      Psolution2[1]=(Ev[1]*tA[1]+Determ[1])/tB[1];
      for(unsigned int i=0; i<2; i++)
	{
	  if(Psolution1[i]>0.0) P1.push_back(Psolution1[i]);
	  if(Psolution2[i]>0.0) P2.push_back(Psolution2[i]);
	}
      if(P1.size()==0 || P2.size()==0) return 0; // no solutions
      //--------- obtaining theta solutions
      double sinT;
      double theta_s[2];
      double check[2];
      double err=0.0001;
      for(unsigned int i=0; i<P1.size(); i++)
	{
	  sinT=R[0]/P1[i];
	  if(fabs(sinT)<=1.0)
	    {
	      theta_s[0]=asin(sinT);
	      theta_s[1]=TMath::Pi()-theta_s[0];
	      check[0]=P1[i]*(Ev[0]-sin(theta_s[0])*e[0]-cos(theta_s[0])*f[0]);
	      check[1]=P1[i]*(Ev[0]-sin(theta_s[1])*e[0]-cos(theta_s[1])*f[0]);
	      if(fabs(check[0]-d[0]/2.0)<err) 
		{
		  dummy.SetPxPyPzE(P1[i]*sin(theta_s[0])*cos(phi1),P1[i]*sin(theta_s[0])*sin(phi1),P1[i]*cos(theta_s[0]),P1[i]);
		  nu_vec1.push_back(dummy);
		}
	      if(fabs(check[1]-d[0]/2.0)<err) 
		{
		  dummy.SetPxPyPzE(P1[i]*sin(theta_s[1])*cos(phi1),P1[i]*sin(theta_s[1])*sin(phi1),P1[i]*cos(theta_s[1]),P1[i]);
		  nu_vec1.push_back(dummy);
		}	      
	    }
	}
      for(unsigned int i=0; i<P2.size(); i++)
	{
	  sinT=R[1]/P2[i];
	  if(fabs(sinT)<=1.0)
	    {
	      theta_s[0]=asin(sinT);
	      theta_s[1]=TMath::Pi()-theta_s[0];
	      check[0]=P2[i]*(Ev[1]-sin(theta_s[0])*e[1]-cos(theta_s[0])*f[1]);
	      check[1]=P2[i]*(Ev[1]-sin(theta_s[1])*e[1]-cos(theta_s[1])*f[1]);
	      if(fabs(check[0]-d[1]/2.0)<err) 
		{
		  dummy.SetPxPyPzE(P2[i]*sin(theta_s[0])*cos(phi2),P2[i]*sin(theta_s[0])*sin(phi2),P2[i]*cos(theta_s[0]),P2[i]);
		  nu_vec2.push_back(dummy);
		}
	      if(fabs(check[1]-d[1]/2.0)<err)
		{
		  dummy.SetPxPyPzE(P2[i]*sin(theta_s[1])*cos(phi2),P2[i]*sin(theta_s[1])*sin(phi2),P2[i]*cos(theta_s[1]),P2[i]);
		  nu_vec2.push_back(dummy);		  
		}
	    }
	}

      if(nu_vec1.size()>0 && nu_vec2.size()>0) solution_code=1;
      else return 0;
    }
  return solution_code;
}



// returns P1, P2, and theta1 & theta2 solutions
int MissingMassCalculator::NuPsolutionV2fast(const TVector2 & met_vec,const TLorentzVector & vec1,const  TLorentzVector & vec2, 
					 const double & mass1, const double & mass2, const double & phi1,const double & phi2, 
						  int & nsol1, int & nsol2) {


  int solution_code=0; // 0== no solution, 1==with solution
  nsol1=0;
  nsol2=0;
  


  //SpeedUp cache sin/cos
  bool updatedphi1=false;
  bool updatedphi2=false;
  bool updatedvec1=false;
  bool updatedvec2=false;
  bool updatedmet=false;
  bool updateddphi=false;

  //assumes vec1 cannot change keeping the same energy

  if (met_vec.Px()!=metvecpxcache || met_vec.Py()!=metvecpycache){
    metvecpxcache=met_vec.Px();
    metvecpycache=met_vec.Py();
    updatedmet=true;
  }

  
  if (vec1.E()!=vec1ecache){
    vec1ecache=vec1.E();
    vec1mcache=vec1.M();
    vec1m2cache=vec1mcache*vec1mcache;
    fcache0=vec1.Pz();
    Evcache[0]=vec1.E();
    tBcache[0]=Evcache[0]*Evcache[0]-fcache0*fcache0;
    if (tBcache[0]<=0.0) return 0;
    updatedvec1=true;
    }

  if (vec2.E()!=vec2ecache){
    vec2ecache=vec2.E();
    vec2mcache=vec2.M();	
    vec2m2cache=vec2mcache*vec2mcache;
    fcache1=vec2.Pz();
    Evcache[1]=vec2.E();
    tBcache[1]=Evcache[1]*Evcache[1]-fcache1*fcache1;
    if (tBcache[1]<=0.0) return 0;
    updatedvec2=true;
  }
  

  if (phi1!=phi1cache){
    phi1cache=phi1;
    sinphi1cache=sin(phi1cache);
    cosphi1cache=cos(phi1cache);
    updatedphi1=true;
  } 

    
  if (phi2!=phi2cache){
    phi2cache=phi2;
    sinphi2cache=sin(phi2cache);
    cosphi2cache=cos(phi2cache);
    updatedphi2=true;
  }



  if (updatedphi1 || updatedphi2){
      dsincache=sin(phi2-phi1);
      updateddphi=true;
  }

  //FIXME could be improve
    if (updatedvec1 || updatedphi1 ){
      ecache[0]=vec1.Px()*cosphi1cache+vec1.Py()*sinphi1cache;
    }

    if (updatedvec2 || updatedphi2 ){
      ecache[1]=vec2.Px()*cosphi2cache+vec2.Py()*sinphi2cache;
    }
    
    if (dsincache==0.) return solution_code;
    
  
  double P1[2];
  double P2[2];
  int nsolp1=0;
  int nsolp2=0;
  

  //std::vector<double> theta1;
  //std::vector<double> theta2;
  // SpeedUp better reuse the same one  TLorentzVector dummy(0.0,0.0,0.0,0.0);

  //------defining/calculating service parameters
  if(SearchMode==0 && (vec1mcache>1.5*mass1 || vec2mcache>1.5*mass2)) return 0; // no solutions if visM>1.5*M(tau)

  double d[2];
  d[0]=fabs(mass1*mass1-vec1m2cache); // !!!! use abs value because mass of vis_tau can be slighty mis-measured 
  d[1]=fabs(mass2*mass2-vec2m2cache); // !!!! use abs value because mass of vis_tau can be slighty mis-measured 
      
  if (updatedmet || updateddphi){
    Rcache[0]=(metvecpxcache*sinphi2cache-metvecpycache*cosphi2cache)/dsincache;
  }
  if (updatedmet || updateddphi){      
    Rcache[1]=(metvecpycache*cosphi1cache-metvecpxcache*sinphi1cache)/dsincache;
  }

  if (updatedmet || updateddphi || updatedvec1){
    ffRRcache[0]=pow(fcache0,2)*pow(Rcache[0],2);
  }
  if (updatedmet || updateddphi || updatedvec2){
    ffRRcache[1]=pow(fcache1,2)*pow(Rcache[1],2);
  }


  
  

  double tA[2];
  double Determ[2];
  //double Psolution1[2];
  //double Psolution2[2];
  double Psol;
  //--------------- obtaining P-solutions
  for(int i=0; i<2; i++)
    {
      tA[i]=Rcache[i]*ecache[i]+d[i]/2.0;
      Determ[i]=pow(Evcache[i],2)*pow(tA[i],2)-(pow(tA[i],2)+ffRRcache[i])*tBcache[i];
      if(Determ[i]<0.0) return 0; // determinant is zero ==> no solutions
      else Determ[i]=sqrt(Determ[i]);
    }

  Psol=(Evcache[0]*tA[0]-Determ[0])/tBcache[0];
  if(Psol>0.0) { 
    P1[nsolp1]=Psol;
    ++nsolp1;
  }
  
  Psol=(Evcache[0]*tA[0]+Determ[0])/tBcache[0];
  if(Psol>0.0){
    P1[nsolp1]=Psol;
    ++nsolp1;
  }
  if (nsolp1==0) return 0;
  
  Psol=(Evcache[1]*tA[1]-Determ[1])/tBcache[1];
  if(Psol>0.0) { 
    P2[nsolp2]=Psol;
    ++nsolp2;
  }
  Psol=(Evcache[1]*tA[1]+Determ[1])/tBcache[1];
  if(Psol>0.0) { 
    P2[nsolp2]=Psol;
    ++nsolp2;
  }
  if (nsolp2==0) return 0;
      
  

//       Psolution1[0]=(Ev[0]*tA[0]-Determ[0])/tB[0];
//       Psolution1[1]=(Ev[0]*tA[0]+Determ[0])/tB[0];
//       Psolution2[0]=(Ev[1]*tA[1]-Determ[1])/tB[1];
//       Psolution2[1]=(Ev[1]*tA[1]+Determ[1])/tB[1];
//       for(int i=0; i<2; i++)
// 	{
// 	  if(Psolution1[i]>0.0) P1.push_back(Psolution1[i]);
// 	  if(Psolution2[i]>0.0) P2.push_back(Psolution2[i]);
// 	}
//       if(P1.size()==0 || P2.size()==0) return 0; // no solutions
      //--------- obtaining theta solutions
  double sinT;
  double err=0.0001;
  for(int i=0; i<nsolp1; ++i)
	{
	  const double & P1i = P1[i];
	  sinT=Rcache[0]/P1i;
	  if(fabs(sinT)<=1.0)
	    {
	      const double costheta0=sqrt(1-sinT*sinT);
	      const double & sintheta0=sinT;
	      const double costheta1=-costheta0;
	      const double & sintheta1=sintheta0;

	      const double check0=P1i*(Evcache[0]-sintheta0*ecache[0]-costheta0*fcache0);
	      const double check1=P1i*(Evcache[0]-sintheta1*ecache[0]-costheta1*fcache0);

	      if(fabs(check0-d[0]/2.0)<err) 
		{
		  nuvecsol1[nsol1].SetPxPyPzE(P1i*sintheta0*cosphi1cache,P1i*sintheta0*sinphi1cache,P1i*costheta0,P1i);
		  ++nsol1;
		}

	      if(fabs(check1-d[0]/2.0)<err) 
		{
		  nuvecsol1[nsol1].SetPxPyPzE(P1i*sintheta1*cosphi1cache,P1i*sintheta1*sinphi1cache,P1i*costheta1,P1i);
		  ++nsol1;
		}	      
	    }
	}

  //SpeedUp no need to calculate nu_vec2 if no solution for nu_vec1
  if (nsol1==0) return 0;
      for(int i=0; i<nsolp2; ++i)
	{
	  const double & P2i = P2[i];
	  sinT=Rcache[1]/P2i;
	  if(fabs(sinT)<=1.0)
	    {
	      const double costheta0=sqrt(1-sinT*sinT);
	      const double & sintheta0=sinT;
	      const double costheta1=-costheta0;
	      const double & sintheta1=sintheta0;

	      const double check0=P2i*(Evcache[1]-sintheta0*ecache[1]-costheta0*fcache1);
	      const double check1=P2i*(Evcache[1]-sintheta1*ecache[1]-costheta1*fcache1);
	      if(fabs(check0-d[1]/2.0)<err) 
		{
		  nuvecsol2[nsol2].SetPxPyPzE(P2i*sintheta0*cosphi2cache,P2i*sintheta0*sinphi2cache,P2i*costheta0,P2i);
		  ++nsol2;
		}
	      if(fabs(check1-d[1]/2.0)<err)
		{
		  nuvecsol2[nsol2].SetPxPyPzE(P2i*sintheta1*cosphi2cache,P2i*sintheta1*sinphi2cache,P2i*costheta1,P2i);
		  ++nsol2;
		}
	    }
	}
      if (nsol2==0) return 0;      
      if (nsol1>4 || nsol2>4)
	{
	  std::cout << "ERROR nsol1 or nsol2 >2 ! should never happen " << std::endl;
	  return 0;
	}

      solution_code=1;

  return solution_code;
}




// dTheta3 parameterization
// This parameterization is obtained for taus with 20<P<220, everything above 220 GeV is pure extrapolation
// Parameterization is based on Z->tautau sample.
// OBSOLETE SHOULD NOT BE USED
double MissingMassCalculator::dTheta3Dfit_parameterization(int tau_type, int dT3dfit_par, int Pfit_par) {
  double fit_para[3][6][4];
  // leptonic tau
  //-par[0]
  //SpeedUp changed name fit_param->fit_para to avoid collision
  fit_para[0][0][0]=-9.82013E-1;
  fit_para[0][0][1]=9.09874E-1;
  fit_para[0][0][2]=0.0;
  fit_para[0][0][3]=0.0;
  //-par[1]
  fit_para[0][1][0]=9.96303E1;
  fit_para[0][1][1]=1.68873E1;
  fit_para[0][1][2]=3.23798E-2;
  fit_para[0][1][3]=0.0;
  //-par[2]
  fit_para[0][2][0]=0.0;
  fit_para[0][2][1]=0.0;
  fit_para[0][2][2]=0.0;
  fit_para[0][2][3]=0.3; // fit value is 2.8898E-1, I use 0.3
  //-par[3]
  fit_para[0][3][0]=9.70055E1; 
  fit_para[0][3][1]=6.46175E1; 
  fit_para[0][3][2]=3.20679E-2;
  fit_para[0][3][3]=0.0;
  //-par[4]
  fit_para[0][4][0]=1.56865;	  
  fit_para[0][4][1]=2.28336E-1;
  fit_para[0][4][2]=1.09396E-1;
  fit_para[0][4][3]=1.99975E-3;
  //-par[5]
  fit_para[0][5][0]=0.0;
  fit_para[0][5][1]=0.0;
  fit_para[0][5][2]=0.0;
  fit_para[0][5][3]=0.66;
  //-----------------------------------------------------------------
  // 1-prong hadronic tau
  //-par[0]
  fit_para[1][0][0]=-2.42674;  
  fit_para[1][0][1]=7.69124E-1;
  fit_para[1][0][2]=0.0;
  fit_para[1][0][3]=0.0;
  //-par[1]
  fit_para[1][1][0]=9.52747E1; 
  fit_para[1][1][1]=1.26319E1; 
  fit_para[1][1][2]=3.09643E-2;
  fit_para[1][1][3]=0.0;
  //-par[2]
  fit_para[1][2][0]=1.71302E1; 
  fit_para[1][2][1]=3.00455E1; 
  fit_para[1][2][2]=7.49445E-2;
  fit_para[1][2][3]=0.0;
  //-par[3]
  fit_para[1][3][0]=1.06137E2; 
  fit_para[1][3][1]=6.01548E1; 
  fit_para[1][3][2]=3.50867E-2;
  fit_para[1][3][3]=0.0;
  //-par[4]
  fit_para[1][4][0]=4.26079E-1;
  fit_para[1][4][1]=1.76978E-1;
  fit_para[1][4][2]=1.43419;   
  fit_para[1][4][3]=0.0;
  //-par[5]
  fit_para[1][5][0]=0.0;
  fit_para[1][5][1]=0.0;
  fit_para[1][5][2]=0.0;
  fit_para[1][5][3]=0.4;
  //-----------------------------------------------------------------
  // 3-prong hadronic tau
  //-par[0]
  fit_para[2][0][0]=-2.43533;  
  fit_para[2][0][1]=6.12947E-1;
  fit_para[2][0][2]=0.0;
  fit_para[2][0][3]=0.0;
  //-par[1]
  fit_para[2][1][0]=9.54202;	  
  fit_para[2][1][1]=2.80011E-1;
  fit_para[2][1][2]=2.49782E-1;
  fit_para[2][1][3]=0.0;
  //-par[2]
  fit_para[2][2][0]=1.61325E1; 
  fit_para[2][2][1]=1.74892E1; 
  fit_para[2][2][2]=7.05797E-2;
  fit_para[2][2][3]=0.0;
  //-par[3]
  fit_para[2][3][0]=1.17093E2; 
  fit_para[2][3][1]=4.70000E1; 
  fit_para[2][3][2]=3.87085E-2;
  fit_para[2][3][3]=0.0;
  //-par[4]
  fit_para[2][4][0]=4.16557E-1;
  fit_para[2][4][1]=1.58902E-1;
  fit_para[2][4][2]=1.53516;   
  fit_para[2][4][3]=0.0;
  //-par[5]
  fit_para[2][5][0]=0.0;
  fit_para[2][5][1]=0.0;
  fit_para[2][5][2]=0.0;
  fit_para[2][5][3]=0.95;

  return fit_para[tau_type][dT3dfit_par][Pfit_par];
}

// returns parameters for dTheta3D pdf
double MissingMassCalculator::dTheta3Dparam(const int & parInd, const double & P_tau, const double *par) {

  if(P_tau<0.0) return 0.0;
  if(parInd==0) {
    return (par[0]+par[1]*P_tau)*0.00125;
  }
  else {
    return par[0]*(exp(-par[1]*P_tau)+par[2]/P_tau)+par[3];
  }
  
}

// returns dTheta3D probability based on ATLAS parameterization
double MissingMassCalculator::dTheta3d_probability(const int & tau_type,const double & dTheta3d,const double & P_tau) {
  double prob=1.0E-10;
  int tau_code=tau_type;
  double cut1=0.0;
  double cut2=0.5*TMath::Pi();
  if(tau_type==3) tau_code=2;
  if(tau_type>3) 
    {
      std::cout<<">>>> WARNING in MissingMassCalculator::dTheta3d_probability() <<<<"<<std::endl;
      std::cout<<"..... wrong tau_type="<<tau_type<<std::endl;
      std::cout<<"..... returning prob="<<prob<<std::endl;
      std::cout<<"____________________________________________________________"<<std::endl;
      return prob;
    }
  TF1* dTheta3Fun=0;
  if(tau_type==0) {
    dTheta3Fun=new TF1("pdf",myDelThetaLepFunc,cut1,cut2,6);
  }
  else {
    dTheta3Fun=new TF1("pdf",myDelThetaHadFunc,cut1,cut2,6);
  }
  
  double param[4];
  param[0]=dTheta3Dfit_parameterization(tau_code,0,0);
  param[1]=dTheta3Dfit_parameterization(tau_code,0,1);
  param[2]=dTheta3Dfit_parameterization(tau_code,0,2);
  param[3]=dTheta3Dfit_parameterization(tau_code,0,3);
  dTheta3Fun->SetParameter(0,dTheta3Dparam(0,P_tau,param));
  param[0]=dTheta3Dfit_parameterization(tau_code,1,0);
  param[1]=dTheta3Dfit_parameterization(tau_code,1,1);
  param[2]=dTheta3Dfit_parameterization(tau_code,1,2);
  param[3]=dTheta3Dfit_parameterization(tau_code,1,3);
  dTheta3Fun->SetParameter(1,dTheta3Dparam(1,P_tau,param));
  param[0]=dTheta3Dfit_parameterization(tau_code,2,0);
  param[1]=dTheta3Dfit_parameterization(tau_code,2,1);
  param[2]=dTheta3Dfit_parameterization(tau_code,2,2);
  param[3]=dTheta3Dfit_parameterization(tau_code,2,3);
  dTheta3Fun->SetParameter(2,dTheta3Dparam(2,P_tau,param));
  param[0]=dTheta3Dfit_parameterization(tau_code,3,0);
  param[1]=dTheta3Dfit_parameterization(tau_code,3,1);
  param[2]=dTheta3Dfit_parameterization(tau_code,3,2);
  param[3]=dTheta3Dfit_parameterization(tau_code,3,3);
  dTheta3Fun->SetParameter(3,dTheta3Dparam(3,P_tau,param));
  param[0]=dTheta3Dfit_parameterization(tau_code,4,0);
  param[1]=dTheta3Dfit_parameterization(tau_code,4,1);
  param[2]=dTheta3Dfit_parameterization(tau_code,4,2);
  param[3]=dTheta3Dfit_parameterization(tau_code,4,3);
  dTheta3Fun->SetParameter(4,dTheta3Dparam(4,P_tau,param));
  param[0]=dTheta3Dfit_parameterization(tau_code,5,0);
  param[1]=dTheta3Dfit_parameterization(tau_code,5,1);
  param[2]=dTheta3Dfit_parameterization(tau_code,5,2);
  param[3]=dTheta3Dfit_parameterization(tau_code,5,3);
  dTheta3Fun->SetParameter(5,dTheta3Dparam(5,P_tau,param));

  prob=dTheta3Fun->Eval(dTheta3d);

  if((prob>1.0 || prob<0.0) && fUseVerbose==1)
    {
      std::cout<<">>>> WARNING in MissingMassCalculator::dTheta3d_probability() <<<<"<<std::endl;
      std::cout<<"..... wrong probability="<<prob<<std::endl;
      std::cout<<"..... debugging: tau_type="<<tau_type<<"dTheta3d="<<dTheta3d<<"  P_tau="<<P_tau<<std::endl;
      std::cout<<"____________________________________________________________"<<std::endl;
      prob=1.0E-10;
    }

  delete dTheta3Fun;
  return prob;
}

//SpeedUp static instantation
double MissingMassCalculator::fit_param[3][6][4];

// returns dTheta3D probability based on ATLAS parameterization
double MissingMassCalculator::dTheta3d_probabilityFast(const int & tau_type,const double & dTheta3d,const  double & P_tau) {
  double prob=1.0E-10;
  int tau_code;
  if(tau_type<3) {
    tau_code=tau_type;
  }
  else if(tau_type==3) {
    tau_code=2;
  }
  else 
    {
      std::cout<<">>>> WARNING in MissingMassCalculator::dTheta3d_probability() <<<<"<<std::endl;
      std::cout<<"..... wrong tau_type="<<tau_type<<std::endl;
      std::cout<<"..... returning prob="<<prob<<std::endl;
      std::cout<<"____________________________________________________________"<<std::endl;
      return prob;
    }


  double myDelThetaParam[6];
  
  for (int i=0;i<6;++i) 
    {
      myDelThetaParam[i]=dTheta3Dparam(i,P_tau,fit_param[tau_code][i]);
    }
  double dTheta3dVal=dTheta3d;
  
  if (tau_type==0) prob=myDelThetaLepFunc(&dTheta3dVal,  myDelThetaParam);
  else prob=myDelThetaHadFunc(&dTheta3dVal,  myDelThetaParam);

  if( fUseVerbose==1 && (prob>1.0 || prob<0.0))
    {
      std::cout<<">>>> WARNING in MissingMassCalculator::dTheta3d_probability() <<<<"<<std::endl;
      std::cout<<"..... wrong probability="<<prob<<std::endl;
      std::cout<<"..... debugging: tau_type="<<tau_type<<"dTheta3d="<<dTheta3d<<"  P_tau="<<P_tau<<std::endl;
      std::cout<<"____________________________________________________________"<<std::endl;
      prob=1.0E-10;
    }

  return prob;
}


// ----- returns dTheta3D lower and upper boundaries:
// limit_code=0: 99% lower limit
// limit_code=1; 99% upper limit
// limit_code=2; 95% upper limit
double MissingMassCalculator::dTheta3DLimit(const int & tau_type,const  int & limit_code,const  double & P_tau) {
  double limit=1.0;
  if(limit_code==0) limit=0.0;
  double par[3]={0.0,0.0,0.0};
  // ---- leptonic tau's
  if(tau_type==0)
    {
      if(limit_code==0) // lower 99% limit
	{
	  par[0]=0.3342;
	  par[1]=-0.3376;
	  par[2]=-0.001377;
	}
      if(limit_code==1) // upper 99% limit
	{
	  par[0]=3.243;
	  par[1]=-12.87;
	  par[2]=0.009656;
	}
      if(limit_code==2) // upper 95% limit
	{
	  par[0]=2.927;
	  par[1]=-7.911;
	  par[2]=0.007783;
	}
    }
  // ---- 1-prong tau's
  if(tau_type==1)
    {
      if(limit_code==0) // lower 99% limit
	{
	  par[0]=0.2673;
	  par[1]=-14.8;
	  par[2]=-0.0004859;
	}
      if(limit_code==1) // upper 99% limit
	{
	  par[0]=9.341;
	  par[1]=-15.88;
	  par[2]=0.0333;
	}
      if(limit_code==2) // upper 95% limit
	{
	  par[0]=6.535;
	  par[1]=-8.649;
	  par[2]=0.00277;
	}
    }
  // ---- 3-prong tau's
  if(tau_type==3)
    {
      if(limit_code==0) // lower 99% limit
	{
	  par[0]=0.2308;
	  par[1]=-15.24;
	  par[2]=-0.0009458;
	}
      if(limit_code==1) // upper 99% limit
	{
	  par[0]=14.58;
	  par[1]=-6.043;
	  par[2]=-0.00928;
	}
      if(limit_code==2) // upper 95% limit
	{
	  par[0]=8.233;
	  par[1]=-0.3018;
	  par[2]=-0.009399;
	}
    }

  if(fabs(P_tau+par[1])>0.0) limit=par[0]/(P_tau+par[1])+par[2];
  if(limit_code==0){
    if (limit<0.0){
      limit=0.0;
    }
    else if ( limit>0.03 ){
      limit=0.03;
    }
  } else {
    if (limit<0.0 || limit>0.5*TMath::Pi()) {
      limit=0.5*TMath::Pi();
    }
    else if (limit<0.05 && limit>0.0) {
      limit=0.05; // parameterization only runs up to P~220 GeV in this regime will set an upper bound of 0.05
    }
  }
    

  return limit;
}

// ----- returns dTheta3D upper boundaries using visP_tau:
// limit_code=0; 99% upper limit
// limit_code=1; 95% upper limit
double MissingMassCalculator::dTheta3DLimitVis(int tau_type, int limit_code, double visP_tau) {
  double limit=0.4;
  double max_limit=0.4;
  if(limit_code==0 && tau_type==0) max_limit=0.2;
  if(limit_code==1 && tau_type==0) max_limit=0.15;
  if(limit_code==0 && tau_type>0) max_limit=0.4;
  if(limit_code==1 && tau_type>0) max_limit=0.19;
  double par[3]={0.0,0.0,0.0};

  // ---- leptonic tau's
  if(tau_type==0)
    {
      if(limit_code==0) // upper 99% limit
	{
	  par[0]=4.289;
	  par[1]=12.67;
	  par[2]=0.01996;
	}
      if(limit_code==1) // upper 95% limit
	{
	  par[0]=3.402;
	  par[1]=16.48;
	  par[2]=0.005309;
	}
    }
  // ---- 1-prong tau's
  if(tau_type==1)
    {
      if(limit_code==0) // upper 99% limit
	{
	  par[0]=7.36;
	  par[1]=-17.8;
	  par[2]=0.1051;
	}
      if(limit_code==1) // upper 95% limit
	{
	  par[0]=3.251;
	  par[1]=-18.02;
	  par[2]=0.04916;
	}
    }
  // ---- 3-prong tau's
  if(tau_type==3)
    {
      if(limit_code==0) // upper 99% limit
	{
	  par[0]=11.88;
	  par[1]=-8.584;
	  par[2]=0.06711;
	}
      if(limit_code==1) // upper 95% limit
	{
	  par[0]=5.048;
	  par[1]=-9.994;
	  par[2]=0.03248;
	}
    }
  
  double limit2=dTheta3DLimit(tau_type,limit_code+1,visP_tau);

  if(fabs(visP_tau+par[1])>0.0) limit=par[0]/(visP_tau+par[1])+par[2];
  if(limit>limit2) limit=limit2;
  if(limit<0.0) limit=0.2;
  if(limit>max_limit) limit=max_limit;
  return limit;
}

// Transfromation from dTheta to dPhi
double MissingMassCalculator::dTheta2dPhi(double etaL, double dTheta)
{
  double dPhi=dTheta;
  if(dTheta==0.0) return dPhi;
  double thetaL=2.0*atan(exp(-etaL));
  if(sin(thetaL)>0.0) dPhi=cos(dTheta)/pow(sin(thetaL),2)-pow(cos(thetaL)/sin(thetaL),2);
  if(fabs(dPhi)<=1.0) dPhi=fabs(acos(dPhi));
  else 
    {
      if(fUseVerbose==1) std::cout<<"Warning: cos(dPhi)>1, cos(dPhi)="<<dPhi<<" dTheta="<<dTheta<<" eta="<<etaL<<", ...returning 3.14"<<std::endl;
      return TMath::Pi();
    }
  return dPhi;
}


// simplistic MetProbability, to be later replaced by MetModel
double MissingMassCalculator::MetProbability(const double & metX,const  double & metY,const  double & MetSigma) {
  double metprob=1.0;
  if(MetSigma>0.0)
    {
      //SpeedUp
      metprob=exp(-0.5*(metX*metX+metY*metY)/(MetSigma*MetSigma));      
    }
  return metprob;
}
// simplistic MetProbability, needs  MetSigma1 and  MetSigma2 from MetModel
double MissingMassCalculator::MetProbability(const double & met1,const  double & met2,const  double & MetSigma1,const  double & MetSigma2) {

  if(MetSigma1>1.0 && MetSigma2>1.0) // it doesn't make sense if MET resolution sigma is <1 GeV
    {
      //SpeedUp
      return exp(-0.5*(met1*met1/(MetSigma1*MetSigma1)+met2*met2/(MetSigma2*MetSigma2)))/(MetSigma1*MetSigma2*2*TMath::Pi());
    }
  else 
    {
      if(fUseVerbose==1) std::cout<<"Warning!!! MissingMassCalculator::MetProbability: either MetSigma1 or MetSigma2 are <1 GeV--- too low, returning prob=1"<<std::endl;
      return 1.;
    }
  
}


// 05/19/11: modified V7. Replaced dR-probability by dTheta3D-probability
//                        Di-lepton scan is the same as for lep-had and had-had channels, there are 2 additional Mnu-scans.
// 06/01/11: includes David's improvements
int MissingMassCalculator::DitauMassCalculatorV9(const TLorentzVector & tau_vec1, const int & tau_type1, 
						 const TLorentzVector & tau_vec2, const int & tau_type2, 
						 const TVector2 & met_vec) {

  int fit_code=0; // 0==bad, 1==good
  ClearDitauStuff(fDitauStuffFit);
  ClearDitauStuff(fDitauStuffHisto);
     
  //------- Settings -------------------------------
  int Niter=Niter_fit1; // number of points for each dR loop
  int NiterMET=Niter_fit2; // number of iterations for each MET scan loop  
  int NiterMnu=Niter_fit3; // number of iterations for Mnu loop
  double Mtau=1.777;
  double Mnu_binSize=Mtau/NiterMnu;

  double METresX=InputInfo.METsigmaL; // MET resolution in direction parallel to leading jet, for MET scan
  double METresY=InputInfo.METsigmaP; // MET resolution in direction perpendicular to leading jet, for MET scan
  double N_METsigma=Nsigma_METscan; // number of sigmas for MET scan
  double METresX_binSize=2*N_METsigma*METresX/NiterMET;
  double METresY_binSize=2*N_METsigma*METresY/NiterMET;

  //-------- end of Settings
  int Leptau1=0; // lepton-tau code
  int Leptau2=0; 
  if(tau_type1==0) Leptau1=1;
  if(tau_type2==0) Leptau2=1;      
  
  double dPhi_nu1[2];
  double dPhi_nu2[2];  
  double dPhi_tmp;
  int phi1_loop=0;
  int phi2_loop=0;
  double phi1=0.0;
  double phi2=0.0;
  int solution=0;
  
  TLorentzVector nu1_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector nu2_tmp(0.0,0.0,0.0,0.0);
  //  TLorentzVector tau1_tmp(0.0,0.0,0.0,0.0);
  // TLorentzVector tau2_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector tautau_tmp(0.0,0.0,0.0,0.0);
  TVector2 metvec_tmp(0.0,0.0);
  
  double tauProb1r=1.0;
  double tauProb2r=1.0;			  
  double metprob=1.0;
  double sign_tmp=0.0; 
  double totalProb=0.0;
  double ditauProb=1.0;

  prob_tmp=0.0;
  
  double met_smear_x=0.0;
  double met_smear_y=0.0;
  double met_smearL=0.0;
  double met_smearP=0.0;
  
  double angle1=0.0;
  double angle2=0.0;
  //--- define histograms for histogram method
  //--- upper limits need to be revisied in the future!!! It may be not enough for some analyses

  fMfit_all->Reset();
  fPXfit_nu1->Reset();
  fPYfit_nu1->Reset();
  fPZfit_nu1->Reset();
  fPXfit_nu2->Reset();
  fPYfit_nu2->Reset();
  fPZfit_nu2->Reset();


  //SpeedUp
    TStopwatch timer;
  if (fSpeedStudy) timer.Start();
  int iter0=0;
  iter1=0;
  iter2=0;
  iter3=0;
  iter4=0;
  const double cosphi_jet=cos(InputInfo.phi_jet);
  const double sinphi_jet=sin(InputInfo.phi_jet);
  const double tau_vec1phi=tau_vec1.Phi();
  const double tau_vec2phi=tau_vec2.Phi();
  const double tau_vec1m=tau_vec1.M();
  const double tau_vec2m=tau_vec2.M();

  int nsol1;
  int nsol2;

  iang1low=0;
  iang1high=0;
  iang2low=0;
  iang2high=0;

  std::vector<TLorentzVector> tauvecsol1(4);
  std::vector<TLorentzVector> tauvecsol2(4);
  std::vector<double> tauvecprob1(4);
  std::vector<double> tauvecprob2(4);
  

  //const double threepi=3*TMath::Pi();
  //const double twopi=2*TMath::Pi();
  //const double onepi=1*TMath::Pi();

  double Mvis=(tau_vec1+tau_vec2).M();
  TLorentzVector met4vec(0.0,0.0,0.0,0.0);
  met4vec.SetPxPyPzE(met_vec.X(),met_vec.Y(),0.0,met_vec.Mod());
  double Meff=(tau_vec1+tau_vec2+met4vec).M();
  
  //---------------------------------------------    
  if((Leptau1+Leptau2)==2) // both tau's are leptonic V9
    {
      //------- Settings -------------------------------	  
      double dPhi_max=dRmax_tau;
      double dPhi_binSize=dPhi_max/Niter;	
      double dMnu_max1=Mtau-tau_vec1m; 
      double dMnu_max2=Mtau-tau_vec2m;
      double Mnu_binSize1=dMnu_max1/NiterMnu;
      double Mnu_binSize2=dMnu_max2/NiterMnu;
      //-------- end of Settings	  	  
      std::vector<TLorentzVector> nuvec1_tmp;
      std::vector<TLorentzVector> nuvec2_tmp;
      double M1=Mtau;
      double M2=Mtau;
      double M_nu1=0.0;
      double M_nu2=0.0;
      //---------------------------------------------
      for(int i1=0; i1<Niter; i1++) //---- loop-1: dPhi for nu1
	{
	  dPhi_tmp=dPhi_binSize*i1;
	  dPhi_nu1[0]=-dPhi_tmp;
	  dPhi_nu1[1]=dPhi_tmp;
	  phi1_loop=i1>0 ? 2 : 1;
	  //---------------------------------------------
	  for(int i2=0; i2<Niter; i2++) //---- loop-2: dPhi for nu2
	    {
	      dPhi_tmp=dPhi_binSize*i2;
	      dPhi_nu2[0]=-dPhi_tmp;
	      dPhi_nu2[1]=dPhi_tmp;
	      phi2_loop=i2>0 ? 2 : 1;		      
	      for(int ij1=0; ij1<phi1_loop; ij1++)
		{
		  phi1=tau_vec1phi+dPhi_nu1[ij1];
		  for(int ij2=0; ij2<phi2_loop; ij2++)
		    {		      
		      phi2=tau_vec2phi+dPhi_nu2[ij2];
		      for(int i3=0; i3<NiterMET+1; i3++) // MET_X scan
			{
			  met_smearL=METresX_binSize*i3-N_METsigma*METresX;
// 			  for(int i4=0; i4<NiterMET+1; i4++) // MET_Y scan
			  for(int i4=0; i4<NiterMET+1 && (i4*i4+i3*i3)<(NiterMET+1)*(NiterMET+1); i4++) // MET_Y scan			    
			    {
			      met_smearP=METresY_binSize*i4-N_METsigma*METresY;
			      //SpeedUp 
			      met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
			      met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
			      metprob=MetProbability(met_smearL,met_smearP,METresX,METresY);
			      if (metprob<=0) continue;
			      
			      metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);

			      for(int i5=0; i5<NiterMnu; i5++) //---- Mnu1 scan 
				{
				  M_nu1=Mnu_binSize1*i5;
				  if(M_nu1>=(Mtau-tau_vec1m)) continue; // checking condition: M_nu < M_tau-M_lep
				  M1=sqrt(Mtau*Mtau-M_nu1*M_nu1);
				  for(int i6=0; i6<NiterMnu; i6++) //---- Mnu2 scan 
				    {
				      M_nu2=Mnu_binSize2*i6;
				      if(M_nu2>=(Mtau-tau_vec2m)) continue; // checking condition: M_nu < M_tau-M_lep
				      M2=sqrt(Mtau*Mtau-M_nu2*M_nu2);
				      //nuvec1_tmp.clear();
				      //nuvec2_tmp.clear();
				      solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nsol1,nsol2);
				      ++iter0;
				      if(solution!=1) continue;

				      ++iter1;				      
				      iter2+=(nsol1*nsol2);
				      //SpeedUp no nested loop to compute individual probability
				      int ngoodsol1=0;
				      for(int j1=0; j1<nsol1; j1++)
					{
					  //SpeedUp reference to avoid picking again in again in array
					  TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
					  TLorentzVector & tauvecsol1j=tauvecsol1[j1];
					  double &tauvecprob1j = tauvecprob1[j1];
					  tauvecprob1j=0.;
					      
					      
					  //SpeedUp sum energy rather than TLV
					  if(tau_vec1.E()+nuvec1_tmpj.E()>=beamEnergy) continue;
					  //SpeedUp use SetXYZM to reassign the mass
					  nuvec1_tmpj.SetXYZM(nuvec1_tmpj.Px(),
							      nuvec1_tmpj.Py(),
							      nuvec1_tmpj.Pz(),
							      M_nu1);

					  tauvecsol1j.SetPxPyPzE(0.,0.,0.,0.);
					  tauvecsol1j+=nuvec1_tmpj;
					  tauvecsol1j+=tau_vec1;

					  angle1=Angle(nuvec1_tmpj,tau_vec1);
					  const double tau1_tmpp=tauvecsol1j.P();
					  //SpeedUp : the low limit is not worth computing (not sure)
					  if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)){++iang1low ; continue;} // lower 99% bound
					  if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)){++iang1high ; continue;} // upper 99% bound
					  tauvecprob1j=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp);  						  
					  ++ngoodsol1;
					  
					  if (fSpeedStudy)
					    {
					      ++iter3;
					      if (iter3 % 10000 == 1) {
						const double aux=dTheta3d_probability(tau_type1,angle1,tau1_tmpp);
						const double ratio=std::abs(aux-tauvecprob1j)/(aux+tauvecprob1j);
						if  (ratio>=1E-3){
						  std::cout << "SpeedUp WARNING mismatch between quick computation and slow "  << iter3 << " ref " << tauProb1r << " check-1 " << aux << std::endl;
							} 
					      } 
					    }
					}

				      if (ngoodsol1==0) continue;
				      int ngoodsol2=0;

				      for(int j2=0; j2<nsol2; j2++)
					{
					  TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
					  TLorentzVector & tauvecsol2j=tauvecsol2[j2];
					  double &tauvecprob2j = tauvecprob2[j2];
					  tauvecprob2j=0.;

					  //SpeedUp
					  nuvec2_tmpj.SetXYZM(nuvec2_tmpj.Px(),
							      nuvec2_tmpj.Py(),
							      nuvec2_tmpj.Pz(),
							      M_nu2);

					  //SpeedUp sum energy rather than TLV
					  if(tau_vec2.E()+nuvec2_tmpj.E()>=beamEnergy) continue;
					  //redundant if(tau_vec1.E()+nuvec1_tmpj.E()+tau_vec2.E()+nuvec2_tmpj.E()>=2.0*beamEnergy) continue;
					  //use += TLV rather than TLV1+TLV2
					  tauvecsol2j.SetPxPyPzE(0.,0.,0.,0.);
					  tauvecsol2j+=nuvec2_tmpj;
					  tauvecsol2j+=tau_vec2;
					  angle2=Angle(nuvec2_tmpj,tau_vec2);
					  const double tau2_tmpp=tauvecsol2j.P();
					  if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)){++iang2low; continue;} // lower 99% bound
					  if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)){++iang2high; continue;} // upper 99% bound
					  tauvecprob2j=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp); 
					  ++ngoodsol2;
					  if (fSpeedStudy)
					    {
					      ++iter3;
					      if (iter3 % 10000 == 1) {
						const double aux2=dTheta3d_probability(tau_type2,angle2,tau2_tmpp);
						const double ratio2=std::abs(aux2-tauvecprob2j)/(aux2+tauvecprob2j);
						if  (ratio2>=1E-3){
						  std::cout << "SpeedUp WARNING mismatch between quick computation and slow "  << iter3 << " ref " << tauProb2r << " check-2 " << aux2 << std::endl;
						} 
					      }
					    }
					}
				      if (ngoodsol2==0) continue;					      

				      // now reloop
				      for (int j1=0; j1<nsol1;++j1)
					{
					  double &tauvecprob1j = tauvecprob1[j1];
					  if (tauvecprob1j==0.) continue;
					  TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
					  TLorentzVector & tauvecsol1j=tauvecsol1[j1];
					  
					  for (int j2=0; j2<nsol2;++j2)
					    {
					      double &tauvecprob2j = tauvecprob2[j2];
					      if (tauvecprob2j==0.) continue;
					      TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
					      TLorentzVector & tauvecsol2j=tauvecsol2[j2];
					      ++iter4;
					      fit_code=1; // at least one solution is found

					      tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
					      tautau_tmp+=tauvecsol1j;
					      tautau_tmp+=tauvecsol2j;
						  
					      const double mtautau=tautau_tmp.M();
					      if(TailCleanUp(tau_type1,tau_vec1,nuvec1_tmpj, 
							     tau_type2,tau_vec2,nuvec2_tmpj, 
							     mtautau,Mvis,Meff,InputInfo.DelPhiTT)==0) continue;	

					      ditauProb=TauProbability(tau_type1,tau_vec1,nuvec1_tmpj, 
								       tau_type2,tau_vec2,nuvec2_tmpj);
					      totalProb=tauvecprob1j*tauvecprob2j*metprob*ditauProb;
					      // totalProb=tauvecprob1j*tauvecprob2j*metprob;
					      // cannot happen if(totalProb<=0.0) continue;

						  
					      fMfit_all->Fill(mtautau,totalProb);
					      fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
					      fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
					      fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
					      fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
					      fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
					      fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
					      if(totalProb>prob_tmp) // fill solution with highest probability                                   
						{
						  sign_tmp=-log10(totalProb);
						  prob_tmp=totalProb;
						  fDitauStuffFit.Mditau_best=mtautau; 
						  fDitauStuffFit.Sign_best=sign_tmp; 
						  fDitauStuffFit.nutau1=nuvec1_tmpj; 
						  fDitauStuffFit.nutau2=nuvec2_tmpj; 
						}
					    }
					}

				    }
				}
			    }
			}
		    }
		}
	    }
	}	  
    }

  int isol=0;

  
  if((Leptau1+Leptau2)==1) // one leptonic tau V9
    {
      //------- Settings -------------------------------	  
      double dPhi_max=dRmax_tau;
      double dPhi_binSize=dPhi_max/Niter;
      double dMnu_max= (Leptau1==1) ? Mtau-tau_vec1m : Mtau-tau_vec2m;
      Mnu_binSize=dMnu_max/NiterMnu;
      //-------- end of Settings	  	  
      std::vector<TLorentzVector> nuvec1_tmp;
      std::vector<TLorentzVector> nuvec2_tmp;
      double M1=Mtau;
      double M2=Mtau;
      double M_nu=0.0;
//       double MnuProb=1.0;

//----- Stuff below are for Winter 2012 lep-had analysis only; it has to be replaced by a more common scheme once over channels are optimized
      double met_det=met_vec.Mod();
      TVector2 mht_vec((tau_vec1+tau_vec2).Px(),(tau_vec1+tau_vec2).Py()); // missing Ht vector for Njet25=0 events
      const double mht=mht_vec.Mod(); 
      double input_metX;
      double input_metY;
      double mht_offset=0.0;
      if(InputInfo.Njet25==0)  // use missing Ht for Njet25=0 events
	{
	  input_metX=-mht_vec.X();
	  input_metY=-mht_vec.Y();
	  // dPhi(l-t) dependence of misHt-trueMET
	  if(met_det<20.0) mht_offset=InputInfo.DelPhiTT>1.8 ? -161.9+235.5*InputInfo.DelPhiTT-101.6*pow(InputInfo.DelPhiTT,2)+13.57*pow(InputInfo.DelPhiTT,3) : 12.0;
	  else mht_offset=115.5-78.1*InputInfo.DelPhiTT+12.83*pow(InputInfo.DelPhiTT,2);
	}
      else // use MET for Njet25>1 events
	{
	  input_metX=met_vec.X();
	  input_metY=met_vec.Y();	  
	}

      //---------------------------------------------
      for(int i1=0; i1<Niter; i1++) //---- loop-1: dPhi for nu1
	{
	  dPhi_tmp=dPhi_binSize*i1;
	  dPhi_nu1[0]=-dPhi_tmp;
	  dPhi_nu1[1]=dPhi_tmp;
	  phi1_loop=i1>0 ? 2 : 1;
	  //---------------------------------------------
	  for(int i2=0; i2<Niter; i2++) //---- loop-2: dPhi for nu2
	    {
	      dPhi_tmp=dPhi_binSize*i2;
	      dPhi_nu2[0]=-dPhi_tmp;
	      dPhi_nu2[1]=dPhi_tmp;
	      phi2_loop=i2>0 ? 2 : 1;
	      for(int i3=0; i3<NiterMnu; i3++) //---- loop-3: virtual neutrino mass 
		{
		  M_nu=Mnu_binSize*i3;
		  if(Leptau1==1){
		    if (M_nu>=(Mtau-tau_vec1m)) continue; // checking condition: M_nu < M_tau-M_lep
		    M1=sqrt(Mtau*Mtau-M_nu*M_nu);
		  }
		  else {
		    if (M_nu>=(Mtau-tau_vec2m)) continue; // checking condition: M_nu < M_tau-M_lep
		    M2=sqrt(Mtau*Mtau-M_nu*M_nu);
		  }
// 		  MnuProb=MnuProbability(M_nu,Mnu_binSize); // Mnu probability
		  for(int ij1=0; ij1<phi1_loop; ij1++)
		    {
		      phi1=tau_vec1phi+dPhi_nu1[ij1];
		      for(int ij2=0; ij2<phi2_loop; ij2++)
			{		      
			  phi2=tau_vec2phi+dPhi_nu2[ij2];
			  for(int i4=0; i4<NiterMET+1; i4++) // MET_X scan
			    {
			      met_smearL=METresX_binSize*i4-N_METsigma*METresX;
			      for(int i5=0; i5<NiterMET+1 && (i4*i4+i5*i5)<(NiterMET+1)*(NiterMET+1); i5++) // MET_Y scan
// 			      for(int i5=0; i5<NiterMET+1; i5++) // MET_Y scan
				{
				  met_smearP=METresY_binSize*i5-N_METsigma*METresY;
				  met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
				  met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
// 				  metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);
				  metvec_tmp.Set(input_metX+met_smear_x,input_metY+met_smear_y);
				  
				  //nuvec1_tmp.clear();
				  //nuvec2_tmp.clear();
				  solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nsol1,nsol2);
				  ++iter0;					  
				  
				  if(solution!=1) continue;
				  ++iter1;					  

				  if (fSpeedStudy)
				    {
				      ++isol;
				      if (isol % 10000==1)
					{
					  const double aux=nuvecsol1[0].E();
					  int solutionaux=NuPsolutionV2(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nuvec1_tmp,nuvec2_tmp);				      
					  if (solutionaux!=solution){
					    std::cout << "SpeedUp WARNING mismatch between quick NuPsolutionV2 computation and slow "  << isol << " ref " << solutionaux << " check " << solution << std::endl;
					  }
					  else {
					    
					    const double ratio=std::abs(aux-nuvec1_tmp[0].E())/(aux+nuvec1_tmp[0].E());
					    if  (ratio>=1E-3){
					      std::cout << "SpeedUp WARNING mismatch between quick NuPsolutionV2 computation and slow "  << isol << " ratio " << ratio << std::endl;
					    } 
					  }
					  
					}
				    }	

				  iter2+=nsol1*nsol2;					  
				  //SpeedUp no nested loop to compute individual probability
				  int ngoodsol1=0;


				  for(int j1=0; j1<nsol1; ++j1)
				    {
				      //SpeedUp reference to avoid picking again in again in array
				      TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
				      TLorentzVector & tauvecsol1j=tauvecsol1[j1];
				      double &tauvecprob1j = tauvecprob1[j1];
				      tauvecprob1j=0.;
				      
				      if(tau_vec1.E()+nuvec1_tmpj.E()>=beamEnergy) continue;
				      if(Leptau1==1)
					{
					  //SpeedUp
					  nuvec1_tmpj.SetXYZM(nuvec1_tmpj.Px(),
							      nuvec1_tmpj.Py(),
							      nuvec1_tmpj.Pz(),
							      M_nu);
					  
					}
				      
				      tauvecsol1j.SetPxPyPzE(0.,0.,0.,0.);
				      tauvecsol1j+=nuvec1_tmpj;
				      tauvecsol1j+=tau_vec1;
				      const double tau1_tmpp=tauvecsol1j.P();
				      
				      angle1=Angle(nuvec1_tmpj,tau_vec1);
				      
				      if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)) {++iang1low ; continue;} // lower 99% bound
				      if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)) {++iang1high ; continue;} // upper 99% bound
				      tauvecprob1j=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp); 
				      ++ngoodsol1;
				      if (fSpeedStudy)
					{
					  ++iter3;
					  if (iter3 % 10000 == 1) {
					    const double aux=dTheta3d_probability(tau_type1,angle1,tau1_tmpp);
					    const double ratio=std::abs(aux-tauvecprob1j)/(aux+tauvecprob1j);
					    if  (ratio>=1E-3){
					      std::cout << "SpeedUp WARNING mismatch between quick dTheta3D_probability computation and slow "  << iter3 << " ref " << tauvecprob1j << " check-1 " << aux << std::endl;
					    } 
					  }
					}
				    }
				  
				  if (ngoodsol1==0) continue;
				  int ngoodsol2=0;				      
				  for(int j2=0; j2<nsol2; ++j2)
				    {
				      TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];	
				      TLorentzVector & tauvecsol2j=tauvecsol2[j2];
				      double &tauvecprob2j = tauvecprob2[j2];
				      tauvecprob2j=0.;
				      
				      
				      if(Leptau2==1)
					{
					  //SpeedUp
					  nuvec2_tmpj.SetXYZM(nuvec2_tmpj.Px(),
							      nuvec2_tmpj.Py(),
							      nuvec2_tmpj.Pz(),
							      M_nu);
					  
					}						      
				      
				      
				      if(tau_vec2.E()+nuvec2_tmpj.E()>=beamEnergy) continue;
				      //redundantif(tau_vec1.E()+nuvec1_tmp[j1].E()+tau_vec2.E()+nuvec2_tmp[j2].E()>=2.0*beamEnergy) continue;
				      // 					      if(CheckSolutions(nuvec1_tmp[j1],tau_vec1,Leptau1)<1) continue;
				      // 					      if(CheckSolutions(nuvec2_tmp[j2],tau_vec2,Leptau2)<1) continue;
				      //SpeedUp use += TLV rather than TLV1+TLV2
				      //tau1_tmp=nuvec1_tmp[j1]+tau_vec1;
				      tauvecsol2j.SetPxPyPzE(0.,0.,0.,0.);
				      tauvecsol2j+=nuvec2_tmpj;
				      tauvecsol2j+=tau_vec2;
				      angle2=Angle(nuvec2_tmpj,tau_vec2);
				      const double tau2_tmpp=tauvecsol2j.P();
				      if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)) { ++iang2low ; continue;} // lower 99% bound
				      if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)) { ++iang2high ; continue;} // upper 99% bound
				      // 					      tauProb1r=dTheta3d_probability(tau_type1,angle1,tau1_tmpp);
				      // 					      tauProb2r=dTheta3d_probability(tau_type2,angle2,tau2_tmpp);
				      tauvecprob2j=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp); 
				      ++ngoodsol2;
				      if (fSpeedStudy)
					{
					  ++iter3;
					  if (iter3 % 10000 == 1) {
					    const double aux2=dTheta3d_probability(tau_type2,angle2,tau2_tmpp);
					    const double ratio2=std::abs(aux2-tauvecprob2j)/(aux2+tauvecprob2j);
					    if  (ratio2>=1E-3){
					      std::cout << "SpeedUp WARNING mismatch between quick dTheta3D_probability computation and slow "  << iter3 << " ref " << tauvecprob2j << " check-2 " << aux2 << std::endl;
					    } 
					  }
					}
				      
				      
				    }
				  if (ngoodsol2==0) continue;
				  if(InputInfo.Njet25==0) metprob=MHtProbability(met_smearL,met_smearP,mht,metvec_tmp.Mod(),mht_offset); // for lep-had Winter 2012 analysis
				  else metprob=MetProbability(met_smearL,met_smearP,METresX,METresY); 
				  if (metprob<=0) continue;
				  
				  // now reloop
				  for (int j1=0; j1<nsol1;++j1)
				    {
				      double &tauvecprob1j = tauvecprob1[j1];
				      if (tauvecprob1j==0.) continue;
				      TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
				      TLorentzVector & tauvecsol1j=tauvecsol1[j1];
				      
				      for (int j2=0; j2<nsol2;++j2)
					{
					  double &tauvecprob2j = tauvecprob2[j2];
					  if (tauvecprob2j==0.) continue;
					  TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
					  TLorentzVector & tauvecsol2j=tauvecsol2[j2];
					  					  
					  ++iter4;
					  					  
					  tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
					  tautau_tmp+=tauvecsol1j;
					  tautau_tmp+=tauvecsol2j;
					  const double mtautau=tautau_tmp.M();

					  if(TailCleanUp(tau_type1,tau_vec1,nuvec1_tmpj, 
							 tau_type2,tau_vec2,nuvec2_tmpj, 
							 mtautau,Mvis,Meff,InputInfo.DelPhiTT)==0) continue;	

					  if(InputInfo.Njet25>0) ditauProb=TauProbability(tau_type1,tau_vec1,nuvec1_tmpj,tau_type2,tau_vec2,nuvec2_tmpj);
					  else ditauProb=TauProbability(tau_type1,tau_vec1,nuvec1_tmpj,tau_type2,tau_vec2,nuvec2_tmpj,met_det); // customized prob for Njet25=0
					  totalProb=tauvecprob1j*tauvecprob2j*metprob*ditauProb;
// 					  //cannot happen if(totalProb<=0.0) continue;				  

					  fit_code=1; // at least one solution is found

					  fMfit_all->Fill(mtautau,totalProb);
					  fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
					  fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
					  fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
					  fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
					  fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
					  fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
					  if(totalProb>prob_tmp) // fill solution with highest probability                                   
					    {
					      sign_tmp=-log10(totalProb);
					      prob_tmp=totalProb;
					      fDitauStuffFit.Mditau_best=mtautau; 
					      fDitauStuffFit.Sign_best=sign_tmp; 
					      fDitauStuffFit.nutau1=nuvec1_tmpj; 
					      fDitauStuffFit.nutau2=nuvec2_tmpj; 
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}	  
    }
  
  
  if((Leptau1+Leptau2)==0) // both tau's are hadronic V9
    {
      //------- Settings -------------------------------	  
      double dPhi_max=dRmax_tau;
      double dPhi_binSize=dPhi_max/Niter;	
      //-------- end of Settings	  	  
      //      std::vector<TLorentzVector> nuvec1_tmp;
      //std::vector<TLorentzVector> nuvec2_tmp;
      
      //---------------------------------------------
      for(int i1=0; i1<Niter; i1++) //---- loop-1: dPhi for nu1
	{
	  dPhi_tmp=dPhi_binSize*i1;
	  dPhi_nu1[0]=-dPhi_tmp;
	  dPhi_nu1[1]=dPhi_tmp;
	  phi1_loop=i1>0 ? 2 : 1;
	  //---------------------------------------------
	  for(int i2=0; i2<Niter; i2++) //---- loop-2: dPhi for nu2
	    {
	      dPhi_tmp=dPhi_binSize*i2;
	      dPhi_nu2[0]=-dPhi_tmp;
	      dPhi_nu2[1]=dPhi_tmp;
	      phi2_loop=i2>0 ? 2 : 1;		      
	      for(int ij1=0; ij1<phi1_loop; ij1++)
		{
		  phi1=tau_vec1phi+dPhi_nu1[ij1];
		  for(int ij2=0; ij2<phi2_loop; ij2++)
		    {		      
		      phi2=tau_vec2phi+dPhi_nu2[ij2];
		      for(int i3=0; i3<NiterMET+1; i3++) // MET_X scan
			{
			  met_smearL=METresX_binSize*i3-N_METsigma*METresX;
			  for(int i4=0; i4<NiterMET+1 && (i4*i4+i3*i3)<(NiterMET+1)*(NiterMET+1); i4++) // MET_Y scan
// 			  for(int i4=0; i4<NiterMET+1; i4++) // MET_Y scan
			    {
			      met_smearP=METresY_binSize*i4-N_METsigma*METresY;
			      met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
			      met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
			      //SpeedUp
			      if (metprob<=0) continue;
			      metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);
			      //nuvec1_tmp.clear();
			      //nuvec2_tmp.clear();
			      solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,Mtau,Mtau,phi1,phi2,nsol1,nsol2);
			      ++iter0;
			      if(solution!=1) continue; // there is solution

			      ++iter1;				      
			      iter2+=(nsol1*nsol2);
			      int ngoodsol1=0;
			      for(int j1=0; j1<nsol1; j1++)
				{
				  //SpeedUp reference to avoid picking again in again in array
				  TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
				  TLorentzVector & tauvecsol1j=tauvecsol1[j1];
				  double &tauvecprob1j = tauvecprob1[j1];
				  tauvecprob1j=0.;
					      

				  if((tau_vec1+nuvec1_tmpj).E()>=beamEnergy) continue;

				  tauvecsol1j.SetPxPyPzE(0.,0.,0.,0.);
				  tauvecsol1j+=nuvec1_tmpj;
				  tauvecsol1j+=tau_vec1;
				  const double tau1_tmpp=tauvecsol1j.P();
				  angle1=Angle(nuvec1_tmpj,tau_vec1);
				  if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)){++iang1low; continue; }// lower 99% bound
				  if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)){++iang1high; continue;} // upper 99% bound
				  tauvecprob1j=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp); 
				  ++ngoodsol1;


				  if (fSpeedStudy)
				    {
				      ++iter3;
				      if (iter3 % 10000 == 1) {
					const double aux=dTheta3d_probability(tau_type1,angle1,tau1_tmpp);
					const double ratio=std::abs(aux-tauvecprob1j)/(aux+tauvecprob1j);
					if  (ratio>=1E-3){
					  std::cout << "SpeedUp WARNING mismatch between quick computation and slow "  << iter3 << " ref-1 " << tauProb1r << " check " << aux << std::endl;
					} 
				      }
				    }
				}
			      

			      if (ngoodsol1==0) continue;			      
			      int ngoodsol2=0;

			      for(int j2=0; j2<nsol2; j2++)
				{
				  TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];	
				  TLorentzVector & tauvecsol2j=tauvecsol2[j2];
				  double &tauvecprob2j = tauvecprob2[j2];
				  tauvecprob2j=0.;

				  tauvecsol2j.SetPxPyPzE(0.,0.,0.,0.);
				  tauvecsol2j+=nuvec2_tmpj;
				  tauvecsol2j+=tau_vec2;
				  if (tauvecsol2j.E()>=beamEnergy) continue;
				  const double tau2_tmpp=tauvecsol2j.P();
				  //redundant if(tau_vec1.E()+nuvec1_tmpj.E()+tau_vec2.E()+nuvec2_tmpj.E()>=2.0*beamEnergy) continue;
				  angle2=Angle(nuvec2_tmpj,tau_vec2);
				  if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)){++iang2low; continue;} // lower 99% bound
				  if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)){++iang2high; continue;} // upper 99% bound
				  tauvecprob2j=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp);
				  ++ngoodsol2;

				  if (fSpeedStudy)
				    {
				      ++iter3;
				      if (iter3 % 10000 == 1) {
					const double aux2=dTheta3d_probability(tau_type2,angle2,tau2_tmpp);
					const double ratio2=std::abs(aux2-tauvecprob2j)/(aux2+tauvecprob2j);
					if  (ratio2>=1E-3){
					  std::cout << "SpeedUp WARNING mismatch between quick computation and slow "  << iter3 << " ref-2 " << tauProb2r << " check " << aux2 << std::endl;
					} 
				      }
				    }
				}
			      if (ngoodsol2==0) continue;			      
			      metprob=MetProbability(met_smearL,met_smearP,METresX,METresY);
			      
			      // now reloop
			      for (int j1=0; j1<nsol1;++j1)
				{
				  double &tauvecprob1j = tauvecprob1[j1];
				  if (tauvecprob1j==0.) continue;
				  TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
				  TLorentzVector & tauvecsol1j=tauvecsol1[j1];
					  
				  for (int j2=0; j2<nsol2;++j2)
				    {
				      double &tauvecprob2j = tauvecprob2[j2];
				      if (tauvecprob2j==0.) continue;
				      TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
				      TLorentzVector & tauvecsol2j=tauvecsol2[j2];

				      ++iter4;
				      
				      tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
				      tautau_tmp+=tauvecsol1j;
				      tautau_tmp+=tauvecsol2j;
				      const double mtautau=tautau_tmp.M();

				      if(TailCleanUp(tau_type1,tau_vec1,nuvec1_tmpj, 
						     tau_type2,tau_vec2,nuvec2_tmpj, 
						     mtautau,Mvis,Meff,InputInfo.DelPhiTT)==0) continue;

				      ditauProb=TauProbability(tau_type1,tau_vec1,nuvec1_tmpj, 
							       tau_type2,tau_vec2,nuvec2_tmpj);
				      totalProb=tauvecprob1j*tauvecprob2j*metprob*ditauProb;
				      //cannot happen if(totalProb<=0.0) continue;

				      fit_code=1; // at least one solution is found
				      sign_tmp=-log10(totalProb);
				      fMfit_all->Fill(mtautau,totalProb);
				      fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
				      fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
				      fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
				      fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
				      fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
				      fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
				      if(totalProb>prob_tmp) // fill solution with highest probability                                   
					{
					  prob_tmp=totalProb;
					  fDitauStuffFit.Mditau_best=mtautau; 
					  fDitauStuffFit.Sign_best=sign_tmp; 
					  fDitauStuffFit.nutau1=nuvec1_tmpj; 
					  fDitauStuffFit.nutau2=nuvec2_tmpj; 
					}
				    }
				}

			    }
			}
		    }
		}
	    }
	}	  
    }

  //SpeedUp
  if (fSpeedStudy) timer.Stop()  ;
  
  if (fSpeedStudy){
    std::cout << "SpeedUp niters=" << iter0 << " " << iter1 
	      << " " << iter2 <<" " << iter3 <<" " << iter4
              << "skip:" << iang1low << " " << iang1high << " " << iang2low << " " << iang2high 
	      << " CPU time = " << timer.CpuTime() << "s real = " << timer.RealTime() << " s " << std::endl;
  }
  

  int max_bin;
  double Px1, Py1, Pz1;
  double Px2, Py2, Pz2;
  if(fMfit_all->GetEntries()>0)
    {
      fMfit_all->Smooth();
      fPXfit_nu1->Smooth();
      fPYfit_nu1->Smooth();
      fPZfit_nu1->Smooth();
      fPXfit_nu2->Smooth();
      fPYfit_nu2->Smooth();
      fPZfit_nu2->Smooth();

      max_bin=fMfit_all->GetMaximumBin();
      fDitauStuffHisto.Mditau_best=fMfit_all->GetBinCenter(max_bin);
      //      if(fDitauStuffHisto.Mditau_best<1.0) return 0; // no solution if mass is too small
      double prob_hist=1.0*(fMfit_all->GetBinContent(max_bin))/(1.0*fMfit_all->GetEntries());
      if(prob_hist>=1.0) prob_hist=1.0;
      if(prob_hist>0.0) fDitauStuffHisto.Sign_best=-log10(prob_hist);
      else fDitauStuffHisto.Sign_best=-1.0;
      if(fDitauStuffHisto.Mditau_best>0.0) fDitauStuffHisto.RMSoverMPV=fMfit_all->GetRMS()/fDitauStuffHisto.Mditau_best;
      //---- getting Nu1
      max_bin=fPXfit_nu1->GetMaximumBin();
      Px1=fPXfit_nu1->GetBinCenter(max_bin);
      max_bin=fPYfit_nu1->GetMaximumBin();
      Py1=fPYfit_nu1->GetBinCenter(max_bin);
      max_bin=fPZfit_nu1->GetMaximumBin();
      Pz1=fPZfit_nu1->GetBinCenter(max_bin);
      //---- getting Nu2
      max_bin=fPXfit_nu2->GetMaximumBin();
      Px2=fPXfit_nu2->GetBinCenter(max_bin);
      max_bin=fPYfit_nu2->GetMaximumBin();
      Py2=fPYfit_nu2->GetBinCenter(max_bin);
      max_bin=fPZfit_nu2->GetMaximumBin();
      Pz2=fPZfit_nu2->GetBinCenter(max_bin);
      //---- setting 4-vecs
      nu1_tmp.SetXYZM(Px1,Py1,Pz1,0.0);
      nu2_tmp.SetXYZM(Px2,Py2,Pz2,0.0);
      fDitauStuffHisto.nutau1=nu1_tmp; 
      fDitauStuffHisto.nutau2=nu2_tmp; 
    }

//   delete fMfit_all;
//   delete fPXfit_nu1;  
//   delete fPYfit_nu1;  
//   delete fPZfit_nu1;  
//   delete fPXfit_nu2;  
//   delete fPYfit_nu2;  
//   delete fPZfit_nu2;  

  if(fUseVerbose==1)
    {
      if(fit_code==0)
        {
          std::cout<<"!!!----> Warning-3 in MissingMassCalculator::DitauMassCalculatorV9() : fit status="<<fit_code<<std::endl;
          std::cout<<"....... No solution is found. Printing input info ......."<<std::endl;
          std::cout<<"  "<<std::endl;

          std::cout<<"  vis Tau-1: Pt="<<tau_vec1.Pt()<<"  M="<<tau_vec1.M()<<" eta="<<tau_vec1.Eta()<<"  phi="<<tau_vec1.Phi()<<"  type="<<tau_type1<<std::endl;
          std::cout<<"  vis Tau-2: Pt="<<tau_vec2.Pt()<<"  M="<<tau_vec2.M()<<" eta="<<tau_vec2.Eta()<<"  phi="<<tau_vec2.Phi()<<"  type="<<tau_type2<<std::endl;
          std::cout<<"  MET="<<met_vec.Mod()<<"  Met_X="<<met_vec.Px()<<"  Met_Y="<<met_vec.Py()<<std::endl;
          std::cout<<" ---------------------------------------------------------- "<<std::endl;
        }
    }	
  return fit_code;
}

												  


// jul/11 completely reshuffled but keep same logic
// results should be identical to V9 but they differ by a fraction of GeV
int MissingMassCalculator::DitauMassCalculatorV9fast(const TLorentzVector & tau_vec1, const int & tau_type1, 
						 const TLorentzVector & tau_vec2, const int & tau_type2, 
						 const TVector2 & met_vec) {

  int fit_code=0; // 0==bad, 1==good
  ClearDitauStuff(fDitauStuffFit);
  ClearDitauStuff(fDitauStuffHisto);
     
  totalProbSum=0;
  totalTautauSum=0;

  fMfit_all->Reset();
  fPXfit_nu1->Reset();
  fPYfit_nu1->Reset();
  fPZfit_nu1->Reset();
  fPXfit_nu2->Reset();
  fPYfit_nu2->Reset();
  fPZfit_nu2->Reset();


  //------- Settings -------------------------------
  int Niter=Niter_fit1; // number of points for each dR loop
  int NiterMET=Niter_fit2; // number of iterations for each MET scan loop  
  int NiterMnu=Niter_fit3; // number of iterations for Mnu loop
  double Mtau=1.777;

  double METresX=InputInfo.METsigmaL; // MET resolution in direction parallel to leading jet, for MET scan
  double METresY=InputInfo.METsigmaP; // MET resolution in direction perpendicular to leading jet, for MET scan
  double N_METsigma=Nsigma_METscan; // number of sigmas for MET scan
  double METresX_binSize=2*N_METsigma*METresX/NiterMET;
  double METresY_binSize=2*N_METsigma*METresY/NiterMET;

  //-------- end of Settings
  int Leptau1=0; // lepton-tau code
  int Leptau2=0; 
  if(tau_type1==0) Leptau1=1;
  if(tau_type2==0) Leptau2=1;      
  
  double phi1=0.0;
  double phi2=0.0;
  int solution=0;
  int isol=0;
  
  
  TLorentzVector nu1_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector nu2_tmp(0.0,0.0,0.0,0.0);
  //  TLorentzVector tau1_tmp(0.0,0.0,0.0,0.0);
  // TLorentzVector tau2_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector tautau_tmp(0.0,0.0,0.0,0.0);
  TVector2 metvec_tmp(0.0,0.0);
  
  double metprob=1.0;
  prob_tmp=0.0;
  
  double met_smear_x=0.0;
  double met_smear_y=0.0;
  double met_smearL=0.0;
  double met_smearP=0.0;
  

  //SpeedUp
    TStopwatch timer;
  if (fSpeedStudy) timer.Start();
  int iter0=0;
  iter1=0;
  iter2=0;
  iter3=0;
  iter4=0;
  const double cosphi_jet=cos(InputInfo.phi_jet);
  const double sinphi_jet=sin(InputInfo.phi_jet);
  const double tau_vec1phi=tau_vec1.Phi();
  const double tau_vec2phi=tau_vec2.Phi();
  const double tau_vec1m=tau_vec1.M();
  const double tau_vec2m=tau_vec2.M();

  int nsol1;
  int nsol2;

  iang1low=0;
  iang1high=0;
  iang2low=0;
  iang2high=0;

  std::vector<TLorentzVector> tauvecsol1(4);
  std::vector<TLorentzVector> tauvecsol2(4);
  std::vector<double> tauvecprob1(4);
  std::vector<double> tauvecprob2(4);
  

  //const double threepi=3*TMath::Pi();
  //const double twopi=2*TMath::Pi();
  //const double onepi=1*TMath::Pi();

  
  //---------------------------------------------    

      //------- Settings -------------------------------	  
      double dPhi_max=dRmax_tau;
      double dPhi_binSize=dPhi_max/Niter;	
      double dMnu_max1=Mtau-tau_vec1m; 
      double dMnu_max2=Mtau-tau_vec2m;
      double Mnu_binSize1=dMnu_max1/NiterMnu;
      double Mnu_binSize2=dMnu_max2/NiterMnu;
      // for non leptonic tau decau number of iteration zero
      int NiterMnu1 = (Leptau1==1) ? NiterMnu : 1;
      int NiterMnu2 = (Leptau2==1) ? NiterMnu : 1;
      //-------- end of Settings	  	  
      std::vector<TLorentzVector> nuvec1_tmp;
      std::vector<TLorentzVector> nuvec2_tmp;
      double M1=Mtau;
      double M2=Mtau;
      double M_nu1=0.0;
      double M_nu2=0.0;

      // same code for ll,lh,hh. Only difference is number of iteration on Mnu
      //---------------------------------------------
      const double phi1min=tau_vec1phi-dPhi_max+dPhi_binSize;
      const double phi1max=tau_vec1phi+dPhi_max-dPhi_binSize/2; // margin to avoid rounding errors
      const double phi2min=tau_vec2phi-dPhi_max+dPhi_binSize;
      const double phi2max=tau_vec2phi+dPhi_max-dPhi_binSize/2; // margin to avoid rounding errors






      // iteration on MET outside
		      for(int i3=0; i3<NiterMET+1; i3++) // MET_X scan
			{
			  met_smearL=METresX_binSize*i3-N_METsigma*METresX;
			  for(int i4=0; i4<NiterMET+1; i4++) // MET_Y scan
			    {
			      met_smearP=METresY_binSize*i4-N_METsigma*METresY;
			      //SpeedUp 
			      met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
			      met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
			      metprob=MetProbability(met_smearL,met_smearP,METresX,METresY);
			      if (metprob<=0) continue;
			      //FIXME could compute metprob only if a solution
			      metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);






      for(phi1=phi1min ; phi1<phi1max; phi1+= dPhi_binSize) //---- loop-1: dPhi for nu1
	{



	  //---------------------------------------------
	  for(phi2=phi2min ; phi2<phi2max; phi2+= dPhi_binSize) //---- loop-1: dPhi for nu2
	    {









			      for(int i5=0; i5<NiterMnu1; i5++) //---- Mnu1 scan 
				{
				  M_nu1=Mnu_binSize1*i5;
				  if(M_nu1>=(Mtau-tau_vec1m)) continue; // checking condition: M_nu < M_tau-M_lep
				  M1=sqrt(Mtau*Mtau-M_nu1*M_nu1);
				  for(int i6=0; i6<NiterMnu2; i6++) //---- Mnu2 scan 
				    {
				      M_nu2=Mnu_binSize2*i6;
				      if(M_nu2>=(Mtau-tau_vec2m)) continue; // checking condition: M_nu < M_tau-M_lep
				      M2=sqrt(Mtau*Mtau-M_nu2*M_nu2);
				      //nuvec1_tmp.clear();
				      //nuvec2_tmp.clear();
				      //std::cout << "phi " << phi1 << " " << phi2 << " tauvec 1 " << tau_vec1.Px() << " " << tau_vec1.Py() << " " << tau_vec1.Pz() << " " << tau_vec1.E() << " " << tau_vec1.M() << " tauvec 2 " << tau_vec2.Px() << " " << tau_vec2.Py() << " " << tau_vec2.Pz() << " " << tau_vec2.E() << " " << tau_vec2.M() << std::endl;
				      //std::cout << " met " << metvec_tmp.X() << " " << metvec_tmp.Y() << " " << metprob << "M " << M1 << " " << M2 << std::endl;
				      solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nsol1,nsol2);
				      ++iter0;
				      //if (iter1<20) std::cout << iter1 << "total " << totalTautauSum << std::endl;



				      //std::cout << " nuvecsol 1 " << nuvecsol1[0].Px() << " " << nuvecsol1[0].Py() << " " << nuvecsol1[0].Pz() << " " << nuvecsol1[0].E() << " " << nuvecsol1[0].M() << " nuvec 2 " << nuvecsol2[0].Px() << " " << nuvecsol2[0].Py() << " " << nuvecsol2[0].Pz() << " " << nuvecsol2[0].E() << " " << nuvecsol2[0].M() << std::endl;


				      if (fSpeedStudy)
					{
					  ++isol;
					  if (isol % 1000==1)
					    {
					      const double aux=nuvecsol2[0].E();
					      int solutionaux=NuPsolutionV2(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nuvec1_tmp,nuvec2_tmp);				      
					      if (solutionaux!=solution){
						std::cout << "SpeedUp WARNING mismatch between quick NuPsolutionV2 computation and slow "  << isol << " ref " << solutionaux << " check " << solution << std::endl;
					      }
					      else if (solution==1) {
					    
						const double ratio=std::abs(aux-nuvec2_tmp[0].E())/(aux+nuvec2_tmp[0].E());
						
						if  (ratio>=1E-3){
						  std::cout << "SpeedUp WARNING mismatch between quick NuPsolutionV2 computation and slow "  << isol << " ratio " << ratio << std::endl;
						} else {
						  // std::cout << "comp with quick NuPSolutionV2 OK" << std::endl;
						}
						
					      }
					  
					    }
					}	

				      if(solution!=1) continue;
				      ++iter1;				      


				      const bool oneSol=refineSolutions (tau_vec1, tau_type1, tau_vec2, tau_type2, 
						       metprob,tauvecsol1,tauvecsol2,
						       nsol1, nsol2 );

				      if (oneSol) fit_code=1;
				      
				      
				    }
				  
				}
			    }
			}
	    }
	}	  



  

  //SpeedUp
  if (fSpeedStudy) timer.Stop()  ;
  
  if (fSpeedStudy){
    std::cout << "SpeedUp niters=" << iter0 << " " << iter1 
	      << " " << iter2 <<" " << iter3 <<" " << iter4
              << "skip:" << iang1low << " " << iang1high << " " << iang2low << " " << iang2high 
	      << " CPU time = " << timer.CpuTime() << "s real = " << timer.RealTime() << " s " << std::endl;
    std::cout << "SpeedUp totalProbSum=" << totalProbSum << " totalTautauSum=" << totalTautauSum << std::endl;
  }
  

  int max_bin;
  double Px1, Py1, Pz1;
  double Px2, Py2, Pz2;
  if(fMfit_all->GetEntries()>0)
    {
      max_bin=fMfit_all->GetMaximumBin();
      fDitauStuffHisto.Mditau_best=fMfit_all->GetBinCenter(max_bin);
      double prob_hist=1.0*(fMfit_all->GetBinContent(max_bin))/(1.0*fMfit_all->GetEntries());
      if(prob_hist>=1.0) prob_hist=1.0;
      if(prob_hist>0.0) fDitauStuffHisto.Sign_best=-log10(prob_hist);
      else fDitauStuffHisto.Sign_best=-1.0;
      if(fDitauStuffHisto.Mditau_best>0.0) fDitauStuffHisto.RMSoverMPV=fMfit_all->GetRMS()/fDitauStuffHisto.Mditau_best;
      //---- getting Nu1
      max_bin=fPXfit_nu1->GetMaximumBin();
      Px1=fPXfit_nu1->GetBinCenter(max_bin);
      max_bin=fPYfit_nu1->GetMaximumBin();
      Py1=fPYfit_nu1->GetBinCenter(max_bin);
      max_bin=fPZfit_nu1->GetMaximumBin();
      Pz1=fPZfit_nu1->GetBinCenter(max_bin);
      //---- getting Nu2
      max_bin=fPXfit_nu2->GetMaximumBin();
      Px2=fPXfit_nu2->GetBinCenter(max_bin);
      max_bin=fPYfit_nu2->GetMaximumBin();
      Py2=fPYfit_nu2->GetBinCenter(max_bin);
      max_bin=fPZfit_nu2->GetMaximumBin();
      Pz2=fPZfit_nu2->GetBinCenter(max_bin);
      //---- setting 4-vecs
      nu1_tmp.SetXYZM(Px1,Py1,Pz1,0.0);
      nu2_tmp.SetXYZM(Px2,Py2,Pz2,0.0);
      fDitauStuffHisto.nutau1=nu1_tmp; 
      fDitauStuffHisto.nutau2=nu2_tmp; 
    }


  if(fUseVerbose==1)
    {
      if(fit_code==0)
        {
          std::cout<<"!!!----> Warning-3 in MissingMassCalculator::DitauMassCalculatorV9fast() : fit status="<<fit_code<<std::endl;
          std::cout<<"....... No solution is found. Printing input info ......."<<std::endl;
          std::cout<<"  "<<std::endl;

          std::cout<<"  vis Tau-1: Pt="<<tau_vec1.Pt()<<"  M="<<tau_vec1.M()<<" eta="<<tau_vec1.Eta()<<"  phi="<<tau_vec1.Phi()<<"  type="<<tau_type1<<std::endl;
          std::cout<<"  vis Tau-2: Pt="<<tau_vec2.Pt()<<"  M="<<tau_vec2.M()<<" eta="<<tau_vec2.Eta()<<"  phi="<<tau_vec2.Phi()<<"  type="<<tau_type2<<std::endl;
          std::cout<<"  MET="<<met_vec.Mod()<<"  Met_X="<<met_vec.Px()<<"  Met_Y="<<met_vec.Py()<<std::endl;
          std::cout<<" ---------------------------------------------------------- "<<std::endl;
        }
    }

  return fit_code;
}


// 05/25/11: modified V9 (became V10). Added adaptive step size for dTheta3d, scan in dTheta3d is transformed into dPhi scan
// 06/01/11: includes David's improvements
int MissingMassCalculator::DitauMassCalculatorV10fast(const TLorentzVector & tau_vec1, const int & tau_type1, 
						 const TLorentzVector & tau_vec2, const int & tau_type2, 
						 const TVector2 & met_vec) {

  int fit_code=0; // 0==bad, 1==good
  ClearDitauStuff(fDitauStuffFit);
  ClearDitauStuff(fDitauStuffHisto);

  totalProbSum=0;

  fMfit_all->Reset();
  fPXfit_nu1->Reset();;
  fPYfit_nu1->Reset();;
  fPZfit_nu1->Reset();;
  fPXfit_nu2->Reset();;
  fPYfit_nu2->Reset();;
  fPZfit_nu2->Reset();;

  int nsol1;
  int nsol2;

     
  //------- Settings -------------------------------
  int Niter=Niter_fit1; // number of points for each dR loop
  int NiterMET=Niter_fit2; // number of iterations for each MET scan loop  
  int NiterMnu=Niter_fit3; // number of iterations for Mnu loop
  double Mtau=1.777;
  double Mnu_binSize=Mtau/NiterMnu;

  double METresX=InputInfo.METsigmaL; // MET resolution in direction parallel to leading jet, for MET scan
  double METresY=InputInfo.METsigmaP; // MET resolution in direction perpendicular to leading jet, for MET scan
  double N_METsigma=Nsigma_METscan; // number of sigmas for MET scan
  double METresX_binSize=2*N_METsigma*METresX/NiterMET;
  double METresY_binSize=2*N_METsigma*METresY/NiterMET;
  if(METresX_binSize<0.1) METresX_binSize=0.1;
  if(METresY_binSize<0.1) METresY_binSize=0.1;
  double METscanX_limit=2*N_METsigma*METresX;
  double METscanY_limit=2*N_METsigma*METresY;

  double dTheta3d_limit1=dTheta3DLimitVis(tau_type1,0,tau_vec1.P()); // using 99% upper limit; setting is hardcoded for now... need to change to use Set..() 
  double dTheta3d_limit2=dTheta3DLimitVis(tau_type2,0,tau_vec2.P()); // using 99% upper limit; setting is hardcoded for now... need to change to use Set..() 
  double scan_dTheta3d_binsize1=dTheta3d_limit1/(1.0*Niter);
  double scan_dTheta3d_binsize2=dTheta3d_limit2/(1.0*Niter);
  if(scan_dTheta3d_binsize1>dTheta3d_binMax) scan_dTheta3d_binsize1=dTheta3d_binMax;
  if(scan_dTheta3d_binsize2>dTheta3d_binMax) scan_dTheta3d_binsize2=dTheta3d_binMax;
  if(scan_dTheta3d_binsize1<dTheta3d_binMin) scan_dTheta3d_binsize1=dTheta3d_binMin;
  if(scan_dTheta3d_binsize2<dTheta3d_binMin) scan_dTheta3d_binsize2=dTheta3d_binMin;

  if(fUseVerbose==1)
    {
      std::cout<<"------------------ Printing scan parameters ---------------------"<<std::endl;
      std::cout<<"  Niter(MET scan)="<<NiterMET<<"  Niter(Mnu scan)="<<NiterMnu<<" Niter(dPhi scan)="<<Niter<<std::endl;
      std::cout<<"  MET_X scan: range="<<2*N_METsigma*METresX<<"  step size="<<METresX_binSize<<std::endl;
      std::cout<<"  MET_Y scan: range="<<2*N_METsigma*METresY<<"  step size="<<METresY_binSize<<std::endl;
      std::cout<<"  Mnu scan: step size="<<Mnu_binSize<<std::endl;
      std::cout<<"  dPhi1 scan: dPhi limit="<<dTheta2dPhi(tau_vec1.Eta(),dTheta3d_limit1)<<std::endl;
      std::cout<<"  dPhi2 scan: dPhi limit="<<dTheta2dPhi(tau_vec2.Eta(),dTheta3d_limit2)<<std::endl;
      std::cout<<"  dTheta3d_limit1="<<dTheta3d_limit1<<" scan_dTheta3d_binsize1="<<scan_dTheta3d_binsize1<<std::endl;
      std::cout<<"  dTheta3d_limit2="<<dTheta3d_limit2<<" scan_dTheta3d_binsize2="<<scan_dTheta3d_binsize2<<std::endl;
      std::cout<<"................................................................."<<std::endl;
    }

  //-------- end of Settings
  int Leptau1=0; // lepton-tau code
  int Leptau2=0; 
  if(tau_type1==0) Leptau1=1;
  if(tau_type2==0) Leptau2=1;      
  
  double dPhi_nu1[2];
  double dPhi_nu2[2];  
  double dPhi_tmp;
  double dTheta3d_tmp1=0.0;
  double dTheta3d_tmp2=0.0;
  int phi1_loop=0;
  int phi2_loop=0;
  double phi1=0.0;
  double phi2=0.0;
  int solution=0;
  
  TLorentzVector nu1_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector nu2_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector tau1_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector tau2_tmp(0.0,0.0,0.0,0.0);
  TLorentzVector tautau_tmp(0.0,0.0,0.0,0.0);
  TVector2 metvec_tmp(0.0,0.0);

  std::vector<TLorentzVector> tauvecsol1(2);
  std::vector<TLorentzVector> tauvecsol2(2);
  
  
  double tauProb1r=1.0;
  double tauProb2r=1.0;			  
  double metprob=1.0;
  double ditauProb=1.0;
  //  double resProb=1.0;
  double sign_tmp=0.0; 
  double totalProb=0.0;
  prob_tmp=0.0;
  
  double met_smear_x=0.0;
  double met_smear_y=0.0;
  double met_smearL=0.0;
  double met_smearP=0.0;
  
  double angle1=0.0;
  double angle2=0.0;
  //--- define histograms for histogram method
  //--- upper limits need to be revisied in the future!!! It may be not enough for some analyses
//   TH1F* fMfit_all=new TH1F("h1","M",2500,0.0,beamEnergy/3.5); // all solutions 
//   TH1F* fPXfit_nu1=new TH1F("h2","Px1",10000,-beamEnergy/3.5,beamEnergy/3.5); // Px for nu1 
//   TH1F* fPYfit_nu1=new TH1F("h3","Py1",10000,-beamEnergy/3.5,beamEnergy/3.5); // Py for nu1 
//   TH1F* fPZfit_nu1=new TH1F("h4","Pz1",10000,-beamEnergy/3.5,beamEnergy/3.5); // Pz for nu1 
//   TH1F* fPXfit_nu2=new TH1F("h5","Px2",10000,-beamEnergy/3.5,beamEnergy/3.5); // Px for nu2 
//   TH1F* fPYfit_nu2=new TH1F("h6","Py2",10000,-beamEnergy/3.5,beamEnergy/3.5); // Py for nu2 
//   TH1F* fPZfit_nu2=new TH1F("h7","Pz2",10000,-beamEnergy/3.5,beamEnergy/3.5); // Pz for nu2 

  //SpeedUp
    TStopwatch timer;
  if (fSpeedStudy) timer.Start();
  const double cosphi_jet=cos(InputInfo.phi_jet);
  const double sinphi_jet=sin(InputInfo.phi_jet);
  const double tau_vec1phi=tau_vec1.Phi();
  const double tau_vec2phi=tau_vec2.Phi();
  const double tau_vec1m=tau_vec1.M();
  const double tau_vec2m=tau_vec2.M();
  //const double threepi=3*TMath::Pi();
  //const double twopi=2*TMath::Pi();
  //const double onepi=1*TMath::Pi();

  double Mvis=(tau_vec1+tau_vec2).M();
  TLorentzVector met4vec(0.0,0.0,0.0,0.0);
  met4vec.SetPxPyPzE(met_vec.X(),met_vec.Y(),0.0,met_vec.Mod());
  double Meff=(tau_vec1+tau_vec2+met4vec).M();
  
  //---------------------------------------------    
  if((Leptau1+Leptau2)==2) // both tau's are leptonic V10
    {
      //------- Settings -------------------------------	  
      double dMnu_max1=Mtau-tau_vec1m; 
      double dMnu_max2=Mtau-tau_vec2m;
      double Mnu_binSize1=dMnu_max1/NiterMnu;
      double Mnu_binSize2=dMnu_max2/NiterMnu;
      //-------- end of Settings	  	  
      std::vector<TLorentzVector> nuvec1_tmp;
      std::vector<TLorentzVector> nuvec2_tmp;
      double M1=Mtau;
      double M2=Mtau;
      double M_nu1=0.0;
      double M_nu2=0.0;
      //---------------------------------------------
      for(int i1=0; i1<Niter && dTheta3d_tmp1<dTheta3d_limit1; i1++) //---- loop-1: dPhi for nu1
	{
	  dTheta3d_tmp1=scan_dTheta3d_binsize1*i1;
	  dPhi_tmp=dTheta2dPhi(tau_vec1.Eta(),dTheta3d_tmp1);
	  if(fabs(dPhi_tmp)>3.1) continue;
	    //	  dPhi_tmp=dPhi_binSize*i1;
	  dPhi_nu1[0]=-dPhi_tmp;
	  dPhi_nu1[1]=dPhi_tmp;
	  phi1_loop=i1>0 ? 2 : 1;
	  //---------------------------------------------
	  for(int i2=0; i2<Niter && dTheta3d_tmp2<dTheta3d_limit2; i2++) //---- loop-2: dPhi for nu2
	    {
	      dTheta3d_tmp2=scan_dTheta3d_binsize2*i2;
	      dPhi_tmp=dTheta2dPhi(tau_vec2.Eta(),dTheta3d_tmp2);
	      if(fabs(dPhi_tmp)>3.1) continue;
	      // 	      dPhi_tmp=dPhi_binSize*i2;
	      dPhi_nu2[0]=-dPhi_tmp;
	      dPhi_nu2[1]=dPhi_tmp;
	      phi2_loop=i2>0 ? 2 : 1;		      
	      for(int ij1=0; ij1<phi1_loop; ij1++)
		{
		  phi1=tau_vec1phi+dPhi_nu1[ij1];
		  for(int ij2=0; ij2<phi2_loop; ij2++)
		    {		      
		      phi2=tau_vec2phi+dPhi_nu2[ij2];
		      for(int i3=0; i3<NiterMET+1 && met_smearL<METscanX_limit; i3++) // MET_X scan, has to be replaced by while() loops
			{
			  met_smearL=METresX_binSize*i3-N_METsigma*METresX;
			  for(int i4=0; i4<NiterMET+1 && met_smearP<METscanY_limit; i4++) // MET_Y scan
			    {
			      met_smearP=METresY_binSize*i4-N_METsigma*METresY;
			      //SpeedUp 
			      met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
			      met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
			      metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);
			      metprob=MetProbability(met_smearL,met_smearP,METresX,METresY);

			      for(int i5=0; i5<NiterMnu; i5++) //---- Mnu1 scan 
				{
				  M_nu1=Mnu_binSize1*i5;
				  if(M_nu1>=(Mtau-tau_vec1m)) continue; // checking condition: M_nu < M_tau-M_lep
				  M1=sqrt(Mtau*Mtau-M_nu1*M_nu1);
				  for(int i6=0; i6<NiterMnu; i6++) //---- Mnu2 scan 
				    {
				      M_nu2=Mnu_binSize2*i6;
				      if(M_nu2>=(Mtau-tau_vec2m)) continue; // checking condition: M_nu < M_tau-M_lep
				      M2=sqrt(Mtau*Mtau-M_nu2*M_nu2);
				      nuvec1_tmp.clear();
				      nuvec2_tmp.clear();
				      solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nsol1,nsol2);
				      if(solution==1) // there is solution
					{
					  for(int j1=0; j1<nsol1; j1++)
					    {
					      //SpeedUp reference to avoid picking again in again in array
					      TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
					      TLorentzVector & tauvecsol1j=tauvecsol1[j1];
					      
					      //SpeedUp sum energy rather than TLV
					      if(tau_vec1.E()+nuvec1_tmpj.E()>=beamEnergy) continue;
					      //SpeedUp use SetXYZM to reassign the mass
					      nuvec1_tmpj.SetXYZM(nuvec1_tmpj.Px(),
								  nuvec1_tmpj.Py(),
								  nuvec1_tmpj.Pz(),
								  M_nu1);


					      tauvecsol1j.SetPxPyPzE(0.,0.,0.,0.);
					      tauvecsol1j+=nuvec1_tmpj;
					      tauvecsol1j+=tau_vec1;
					      angle1=Angle(nuvec1_tmpj,tau_vec1);
					      const double tau1_tmpp=tauvecsol1j.P();
					      if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)) continue; // lower 99% bound
					      if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)) continue; // upper 99% bound
					      
					      tauProb1r=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp);

					      for(int j2=0; j2<nsol2; j2++)
						{
						  TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
						  TLorentzVector & tauvecsol2j=tauvecsol2[j2];
						  //SpeedUp
						  nuvec2_tmpj.SetXYZM(nuvec2_tmpj.Px(),
								      nuvec2_tmpj.Py(),
								      nuvec2_tmpj.Pz(),
								      M_nu2);

						  //SpeedUp sum energy rather than TLV
						  if(tau_vec2.E()+nuvec2_tmpj.E()>=beamEnergy) continue;
						  //if(tau_vec1.E()+nuvec1_tmpj.E()+tau_vec2.E()+nuvec2_tmpj.E()>=2.0*beamEnergy) continue;
// 					  if(CheckSolutions(nuvec1_tmp[j1],tau_vec1,Leptau1)<1) continue;
// 					  if(CheckSolutions(nuvec2_tmp[j2],tau_vec2,Leptau2)<1) continue;
						  //use += TLV rather than TLV1+TLV2
						  tauvecsol2j.SetPxPyPzE(0.,0.,0.,0.);
						  tauvecsol2j+=nuvec2_tmpj;
						  tauvecsol2j+=tau_vec2;
						  angle2=Angle(nuvec2_tmpj,tau_vec2);
						  const double tau2_tmpp=tauvecsol2j.P();
						  if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)) continue; // lower 99% bound
						  if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)) continue; // upper 99% bound
						  tauProb2r=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp);

						  //						  resProb=ResonanceProbability(tau1_tmp,tau2_tmp);
// 						  totalProb=tauProb1r*tauProb2r*metprob*resProb;
						  totalProb=tauProb1r*tauProb2r*metprob;
						  totalProbSum+=totalProb;				      
						  if(totalProb<=0.0) continue;
						  fit_code=1; // at least one solution is found

						  tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
						  tautau_tmp+=tauvecsol1j;
						  tautau_tmp+=tauvecsol2j;
						  
						  const double mtautau=tautau_tmp.M();
						  
						  fMfit_all->Fill(mtautau,totalProb);
						  fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
						  fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
						  fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
						  fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
						  fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
						  fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
						  if(totalProb>prob_tmp) // fill solution with highest probability                                   
						    {
						      prob_tmp=totalProb;
						      fDitauStuffFit.Mditau_best=mtautau; 
						      fDitauStuffFit.Sign_best=sign_tmp; 
						      fDitauStuffFit.nutau1=nuvec1_tmpj; 
						      fDitauStuffFit.nutau2=nuvec2_tmpj; 
						    }
						}
					    }
					}
				      else continue; // no solution
				    }
				}
			    }
			}
		    }
		}
	    }
	}	  
    }

  if((Leptau1+Leptau2)==1) // one leptonic tau V10
    {
      //------- Settings -------------------------------	  
      double dMnu_max= (Leptau1==1) ? Mtau-tau_vec1m : Mtau-tau_vec2m;
      Mnu_binSize=dMnu_max/NiterMnu;
      //-------- end of Settings	  	  
      std::vector<TLorentzVector> nuvec1_tmp;
      std::vector<TLorentzVector> nuvec2_tmp;
      double M1=Mtau;
      double M2=Mtau;
      double M_nu=0.0;
      
      //---------------------------------------------
      for(int i1=0; i1<Niter && dTheta3d_tmp1<dTheta3d_limit1; i1++) //---- loop-1: dPhi for nu1
	{
	  dTheta3d_tmp1=scan_dTheta3d_binsize1*i1;
	  dPhi_tmp=dTheta2dPhi(tau_vec1.Eta(),dTheta3d_tmp1);
	  if(fabs(dPhi_tmp)>3.1) continue;
	  dPhi_nu1[0]=-dPhi_tmp;
	  dPhi_nu1[1]=dPhi_tmp;
	  phi1_loop=i1>0 ? 2 : 1;
	  //---------------------------------------------
	  for(int i2=0; i2<Niter && dTheta3d_tmp2<dTheta3d_limit2; i2++) //---- loop-2: dPhi for nu2
	    {
	      dTheta3d_tmp2=scan_dTheta3d_binsize2*i2;
	      dPhi_tmp=dTheta2dPhi(tau_vec2.Eta(),dTheta3d_tmp2);
	      if(fabs(dPhi_tmp)>3.1) continue;
	      dPhi_nu2[0]=-dPhi_tmp;
	      dPhi_nu2[1]=dPhi_tmp;
	      phi2_loop=i2>0 ? 2 : 1;
	      for(int i3=0; i3<NiterMnu; i3++) //---- loop-3: virtual neutrino mass 
		{
		  M_nu=Mnu_binSize*i3;
		  if(Leptau1==1 && M_nu>=(Mtau-tau_vec1m)) continue; // checking condition: M_nu < M_tau-M_lep
		  if(Leptau2==1 && M_nu>=(Mtau-tau_vec2m)) continue; // checking condition: M_nu < M_tau-M_lep
		  if(Leptau1==1) M1=sqrt(Mtau*Mtau-M_nu*M_nu);
		  else M2=sqrt(Mtau*Mtau-M_nu*M_nu);
		  
		  for(int ij1=0; ij1<phi1_loop; ij1++)
		    {
		      phi1=tau_vec1phi+dPhi_nu1[ij1];
		      for(int ij2=0; ij2<phi2_loop; ij2++)
			{		      
			  phi2=tau_vec2phi+dPhi_nu2[ij2];
			  for(int i4=0; i4<NiterMET+1 && met_smearL<METscanX_limit; i4++) // MET_X scan
			    {
			      met_smearL=METresX_binSize*i4-N_METsigma*METresX;

			      for(int i5=0; i5<NiterMET+1 && met_smearP<METscanY_limit; i5++) // MET_Y scan
				{
				  met_smearP=METresY_binSize*i5-N_METsigma*METresY;
				  met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
				  met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
				  metprob=MetProbability(met_smearL,met_smearP,METresX,METresY);

				  metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);

				  nuvec1_tmp.clear();
				  nuvec2_tmp.clear();
				  solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,M1,M2,phi1,phi2,nsol1,nsol2);
				  if(solution==1) // there is solution
				    {
				      for(int j1=0; j1<nsol1; j1++)
					{
					  //SpeedUp reference to avoid picking again in again in array
					  TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
					  if(Leptau1==1)
					    {
					      //SpeedUp
					      nuvec1_tmpj.SetXYZM(nuvec1_tmpj.Px(),
								  nuvec1_tmpj.Py(),
								  nuvec1_tmpj.Pz(),
								  M_nu);
					    }

					  if(tau_vec1.E()+nuvec1_tmpj.E()>=beamEnergy) continue;
					  tau1_tmp.SetPxPyPzE(0.,0.,0.,0.);
					  tau1_tmp+=nuvec1_tmpj;
					  tau1_tmp+=tau_vec1;
					  const double tau1_tmpp=tau1_tmp.P();

					  angle1=Angle(nuvec1_tmpj,tau_vec1);
					  if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)) continue; // lower 99% bound
					  if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)) continue; // upper 99% bound
					  tauProb1r=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp);

					  for(int j2=0; j2<nsol2; j2++)
					    {
					      TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];	
					      if(Leptau2==1)
						{
						  //SpeedUp
						  nuvec2_tmpj.SetXYZM(nuvec2_tmpj.Px(),
								      nuvec2_tmpj.Py(),
								      nuvec2_tmpj.Pz(),
								      M_nu);

						}						      
					      if(tau_vec2.E()+nuvec2_tmpj.E()>=beamEnergy) continue;
					      //SpeedUp uselessif(tau_vec1.E()+nuvec1_tmp[j1].E()+tau_vec2.E()+nuvec2_tmp[j2].E()>=2.0*beamEnergy) continue;
// 					      if(CheckSolutions(nuvec1_tmp[j1],tau_vec1,Leptau1)<1) continue;
// 					      if(CheckSolutions(nuvec2_tmp[j2],tau_vec2,Leptau2)<1) continue;
					      //SpeedUp use += TLV rather than TLV1+TLV2
					      //tau1_tmp=nuvec1_tmp[j1]+tau_vec1;
					      tau2_tmp.SetPxPyPzE(0.,0.,0.,0.);
					      tau2_tmp+=nuvec2_tmpj;
					      tau2_tmp+=tau_vec2;
					      angle2=Angle(nuvec2_tmpj,tau_vec2);
					      const double tau2_tmpp=tau2_tmp.P();
					      if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)) continue; // lower 99% bound
					      if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)) continue; // upper 99% bound

					      tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
					      tautau_tmp+=tau1_tmp;
					      tautau_tmp+=tau2_tmp;
					      const double mtautau=tautau_tmp.M();
					  
					      if(TailCleanUp(tau_type1,tau_vec1,nuvec1_tmpj, 
							     tau_type2,tau_vec2,nuvec2_tmpj, 
							     mtautau,Mvis,Meff,InputInfo.DelPhiTT)==0) continue;
			      
					      tauProb2r=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp);
					      ditauProb=TauProbability(tau_type1,tau_vec1,nuvec1_tmpj, 
								       tau_type2,tau_vec2,nuvec2_tmpj);
					      totalProb=tauProb1r*tauProb2r*metprob*ditauProb;
// 					      totalProb=tauProb1r*tauProb2r*metprob;
					      totalProbSum+=totalProb;	

					      if(totalProb<=0.0) continue;
					      fit_code=1; // at least one solution is found


// 					      tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
// 					      tautau_tmp+=tau1_tmp;
// 					      tautau_tmp+=tau2_tmp;
// 					      const double mtautau=tautau_tmp.M();

					      fMfit_all->Fill(mtautau,totalProb);
					      fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
					      fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
					      fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
					      fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
					      fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
					      fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
					      if(totalProb>prob_tmp) // fill solution with highest probability                                   
						{
						  prob_tmp=totalProb;
						  fDitauStuffFit.Mditau_best=mtautau; 
						  fDitauStuffFit.Sign_best=sign_tmp; 
						  fDitauStuffFit.nutau1=nuvec1_tmpj; 
						  fDitauStuffFit.nutau2=nuvec2_tmpj; 
						}
					    }
					}
				    }
				  else continue; // no solution
				}
			    }
			}
		    }
		}
	    }
	}	  
    }
  
  
  if((Leptau1+Leptau2)==0) // both tau's are hadronic V10
    {
      std::vector<TLorentzVector> nuvec1_tmp;
      std::vector<TLorentzVector> nuvec2_tmp;
      
      //---------------------------------------------
      for(int i1=0; i1<Niter && dTheta3d_tmp1<dTheta3d_limit1; i1++) //---- loop-1: dPhi for nu1
	{
	  dTheta3d_tmp1=scan_dTheta3d_binsize1*i1;
	  dPhi_tmp=dTheta2dPhi(tau_vec1.Eta(),dTheta3d_tmp1);
	  if(fabs(dPhi_tmp)>3.1) continue;
	  dPhi_nu1[0]=-dPhi_tmp;
	  dPhi_nu1[1]=dPhi_tmp;
	  phi1_loop=i1>0 ? 2 : 1;
	  //---------------------------------------------
	  for(int i2=0; i2<Niter && dTheta3d_tmp2<dTheta3d_limit2; i2++) //---- loop-2: dPhi for nu2
	    {
	      dTheta3d_tmp2=scan_dTheta3d_binsize2*i2;
	      dPhi_tmp=dTheta2dPhi(tau_vec2.Eta(),dTheta3d_tmp2);
	      if(fabs(dPhi_tmp)>3.1) continue;
	      dPhi_nu2[0]=-dPhi_tmp;
	      dPhi_nu2[1]=dPhi_tmp;
	      phi2_loop=i2>0 ? 2 : 1;		      
	      for(int ij1=0; ij1<phi1_loop; ij1++)
		{
		  phi1=tau_vec1phi+dPhi_nu1[ij1];
		  for(int ij2=0; ij2<phi2_loop; ij2++)
		    {		      
		      phi2=tau_vec2phi+dPhi_nu2[ij2];
		      for(int i3=0; i3<NiterMET+1 && met_smearL<METscanX_limit; i3++) // MET_X scan
			{
			  met_smearL=METresX_binSize*i3-N_METsigma*METresX;
			  for(int i4=0; i4<NiterMET+1 && met_smearP<METscanY_limit; i4++) // MET_Y scan
			    {
			      met_smearP=METresY_binSize*i4-N_METsigma*METresY;
			      metprob=MetProbability(met_smearL,met_smearP,METresX,METresY);
			      met_smear_x=met_smearL*cosphi_jet-met_smearP*sinphi_jet;
			      met_smear_y=met_smearL*sinphi_jet+met_smearP*cosphi_jet;
			      metvec_tmp.Set(met_vec.X()+met_smear_x,met_vec.Y()+met_smear_y);

			      nuvec1_tmp.clear();
			      nuvec2_tmp.clear();
			      solution=NuPsolutionV2fast(metvec_tmp,tau_vec1,tau_vec2,Mtau,Mtau,phi1,phi2,nsol1,nsol2);
			      if(solution==1) // there is solution
				{
				  for(int j1=0; j1<nsol1; j1++)
				    {
				      TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
				      if(tau_vec1.E()+nuvec1_tmpj.E()>=beamEnergy) continue;
				      tau1_tmp=nuvec1_tmpj+tau_vec1;
				      const double tau1_tmpp=tau1_tmp.P();
				      angle1=Angle(nuvec1_tmpj,tau_vec1);
				      if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)) continue; // lower 99% bound
				      if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)) continue; // upper 99% bound
				      tauProb1r=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp); // returns ZERO's ???????


				      for(int j2=0; j2<nsol2; j2++)
					{
					  TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
					  if(tau_vec2.E()+nuvec2_tmpj.E()>=beamEnergy) continue;
				          //SpeedUp useless if(tau_vec1.E()+nuvec1_tmpj.E()+tau_vec2.E()+nuvec2_tmpj.E()>=2.0*beamEnergy) continue;
// 					  if(CheckSolutions(nuvec1_tmp[j1],tau_vec1,Leptau1)<1) continue;
// 					  if(CheckSolutions(nuvec2_tmp[j2],tau_vec2,Leptau2)<1) continue;
					  tau2_tmp=nuvec2_tmpj+tau_vec2;
					  const double tau2_tmpp=tau2_tmp.P();
					  angle2=Angle(nuvec2_tmpj,tau_vec2);

					  if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)) continue; // lower 99% bound
					  if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)) continue; // upper 99% bound

					  tauProb2r=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp);

// 					  resProb=ResonanceProbability(tau1_tmp,tau2_tmp);
// 					  totalProb=tauProb1r*tauProb2r*metprob*resProb;
					  totalProb=tauProb1r*tauProb2r*metprob;

					  totalProbSum+=totalProb;				      
					  if(totalProb<=0.0) continue;


					  tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
					  tautau_tmp+=tau1_tmp;
					  tautau_tmp+=tau2_tmp;
					  const double mtautau=tautau_tmp.M();




					  fit_code=1; // at least one solution is found
					  sign_tmp=-log10(totalProb);
					  fMfit_all->Fill(mtautau,totalProb);
					  fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
					  fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
					  fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
					  fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
					  fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
					  fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
					  if(totalProb>prob_tmp) // fill solution with highest probability                                   
					    {
					      prob_tmp=totalProb;
					      fDitauStuffFit.Mditau_best=mtautau; 
					      fDitauStuffFit.Sign_best=sign_tmp; 
					      fDitauStuffFit.nutau1=nuvec1_tmpj; 
					      fDitauStuffFit.nutau2=nuvec2_tmpj; 
					    }
					}
				    }
				}
			      else continue; // no solution
			    }
			}
		    }
		}
	    }
	}	  
    }

  //SpeedUp
  if (fSpeedStudy) timer.Stop()  ;
  
  if (fSpeedStudy){
    std::cout<< " CPU time = " << timer.CpuTime() << "s real = " << timer.RealTime() << " s " << std::endl;
  }
  

  int max_bin;
  double Px1, Py1, Pz1;
  double Px2, Py2, Pz2;
  if(fMfit_all->GetEntries()>0)
    {
      max_bin=fMfit_all->GetMaximumBin();
      fDitauStuffHisto.Mditau_best=fMfit_all->GetBinCenter(max_bin);
      double prob_hist=1.0*(fMfit_all->GetBinContent(max_bin))/(1.0*fMfit_all->GetEntries());
      if(prob_hist>=1.0) prob_hist=1.0;
      if(prob_hist>0.0) fDitauStuffHisto.Sign_best=-log10(prob_hist);
      else fDitauStuffHisto.Sign_best=-1.0;
      if(fDitauStuffHisto.Mditau_best>0.0) fDitauStuffHisto.RMSoverMPV=fMfit_all->GetRMS()/fDitauStuffHisto.Mditau_best;
      //---- getting Nu1
      max_bin=fPXfit_nu1->GetMaximumBin();
      Px1=fPXfit_nu1->GetBinCenter(max_bin);
      max_bin=fPYfit_nu1->GetMaximumBin();
      Py1=fPYfit_nu1->GetBinCenter(max_bin);
      max_bin=fPZfit_nu1->GetMaximumBin();
      Pz1=fPZfit_nu1->GetBinCenter(max_bin);
      //---- getting Nu2
      max_bin=fPXfit_nu2->GetMaximumBin();
      Px2=fPXfit_nu2->GetBinCenter(max_bin);
      max_bin=fPYfit_nu2->GetMaximumBin();
      Py2=fPYfit_nu2->GetBinCenter(max_bin);
      max_bin=fPZfit_nu2->GetMaximumBin();
      Pz2=fPZfit_nu2->GetBinCenter(max_bin);
      //---- setting 4-vecs
      nu1_tmp.SetXYZM(Px1,Py1,Pz1,0.0);
      nu2_tmp.SetXYZM(Px2,Py2,Pz2,0.0);
      fDitauStuffHisto.nutau1=nu1_tmp; 
      fDitauStuffHisto.nutau2=nu2_tmp; 
    }


  if(fUseVerbose==1)
    {
      if(fit_code==0)
        {
          std::cout<<"!!!----> Warning-3 in MissingMassCalculator::DitauMassCalculatorV10fast() : fit status="<<fit_code<<std::endl;
          std::cout<<"....... No solution is found. Printing input info ......."<<std::endl;
          std::cout<<"  "<<std::endl;

          std::cout<<"  vis Tau-1: Pt="<<tau_vec1.Pt()<<"  M="<<tau_vec1.M()<<" eta="<<tau_vec1.Eta()<<"  phi="<<tau_vec1.Phi()<<"  type="<<tau_type1<<std::endl;
          std::cout<<"  vis Tau-2: Pt="<<tau_vec2.Pt()<<"  M="<<tau_vec2.M()<<" eta="<<tau_vec2.Eta()<<"  phi="<<tau_vec2.Phi()<<"  type="<<tau_type2<<std::endl;
          std::cout<<"  MET="<<met_vec.Mod()<<"  Met_X="<<met_vec.Px()<<"  Met_Y="<<met_vec.Py()<<std::endl;
          std::cout<<" ---------------------------------------------------------- "<<std::endl;
        }
    }

  return fit_code;
}



// standard collinear approximation
// it returns code=0 if collinear approximation can't be applied
// and code=1 and Mrec if collinear approximation was applied
int MissingMassCalculator::StandardCollApprox(TLorentzVector tau_vec1, TLorentzVector tau_vec2, TVector2 met_vec, double &Mrec) {
  int code=0;
  Mrec=0.0;
  double P_nu1=0.0;
  double P_nu2=0.0;
  int coll_code=NuPsolution(met_vec,tau_vec1.Theta(),tau_vec1.Phi(),tau_vec2.Theta(),tau_vec2.Phi(),P_nu1,P_nu2);
  if(coll_code==1)
    {
      code=1;
      TLorentzVector nu1(P_nu1*sin(tau_vec1.Theta())*cos(tau_vec1.Phi()),P_nu1*sin(tau_vec1.Theta())*sin(tau_vec1.Phi()),P_nu1*cos(tau_vec1.Theta()),P_nu1);
      TLorentzVector nu2(P_nu2*sin(tau_vec2.Theta())*cos(tau_vec2.Phi()),P_nu2*sin(tau_vec2.Theta())*sin(tau_vec2.Phi()),P_nu2*cos(tau_vec2.Theta()),P_nu2);
      Mrec=(nu1+nu2+tau_vec1+tau_vec2).M();
    }
  return code;
}

// Sasha: keep this for now, may need in the future
// returns analytical P(theta)-probability for given tau-topology 
// decayType==1 for leptonic decays and 0 for hadronic decays 
// uses product of probabilities
double MissingMassCalculator::AngularProbability(TLorentzVector nu_vec, TLorentzVector vis_vec, int decayType) {
  double prob=0.0;
  double M=1.777;
  double angl=0.0;
  double P=(vis_vec+nu_vec).P();
  double V=P/sqrt(P*P+M*M); // tau speed
  double dA=dRmax_tau/(2.0*Niter_fit1);

  if(decayType==1) // leptonic tau for now
    {
      // exact formular for energy probability is sqrt(1-V^2)/(2*V*p_0)*dE
      double m_1=nu_vec.M();
      double m_2=vis_vec.M();
      double E_nu=(M*M+m_1*m_1-m_2*m_2)/(2.0*M);
      if(E_nu<=m_1) return 0.0;
      double P_nu=sqrt(pow(E_nu,2)-pow(m_1,2));
      double prob1=0.5*sqrt(1-V*V)/(P_nu*V); // energy probability

      angl=Angle(vis_vec,vis_vec+nu_vec); // using lepton direction
      double det1=1.0-V*cos(angl+dA);
      double det2= (angl-dA)>0.0 ? 1.0-V*cos(angl-dA) : 1.0-V;
      double prob2=fabs(1.0/det1-1.0/det2)*(1.0-V*V)/(2.0*V); // using massless approximation for leptons
      prob=prob1*prob2;
    }
  if(decayType==0) // hadronic tau 
    {
      // exact formula for energy probability is sqrt(1-V^2)/(2*V*p_0)*dE
      // drop p_0 because it's a contstant for a given hadronic tau
      double prob1=0.5*sqrt(1-V*V)/V; // "energy" probability

      angl=Angle(nu_vec,vis_vec+nu_vec); // using neutrino direction
      double det1=1.0-V*cos(angl+dA);
      double det2= (angl-dA)>0.0 ? 1.0-V*cos(angl-dA) : 1.0-V;
      double prob2=fabs(1.0/det1-1.0/det2)*(1.0-V*V)/(2.0*V);
      prob=prob1*prob2;
    }
  return prob;
}

//----- checks kinematics of solutions
// returns 1 if kinematics consistent with tau, 0 otherwise
int MissingMassCalculator::CheckSolutions(TLorentzVector nu_vec, TLorentzVector vis_vec, int decayType)
{
  int passcode=1;
  double M=1.777;
  double P=(vis_vec+nu_vec).P();
  double V=P/sqrt(P*P+M*M); // tau speed
  if(decayType==0) // hadronic tau
    {
      double m=vis_vec.M();
      if(m>M) m=2.0*M-m; // artificial hack: visible mass can be larger than M due to smearing
      if(m<=0.0) m=0.0; // this should not happen, but if it does things are already screwed up 
      double E_nu=(M*M-m*m)/(2.0*M);
      double E_vis=(M*M+m*m)/(2.0*M);
      double P_vis=sqrt(pow(E_vis,2)-m*m);
      double Enu_lim[2];
      Enu_lim[0]=E_nu*sqrt((1.0-V)/(1.0+V));
      Enu_lim[1]=E_nu*sqrt((1.0+V)/(1.0-V));
      if(nu_vec.E()<Enu_lim[0] || nu_vec.E()>Enu_lim[1])
	{
 	  if(fUseVerbose==1) std::cout<<" DEBUG: MissingMassCalculator::CheckSolutions --- neutrino energy outside range "<<std::endl;
	  return 0;
	}
      double Evis_lim[2];
      Evis_lim[0]=(E_vis-V*P_vis)/sqrt(1.0-V*V);
      Evis_lim[1]=(E_vis+V*P_vis)/sqrt(1.0-V*V);
      if(vis_vec.E()<Evis_lim[0] || vis_vec.E()>Evis_lim[1])
	{
	  if(fUseVerbose==1) std::cout<<" DEBUG: MissingMassCalculator::CheckSolutions --- visTau energy outside range "<<std::endl;
	  return 0;
	}         
    }
  if(decayType==1) // leptonic tau
    {
      if((nu_vec+vis_vec).M()>M) return 0;
      if((nu_vec+vis_vec).M()<0.9*M) return 0;
      double m1=nu_vec.M();
      double m2=vis_vec.M();
      if(m2>M) m2=2.0*M-m2; // artificial hack: visible mass can be larger than M due to smearing
      if(m2<=0.0) m2=0.0; // this should not happen, but if it does things are already screwed up 
      double E_nu=(M*M+m1*m1-m2*m2)/(2.0*M);
      if(E_nu<=m1) return 0;
      double P_nu=sqrt(pow(E_nu,2)-pow(m1,2));
      double E_vis=(M*M+m2*m2-m1*m1)/(2.0*M);
      if(E_vis<=m2) return 0;
      double P_vis=sqrt(pow(E_vis,2)-pow(m2,2));
      double Enu_lim[2];
      Enu_lim[0]=(E_nu-V*P_nu)/sqrt(1.0-V*V);
      Enu_lim[1]=(E_nu+V*P_nu)/sqrt(1.0-V*V);
      if(nu_vec.E()<Enu_lim[0] || nu_vec.E()>Enu_lim[1])
	{
 	  if(fUseVerbose==1) std::cout<<" DEBUG: MissingMassCalculator::CheckSolutions --- neutrino energy outside range "<<std::endl;
	  return 0;
	}
      double Evis_lim[2];
      Evis_lim[0]=(E_vis-V*P_vis)/sqrt(1.0-V*V);
      Evis_lim[1]=(E_vis+V*P_vis)/sqrt(1.0-V*V);
      if(vis_vec.E()<Evis_lim[0] || vis_vec.E()>Evis_lim[1])
	{
 	  if(fUseVerbose==1) std::cout<<" DEBUG: MissingMassCalculator::CheckSolutions --- visTau energy outside range "<<std::endl;
	  return 0;
	}     
    }

  return passcode;
}

// returns Mnu probability according pol6 parameterization
double MissingMassCalculator::MnuProbability(double mnu, double binsize)
{
  double prob=1.0;
  double norm=4851900.0;
  double p[7];
  p[0]=-288.6/norm; p[1]=6.217E4/(2.0*norm); p[2]=2.122E4/(3.0*norm); p[3]=-9.067E4/(4.0*norm); 
  p[4]=1.433E5/(5.0*norm); p[5]=-1.229E5/(6.0*norm); p[6]=3.434E4/(7.0*norm);
  double int1=0.0;
  double int2=0.0;
  double x1= mnu+0.5*binsize < 1.777-0.113 ? mnu+0.5*binsize : 1.777-0.113;
  double x2= mnu-0.5*binsize > 0.0 ? mnu-0.5*binsize : 0.0;
  for(int i=0; i<7; i++)
    {
      int1=p[i]*pow(x1,i+1)+int1;
      int2=p[i]*pow(x2,i+1)+int2;
    }
  prob=int1-int2;
  if(prob<0.0) 
    {
      if(fUseVerbose==1) std::cout<<" Warning in MissingMassCalculator::MnuProbability: negative probability!!! "<<std::endl;
      return 0.0;
    }
  if(prob>1.0)
    {
      if(fUseVerbose==1) std::cout<<" Warning in MissingMassCalculator::MnuProbability: probability > 1!!! "<<std::endl;
      return 1.0;
    }
  return prob;
}

double MissingMassCalculator::Angle(const TLorentzVector & vec1, const TLorentzVector & vec2) {
  //SpeedUp (both are equivalent in fact)
  return acos((vec1.Px()*vec2.Px()+vec1.Py()*vec2.Py()+vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P()));
  //return vec1.Angle(vec2.Vect());
}



// returns probability of angle between two tau's
// assuming massless tau's for now, should be small effect for M>M_Z
double MissingMassCalculator::ResonanceProbability(TLorentzVector vec1, TLorentzVector vec2) {

  double prob=1.0;
  double boson_P=(vec1+vec2).P();
  if(boson_P==0.0) return 1.0;
  double boson_E=(vec1+vec2).E();
  double boson_V=0.0;
  if(boson_E>0.0) boson_V=boson_P/boson_E;
  else return 1.0E-10;

  double testMass=(vec1+vec2).M();
  double m=1.777; // tau mass
  double E_tau=testMass/2.0;
  double P_tau=sqrt(pow(E_tau,2)-m*m);
  double Evis_lim[2];
  Evis_lim[0]=(E_tau-boson_V*P_tau)/sqrt(1.0-boson_V*boson_V);
  Evis_lim[1]=(E_tau+boson_V*P_tau)/sqrt(1.0-boson_V*boson_V);
  if(vec1.E()<Evis_lim[0] || vec1.E()>Evis_lim[1]) return 1.0E-20;
  if(vec2.E()<Evis_lim[0] || vec2.E()>Evis_lim[1]) return 1.0E-20;

  double prob1=0.5*sqrt(1-boson_V*boson_V)/(P_tau*boson_V);

  if(vec1.P()*vec2.P()>0)
    {
      double theta=acos((vec1.Px()*vec2.Px()+vec1.Py()*vec2.Py()+vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P()));
      if(boson_V>0.0 && boson_V<1.0)
	{
	  if(boson_V<cos(theta/2)) return 1.0E-10;
	  double num=(1.0-boson_V*boson_V)*cos(theta/2);
	  double denom=4*boson_V*sin(theta/2)*sin(theta/2)*sqrt(boson_V*boson_V-cos(theta/2)*cos(theta/2));
	  prob=num/denom; 
	}
      else
	{
	  if(fabs(theta-TMath::Pi())>0.0001) return 1.0E-10;
	}
    }
  else return 1.0E-10;
  prob=prob*prob1;
  if(prob<1.0E-20) prob=1.0E-20;
  return prob;
}


bool MissingMassCalculator::refineSolutions (const TLorentzVector & tau_vec1, const int & tau_type1, 
					     const TLorentzVector & tau_vec2, const int & tau_type2, 
					     const double & probScale,
					     std::vector<TLorentzVector> & tauvecsol1,
					     std::vector<TLorentzVector> & tauvecsol2,
					     int nsol1, int nsol2 ) {


  
  bool oneSol=false;
  
  std::vector<double> tauvecprob1(4);
  std::vector<double> tauvecprob2(4);
  TLorentzVector tautau_tmp(0.0,0.0,0.0,0.0);
 

  iter2+=(nsol1*nsol2);
  int ngoodsol1=0;
  for(int j1=0; j1<nsol1; j1++)
    {
      //SpeedUp reference to avoid picking again in again in array
      TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
      TLorentzVector & tauvecsol1j=tauvecsol1[j1];
      double &tauvecprob1j = tauvecprob1[j1];
      tauvecprob1j=0.;
					      

      if((tau_vec1+nuvec1_tmpj).E()>=beamEnergy) continue;

      tauvecsol1j.SetPxPyPzE(0.,0.,0.,0.);
      tauvecsol1j+=nuvec1_tmpj;
      tauvecsol1j+=tau_vec1;
      const double tau1_tmpp=tauvecsol1j.P();
      const double angle1=Angle(nuvec1_tmpj,tau_vec1);
      if(angle1<dTheta3DLimit(tau_type1,0,tau1_tmpp)){++iang1low; continue; }// lower 99% bound
      if(angle1>dTheta3DLimit(tau_type1,1,tau1_tmpp)){++iang1high; continue;} // upper 99% bound
      tauvecprob1j=dTheta3d_probabilityFast(tau_type1,angle1,tau1_tmpp); 
      ++ngoodsol1;


      if (fSpeedStudy)
	{
	  ++iter3;
	  if (iter3 % 10000 == 1) {
	    const double aux=dTheta3d_probability(tau_type1,angle1,tau1_tmpp);
	    const double ratio=std::abs(aux-tauvecprob1j)/(aux+tauvecprob1j);
	    if  (ratio>=1E-3){
	      std::cout << "SpeedUp WARNING mismatch between quick computation and slow "  << iter3 << " ref-1 " << tauvecprob1j << " check " << aux << std::endl;
	    } 
	  }
	}
    }
  
 
  if (ngoodsol1==0) return false;			      
  int ngoodsol2=0;
  
  for(int j2=0; j2<nsol2; j2++)
    {
      TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];	
      TLorentzVector & tauvecsol2j=tauvecsol2[j2];
      double &tauvecprob2j = tauvecprob2[j2];
      tauvecprob2j=0.;

      tauvecsol2j.SetPxPyPzE(0.,0.,0.,0.);
      tauvecsol2j+=nuvec2_tmpj;
      tauvecsol2j+=tau_vec2;
      if (tauvecsol2j.E()>=beamEnergy) continue;
      const double tau2_tmpp=tauvecsol2j.P();
				  //redundant if(tau_vec1.E()+nuvec1_tmpj.E()+tau_vec2.E()+nuvec2_tmpj.E()>=2.0*beamEnergy) continue;
      const double angle2=Angle(nuvec2_tmpj,tau_vec2);
      if(angle2<dTheta3DLimit(tau_type2,0,tau2_tmpp)){++iang2low; continue;} // lower 99% bound
      if(angle2>dTheta3DLimit(tau_type2,1,tau2_tmpp)){++iang2high; continue;} // upper 99% bound
      tauvecprob2j=dTheta3d_probabilityFast(tau_type2,angle2,tau2_tmpp);
      ++ngoodsol2;
      
      if (fSpeedStudy)
	{
	  ++iter3;
	  if (iter3 % 10000 == 1) {
	    const double aux2=dTheta3d_probability(tau_type2,angle2,tau2_tmpp);
	    const double ratio2=std::abs(aux2-tauvecprob2j)/(aux2+tauvecprob2j);
	    if  (ratio2>=1E-3){
	      std::cout << "SpeedUp WARNING mismatch between quick computation and slow "  << iter3 << " ref-2 " << tauvecprob2j << " check " << aux2 << std::endl;
	    } 
	  }
	}
    }

 
  if (ngoodsol2==0) return false;			      
			      
  // now reloop
  for (int j1=0; j1<nsol1;++j1)
    {
      double &tauvecprob1j = tauvecprob1[j1];
      if (tauvecprob1j==0.) continue;
      TLorentzVector & nuvec1_tmpj=nuvecsol1[j1];
      TLorentzVector & tauvecsol1j=tauvecsol1[j1];
					  
      for (int j2=0; j2<nsol2;++j2)
	{
	  double &tauvecprob2j = tauvecprob2[j2];
	  if (tauvecprob2j==0.) continue;
	  TLorentzVector & nuvec2_tmpj=nuvecsol2[j2];
	  TLorentzVector & tauvecsol2j=tauvecsol2[j2];

	  const double totalProb=tauvecprob1j*tauvecprob2j*probScale;
	  //cannot happen if(totalProb<=0.0) continue;

	  ++iter4;

	  
	  totalProbSum+=totalProb;				      
	  tautau_tmp.SetPxPyPzE(0.,0.,0.,0.);
	  tautau_tmp+=tauvecsol1j;
	  tautau_tmp+=tauvecsol2j;
	  const double mtautau=tautau_tmp.M();
	  totalTautauSum+=mtautau;
	  //std::cout << iter1 << " " << mtautau << " totalprob " << totalProb << " tauvecsolE " << tauvecsol1[0].E() << " " << tauvecsol2[0].E()<< " tauvecsolM " << tauvecsol1[0].M() << " " << tauvecsol2[0].M() << std::endl;

	  const double sign_tmp=-log10(totalProb);
	  oneSol=true;
	  fMfit_all->Fill(mtautau,totalProb);
	  fPXfit_nu1->Fill(nuvec1_tmpj.Px(),totalProb); 
	  fPYfit_nu1->Fill(nuvec1_tmpj.Py(),totalProb); 
	  fPZfit_nu1->Fill(nuvec1_tmpj.Pz(),totalProb); 
	  fPXfit_nu2->Fill(nuvec2_tmpj.Px(),totalProb); 
	  fPYfit_nu2->Fill(nuvec2_tmpj.Py(),totalProb); 
	  fPZfit_nu2->Fill(nuvec2_tmpj.Pz(),totalProb); 
	  if(totalProb>prob_tmp) // fill solution with highest probability                                   
	    {
	      prob_tmp=totalProb;
	      fDitauStuffFit.Mditau_best=mtautau; 
	      fDitauStuffFit.Sign_best=sign_tmp; 
	      fDitauStuffFit.nutau1=nuvec1_tmpj; 
	      fDitauStuffFit.nutau2=nuvec2_tmpj; 
	    }

	}
    }


 
  return oneSol;
  
}
	
int MissingMassCalculator::TailCleanUp(const int & type1, const TLorentzVector & vis1, const TLorentzVector & nu1, 
				       const int & type2, const TLorentzVector & vis2, const TLorentzVector & nu2, 
				       const double & mmc_mass, const double & vis_mass, const double & eff_mass, 
				       const double & dphiTT) {
  int pass_code=1;
  if(fUseTailCleanup==0) return pass_code;  
  //the Clean-up cuts are specifically for rel16 analyses.
  // the will change in rel17 analyses and after the MMC is updated

  if(type1==0 && type2==0) // lepton-lepton channel
    {
      const double MrecoMvis=mmc_mass/vis_mass;
      if(MrecoMvis>2.6) return 0; 
      const double MrecoMeff=mmc_mass/eff_mass;
      if(MrecoMeff>1.9) return 0;
      const double e1p1=nu1.E()/vis1.P();
      const double e2p2=nu2.E()/vis2.P();
      if((e1p1+e2p2)>4.5) return 0; 
      if(e2p2>4.0) return 0; 
      if(e1p1>3.0) return 0; 
    }

//   // These are old cuts for lep-had, commented out
//   if((type1==0 || type2==0) && (type1+type2>0)) // lepton-hadron channel
//     {
//       double MrecoMvis=mmc_mass/vis_mass;
//       if((type1==1 || type2==1) && MrecoMvis>2.2) return 0; // for lep + 1-prong, changed on 09/23/11, based on optimization studies
//       if((type1==3 || type2==3) && MrecoMvis>2.0) return 0; // for lep + 3-prong
//       double MrecoMeff=mmc_mass/eff_mass;
//       if((type1==1 || type2==1) && MrecoMeff>1.4) return 0; // for lep + 1-prong, added on 09/23/11
//       if((type1==3 || type2==3) && MrecoMeff>1.3) return 0; // for lep + 3-prong, added on 09/23/11
//       if(type1==0)
//         {
//           double e1p1=nu1.E()/vis1.P();
//           if(e1p1>2.5) return 0; // lepton, same for 1- and 3-prongs
//           double e2p2=nu2.E()/vis2.P();
//           if(type2==1 && e2p2>2.0) return 0; // hadron, for 1-prong only
//           if((e1p1+e2p2)>3.0) return 0; // same for 1- and 3-prongs
//         }
//       if(type2==0)
//         {
//           double e1p1=nu1.E()/vis1.P();
//           if(type1==1 && e1p1>2.0) return 0; // hadron, for 1-prong only
//           double e2p2=nu2.E()/vis2.P();   
//           if(e2p2>2.5) return 0; // lepton, same for 1- and 3-prongs
//           if((e1p1+e2p2)>3.0) return 0;  // same for 1- and 3-prongs
//         }
//     }


//-------- these are new cuts for lep-had analysis for Moriond
  if((type1==0 || type2==0) && (type1+type2>0)) // lepton-hadron channel
    {
      if(InputInfo.Njet25>0)
	{
	  const double MrecoMvis=mmc_mass/vis_mass;
	  const double MrecoMeff=mmc_mass/eff_mass;
	  const double x=dphiTT>1.5 ? dphiTT : 1.5;
	  if((MrecoMeff+MrecoMvis)>5.908-1.881*x+0.2995*x*x) return 0;
	}
    }

  if(type1>0 && type2>0) // hadron-hadron channel
    {
      const double MrecoMvis=mmc_mass/vis_mass;
      if(MrecoMvis>1.7) return 0; 
      const double e1p1=nu1.E()/vis1.P();
      const double e2p2=nu2.E()/vis2.P();
      if((e1p1+e2p2)>1.5) return 0; 
      if(e2p2>1.2) return 0; 
    }

  return pass_code;
}

double MissingMassCalculator::TauProbability(const int & type1, const TLorentzVector & vis1, const TLorentzVector & nu1, 
					     const int & type2, const TLorentzVector & vis2, const TLorentzVector & nu2) {
  double prob=1.0;
  double prob1=1.0;
  double prob2=1.0;
  const double mtau=1.777;
  const double R1=nu1.E()/vis1.E();
  const double R2=nu2.E()/vis2.E();
  //--- dealing with 1st tau
  double m1=nu1.M();
  double m2=vis1.M();
  double E1=0.5*(mtau*mtau+m1*m1-m2*m2)/mtau;
  double E2=mtau-E1;
  if(E1<=m1 || E1>=mtau) 
    {
      if(fUseVerbose==1) std::cout<<" Warning in MissingMassCalculator::TauProbability: bad E1, returning 0 "<<std::endl;
      return 0.0;
    }
  if(E2<=m2 || E2>=mtau)
    {
      if(fUseVerbose==1) std::cout<<" Warning in MissingMassCalculator::TauProbability: bad E2, returning 0 "<<std::endl;
      return 0.0;
    }
  double p=(nu1+vis1).P();
  double V=p/sqrt(p*p+mtau*mtau);
  double p0;
  if(type1==0) p0=sqrt(E2*E2-m2*m2); // leptonic tau
  else p0=E1; // hadronic tau
  prob1=0.5*mtau/(p0*V*pow(R1+1.0,2));

  //--- dealing with 2nd tau
  m1=nu2.M();
  m2=vis2.M();
  E1=0.5*(mtau*mtau+m1*m1-m2*m2)/mtau;
  E2=mtau-E1;
  if(E1<=m1 || E1>=mtau) 
    {
      if(fUseVerbose==1) std::cout<<" Warning in MissingMassCalculator::TauProbability: bad E1, returning 0 "<<std::endl;
      return 0.0;
    }
  if(E2<=m2 || E2>=mtau)
    {
      if(fUseVerbose==1) std::cout<<" Warning in MissingMassCalculator::TauProbability: bad E2, returning 0 "<<std::endl;
      return 0.0;
    }
  p=(nu2+vis2).P();
  V=p/sqrt(p*p+mtau*mtau);
  if(type2==0) p0=sqrt(E2*E2-m2*m2); // leptonic tau
  else p0=E1; // hadronic tau
  prob2=0.5*mtau/(p0*V*pow(R2+1.0,2));
  prob=prob1*prob2;
  return prob;
}

// --------- Updated version of TauProbability for lep-had events with Njet25=0, takes into account Winter-2012 analysis cuts
double MissingMassCalculator::TauProbability(const int & type1, const TLorentzVector & vis1, const TLorentzVector & nu1, 
					     const int & type2, const TLorentzVector & vis2, const TLorentzVector & nu2, const double & detmet) {
  double prob=1.0;
  if(InputInfo.Njet25==0)
    {
      if(detmet<20.0) // low MET Njet=0 category 
	{
	  const double R1=nu1.P()/vis1.P();
	  const double R2=nu2.P()/vis2.P();
	  const double lep_p1[4]={0.417,0.64,0.52,0.678};
	  const double lep_p2[4]={0.23,0.17,0.315,0.319};
	  const double lep_p3[4]={0.18,0.33,0.41,0.299};
	  const double lep_p4[4]={0.033,0.109,0.129,0.096};
	  const double lep_p5[4]={0.145,0.107,0.259,0.304};
	  int ind=3;
	  int indT=3;
	  const double m_1pr[4]={-0.15,-0.13,-0.25,-0.114};
	  const double s_1pr[4]={0.40,0.54,0.62,0.57};
	  const double m_3pr[4]={-1.08,-1.57,-0.46,-0.39};
	  const double s_3pr[4]={0.53,0.85,0.61,0.53};
	  double Ptau=0.0;
	  double Plep=0.0;
	  if(type1==1 || type1==3) 
	    {
	      Ptau=(nu1+vis1).P();
	      Plep=(nu2+vis2).P();
	    }
	  if(type2==1 || type2==3) 
	    {
	      Ptau=(nu2+vis2).P();
	      Plep=(nu1+vis1).P();
	    }
	  if(Plep<50.0 && Plep>=45.0) ind=2;
	  if(Plep<45.0 && Plep>=40.0) ind=1;
	  if(Plep<40.0) ind=0;
	  if(Ptau<50.0 && Ptau>=45.0) indT=2;
	  if(Ptau<45.0 && Ptau>=40.0) indT=1;
	  if(Ptau<40.0) indT=0;
	  if(type1==0) prob=prob*(lep_p5[ind]*TMath::Gaus(R1,lep_p1[ind],lep_p2[ind])+TMath::Landau(R1,lep_p3[ind],lep_p4[ind]))/(1+lep_p5[ind]);
	  if(type2==0) prob=prob*(lep_p5[ind]*TMath::Gaus(R2,lep_p1[ind],lep_p2[ind])+TMath::Landau(R2,lep_p3[ind],lep_p4[ind]))/(1+lep_p5[ind]);
	  if(type1==1) prob=prob*TMath::Gaus(R1,m_1pr[indT],s_1pr[indT]);
	  if(type2==1) prob=prob*TMath::Gaus(R2,m_1pr[indT],s_1pr[indT]);
	  if(type1==3) prob=prob*TMath::Gaus(R1,m_3pr[indT],s_3pr[indT]);
	  if(type2==3) prob=prob*TMath::Gaus(R2,m_3pr[indT],s_3pr[indT]); 
	}
      else // high MET Njet=0 category
	{
	  const double R1=nu1.P()/vis1.P();
	  const double R2=nu2.P()/vis2.P();
	  const double lep_p1[4]={0.441,0.64,0.79,0.8692};
	  const double lep_p2[4]={0.218,0.28,0.29,0.3304};
	  const double lep_p3[4]={0.256,0.33,0.395,0.4105};
	  const double lep_p4[4]={0.048,0.072,0.148,0.1335};
	  const double lep_p5[4]={0.25,0.68,0.10,0.2872};
	  int ind=3;
	  const double p_1prong=-3.706;
	  const double p_3prong=-5.845;
	  double Ptau=0.0;
	  double Plep=0.0;
	  if(type1==1 || type1==3) 
	    {
	      Ptau=(nu1+vis1).P();
	      Plep=(nu2+vis2).P();
	    }
	  if(type2==1 || type2==3) 
	    {
	      Ptau=(nu2+vis2).P();
	      Plep=(nu1+vis1).P();
	    }
	  if(Plep<50.0 && Plep>=45.0) ind=2;
	  if(Plep<45.0 && Plep>=40.0) ind=1;
	  if(Plep<40.0) ind=0;
	  const double scale1prong=Ptau>45.0 ? 1.0 : -1.019/((Ptau*0.0074-0.059)*p_1prong);
	  const double scale3prong=Ptau>40.0 ? 1.0 : -1.24/((Ptau*0.0062-0.033)*p_3prong);
	  if(type1==0) prob=prob*(lep_p5[ind]*TMath::Gaus(R1,lep_p1[ind],lep_p2[ind])+TMath::Landau(R1,lep_p3[ind],lep_p4[ind]))/(1+lep_p5[ind]);
	  if(type2==0) prob=prob*(lep_p5[ind]*TMath::Gaus(R2,lep_p1[ind],lep_p2[ind])+TMath::Landau(R2,lep_p3[ind],lep_p4[ind]))/(1+lep_p5[ind]);
	  if(type1==1) prob=prob*exp(p_1prong*R1*scale1prong);
	  if(type2==1) prob=prob*exp(p_1prong*R2*scale1prong);
	  if(type1==3) prob=prob*exp(p_3prong*R1*scale3prong);
	  if(type2==3) prob=prob*exp(p_3prong*R2*scale3prong);
	}
    }
  return prob;
}


//-------- This function applies correction to compensate for the off-set
double MissingMassCalculator::MassScale(int method, double mass, const int & tau_type1, const int & tau_type2) {
  double Fscale=1.0;
  // calibration for rel16 lep-had analysis only
  if(fApplyMassScale==1) 
    {
      if((tau_type1+tau_type2)>0 && (tau_type1==0 || tau_type2==0))
	{
	  if(method!=1) return 1.0;
// 	  float p0, p1, p2, p3;
// 	  if(tau_type1==1 || tau_type2==1) // 1-prong tau's
// 	    {
// 	      p0=3.014; p1=-71.86; p2=1.018; p3=0.8912;
// 	      if(mass>91.2) Fscale=p0/(p1+p2*mass)+p3;
// 	      else Fscale=p0/(p1+p2*91.2)+p3;
// 	    }
// 	  if(tau_type1==3 || tau_type2==3) // 3-prong tau's
// 	    {
// 	      p0=0.4576; p1=-84.22; p2=0.9783; p3=0.9136;
// 	      if(mass>91.2) Fscale=p0/(p1+p2*mass)+p3;
// 	      else Fscale=p0/(p1+p2*91.2)+p3;
// 	    }
// 	  if(Fscale>1.0) Fscale=1.0;
// 	  if(Fscale<0.89) Fscale=0.89;

	  float p0, p1, p2, p3, p4, p5, p6, p7;
	  if(tau_type1==1 || tau_type2==1) return 1.0; // 1-prong tau's
	  if(tau_type1==3 || tau_type2==3) // 3-prong tau's
	    {
	      p0=3.014; p1=-71.86; p2=1.018; p3=0.8912;
	      p4=0.4576; p5=-84.22; p6=0.9783; p7=0.9136;
	      double scale1=p0/(p1+p2*mass)+p3;
	      double scale3=p4/(p5+p6*mass)+p7;
	      if(mass>91.2) Fscale=scale3/scale1;
	      else 
		{
		  scale1=p0/(p1+p2*91.2)+p3;
		  scale3=p4/(p5+p6*91.2)+p7;
		  Fscale=scale3/scale1;
		}
	    }
	  if(Fscale>1.0) Fscale=1.0;
	  if(Fscale<0.95) Fscale=0.95;
	}
    }
  return 1.0/Fscale;
}

double MissingMassCalculator::MHtProbability(const double & d_mhtX, const double & d_mhtY, const double & mht, 
					     const double & trueMetGuess, const double & mht_offset) {
  double prob=1.0;
  if(MHtSigma1>0.0 && MHtSigma2>0.0 && MHtGaussFr>0.0)
    {
      prob=(exp(-0.5*pow(d_mhtX/MHtSigma1,2))+MHtGaussFr*exp(-0.5*pow(d_mhtX/MHtSigma2,2)));
      prob=prob*(exp(-0.5*pow(d_mhtY/MHtSigma1,2))+MHtGaussFr*exp(-0.5*pow(d_mhtY/MHtSigma2,2)));
      const double _arg=(mht-trueMetGuess-mht_offset)/MHtSigma1;
      prob=prob*exp(-0.25*pow(_arg,2)); // assuming sqrt(2)*sigma
    }
  return prob;
}

}
}
