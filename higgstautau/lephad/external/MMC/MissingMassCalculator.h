/*
 MissingMassCalculator
 Author: Aliaksandr Pranko (appranko@lbl.gov)
 MissingMassCalculator is designed to reconstruct mass in
 events where two particles decay into states with missing ET.
*/
#ifndef MyMissingMassCalculator_h
#define MyMissingMassCalculator_h

// #if !defined (__CINT__) || defined (__MAKECINT__)

#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <vector>
#include <TObject.h>
/* #include "JetResolution/JERProvider.h" */
/* #include "JetResolution/JERProviderAthena.h" */
#include "JERProviderMMC.h"

class TH1F;


// #endif

namespace AnalysisFramework
{
namespace External
{

class MissingMassCalculator: public TObject {
 public:
    
//    #if !defined(__CINT__) | defined(__MAKECINT__)
//        // Define the class for the cint dictionary
        ClassDef(MissingMassCalculator,1);
//    #endif

protected:

  //---------------- structutures
  struct DitauStuff {
    double Mditau_best; // best fitted M(ditau)
    double Sign_best;   // best significance of M(ditau) fit
    TLorentzVector nutau1;  // fitted 4-vec for neutrino from tau-1
    TLorentzVector nutau2;  // fitted 4-vec for neutrino from tau-2
    double RMSoverMPV;
  };

  struct InputInfoStuff {
    TVector2 MetVec;
    TLorentzVector vistau1;
    TLorentzVector vistau2;
    int type_visTau1;
    int type_visTau2;
    int Nprong_tau1;
    int Nprong_tau2;
    int dataType;
    double phi_jet;
    double METsigmaP;
    double METsigmaL;
    double SumEt;
    std::vector<TLorentzVector> jet4vecs;
    int Njet25;
    double DelPhiTT;
  };

  struct OutputInfoStuff {
    int FitStatus;
    double FitSignificance[2];
    double FittedMass[2];
    TLorentzVector nuvec1[2];
    TLorentzVector objvec1[2];
    TLorentzVector nuvec2[2];
    TLorentzVector objvec2[2];
    TLorentzVector totalvec[2];
    TVector2 FittedMetVec[2];
    double RMS2MPV;
  };
  
  //SpeedUp static array for efficient access
  static double fit_param[3][6][4];
  //cache quantities for efficient NuPSolution calculation
  double phi1cache,phi2cache,cosphi1cache,cosphi2cache,sinphi1cache, sinphi2cache,dsincache;
  // for superfast version
  double metvecpxcache,metvecpycache,vec1ecache,vec1mcache,vec1m2cache,vec2ecache,vec2mcache,vec2m2cache;
  double fcache0,fcache1,Evcache[2],tBcache[2],ecache[2],Rcache[2],ffRRcache[2];
  std::vector<TLorentzVector> nuvecsol1;
  std::vector<TLorentzVector> nuvecsol2;

  int iter1,iter2,iter3,iter4,iang1low,iang1high,iang2low,iang2high;
  double prob_tmp;
  
  double totalProbSum;
  double totalTautauSum;
  
  //--- define histograms for histogram method
  //--- upper limits need to be revisied in the future!!! It may be not enough for some analyses
  TH1F* fMfit_all;
  TH1F* fPXfit_nu1;
  TH1F* fPYfit_nu1;
  TH1F* fPZfit_nu1;
  TH1F* fPXfit_nu2;
  TH1F* fPYfit_nu2;
  TH1F* fPZfit_nu2;

  JERProviderMMC *myJER;

  // for intermediate calc
  TLorentzVector TLVdummy;

  //---------------- protected variables
  DitauStuff fDitauStuffFit; // results based on fit method
  DitauStuff fDitauStuffHisto; // results based on histo method
  InputInfoStuff InputInfo; // input info
  OutputInfoStuff OutputInfo; // output info

  int fUseVerbose; // code to turn ON printouts for debugging
  bool fSpeedStudy; // code to turn ON speed study
  int AlgorithmVersion; // version of the algorithm
  int SearchMode; // search Mode: 0=di-tau, 1=WW, 2=W->tau+nu
  int fApplyMassScale; // switch to apply mass scale correction
  int fUseTailCleanup; // switch to apply tail clean-up

  int Niter_fit1; // number of iterations for dR-dPhi scan 
  int Niter_fit2; // number of iterations for MET-scan 
  int Niter_fit3; // number of iterations for Mnu-scan 
  int fJERsyst; // switch for JER systematics
  int METresSyst; // switch to turn on/off MET resolution systematics
  double dTheta3d_binMin; // minimal step size for dTheta3D
  double dTheta3d_binMax; // maximum step size for dTheta3D
  double dRmax_tau; // maximum dR(nu-visTau)
  double Nsigma_METscan; // number of sigmas for MET-scan
  double beamEnergy; // beam energy (Tevatron=980, LHC-1=3500.0) 
  int METScanScheme; // MET-scan scheme: 0- use JER; 1- use simple sumEt & missingHt for Njet=0 events in (lep-had) 
  double MHtSigma1; // sigma of 1st Gaussian in missing Ht resolution
  double MHtSigma2; // sigma of 2nd Gaussian in missing Ht resolution
  double MHtGaussFr; // relative fraction of 2nd Gaussian

  //---------------- protected functions
  void ClearDitauStuff(DitauStuff &fStuff);
  void ClearInputStuff(InputInfoStuff &fStuff);
  void ClearOutputStuff(OutputInfoStuff &fStuff);
  void DoOutputInfo(OutputInfoStuff &fStuff);
  void DoMetResolution(InputInfoStuff &fStuff);
  void PrintResults(OutputInfoStuff fStuff);
  void FinalizeInputStuff(InputInfoStuff &fStuff);
  int NuPsolution(TVector2 met_vec, double theta1, double phi1, 
		  double theta2, double phi2, double &P1, double &P2); // keep this version for simple tests

  //SpeedUp disable inlining (to get correct call graph from callgrind)
  #define inline
  inline int NuPsolutionV2(const TVector2 & met_vec, const TLorentzVector & vec1,const TLorentzVector & vec2, 
		    const double & mass1,const  double & mass2, const double & phi1,const  double & phi2, 
		    std::vector<TLorentzVector> &nu_vec1, std::vector<TLorentzVector> &nu_vec2);

  inline int NuPsolutionV2fast(const TVector2 & met_vec, const TLorentzVector & vec1,const TLorentzVector & vec2, 
		    const double & mass1,const  double & mass2, const double & phi1,const  double & phi2, 
				    int & nsol1,int & nsol2);
  
  double dTheta3Dfit_parameterization(int tau_type, int dT3dfit_par, int Pfit_par);
  inline double dTheta3Dparam(const int & parInd,const double & P_tau,const double *par); 
  inline double dTheta3d_probability(const int & tau_type,const double & dTheta3d,const  double & P_tau);
  inline double dTheta3d_probabilityFast(const int & tau_type,const double & dTheta3d,const double & P_tau);
  inline double MetProbability(const double & metX,const  double & metY, const double & MetSigma);
  inline double MetProbability(const double & met1,const  double & met2,const  double & MetSigma1, const double & MetSigma2);
  inline double dTheta3DLimit(const int & tau_type, const int & limit_code,const double & P_tau);
  inline double dTheta3DLimitVis(int tau_type, int limit_code, double visP_tau); 
  inline double dTheta2dPhi(double etaL, double dTheta);
  inline double MHtProbability(const double & d_mhtX, const double & d_mhtY, const double & mht, 
			       const double & trueMetGuess, const double & mht_offset);

  inline double Angle(const TLorentzVector & vec1, const TLorentzVector & vec2);
  double AngularProbability(TLorentzVector nu_vec, TLorentzVector vis_vec, int decayType);
  int CheckSolutions(TLorentzVector nu_vec, TLorentzVector vis_vec, int decayType);
  double MnuProbability(double mnu, double binsize);
  double ResonanceProbability(TLorentzVector vec1, TLorentzVector vec2);
  inline int TailCleanUp(const int & type1, const TLorentzVector & vis1, const TLorentzVector & nu1, 
			 const int & type2, const TLorentzVector & vis2, const TLorentzVector & nu2, 
			 const double & mmc_mass, const double & vis_mass, const double & eff_mass, const double & dphiTT);

  inline double TauProbability(const int & type1, const TLorentzVector & vis1, const TLorentzVector & nu1, 
			       const int & type2, const TLorentzVector & vis2, const TLorentzVector & nu2);
  inline double TauProbability(const int & type1, const TLorentzVector & vis1, const TLorentzVector & nu1, 
			       const int & type2, const TLorentzVector & vis2, const TLorentzVector & nu2, const double & detmet);


  bool refineSolutions (const TLorentzVector & tau_vec1, const int & tau_type1, 
						 const TLorentzVector & tau_vec2, const int & tau_type2, 
						      const double & probScale,
						      std::vector<TLorentzVector> & tauvecsol1,
						      std::vector<TLorentzVector> & tauvecsol2,
					       int nsol1, int nsol2 );
  double MassScale(int method, double mass, const int & tau_type1, const int & tau_type2);

  // V9: based on V7, dR-probability functions replaced by dTheta3D-probability functions
  // lep-lep channel is modified to have two dPhi scans, MET_x,y scans, and two Mnu scans. 
  int DitauMassCalculatorV9(const TLorentzVector & tau_vec1, const int & tau_type1, 
			    const TLorentzVector & tau_vec2, const int & tau_type2, 
			    const TVector2 & met_vec);  

  //
  // rewritten with same logic
  int DitauMassCalculatorV9fast(const TLorentzVector & tau_vec1, const int & tau_type1, 
			    const TLorentzVector & tau_vec2, const int & tau_type2, 
			    const TVector2 & met_vec);  

  // V10: based on V9, adaptive step size in dTheta3D transformed into dPhi scan.
  // adaptive step size for MET scan, also includes David's improvements
  int DitauMassCalculatorV10fast(const TLorentzVector & tau_vec1, const int & tau_type1, 
				 const TLorentzVector & tau_vec2, const int & tau_type2, 
				 const TVector2 & met_vec);  


 //----------------------------------------------
  //
  //   >>>>>>>>>>>>>   Public methods <<<<<<<<<<<< 
  //
  //______________________________________________

public:
  MissingMassCalculator(std::string JERProviderFile="JERProviderPlots.root");
  ~MissingMassCalculator();

  int RunMissingMassCalculator();
  
  //-------- Set Input Parameters
  void SetSearchMode(int val) { SearchMode=val; }
  void SetUseVerbose(int val) { fUseVerbose=val; }
  void SetSpeedStudy(int val) { fSpeedStudy=val; }
  void SetAlgorithmVersion(int val) { AlgorithmVersion=val; }

  void SetNiterFit1(int val) { Niter_fit1=val; } // number of iterations per loop in dPhi loop
  void SetNiterFit2(int val) { Niter_fit2=val; } // number of iterations per loop in MET loop
  void SetNiterFit3(int val) { Niter_fit3=val; } // number of iterations per loop in Mnu loop
  void SetMaxDRtau(double val) { dRmax_tau=val; } // max value of dR
  void SetNsigmaMETscan(double val) { Nsigma_METscan=val; } // number of sigma's for MET-scan
  void SetBeamEnergy(double val) { beamEnergy=val; } // beam energy
  void SetdTheta3d_binMax(double val) { dTheta3d_binMax=val; } // maximum step size for dTheta3D
  void SetdTheta3d_binMin(double val) { dTheta3d_binMin=val; } // minimal step size for dTheta3D

  void SetMetVec(TVector2 vec);
  void SetVisTauVec(int i, TLorentzVector vec);
  void SetVisTauType(int i, int tautype);
  void SetNprong(int i, int nprong);
  void SetSumEt(double sumEt);
  void SetIsData(int val);
  void SetMetScanParams(double phi, double sigmaP, double sigmaL);
  void SetMetScanParamsUE(double sumEt, double phi_scan=0.0, int data_code=0);
  void SetMetScanParamsJets(std::vector<TLorentzVector> jets);
  void SetJERsyst(int val) { fJERsyst=val; }
  void SetMETresSyst(int val) { METresSyst=val; } // MET resolution systematics: +/-1 sigma
  void SetApplyMassScale(int val) { fApplyMassScale=val; }
  void SetUseTailCleanup(int val) { fUseTailCleanup=val; }
  void SetNjet25(int val);
  void SetMETScanScheme(int val) { METScanScheme=val; }

  void SetMHtSigma1(double val) { MHtSigma1=val; }
  void SetMHtSigma2(double val) { MHtSigma2=val; } 
  void SetMHtGaussFr(double val) { MHtGaussFr=val; }

  //-------- Get results;
  int GetFitStatus(); // return fit status
  double GetFittedMass(int fitcode); // returns fitted Mass
  double GetRms2Mpv(); // returns RMS/MPV according to histogram method
  TLorentzVector GetNeutrino4vec(int fitcode, int ind); // returns neutrino 4-vec
  double GetFitSignificance(int fitcode); // returns fit significance
  TLorentzVector GetTau4vec(int fitcode, int ind); // returns full tau 4-vec
  TLorentzVector GetResonanceVec(int fitcode); // returns 4-vec for resonance 
  TVector2 GetFittedMetVec(int fitcode); // returns 2-vec for fitted MET
  int StandardCollApprox(TLorentzVector tau_vec1, TLorentzVector tau_vec2, TVector2 met_vec, double &Mrec); // standard collinear approximation


};

}
}

#endif
