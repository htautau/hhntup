#include "TMath.h"
#include "TMatrixD.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TVector2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TStyle.h"
#include "TROOT.h"
#include <TStopwatch.h>
#include <TDatime.h>
//#include "MissingMassCalculator.h"

namespace TER {


int myVBFcut(const TLorentzVector &tau1, const TLorentzVector &tau2, 
        const TLorentzVector &jet1, const TLorentzVector &jet2) {
    int passcode=1;
    if(jet1.Pt()<50.0) return 0;
    if(jet2.Pt()<30.0) return 0;
    if(fabs(jet1.Eta()-jet2.Eta())<3.0) return 0;
    if((tau1+tau2).M()<40.0) return 0;
    return passcode;
}


bool CollinearApproximationMatrix(const TLorentzVector &k1, const TLorentzVector &k2, 
        double METpx, double METpy, 
        double &mass, double &xp1, double &xp2){
    TMatrixD K(2, 2);
    K(0, 0) = k1.Px();      K(0, 1) = k2.Px();
    K(1, 0) = k1.Py();      K(1, 1) = k2.Py();
    if(K.Determinant()==0)
        return false;
    TMatrixD M(2, 1);
    M(0, 0) = METpx;
    M(1, 0) = METpy;
    TMatrixD Kinv = K.Invert();
    TMatrixD X(2, 1);      
    X = Kinv*M;     
    double X1 = X(0, 0);    double X2 = X(1, 0);
    double x1 = 1./(1.+X1); double x2 = 1./(1.+X2);               
    TLorentzVector p1 = k1*(1/x1);
    TLorentzVector p2 = k2*(1/x2);
    double m = (p1+p2).M(); 
    mass = m; xp1 = x1; xp2 = x2;                   
    return true;
}


int MyLepHadCuts(int useLHcut, int lepbit, int SLTbit, TLorentzVector tau0, TLorentzVector tau1, TVector2 met_vect, int Njet30, std::vector<TLorentzVector> jetvec)
{
    int pass=1;

    const double tau1Pt_slt_cut=26.0; // 1st tau is always a lepton, SLT cuts
    const double tau2Pt_slt_cut=20.0; // 1nd tau is always a hadronic tau, SLT cuts
    double tau1Pt_ltt_cut; // lepton, LTT cuts
    if(lepbit==0) tau1Pt_ltt_cut=20.0;
    if(lepbit==1) tau1Pt_ltt_cut=17.0;
    const double tau2Pt_ltt_cut=25.0; // hadronic tau, LLT cuts
    const double Mt_cut=70.0;
    const double PtH_cut=100.0;

    TLorentzVector met4vec(met_vect.X(),met_vect.Y(),0.0,met_vect.Mod());

    //---------------- MVA categories
    if(useLHcut==0) // preselection
    {
        if(SLTbit==0)
        {
            if(tau0.Pt()<tau1Pt_ltt_cut && tau0.Pt()>=tau1Pt_slt_cut) return 0;
            if(tau1.Pt()<tau2Pt_ltt_cut) return 0;
        }
        if(SLTbit==1)
        {
            if(tau0.Pt()<tau1Pt_slt_cut) return 0;
            if(tau1.Pt()<tau2Pt_slt_cut) return 0;
        }
        if((met4vec+tau0).M()>Mt_cut) return 0;
    }
    if(useLHcut==1) //  MVA Rest
    {
        if(SLTbit==0)
        {
            if(tau0.Pt()<tau1Pt_ltt_cut && tau0.Pt()>=tau1Pt_slt_cut) return 0;
            if(tau1.Pt()<tau2Pt_ltt_cut) return 0;
        }
        if(SLTbit==1)
        {
            if(tau0.Pt()<tau1Pt_slt_cut) return 0;
            if(tau1.Pt()<tau2Pt_slt_cut) return 0;
        }
        if((met4vec+tau0).M()>Mt_cut) return 0;
        if(Njet30>0) return 0;
        if(SLTbit==1)
        {
            if((met4vec+tau0+tau1).Pt()>=PtH_cut) return 0;
            if(Njet30>=2)
            {
                if(myVBFcut(tau0,tau1,jetvec[0],jetvec[1])==1) return 0;
            }      
        }
    }
    if(useLHcut==2) //  MVA 1-Jet
    {
        if(SLTbit==0)
        {
            if(tau0.Pt()<tau1Pt_ltt_cut && tau0.Pt()>=tau1Pt_slt_cut) return 0;
            if(tau1.Pt()<tau2Pt_ltt_cut) return 0;
        }
        if(SLTbit==1)
        {
            if(tau0.Pt()<tau1Pt_slt_cut) return 0;
            if(tau1.Pt()<tau2Pt_slt_cut) return 0;
        }
        if((met4vec+tau0).M()>Mt_cut) return 0;
        if(Njet30==0) return 0;
        if(SLTbit==1)
        {
            if((met4vec+tau0+tau1).Pt()>=PtH_cut) return 0;
            if(Njet30>=2)
            {
                if(myVBFcut(tau0,tau1,jetvec[0],jetvec[1])==1) return 0;
            }      
        }
    }

    if(useLHcut==3) //  MVA Boosted 
    {
        if(SLTbit==0) return 0;
        if(tau0.Pt()<tau1Pt_slt_cut) return 0;
        if(tau1.Pt()<tau2Pt_slt_cut) return 0;
        if((met4vec+tau0).M()>Mt_cut) return 0;
        if((met4vec+tau0+tau1).Pt()<PtH_cut) return 0;
        if(Njet30>=2)
        {
            if(myVBFcut(tau0,tau1,jetvec[0],jetvec[1])==1) return 0;
        }
    }

    if(useLHcut==4) //  MVA VBF 
    {
        if(SLTbit==0) return 0;
        if(tau0.Pt()<tau1Pt_slt_cut) return 0;
        if(tau1.Pt()<tau2Pt_slt_cut) return 0;
        if((met4vec+tau0).M()>Mt_cut) return 0;
        if(Njet30<2) return 0;
        if(myVBFcut(tau0,tau1,jetvec[0],jetvec[1])==0) return 0;
    }


    return pass;
}


double MyCrystallBall(double *x, double *par) {
    double fit=0.0;
    double arg=-(x[0]-par[0])/par[1];
    double n=par[3];
    double termA=pow(n/fabs(par[2]),n)*exp(-pow(par[2],2)/2.0);
    double termB=n/fabs(par[2])-fabs(par[2]);
    double termC=exp(-pow(par[2],2)/2.0)*n/(fabs(par[2])*(n-1));
    double termD=sqrt(TMath::Pi()/2.0)*(1.0+erf(fabs(par[2])/sqrt(2.0)));
    double termN=1.0/(par[1]*(termC+termD));
    double argP=termB-arg;
    if(arg>-par[2]) fit=termN*exp(-0.5*pow(arg,2.0));
    else fit=termN*termA/pow(fabs(argP),n);
    return fit;
}

double TERafla(int tau_type, int eta_bin, double Pt_scan)
{
    int type_ind=tau_type/3;
    double ter_alfa_par[2][10][6];

    ter_alfa_par[0][0][0]=1.42121;
    ter_alfa_par[0][0][1]=0.222565;
    ter_alfa_par[0][0][2]=42.8729;
    ter_alfa_par[0][0][3]=15;
    ter_alfa_par[0][0][4]=0.0;
    ter_alfa_par[0][0][5]=0.0;
    ter_alfa_par[0][1][0]=1.49137;
    ter_alfa_par[0][1][1]=0.186306;
    ter_alfa_par[0][1][2]=42.5663;
    ter_alfa_par[0][1][3]=15;
    ter_alfa_par[0][1][4]=0.00482475;
    ter_alfa_par[0][1][5]=0.0;  
    ter_alfa_par[0][2][0]=1.44915;
    ter_alfa_par[0][2][1]=0.284083;
    ter_alfa_par[0][2][2]=46.8241;
    ter_alfa_par[0][2][3]=15;
    ter_alfa_par[0][2][4]=0.00277167;
    ter_alfa_par[0][2][5]=0.0;
    ter_alfa_par[0][3][0]=1.30058;
    ter_alfa_par[0][3][1]=0.215696;
    ter_alfa_par[0][3][2]=42.5161;
    ter_alfa_par[0][3][3]=10.0432;
    ter_alfa_par[0][3][4]=0.00886606;
    ter_alfa_par[0][3][5]=0.0;
    ter_alfa_par[0][4][0]=1.30862;
    ter_alfa_par[0][4][1]=0.291135;
    ter_alfa_par[0][4][2]=41.935;
    ter_alfa_par[0][4][3]=8.89068;
    ter_alfa_par[0][4][4]=0.00213514;
    ter_alfa_par[0][4][5]=0.0;
    ter_alfa_par[0][5][0]=1.10272;
    ter_alfa_par[0][5][1]=0.321954;
    ter_alfa_par[0][5][2]=46.8587;
    ter_alfa_par[0][5][3]=15;
    ter_alfa_par[0][5][4]=0.0;
    ter_alfa_par[0][5][5]=0.0;
    ter_alfa_par[0][6][0]=0.834864;
    ter_alfa_par[0][6][1]=0.303157;
    ter_alfa_par[0][6][2]=43.2212;
    ter_alfa_par[0][6][3]=15;
    ter_alfa_par[0][6][4]=0.0366042;
    ter_alfa_par[0][6][5]=0.0;
    ter_alfa_par[0][7][0]=1.20809;
    ter_alfa_par[0][7][1]=0.342144;
    ter_alfa_par[0][7][2]=45.2857;
    ter_alfa_par[0][7][3]=15;
    ter_alfa_par[0][7][4]=0.0;
    ter_alfa_par[0][7][5]=0.0;  
    ter_alfa_par[0][8][0]=1.14769;
    ter_alfa_par[0][8][1]=0.194716;
    ter_alfa_par[0][8][2]=41.244;
    ter_alfa_par[0][8][3]=8.1812;
    ter_alfa_par[0][8][4]=0.00870686;
    ter_alfa_par[0][8][5]=0.0;
    ter_alfa_par[0][9][0]=1.10921;
    ter_alfa_par[0][9][1]=0.244763;
    ter_alfa_par[0][9][2]=41.4174;
    ter_alfa_par[0][9][3]=8.00588;
    ter_alfa_par[0][9][4]=0.00956358;
    ter_alfa_par[0][9][5]=0.0;  

    ter_alfa_par[1][0][0]=1.39132;
    ter_alfa_par[1][0][1]=0.272645;
    ter_alfa_par[1][0][2]=37.1485;
    ter_alfa_par[1][0][3]=7.10422;
    ter_alfa_par[1][0][4]=0.0;
    ter_alfa_par[1][0][5]=0.000594665;
    ter_alfa_par[1][1][0]=1.01585;
    ter_alfa_par[1][1][1]=0.429444;
    ter_alfa_par[1][1][2]=45.7367;
    ter_alfa_par[1][1][3]=40.0045;
    ter_alfa_par[1][1][4]=0.0346375;
    ter_alfa_par[1][1][5]=0.0;
    ter_alfa_par[1][2][0]=0.810001;
    ter_alfa_par[1][2][1]=0.508107;
    ter_alfa_par[1][2][2]=58.0322;
    ter_alfa_par[1][2][3]=29.7679;
    ter_alfa_par[1][2][4]=0.048485;
    ter_alfa_par[1][2][5]=0.0;
    ter_alfa_par[1][3][0]=0.878782;
    ter_alfa_par[1][3][1]=0.364332;
    ter_alfa_par[1][3][2]=57.3821;
    ter_alfa_par[1][3][3]=34.455;
    ter_alfa_par[1][3][4]=0.0355081;
    ter_alfa_par[1][3][5]=0.0;
    ter_alfa_par[1][4][0]=1.06262;
    ter_alfa_par[1][4][1]=0.22186;
    ter_alfa_par[1][4][2]=38.7808;
    ter_alfa_par[1][4][3]=5;
    ter_alfa_par[1][4][4]=0.0376183;
    ter_alfa_par[1][4][5]=0.0;
    ter_alfa_par[1][5][0]=1.08259;
    ter_alfa_par[1][5][1]=0.114432;
    ter_alfa_par[1][5][2]=38.6835;
    ter_alfa_par[1][5][3]=5;
    ter_alfa_par[1][5][4]=0.0202652;
    ter_alfa_par[1][5][5]=0.000254057;
    ter_alfa_par[1][6][0]=1.07754;
    ter_alfa_par[1][6][1]=0.0973914;
    ter_alfa_par[1][6][2]=42.8723;
    ter_alfa_par[1][6][3]=50;
    ter_alfa_par[1][6][4]=0.0;
    ter_alfa_par[1][6][5]=0.00111601;
    ter_alfa_par[1][7][0]=1.15652;
    ter_alfa_par[1][7][1]=0.180654;
    ter_alfa_par[1][7][2]=42.3473;
    ter_alfa_par[1][7][3]=5;
    ter_alfa_par[1][7][4]=0.0173592;
    ter_alfa_par[1][7][5]=0.0;
    ter_alfa_par[1][8][0]=1.16701;
    ter_alfa_par[1][8][1]=0.154514;
    ter_alfa_par[1][8][2]=38.1972;
    ter_alfa_par[1][8][3]=5;
    ter_alfa_par[1][8][4]=0.0;
    ter_alfa_par[1][8][5]=0.00103872;
    ter_alfa_par[1][9][0]=1.00791;
    ter_alfa_par[1][9][1]=0.225649;
    ter_alfa_par[1][9][2]=41.6963;
    ter_alfa_par[1][9][3]=5;
    ter_alfa_par[1][9][4]=0.0373697;
    ter_alfa_par[1][9][5]=0.0;  

    double x=Pt_scan;
    double arg=(x-ter_alfa_par[type_ind][eta_bin][2])/ter_alfa_par[type_ind][eta_bin][3];
    double alfa=ter_alfa_par[type_ind][eta_bin][0]+ter_alfa_par[type_ind][eta_bin][1]*exp(-0.5*arg*arg)
        +ter_alfa_par[type_ind][eta_bin][4]*sqrt(x)+ter_alfa_par[type_ind][eta_bin][5]*x;
    if(type_ind==1 && alfa<1.202) alfa=1.202;
    return alfa;
}

//-------------- simple Gaussian resolution, assuming perfect TES
// input tau_type=1 or 3;
// vec-- 4-vec of input detector tau
// Pt_scan value of tau Pt from TER scan
double TERSigma(int tau_type, int eta_bin, double Pt_scan)
{
    int type_ind=tau_type/3;

    // ter_sigma_par should be declared globally and intialized in constructer
    // but for now, I put it here

    double ter_sigma_par[2][10][3];

    ter_sigma_par[0][0][0]=0.412255; 
    ter_sigma_par[0][0][1]=0.0172081;
    ter_sigma_par[0][0][2]=0.598593; 
    ter_sigma_par[0][1][0]=0.345992;  
    ter_sigma_par[0][1][1]=0.0190559; 
    ter_sigma_par[0][1][2]=0.830882; 
    ter_sigma_par[0][2][0]=0.378877; 
    ter_sigma_par[0][2][1]=0.0182099;
    ter_sigma_par[0][2][2]=0.611812; 
    ter_sigma_par[0][3][0]=0.336333; 
    ter_sigma_par[0][3][1]=0.0276647;
    ter_sigma_par[0][3][2]=0.658391; 
    ter_sigma_par[0][4][0]=0.542694;  
    ter_sigma_par[0][4][1]=0.0154589; 
    ter_sigma_par[0][4][2]=0.0184377; 
    ter_sigma_par[0][5][0]=0.543344;  
    ter_sigma_par[0][5][1]=0.0147195; 
    ter_sigma_par[0][5][2]=0.180618;  
    ter_sigma_par[0][6][0]=0.181561;  
    ter_sigma_par[0][6][1]=0.0342868; 
    ter_sigma_par[0][6][2]=1.13322; 
    ter_sigma_par[0][7][0]=0.166558; 
    ter_sigma_par[0][7][1]=0.0148713;
    ter_sigma_par[0][7][2]=1.20254; 
    ter_sigma_par[0][8][0]=0.102669; 
    ter_sigma_par[0][8][1]=0.027544; 
    ter_sigma_par[0][8][2]=1.11598; 
    ter_sigma_par[0][9][0]=0.0751179; 
    ter_sigma_par[0][9][1]=0.0282496; 
    ter_sigma_par[0][9][2]=1.19985; 

    ter_sigma_par[1][0][0]=0.132539;
    ter_sigma_par[1][0][1]=0.041477;
    ter_sigma_par[1][0][2]=1.80224;
    ter_sigma_par[1][1][0]=0.285899;
    ter_sigma_par[1][1][1]=0.0285415;
    ter_sigma_par[1][1][2]=1.26812;
    ter_sigma_par[1][2][0]=0.389067;
    ter_sigma_par[1][2][1]=0.0255133;
    ter_sigma_par[1][2][2]=0.625702;
    ter_sigma_par[1][3][0]=0.480427;
    ter_sigma_par[1][3][1]=0.0258971;
    ter_sigma_par[1][3][2]=0.24598;
    ter_sigma_par[1][4][0]=0.29975;
    ter_sigma_par[1][4][1]=0.0397503;
    ter_sigma_par[1][4][2]=0.780916;
    ter_sigma_par[1][5][0]=0.11815;
    ter_sigma_par[1][5][1]=0.0542567;
    ter_sigma_par[1][5][2]=1.38387;   
    ter_sigma_par[1][6][0]=0.00439531;
    ter_sigma_par[1][6][1]=0.0475915;
    ter_sigma_par[1][6][2]=1.95818;
    ter_sigma_par[1][7][0]=0.266689;
    ter_sigma_par[1][7][1]=0.0167383;
    ter_sigma_par[1][7][2]=0.912125;
    ter_sigma_par[1][8][0]=0.200842;
    ter_sigma_par[1][8][1]=0.0206279;
    ter_sigma_par[1][8][2]=0.946193;
    ter_sigma_par[1][9][0]=0.144839;
    ter_sigma_par[1][9][1]=0.0258365;
    ter_sigma_par[1][9][2]=1.06918;

    double ter_sigma=ter_sigma_par[type_ind][eta_bin][0]/sqrt(Pt_scan)+
        ter_sigma_par[type_ind][eta_bin][1]+
        ter_sigma_par[type_ind][eta_bin][2]/Pt_scan;

    return ter_sigma;
}

int EtaBin(double t_eta)
{
    int eta_bin=(int)(fabs(t_eta)/0.25);
    if(eta_bin>9) eta_bin=9;
    return eta_bin;
}

}

/*
//------- Systematics study for lep-had channel, rel16 analysis
//------- slightly modified MMC code
void runMMC_lephad_tag15_8TeV_TERtest_ntup(int useLHcut=0, int events=-1){

    //   TF1 *ter_CB=new TF1("CB",MyCrystallBall,0.0,3.0,5);
    //   TF1 *ter_G=new TF1("G",MyCrystallBall,0.0,3.0,5);

    //define variables
    TString sample_name,channelType;
    MissingMassCalculator fMMC;
    //   fMMC.SetNsigmaMETscan(4.0);
    //   fMMC.SetUseTailCleanup(0);
    fMMC.SetCalibrationSet(MMCCalibrationSet::MMC2012);

    double MMC_time,cpu_time,real_time;

    TH1F* h_mmc_meth0=new TH1F("h_1","MMC method-0",100,0.0,500.0); 
    TH1F* h_mmc_meth1=new TH1F("h_2","MMC method-1",100,0.0,500.0); 
    TH1F* h_mmc_meth0_cb=new TH1F("h_3","MMC method-0, CB smear",100,0.0,500.0); 
    TH1F* h_mmc_meth1_cb=new TH1F("h_4","MMC method-1, CB smear",100,0.0,500.0); 
    TH1F* h_mmc_meth0_g=new TH1F("h_5","MMC method-0, G smear",100,0.0,500.0); 
    TH1F* h_mmc_meth1_g=new TH1F("h_6","MMC method-1, G smear",100,0.0,500.0); 

    TH1F* h_mmc_ratio_meth0=new TH1F("h_7","MMC method-0: (G-smear)/(CB smear)",100,0.0,500.0);
    TH1F* h_mmc_ratio_meth1=new TH1F("h_8","MMC method-1: (G-smear)/(CB smear)",100,0.0,500.0);

    TH1F* h_tau1res_cb=new TH1F("h_10","Tau-1 resolution: Crystal Ball",200,0.0,2.0);
    TH1F* h_tau1res_g=new TH1F("h_10","Tau-1 resolution: Gaussian",200,0.0,2.0);

    h_mmc_meth0->Sumw2(); 
    h_mmc_meth1->Sumw2(); 
    h_mmc_meth0_cb->Sumw2(); 
    h_mmc_meth1_cb->Sumw2(); 
    h_mmc_meth0_g->Sumw2(); 
    h_mmc_meth1_g->Sumw2(); 
    h_mmc_ratio_meth0->Sumw2();
    h_mmc_ratio_meth1->Sumw2();
    h_tau1res_cb->Sumw2();
    h_tau1res_g->Sumw2();


    //---------------------------------------------------------------------------
    // ------------- output variables
    int cat_def=-10; // event category
    int mmc_stat=0;  // MMC status
    int col_stat=0;  // collinear approximation status
    double cpu=0.0;  // CPU time
    double mmc_mass0=0.0; // method-0
    double mmc_mass1=0.0; // method-1 (histo)
    double mmc_mass2=0.0; // method-2 (4-vec)
    double mmc_mass0_cb=0.0; // method-0, Crystal Ball resolution
    double mmc_mass1_cb=0.0; // method-1, Crystal Ball resolution
    double mmc_mass0_g=0.0; // method-0, Gaussian resolution
    double mmc_mass1_g=0.0; // method-1, Gaussian resolution
    double res_cb=0.0; // tau, Crystal Ball resolution (applied to truth vis-tau)
    double res_g=0.0; // tau, Gaussian resolution (applied to truth vis-tau)
    int cutpass_def=-1;
    int cutpass_cb=-1;
    int cutpass_g=-1;

    double mvis=0.0;
    double meff=0.0;
    double mcoll=0.0; // official coll app mass
    //------------
    double coll_x1=-10.0; // x1 in collinear approximation
    double coll_x2=-10.0; // x2 in collinear approximation
    double evnt_weight=1.0; // event weight

    int mmc_type0=-1;
    int mmc_type1=-1;
    int njet30=-1;
    double sumet=0.0;
    double det_met_x;
    double det_met_y;
    double true_met_x;
    double true_met_y;
    double sumjet_pt;
    double sumjet_eta;
    double sumjet_phi;
    double sumjet_m;
    double true_mass;
    // ---------- true original tau's
    double true_tau0_pt;
    double true_tau0_eta;
    double true_tau0_phi;
    double true_tau0_m;
    double true_tau1_pt;
    double true_tau1_eta;
    double true_tau1_phi;
    double true_tau1_m;
    // --------- true original neutrino's
    double true_nu0_pt;
    double true_nu0_eta;
    double true_nu0_phi;
    double true_nu0_m;
    double true_nu1_pt;
    double true_nu1_eta;
    double true_nu1_phi;
    double true_nu1_m;
    //----- visible tau's at detector level
    double det_tau0_pt;
    double det_tau0_eta;
    double det_tau0_phi;
    double det_tau0_m;
    double det_tau1_pt;
    double det_tau1_eta;
    double det_tau1_phi;
    double det_tau1_m;
    // fully reconstructed tau's from MMC, method-0 (highest probability solution)
    double meth0_tau0_pt;
    double meth0_tau0_eta;
    double meth0_tau0_phi;
    double meth0_tau0_m;
    double meth0_tau1_pt;
    double meth0_tau1_eta;
    double meth0_tau1_phi;
    double meth0_tau1_m;
    // fully reconstructed tau's from MMC, method-0 (4-vec method)
    double meth2_tau0_pt;
    double meth2_tau0_eta;
    double meth2_tau0_phi;
    double meth2_tau0_m;
    double meth2_tau1_pt;
    double meth2_tau1_eta;
    double meth2_tau1_phi;
    double meth2_tau1_m;
    // fully reconstructed nu's from MMC, method-0 (highest probability solution)
    double meth0_nu0_pt;
    double meth0_nu0_eta;
    double meth0_nu0_phi;
    double meth0_nu0_m;
    double meth0_nu1_pt;
    double meth0_nu1_eta;
    double meth0_nu1_phi;
    double meth0_nu1_m;
    // fully reconstructed nu's from MMC, method-0 (4-vec method)
    double meth2_nu0_pt;
    double meth2_nu0_eta;
    double meth2_nu0_phi;
    double meth2_nu0_m;
    double meth2_nu1_pt;
    double meth2_nu1_eta;
    double meth2_nu1_phi;
    double meth2_nu1_m;
    //__________________________________ end of output variables

    TString outName; 

    if(useLHcut==0) outName= "ana_embedding_lephad_tag15_MVAcut0_preselection_TERstudy_091413.root";
    if(useLHcut==1) outName= "ana_embedding_lephad_tag15_MVAcut1_rest_TERstudy_091413.root";
    if(useLHcut==2) outName= "ana_embedding_lephad_tag15_MVAcut2_1jet_TERstudy_091413.root";
    if(useLHcut==3) outName= "ana_embedding_lephad_tag15_MVAcut3_boosted_TERstudy_091413.root";
    if(useLHcut==4) outName= "ana_embedding_lephad_tag15_MVAcut4_vbf_TERstudy_091413.root";  

    TFile *fout  = new TFile(outName,"RECREATE");
    TTree *mytree = new TTree("Default","Tree");

    mytree->Branch("category",&cat_def,"cat_def/I");
    mytree->Branch("mmc_stat", &mmc_stat, "mmc_stat/I");
    mytree->Branch("col_stat", &col_stat, "col_stat/I");
    mytree->Branch("time", &cpu, "cpu/D");
    mytree->Branch("mmc_mass0", &mmc_mass0, "mmc_mass0/D");
    mytree->Branch("mmc_mass1", &mmc_mass1, "mmc_mass1/D");
    mytree->Branch("mmc_mass2", &mmc_mass2, "mmc_mass2/D");
    mytree->Branch("mvis", &mvis, "mvis/D");
    mytree->Branch("meff", &meff, "meff/D");
    mytree->Branch("mcoll", &mcoll, "mcoll/D"); 

    mytree->Branch("coll_x1", &coll_x1, "coll_x1/D");
    mytree->Branch("coll_x2", &coll_x2, "coll_x2/D");
    mytree->Branch("evnt_weight", &evnt_weight, "evnt_weight/D"); 

    mytree->Branch("mmc_type0", &mmc_type0, "mmc_type0/I");
    mytree->Branch("mmc_type1", &mmc_type1, "mmc_type1/I");
    mytree->Branch("njet30", &njet30, "njet30/I");
    mytree->Branch("sumet", &sumet, "sumet/D");
    mytree->Branch("det_met_x", &det_met_x, "det_met_x/D");
    mytree->Branch("det_met_y", &det_met_y, "det_met_y/D");
    mytree->Branch("true_met_x", &true_met_x, "true_met_x/D");
    mytree->Branch("true_met_y", &true_met_y, "true_met_y/D");
    mytree->Branch("sumjet_pt", &sumjet_pt, "sumjet_pt/D");
    mytree->Branch("sumjet_eta", &sumjet_eta, "sumjet_eta/D");
    mytree->Branch("sumjet_phi", &sumjet_phi, "sumjet_phi/D");
    mytree->Branch("sumjet_m", &sumjet_m, "sumjet_m/D");
    mytree->Branch("true_mass", &true_mass, "true_mass/D");

    mytree->Branch("mmc_mass0_cb",&mmc_mass0_cb,"mmc_mass0_cb/D");
    mytree->Branch("mmc_mass1_cb",&mmc_mass1_cb,"mmc_mass1_cb/D");
    mytree->Branch("mmc_mass0_g",&mmc_mass0_g,"mmc_mass0_g/D");
    mytree->Branch("mmc_mass1_g",&mmc_mass1_g,"mmc_mass1_g/D");
    mytree->Branch("res_cb",&res_cb,"res_cb/D");
    mytree->Branch("res_g",&res_g,"res_g/D");
    mytree->Branch("cutpass_def",&cutpass_def,"cutpass_def/I");
    mytree->Branch("cutpass_cb",&cutpass_cb,"cutpass_cb/I");
    mytree->Branch("cutpass_g",&cutpass_g,"cutpass_g/I");
    // ---------- true original tau's
    mytree->Branch("true_tau0_pt", &true_tau0_pt, "true_tau0_pt/D");
    mytree->Branch("true_tau0_eta", &true_tau0_eta, "true_tau0_eta/D");
    mytree->Branch("true_tau0_phi", &true_tau0_phi, "true_tau0_phi/D");
    mytree->Branch("true_tau0_m", &true_tau0_m, "true_tau0_m/D");
    mytree->Branch("true_tau1_pt", &true_tau1_pt, "true_tau1_pt/D");
    mytree->Branch("true_tau1_eta", &true_tau1_eta, "true_tau1_eta/D");
    mytree->Branch("true_tau1_phi", &true_tau1_phi, "true_tau1_phi/D");
    mytree->Branch("true_tau1_m", &true_tau1_m, "true_tau1_m/D");
    // --------- true original neutrino's
    mytree->Branch("true_nu0_pt", &true_nu0_pt, "true_nu0_pt/D");
    mytree->Branch("true_nu0_eta", &true_nu0_eta, "true_nu0_eta/D");
    mytree->Branch("true_nu0_phi", &true_nu0_phi, "true_nu0_phi/D");
    mytree->Branch("true_nu0_m", &true_nu0_m, "true_nu0_m/D");
    mytree->Branch("true_nu1_pt", &true_nu1_pt, "true_nu1_pt/D");
    mytree->Branch("true_nu1_eta", &true_nu1_eta, "true_nu1_eta/D");
    mytree->Branch("true_nu1_phi", &true_nu1_phi, "true_nu1_phi/D");
    mytree->Branch("true_nu1_m", &true_nu1_m, "true_nu1_m/D");
    //----- visible tau's at detector level
    mytree->Branch("det_tau0_pt", &det_tau0_pt, "det_tau0_pt/D");
    mytree->Branch("det_tau0_eta", &det_tau0_eta, "det_tau0_eta/D");
    mytree->Branch("det_tau0_phi", &det_tau0_phi, "det_tau0_phi/D");
    mytree->Branch("det_tau0_m", &det_tau0_m, "det_tau0_m/D");
    mytree->Branch("det_tau1_pt", &det_tau1_pt, "det_tau1_pt/D");
    mytree->Branch("det_tau1_eta", &det_tau1_eta, "det_tau1_eta/D");
    mytree->Branch("det_tau1_phi", &det_tau1_phi, "det_tau1_phi/D");
    mytree->Branch("det_tau1_m", &det_tau1_m, "det_tau1_m/D");
    // fully reconstructed tau's from MMC, method-0 (highest probability solution)
    mytree->Branch("meth0_tau0_pt", &meth0_tau0_pt, "meth0_tau0_pt/D");
    mytree->Branch("meth0_tau0_eta", &meth0_tau0_eta, "meth0_tau0_eta/D");
    mytree->Branch("meth0_tau0_phi", &meth0_tau0_phi, "meth0_tau0_phi/D");
    mytree->Branch("meth0_tau0_m", &meth0_tau0_m, "meth0_tau0_m/D");
    mytree->Branch("meth0_tau1_pt", &meth0_tau1_pt, "meth0_tau1_pt/D");
    mytree->Branch("meth0_tau1_eta", &meth0_tau1_eta, "meth0_tau1_eta/D");
    mytree->Branch("meth0_tau1_phi", &meth0_tau1_phi, "meth0_tau1_phi/D");
    mytree->Branch("meth0_tau1_m", &meth0_tau1_m, "meth0_tau1_m/D");
    // fully reconstructed tau's from MMC, method-0 (4-vec method)
    mytree->Branch("meth2_tau0_pt", &meth2_tau0_pt, "meth2_tau0_pt/D");
    mytree->Branch("meth2_tau0_eta", &meth2_tau0_eta, "meth2_tau0_eta/D");
    mytree->Branch("meth2_tau0_phi", &meth2_tau0_phi, "meth2_tau0_phi/D");
    mytree->Branch("meth2_tau0_m", &meth2_tau0_m, "meth2_tau0_m/D");
    mytree->Branch("meth2_tau1_pt", &meth2_tau1_pt, "meth2_tau1_pt/D");
    mytree->Branch("meth2_tau1_eta", &meth2_tau1_eta, "meth2_tau1_eta/D");
    mytree->Branch("meth2_tau1_phi", &meth2_tau1_phi, "meth2_tau1_phi/D");
    mytree->Branch("meth2_tau1_m", &meth2_tau1_m, "meth2_tau1_m/D");
    // fully reconstructed nu's from MMC, method-0 (highest probability solution)
    mytree->Branch("meth0_nu0_pt", &meth0_nu0_pt, "meth0_nu0_pt/D");
    mytree->Branch("meth0_nu0_eta", &meth0_nu0_eta, "meth0_nu0_eta/D");
    mytree->Branch("meth0_nu0_phi", &meth0_nu0_phi, "meth0_nu0_phi/D");
    mytree->Branch("meth0_nu0_m", &meth0_nu0_m, "meth0_nu0_m/D");
    mytree->Branch("meth0_nu1_pt", &meth0_nu1_pt, "meth0_nu1_pt/D");
    mytree->Branch("meth0_nu1_eta", &meth0_nu1_eta, "meth0_nu1_eta/D");
    mytree->Branch("meth0_nu1_phi", &meth0_nu1_phi, "meth0_nu1_phi/D");
    mytree->Branch("meth0_nu1_m", &meth0_nu1_m, "meth0_nu1_m/D");
    // fully reconstructed nu's from MMC, method-0 (4-vec method)
    mytree->Branch("meth2_nu0_pt", &meth2_nu0_pt, "meth2_nu0_pt/D");
    mytree->Branch("meth2_nu0_eta", &meth2_nu0_eta, "meth2_nu0_eta/D");
    mytree->Branch("meth2_nu0_phi", &meth2_nu0_phi, "meth2_nu0_phi/D");
    mytree->Branch("meth2_nu0_m", &meth2_nu0_m, "meth2_nu0_m/D");
    mytree->Branch("meth2_nu1_pt", &meth2_nu1_pt, "meth2_nu1_pt/D");
    mytree->Branch("meth2_nu1_eta", &meth2_nu1_eta, "meth2_nu1_eta/D");
    mytree->Branch("meth2_nu1_phi", &meth2_nu1_phi, "meth2_nu1_phi/D");
    mytree->Branch("meth2_nu1_m", &meth2_nu1_m, "meth2_nu1_m/D");
    //__________________________________ end of output variables


    TString fFilename;
    TString TreeName="tau";
    TString inputPath = "/Users/pronko/work/ATLAS/atlas_analysis/Htautau/results/TER_studies/lephad_091413/";

    int Nsamples=9;

    int N_mmcOK=0;
    int N_mmcBad=0;

    TFile *fin;

    for(int j=0; j<Nsamples; j++)
    {

        if(j==0) fFilename ="group.phys-higgs.364585_030337._00026.merge.ntuple.root";
        if(j==1) fFilename ="group.phys-higgs.364585_030337._00037.merge.ntuple.root";
        if(j==2) fFilename ="group.phys-higgs.364585_030337._00041.merge.ntuple.root";
        if(j==3) fFilename ="group.phys-higgs.364585_030337._00066.merge.ntuple.root";
        if(j==4) fFilename ="group.phys-higgs.364585_030337._00266.merge.ntuple.root";
        if(j==5) fFilename ="group.phys-higgs.364585_030337._00466.merge.ntuple.root";
        if(j==6) fFilename ="group.phys-higgs.364585_030337._00493.merge.ntuple.root";
        if(j==7) fFilename ="group.phys-higgs.364585_030337._00502.merge.ntuple.root";
        if(j==8) fFilename ="group.phys-higgs.364585_030337._00507.merge.ntuple.root";

        std::cout<<"Open file "<<inputPath+fFilename<<std::endl;
        fin = new TFile(inputPath+fFilename,"READ");      

        TTree* tree = (TTree*) fin->Get(TreeName);
        int TotalEntries = (int) tree->GetEntries();

        std::cout<<"Total number of events to be processed: "<<TotalEntries<<std::endl;

        //--------- global
        int sltbit;
        int osbit;
        double eventWeight;
        double met_x;
        double met_y;
        double sum_et;
        int lepton_type_det;
        int tau_type_true;
        //--- 1st tau is always a lepton
        double tau1_pt;
        double tau1_eta;
        double tau1_phi;
        //--- 2nd tau is always hadronic tau
        std::vector<double> *tau2_pt=0;
        std::vector<double> *tau2_eta=0;
        std::vector<double> *tau2_phi=0;
        std::vector<int> *tau2_ntrack=0;
        std::vector<int> *tau2_eveto=0;
        std::vector<int> *tau2_muveto=0;
        //------ jets
        std::vector<double> *jet_Pt=0;
        std::vector<double> *jet_Phi=0;
        std::vector<double> *jet_Eta=0;
        std::vector<double> *jet_M=0;
        std::vector<double> *jet_jvf=0;
        //------ true tau's
        double true_tau1_Pt; // tau-1 is always a lepton
        double true_tau1_Phi;
        double true_tau1_Eta;
        double true_tau1_M;
        double true_tau2_Pt; // tau-2 is always a hadronic tau 
        double true_tau2_Phi;
        double true_tau2_Eta;
        double true_tau2_M;
        //------ true nu's
        double true_nu1_Pt;
        double true_nu1_Phi;
        double true_nu1_Eta;
        double true_nu1_M;
        double true_nu2_Pt;
        double true_nu2_Phi;
        double true_nu2_Eta;
        double true_nu2_M;

        //--------- getting detector tau's
        tree->SetBranchAddress("isSLT", &sltbit);
        tree->SetBranchAddress("isOS", &osbit);
        tree->SetBranchAddress("true_tau_lep_decay_mode", &lepton_type_det);
        tree->SetBranchAddress("true_tau_had_decay_mode", &tau_type_true);

        tree->SetBranchAddress("lep_pt", &tau1_pt);
        tree->SetBranchAddress("lep_eta", &tau1_eta);
        tree->SetBranchAddress("lep_phi", &tau1_phi);

        tree->SetBranchAddress("tau_pt", &tau2_pt);
        tree->SetBranchAddress("tau_eta", &tau2_eta);
        tree->SetBranchAddress("tau_phi", &tau2_phi);
        tree->SetBranchAddress("tau_numTrack", &tau2_ntrack);
        tree->SetBranchAddress("tau_EleBADTMedium", &tau2_eveto);
        tree->SetBranchAddress("tau_MuonVeto", &tau2_muveto);

        tree->SetBranchAddress("met_Etx", &met_x);
        tree->SetBranchAddress("met_Ety", &met_y);
        tree->SetBranchAddress("met_sumEt", &sum_et);
        tree->SetBranchAddress("evt_tot_weight", &eventWeight);

        tree->SetBranchAddress("jet_pt", &jet_Pt);
        tree->SetBranchAddress("jet_eta", &jet_Eta);
        tree->SetBranchAddress("jet_phi", &jet_Phi);
        tree->SetBranchAddress("jet_m", &jet_M);
        tree->SetBranchAddress("jet_JVF", &jet_jvf);

        tree->SetBranchAddress("true_vis_lep_pt", &true_tau1_Pt);
        tree->SetBranchAddress("true_vis_lep_eta", &true_tau1_Eta);
        tree->SetBranchAddress("true_vis_lep_phi", &true_tau1_Phi);
        tree->SetBranchAddress("true_vis_lep_m", &true_tau1_M);
        tree->SetBranchAddress("true_vis_tau_pt", &true_tau2_Pt);
        tree->SetBranchAddress("true_vis_tau_eta", &true_tau2_Eta);
        tree->SetBranchAddress("true_vis_tau_phi", &true_tau2_Phi);
        tree->SetBranchAddress("true_vis_tau_m", &true_tau2_M);

        tree->SetBranchAddress("true_nu_lep_pt", &true_nu1_Pt);
        tree->SetBranchAddress("true_nu_lep_eta", &true_nu1_Eta);
        tree->SetBranchAddress("true_nu_lep_phi", &true_nu1_Phi);
        tree->SetBranchAddress("true_nu_lep_m", &true_nu1_M);
        tree->SetBranchAddress("true_nu_had_pt", &true_nu2_Pt);
        tree->SetBranchAddress("true_nu_had_eta", &true_nu2_Eta);
        tree->SetBranchAddress("true_nu_had_phi", &true_nu2_Phi);
        tree->SetBranchAddress("true_nu_had_m", &true_nu2_M);

        int tau0Type,tau1Type;      

        std::vector<TLorentzVector> det_tau;
        std::vector<TLorentzVector> true_nu;
        std::vector<TLorentzVector> true_tau;
        std::vector<int> det_tau_Nprong;

        int Nloop= events>0 ? events : TotalEntries;

        double m_1prong = 0.8;
        double m_3prong = 1.2;

        for(int iEntry=0; iEntry<Nloop; iEntry++)
        {	

            true_nu.clear();
            true_tau.clear();
            det_tau.clear();
            det_tau_Nprong.clear();

            if ( iEntry % 10  == 0 ) std::cout<<"Processed events : "<<iEntry<<std::endl;  
            //get variables needed for MMC  

            tree->GetEntry(iEntry);

            if(osbit==0) continue;

            // ---------- filling true tau's
            TLorentzVector tmp(0.0,0.0,0.0,0.0);
            tmp.SetPtEtaPhiM(true_tau1_Pt,true_tau1_Eta,true_tau1_Phi,true_tau1_M); // 1st tau is always a lepton
            true_tau.push_back(tmp);
            tmp.SetPtEtaPhiM(true_tau2_Pt,true_tau2_Eta,true_tau2_Phi,true_tau2_M);
            true_tau.push_back(tmp);
            tmp.SetPtEtaPhiM(true_nu1_Pt,true_nu1_Eta,true_nu1_Phi,true_nu1_M); // 1st neutrino is always a combined lepton neutrino
            true_nu.push_back(tmp);
            tmp.SetPtEtaPhiM(true_nu2_Pt,true_nu2_Eta,true_nu2_Phi,0.0);
            true_nu.push_back(tmp);
            // ------------------- filling detector leptons & tau's
            TLorentzVector tau0,tau1;
            tau0.SetPtEtaPhiM(tau1_pt,tau1_eta,tau1_phi,true_tau1_M); // 1st tau is always a lepton 

            int nTau=tau2_pt->size();

            if(nTau<1) continue;

            for(int i=0; i<nTau; i++)
            {
                bool goodtau=true;
                double tau_m = tau2_ntrack->at(i)==1 ? m_1prong : m_3prong;
                tmp.SetPtEtaPhiM(tau2_pt->at(i),tau2_eta->at(i),tau2_phi->at(i),tau_m);
                if(tau2_muveto->at(i)==0) goodtau=false;
                if(tmp.Pt()>15.0 && tau2_eveto->at(i)==0) goodtau=false;
                if(tmp.DeltaR(true_tau[0])<0.2) goodtau=false;
                if(tmp.DeltaR(true_tau[1])>0.2) goodtau=false;
                if(tmp.DeltaR(tau0)<0.2) goodtau=false;
                if(goodtau==true)
                {
                    det_tau.push_back(tmp);
                    det_tau_Nprong.push_back(tau2_ntrack->at(i));
                }
            }

            if(det_tau.size()<1) continue;

            tau1=det_tau[0]; // if there is more than one detector tau, take the first one
            evnt_weight=eventWeight;
            tau0Type = lepton_type_det; 
            tau1Type = det_tau_Nprong[0]*10; 
            sumet = sum_et;

            //prepare quantities to run MMC
            TVector2 met_vect(met_x,met_y);
            TVector2 met_true((true_nu[0]+true_nu[1]).X(),(true_nu[0]+true_nu[1]).Y());
            TLorentzVector met4vec(met_vect.Px(),met_vect.Py(),0.0,met_vect.Mod());

            int nJet = jet_Pt->size();

            int Njet30=0;

            //make jet TLorentzVector    
            std::vector<TLorentzVector> jetvec;
            jetvec.clear();
            TLorentzVector jetP4(0. , 0. , 0. , 0.);
            TLorentzVector total_jet30(0. , 0. , 0. , 0.);

            //---- getting jets
            for(int i=0; i<nJet; i++)
            {
                jetP4.SetPtEtaPhiM(jet_Pt->at(i),jet_Eta->at(i),jet_Phi->at(i),jet_M->at(i));
                if(jetP4.DeltaR(tau0)<0.4 || jetP4.DeltaR(tau1)<0.4) continue;
                bool good_jet=true;
                if(fabs(jetP4.Eta())<2.4 && fabs(jet_jvf->at(i))<0.5 && jetP4.Pt()<50.0) good_jet=false;
                if(good_jet==true) jetvec.push_back(jetP4);
                if(jetP4.Pt()>30.0 && good_jet==true) 
                {
                    Njet30++;
                    total_jet30=total_jet30+jetP4;
                }
            }

            //====================================================================================
            //---------------    Very quick cuts -------------------------------------------------
            if(useLHcut==2 && Njet30==0) continue; // 1-jet
            if(useLHcut==3 && sltbit==0) continue; // Boosted
            if(useLHcut==4) // VBF
            {
                if(sltbit==0) continue;
                if(Njet30<2) continue; // VBF
                if(jetvec[0].Pt()<50.0) continue;
                if(jetvec[1].Pt()<30.0) continue;
                if(fabs(jetvec[0].Eta()-jetvec[1].Eta())<3.0) continue;

            }	  

            //====================================================================================
            //------- TER smearing
            //
            float ter_sigma=TERSigma(det_tau_Nprong[0],EtaBin(tau1.Eta()),tau1.Pt());
            float ter_alfa=TERafla(det_tau_Nprong[0],EtaBin(tau1.Eta()),tau1.Pt());
            //------ Crystal Ball TER
            TF1 *ter_CBf=new TF1("CB",MyCrystallBall,0.0,3.0,4);
            ter_CBf->SetParameter(0,1.0); // mean 
            ter_CBf->SetParameter(1,ter_sigma); // sigma
            ter_CBf->SetParameter(2,ter_alfa); // alfa
            if(det_tau_Nprong[0]==1) ter_CBf->SetParameter(3,100.0); // CB power for power series
            if(det_tau_Nprong[0]==3) ter_CBf->SetParameter(3,147.0); // CB power for power series
            double terCB_x=ter_CBf->GetRandom(0.0,3.0);
            //------ Gaussian TER
            TF1 *ter_Gf=new TF1("G",MyCrystallBall,0.0,3.0,4);
            ter_Gf->SetParameter(0,1.0); // mean 
            ter_Gf->SetParameter(1,ter_sigma); // sigma
            ter_Gf->SetParameter(2,5.0); // alfa
            ter_Gf->SetParameter(3,0.0); // CB power for power series
            double terG_x=ter_Gf->GetRandom(0.0,3.0);

            std::cout<<" TER_CB_x="<<terCB_x<<std::endl;
            std::cout<<" TER_G_x="<<terG_x<<std::endl;

            TLorentzVector tau1_cb(0.0,0.0,0.0,0.0);
            tau1_cb.SetPtEtaPhiM(true_tau[1].Pt()/terCB_x,tau1.Eta(),tau1.Phi(),tau1.M());
            TLorentzVector tau1_g(0.0,0.0,0.0,0.0);
            tau1_g.SetPtEtaPhiM(true_tau[1].Pt()/terG_x,tau1.Eta(),tau1.Phi(),tau1.M());

            // 	  double dMet_x=tau1.Px()-true_tau[1].Px();
            // 	  double dMet_y=tau1.Py()-true_tau[1].Py();
            // 	  double dMetCB_x=tau1_cb.Px()-true_tau[1].Px();
            // 	  double dMetCB_y=tau1_cb.Py()-true_tau[1].Py();
            // 	  double dMetG_x=tau1_g.Px()-true_tau[1].Px();
            // 	  double dMetG_y=tau1_g.Py()-true_tau[1].Py();

            TVector2 met_vect_cb(met_vect.X()-tau1.Px()+tau1_cb.Px(),met_vect.Y()-tau1.Py()+tau1_cb.Py());
            TVector2 met_vect_g(met_vect.X()-tau1.Px()+tau1_g.Px(),met_vect.Y()-tau1.Py()+tau1_g.Py());
            // 	  TVector2 met_vect_cb(met_vect.X()-dMet_x+dMetCB_x,met_vect.Y()-dMet_y+dMetCB_y);
            // 	  TVector2 met_vect_g(met_vect.X()-dMet_x+dMetG_x,met_vect.Y()-dMet_y+dMetG_y);


            TLorentzVector met4vec_cb(met_vect_cb.Px(),met_vect_cb.Py(),0.0,met_vect_cb.Mod());
            TLorentzVector met4vec_g(met_vect_g.Px(),met_vect_g.Py(),0.0,met_vect_g.Mod());

            if(useLHcut==3) // quick filter for Boosted events
            {
                if((tau0+tau1+met4vec).Pt()<100.0 && (tau0+tau1_cb+met4vec_cb).Pt()<100.0 && (tau0+tau1_g+met4vec_g).Pt()<100.0) continue;
            }

            cat_def=useLHcut;

            //prepare quantities to run MMC
            TLorentzVector MMC_rec4vec, nu0, nu1;

            double xp1, xp2;
            if(CollinearApproximationMatrix(tau0,tau1,met_vect.X(),met_vect.Y(), 
                        mcoll,xp1,xp2)==true)
            {
                coll_x1=xp1;
                coll_x2=xp2; 
            }

            //----------------- initialize output vairables
            mmc_stat=0;
            cpu=0.0;  // CPU time
            mmc_mass0=0.0; // method-0
            mmc_mass1=0.0; // method-1 (histo)
            mmc_mass2=0.0; // method-2 (4-vec)

            mmc_mass0_cb=0.0; // method-0, Crystal Ball resolution
            mmc_mass1_cb=0.0; // method-1, Crystal Ball resolution
            mmc_mass0_g=0.0; // method-0, Gaussian resolution
            mmc_mass1_g=0.0; // method-1, Gaussian resolution
            res_cb=1.0/terCB_x; 
            res_g=1.0/terG_x; 
            cutpass_def=0;
            cutpass_cb=0;
            cutpass_g=0;

            mvis=(tau0+tau1).M();
            meff=(tau0+tau1+met4vec).M();
            mmc_type0=tau0Type;
            mmc_type1=tau1Type;
            njet30=Njet30;
            det_met_x=met_vect.X();
            det_met_y=met_vect.Y();
            true_met_x=met_true.X();
            true_met_y=met_true.Y();
            if(njet30>0)
            {
                sumjet_pt=total_jet30.Pt();
                sumjet_eta=total_jet30.Eta();
                sumjet_phi=total_jet30.Phi();
                sumjet_m=total_jet30.M();
            }
            else 
            {
                sumjet_pt=0.0;
                sumjet_eta=-10.0;
                sumjet_phi=-10.0;
                sumjet_m=0.0;
            }
            // ---------- true original tau's
            true_tau0_pt=true_tau[0].Pt();
            true_tau0_eta=true_tau[0].Eta();
            true_tau0_phi=true_tau[0].Phi();
            true_tau0_m=true_tau[0].M();
            true_tau1_pt=true_tau[1].Pt();
            true_tau1_eta=true_tau[1].Eta();
            true_tau1_phi=true_tau[1].Phi();
            true_tau1_m=true_tau[1].M();
            // --------- true original neutrino's
            true_nu0_pt=true_nu[0].Pt();
            true_nu0_eta=true_nu[0].Eta();
            true_nu0_phi=true_nu[0].Phi();
            true_nu0_m=true_nu[0].M();
            true_nu1_pt=true_nu[1].Pt();
            true_nu1_eta=true_nu[1].Eta();
            true_nu1_phi=true_nu[1].Phi();
            true_nu1_m=true_nu[1].M();

            true_mass=(true_tau[0]+true_tau[1]+true_nu[0]+true_nu[1]).M();
            //----- visible tau's at detector level
            det_tau0_pt=tau0.Pt(); 
            det_tau0_eta=tau0.Eta();
            det_tau0_phi=tau0.Phi();
            det_tau0_m=tau0.M();	 
            det_tau1_pt=tau1.Pt(); 
            det_tau1_eta=tau1.Eta();
            det_tau1_phi=tau1.Phi();
            det_tau1_m=tau1.M();  
            // fully reconstructed tau's from MMC, method-0 (highest probability solution)
            meth0_tau0_pt=0.0;
            meth0_tau0_eta=0.0;
            meth0_tau0_phi=0.0;
            meth0_tau0_m=0.0;
            meth0_tau1_pt=0.0;
            meth0_tau1_eta=0.0;
            meth0_tau1_phi=0.0;
            meth0_tau1_m=0.0;
            // fully reconstructed tau's from MMC, method-0 (4-vec method)
            meth2_tau0_pt=0.0;
            meth2_tau0_eta=0.0;
            meth2_tau0_phi=0.0;
            meth2_tau0_m=0.0;
            meth2_tau1_pt=0.0;
            meth2_tau1_eta=0.0;
            meth2_tau1_phi=0.0;
            meth2_tau1_m=0.0;
            // fully reconstructed nu's from MMC, method-0 (highest probability solution)
            meth0_nu0_pt=0.0;
            meth0_nu0_eta=0.0;
            meth0_nu0_phi=0.0;
            meth0_nu0_m=0.0;
            meth0_nu1_pt=0.0;
            meth0_nu1_eta=0.0;
            meth0_nu1_phi=0.0;
            meth0_nu1_m=0.0;
            // fully reconstructed nu's from MMC, method-0 (4-vec method)
            meth2_nu0_pt=0.0;
            meth2_nu0_eta=0.0;
            meth2_nu0_phi=0.0;
            meth2_nu0_m=0.0;
            meth2_nu1_pt=0.0;
            meth2_nu1_eta=0.0;
            meth2_nu1_phi=0.0;
            meth2_nu1_m=0.0;

            TDatime startTime; 
            TStopwatch timer;
            timer.Start();

            //---------------------- running MMC
            if(cat_def>=0)
            {

                cutpass_def=MyLepHadCuts(useLHcut,tau0Type,sltbit,tau0,tau1,met_vect,Njet30,jetvec);
                cutpass_cb=MyLepHadCuts(useLHcut,tau0Type,sltbit,tau0,tau1_cb,met_vect_cb,Njet30,jetvec);
                cutpass_g=MyLepHadCuts(useLHcut,tau0Type,sltbit,tau0,tau1_g,met_vect_g,Njet30,jetvec);

                fMMC.SetMetVec(met_vect);
                fMMC.SetNjet25(Njet30);    
                fMMC.SetVisTauVec (0,  tau0);  
                fMMC.SetVisTauVec (1,  tau1);  
                fMMC.SetVisTauType(0,  tau0Type);  
                fMMC.SetVisTauType(1,  tau1Type);
                fMMC.SetSumEt(sumet);

                mmc_stat=fMMC.RunMissingMassCalculator(); // run MMC 
                mmc_stat=fMMC.GetFitStatus();
                TDatime endTime; 
                MMC_time = endTime.GetTime()- startTime.GetTime(); 
                timer.Stop();
                cpu_time = timer.CpuTime();
                real_time =timer.RealTime();
                cpu=real_time;
                //DEBUG messages
                std::cout <<"DEBUG MMC done. Time : "<< endTime.GetTime()-startTime.GetTime() <<" fit status="<< mmc_stat << std::endl;
                std::cout <<"DEBUG CPU time = " << timer.CpuTime() << "s real = " << timer.RealTime() << " s " << std::endl;
                std::cout <<"    "<<std::endl;

                if(fMMC.GetFitStatus()==1)
                {  
                    N_mmcOK++;
                    mmc_mass0=fMMC.GetFittedMass(0); 
                    mmc_mass1=fMMC.GetFittedMass(1); 
                    std::cout <<"Default MMC done. category="<<cat_def<<" Mass = "<<mmc_mass1<<std::endl;
                    std::cout <<"    "<<std::endl;      
                    MMC_rec4vec=fMMC.GetResonanceVec(1);
                    mmc_mass2=MMC_rec4vec.M();
                    TLorentzVector meth0_tau0=fMMC.GetTau4vec(0,0);
                    TLorentzVector meth0_tau1=fMMC.GetTau4vec(0,1);
                    TLorentzVector meth0_nu0=fMMC.GetNeutrino4vec(0,0);
                    TLorentzVector meth0_nu1=fMMC.GetNeutrino4vec(0,1);
                    TLorentzVector meth2_tau0=fMMC.GetTau4vec(1,0);
                    TLorentzVector meth2_tau1=fMMC.GetTau4vec(1,1);
                    TLorentzVector meth2_nu0=fMMC.GetNeutrino4vec(1,0);
                    TLorentzVector meth2_nu1=fMMC.GetNeutrino4vec(1,1);

                    h_mmc_meth0->Fill(mmc_mass0,eventWeight); 
                    h_mmc_meth1->Fill(mmc_mass1,eventWeight);

                    // fully reconstructed tau's from MMC, method-0 (highest probability solution)
                    meth0_tau0_pt=meth0_tau0.Pt();
                    meth0_tau0_eta=meth0_tau0.Eta();
                    meth0_tau0_phi=meth0_tau0.Phi();
                    meth0_tau0_m=meth0_tau0.M();
                    meth0_tau1_pt=meth0_tau1.Pt();
                    meth0_tau1_eta=meth0_tau1.Eta();
                    meth0_tau1_phi=meth0_tau1.Phi();
                    meth0_tau1_m=meth0_tau1.M();
                    // fully reconstructed tau's from MMC, method-0 (4-vec method)
                    meth2_tau0_pt=meth2_tau0.Pt();
                    meth2_tau0_eta=meth2_tau0.Eta();
                    meth2_tau0_phi=meth2_tau0.Phi();
                    meth2_tau0_m=meth2_tau0.M();
                    meth2_tau1_pt=meth2_tau1.Pt();
                    meth2_tau1_eta=meth2_tau1.Eta();
                    meth2_tau1_phi=meth2_tau1.Phi();
                    meth2_tau1_m=meth2_tau1.M();
                    // fully reconstructed nu's from MMC, method-0 (highest probability solution)
                    meth0_nu0_pt=meth0_nu0.Pt();
                    meth0_nu0_eta=meth0_nu0.Eta();
                    meth0_nu0_phi=meth0_nu0.Phi();
                    meth0_nu0_m=meth0_nu0.M();
                    meth0_nu1_pt=meth0_nu1.Pt();
                    meth0_nu1_eta=meth0_nu1.Eta();
                    meth0_nu1_phi=meth0_nu1.Phi();
                    meth0_nu1_m=meth0_nu1.M();
                    // fully reconstructed nu's from MMC, method-0 (4-vec method)
                    meth2_nu0_pt=meth2_nu0.Pt();
                    meth2_nu0_eta=meth2_nu0.Eta();
                    meth2_nu0_phi=meth2_nu0.Phi();
                    meth2_nu0_m=meth2_nu0.M();
                    meth2_nu1_pt=meth2_nu1.Pt();
                    meth2_nu1_eta=meth2_nu1.Eta();
                    meth2_nu1_phi=meth2_nu1.Phi();
                    meth2_nu1_m=meth2_nu1.M();
                } 
                else  
                { 
                    N_mmcBad++;
                    std::cout <<"    "<<std::endl;
                } 

                //----------------- running MMC for Crystal Ball smearing 
                fMMC.SetMetVec(met_vect_cb);
                fMMC.SetNjet25(Njet30);    
                fMMC.SetVisTauVec (0,  tau0);  
                fMMC.SetVisTauVec (1,  tau1_cb);  
                fMMC.SetVisTauType(0,  tau0Type);  
                fMMC.SetVisTauType(1,  tau1Type);
                fMMC.SetSumEt(sumet);
                int mmc_stat_cb=fMMC.RunMissingMassCalculator(); // run MMC	      
                if(fMMC.GetFitStatus()==1)
                {  
                    mmc_mass0_cb=fMMC.GetFittedMass(0); 
                    mmc_mass1_cb=fMMC.GetFittedMass(1); 
                    //DEBUG messages
                    std::cout <<"CB smearing: MMC-1="<<mmc_mass1_cb<<std::endl;
                    std::cout <<"    "<<std::endl;
                }

                //----------------- running MMC for Gaussian smearing 
                fMMC.SetMetVec(met_vect_g);
                fMMC.SetNjet25(Njet30);    
                fMMC.SetVisTauVec (0,  tau0);  
                fMMC.SetVisTauVec (1,  tau1_g);  
                fMMC.SetVisTauType(0,  tau0Type);  
                fMMC.SetVisTauType(1,  tau1Type);
                fMMC.SetSumEt(sumet);
                int mmc_stat_g=fMMC.RunMissingMassCalculator(); // run MMC	      	      
                if(fMMC.GetFitStatus()==1)
                {  
                    mmc_mass0_g=fMMC.GetFittedMass(0); 
                    mmc_mass1_g=fMMC.GetFittedMass(1); 
                    //DEBUG messages
                    std::cout <<"G smearing: MMC-1="<<mmc_mass1_g<<std::endl;
                    std::cout <<"    "<<std::endl;
                }

                h_mmc_meth0_cb->Fill(mmc_mass0_cb,eventWeight); 
                h_mmc_meth1_cb->Fill(mmc_mass1_cb,eventWeight); 
                h_mmc_meth0_g->Fill(mmc_mass0_g,eventWeight); 
                h_mmc_meth1_g->Fill(mmc_mass1_g,eventWeight); 
                h_tau1res_cb->Fill(terCB_x,eventWeight);
                h_tau1res_g->Fill(terG_x,eventWeight);

            }
            mytree->Fill();
        }//entries-loop
    }

    h_mmc_ratio_meth0->Divide(h_mmc_meth0_g,h_mmc_meth0_cb,1.0,1.0);
    h_mmc_ratio_meth1->Divide(h_mmc_meth1_g,h_mmc_meth1_cb,1.0,1.0);

    //   fout->cd();
    //   mytree->Write();
    //   tes_up->Write();
    //   tes_down->Write();
    fout->Write();
    fout->Close();
    delete fout;

    std::cout<<"===> Default MMC efficiency: "<<1.0*N_mmcOK/(1.0*(N_mmcOK+N_mmcBad))<<std::endl;

    TCanvas *c3 = new TCanvas("c3","MMC method-0: black=default, red=Crystal Ball res, blue=Gaussian res",200,10,700,500);
    c3->SetBorderMode(0);
    c3->SetLeftMargin(0.13);
    c3->SetBottomMargin(0.12);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(1);
    gPad->SetTicky();
    gPad->SetTickx();
    gPad->SetFillColor(0);
    gPad->SetFillStyle(0);
    gPad->SetBorderSize(0);
    h_mmc_meth0->DrawNormalized("hist");
    h_mmc_meth0_cb->SetLineColor(2);
    h_mmc_meth0_cb->DrawNormalized("hist same");
    h_mmc_meth0_g->SetLineColor(4);
    h_mmc_meth0_g->DrawNormalized("hist same");

    TCanvas *c4 = new TCanvas("c4","MMC method-1: black=default, red=Crystal Ball res, blue=Gaussian res",200,10,700,500);
    c4->SetBorderMode(0);
    c4->SetLeftMargin(0.13);
    c4->SetBottomMargin(0.12);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(1);
    gPad->SetTicky();
    gPad->SetTickx();
    gPad->SetFillColor(0);
    gPad->SetFillStyle(0);
    gPad->SetBorderSize(0);
    h_mmc_meth1->DrawNormalized("hist");
    h_mmc_meth1_cb->SetLineColor(2);
    h_mmc_meth1_cb->DrawNormalized("hist same");
    h_mmc_meth1_g->SetLineColor(4);
    h_mmc_meth1_g->DrawNormalized("hist same");


    TCanvas *c5 = new TCanvas("c5","MMC method-0: (Gaussian res)/(Crystal Ball res)",200,10,700,500);
    c5->SetBorderMode(0);
    c5->SetLeftMargin(0.13);
    c5->SetBottomMargin(0.12);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(1);
    gPad->SetTicky();
    gPad->SetTickx();
    gPad->SetFillColor(0);
    gPad->SetFillStyle(0);
    gPad->SetBorderSize(0);
    h_mmc_ratio_meth0->Draw();

    TCanvas *c6 = new TCanvas("c6","MMC method-1: (Gaussian res)/(Crystal Ball res)",200,10,700,500);
    c6->SetBorderMode(0);
    c6->SetLeftMargin(0.13);
    c6->SetBottomMargin(0.12);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(1);
    gPad->SetTicky();
    gPad->SetTickx();
    gPad->SetFillColor(0);
    gPad->SetFillStyle(0);
    gPad->SetBorderSize(0);
    h_mmc_ratio_meth1->Draw();

    TCanvas *c2 = new TCanvas("c2","Crystal Ball: tau-1 resolution)",200,10,700,500);
    c2->SetBorderMode(0);
    c2->SetLeftMargin(0.13);
    c2->SetBottomMargin(0.12);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(1);
    gPad->SetTicky();
    gPad->SetTickx();
    gPad->SetFillColor(0);
    gPad->SetFillStyle(0);
    gPad->SetBorderSize(0);
    h_tau1res_cb->DrawNormalized("hist");

    TCanvas *c_2 = new TCanvas("c_2","Gaussian: tau-1 resolution)",200,10,700,500);
    c_2->SetBorderMode(0);
    c_2->SetLeftMargin(0.13);
    c_2->SetBottomMargin(0.12);
    gStyle->SetOptStat(1111);
    gStyle->SetOptTitle(1);
    gPad->SetTicky();
    gPad->SetTickx();
    gPad->SetFillColor(0);
    gPad->SetFillStyle(0);
    gPad->SetBorderSize(0);
    h_tau1res_g->DrawNormalized("hist");

} 

*/
