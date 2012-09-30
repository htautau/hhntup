/*
Author: Aaron Armbruster
Date:   2012-05-25
Email:  armbrusa@umich.edu
Description: Script to run asymptotic CLs.



/////////////////////
//////PREAMBLE///////
/////////////////////

The script uses an iterative process to find the crossing of qmu with the qmu95(mu/sigma) curve, 
where qmu95(mu/sigma) is found assuming asymptotic formula for the distribution of the 
test statistic f(qmu|mu') (arxiv 1007.1727) and of the test statistic qmu (or tilde)

The sequence is

mu_i+1 = mu_i - gamma_i*(mu_i - mu'_i)

where gamma_i is a dynamic damping factor used for convergence (nominal gamma_i = 1), and mu'_i is 
determined by extrapolating the test statistic to the qmu95 curve assuming qmu is parabolic:

qmu'_i = (mu'_i - muhat)^2 / sigma_i^2 = qmu95(mu'_i / sigma_i)

where sigma_i is determined by computing qmu_i (not '):

sigma_i = (mu_i - muhat) / sqrt(qmu_i)

At the crossing qmu_N = qmu95 the assumption that qmu is a parabola goes away, 
so we're not ultimately dependent on this assumption beyond its use in the asymptotic formula.

The sequence ends when the relative correction factor gamma*(mu_i - mu'_i) / mu_i is less than some
specified precision (0.005 by default)




///////////////////////////
//////AFTER RUNNING////////
///////////////////////////


The results will be printed as well as stored in a root file in the folder 'root-files/<folder>', where <folder>
is specified by you (default 'test')

The root file has a 7-bin TH1D named 'limit', where each bin is filled with the upper limit values in this order:

1: Observed
2: Median
3: +2 sigma
4: +1 sigma
5: -1 sigma
6: -2 sigma
7: mu=0 fit status (only meaningful if asimov data is generated within the macro)

It will also store the result of the old bands procedure in a TH1D named 'limit_old'. 




//////////////////////////


This version is functionally fully consistent with the previous tag.

NOTE: The script runs significantly faster when compiled
*/




#include "TMath.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/QuantFuncMathCore.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"

#include "RooWorkspace.h"
#include "RooNLLVar.h"
#include "RooStats/ModelConfig.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "Math/MinimizerOptions.h"
#include "TStopwatch.h"
#include "RooMinimizerFcn.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooProduct.h"

#include <map>
#include <iostream>
#include <sstream>

using namespace std;
using namespace RooFit;
using namespace RooStats;

//band configuration
bool betterBands           = 1; // (recommendation = 1) improve bands by using a more appropriate asimov dataset for those points
bool betterNegativeBands   = 0; // (recommendation = 0) improve also the negative bands
bool profileNegativeAtZero = 0; // (recommendation = 0) profile asimov for negative bands at zero


//other configuration
string defaultMinimizer    = "Minuit2";     // or "Minuit"
int defaultPrintLevel      = -1;            // Minuit print level
int defaultStrategy        = 1;             // Minimization strategy. 0-2. 0 = fastest, least robust. 2 = slowest, most robust
bool killBelowFatal        = 1;             // In case you want to suppress RooFit warnings further, set to 1
bool doBlind               = 0;             // in case your analysis is blinded
bool conditionalExpected   = 1 && !doBlind; // Profiling mode for Asimov data: 0 = conditional MLEs, 1 = nominal MLEs
bool doTilde               = 1;             // bound mu at zero if true and do the \tilde{q}_{mu} asymptotics
bool doExp                 = 1;             // compute expected limit
bool doObs                 = 1 && !doBlind; // compute observed limit
double precision           = 0.005;         // precision in mu that defines iterative cutoff
bool verbose               = 0;             // 1 = very spammy








//don't touch!
map<RooNLLVar*, double> map_nll_muhat;
map<RooNLLVar*, double> map_muhat;
map<RooDataSet*, RooNLLVar*> map_data_nll;
map<RooNLLVar*, string> map_snapshots;
RooWorkspace* w = NULL;
ModelConfig* mc = NULL;
RooDataSet* data = NULL;
RooRealVar* firstPOI = NULL;
RooNLLVar* asimov_0_nll = NULL;
RooNLLVar* obs_nll = NULL;
int nrMinimize=0;
int direction=1;
int global_status=0;
double target_CLs=0.05;



//main
void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       string folder,
		       string mass,
		       double CL);

 //for backwards compatibility
void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName = "combWS",
		       const char* modelConfigName = "ModelConfig",
		       const char* dataName = "combData",
		       const char* asimovDataName = "asimovData_0",
		       const char* conditionalSnapshot = "conditionalGlobs_0",
		       const char* nominalSnapshot = "nominalGlobs",
		       string folder = "test",
		       double mass = 130,
		       double CL = 0.95);

double getLimit(RooNLLVar* nll, double initial_guess = 0);
double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu);
double getQmu(RooNLLVar* nll, double mu);
void saveSnapshot(RooNLLVar* nll, double mu);
void loadSnapshot(RooNLLVar* nll, double mu);
double getNLL(RooNLLVar* nll);
double findCrossing(double sigma_obs, double sigma, double muhat);
void setMu(double mu);
double getQmu95_brute(double sigma, double mu);
double getQmu95(double sigma, double mu);
double calcCLs(double qmu_tilde, double sigma, double mu);
double calcPmu(double qmu_tilde, double sigma, double mu);
double calcPb(double qmu_tilde, double sigma, double mu);
double calcDerCLs(double qmu, double sigma, double mu);
int minimize(RooNLLVar* nll);
int minimize(RooAbsReal* nll);
RooDataSet* makeAsimovData2(RooDataSet* conditioningData, double mu_true, double mu_prof = -999, string* mu_str = NULL, string* mu_prof_str = NULL);
RooDataSet* makeAsimovData2(RooNLLVar* conditioningNLL, double mu_true, double mu_prof = -999, string* mu_str = NULL, string* mu_prof_str = NULL);


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////



void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       const char* conditionalSnapshot,
		       const char* nominalSnapshot,
		       string folder,
		       double mass,
		       double CL)
{
  stringstream smass;
  smass << mass;

  conditionalSnapshot = ""; // warningless compile
  nominalSnapshot = "";     // warningless compile

  runAsymptoticsCLs(infile, workspaceName, modelConfigName, dataName, asimovDataName, folder, smass.str(), CL);
}



void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       string folder,
		       string mass,
		       double CL)
{
  TStopwatch timer;
  timer.Start();

  if (killBelowFatal) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(defaultMinimizer.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(defaultStrategy);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(defaultPrintLevel);
  //RooNLLVar::SetIgnoreZeroEntries(1);

//check inputs
  TFile f(infile);
  w = (RooWorkspace*)f.Get(workspaceName);
  if (!w)
  {
    cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
    return;
  }

  mc = (ModelConfig*)w->obj(modelConfigName);
  if (!mc)
  {
    cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
    return;
  }
  firstPOI = (RooRealVar*)mc->GetParametersOfInterest()->first();

  data = (RooDataSet*)w->data(dataName);
  if (!data)
  {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    return;
  }

  RooAbsPdf* pdf = mc->GetPdf();
  obs_nll = (RooNLLVar*)pdf->createNLL(*data);
  map_snapshots[obs_nll] = "nominalGlobs";
  map_data_nll[data] = obs_nll;
  w->saveSnapshot("nominalGlobs",*mc->GetGlobalObservables());
  w->saveSnapshot("nominalNuis",*mc->GetNuisanceParameters());

  RooDataSet* asimovData_0 = (RooDataSet*)w->data(asimovDataName);
  if (!asimovData_0)
  {
    asimovData_0 = makeAsimovData2((conditionalExpected ? obs_nll : (RooNLLVar*)NULL), 0., 0.);
  }
  
  asimov_0_nll = (RooNLLVar*)pdf->createNLL(*asimovData_0);
  map_snapshots[asimov_0_nll] = "conditionalGlobs_0";
  map_data_nll[asimovData_0] = asimov_0_nll;
  setMu(0);
  map_muhat[asimov_0_nll] = 0;
  saveSnapshot(asimov_0_nll, 0);
  w->loadSnapshot("conditionalNuis_0");
  w->loadSnapshot("conditionalGlobs_0");
  map_nll_muhat[asimov_0_nll] = asimov_0_nll->getVal();

  
  target_CLs=1-CL;

  
  double med_limit = doExp ? getLimit(asimov_0_nll, 1.0) : 1.0;

  double sigma = med_limit/sqrt(3.84); // pretty close
  double mu_up_p2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 2), 1) + 2);
  double mu_up_p1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf( 1), 1) + 1);
  double mu_up_n1_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-1), 1) - 1);
  double mu_up_n2_approx = sigma*(ROOT::Math::gaussian_quantile(1 - target_CLs*ROOT::Math::gaussian_cdf(-2), 1) - 2);

  double mu_up_p2 = mu_up_p2_approx;
  double mu_up_p1 = mu_up_p1_approx;
  double mu_up_n1 = mu_up_n1_approx;
  double mu_up_n2 = mu_up_n2_approx;

  firstPOI->setRange(-5*sigma, 5*sigma);
  if (betterBands && doExp) // no better time than now to do this
  {
    //find quantiles, starting with +2, since we should be at +1.96 right now

    double init_targetCLs = target_CLs;
    firstPOI->setRange(-5*sigma, 5*sigma);
    for (int N=2;N>=-2;N--)
    {
      if (N < 0 && !betterNegativeBands) continue;
      if (N == 0) continue;
      target_CLs=2*(1-ROOT::Math::gaussian_cdf(fabs(N))); // change this so findCrossing looks for sqrt(qmu95)=2
      if (N < 0) direction = -1;

      //get the acual value
      double NtimesSigma = getLimit(asimov_0_nll, N*med_limit/sqrt(3.84)); // use N * sigma(0) as an initial guess

      sigma = NtimesSigma/N;
      cout << endl;
      cout << "Found N * sigma = " << N << " * " << sigma << endl;

      string muStr,muStrPr;
      w->loadSnapshot("conditionalGlobs_0");
      double pr_val = NtimesSigma;
      if (N < 0 && profileNegativeAtZero) pr_val = 0;
      RooDataSet* asimovData_N = makeAsimovData2(asimov_0_nll, NtimesSigma, pr_val, &muStr, &muStrPr);


      RooNLLVar* asimov_N_nll = (RooNLLVar*)pdf->createNLL(*asimovData_N);
      map_data_nll[asimovData_N] = asimov_N_nll;
      map_snapshots[asimov_N_nll] = "conditionalGlobs"+muStrPr;
      w->loadSnapshot(map_snapshots[asimov_N_nll].c_str());
      w->loadSnapshot(("conditionalNuis"+muStrPr).c_str());
      setMu(NtimesSigma);

      double nll_val = asimov_N_nll->getVal();
      saveSnapshot(asimov_N_nll, NtimesSigma);
      map_muhat[asimov_N_nll] = NtimesSigma;
      if (N < 0 && doTilde)
      {
	setMu(0);
	firstPOI->setConstant(1);
	nll_val = getNLL(asimov_N_nll);
      }
      map_nll_muhat[asimov_N_nll] = nll_val;

      target_CLs = init_targetCLs;
      direction=1;
      double initial_guess = findCrossing(NtimesSigma/N, NtimesSigma/N, NtimesSigma);
      double limit = getLimit(asimov_N_nll, initial_guess);
      if (N == 2) mu_up_p2 = limit;
      else if (N == 1) mu_up_p1 = limit;
      else if (N ==-1) mu_up_n1 = limit;
      else if (N ==-2) mu_up_n2 = limit;
      //return;
    }
    direction = 1;
    target_CLs = init_targetCLs;

  }

  w->loadSnapshot("conditionalNuis_0");
  double obs_limit = doObs ? getLimit(obs_nll, med_limit) : 0;

  if (betterBands) cout << "Guess for bands" << endl;
  cout << "+2sigma:  " << mu_up_p2_approx << endl;
  cout << "+1sigma:  " << mu_up_p1_approx << endl;
  cout << "-1sigma:  " << mu_up_n1_approx << endl;
  cout << "-2sigma:  " << mu_up_n2_approx << endl;
  if (betterBands)
  {
    cout << endl;
    cout << "Correct bands" << endl;
    cout << "+2sigma:  " << mu_up_p2 << endl;
    cout << "+1sigma:  " << mu_up_p1 << endl;
    cout << "-1sigma:  " << mu_up_n1 << endl;
    cout << "-2sigma:  " << mu_up_n2 << endl;

  }

  cout << "Median:   " << med_limit << endl;
  cout << "Observed: " << obs_limit << endl;
  cout << endl;


  system(("mkdir -vp root-files/" + folder).c_str());
  stringstream fileName;
  fileName << "root-files/" << folder << "/" << mass << ".root";
  TFile fout(fileName.str().c_str(),"recreate");

  TH1D* h_lim = new TH1D("limit","limit",7,0,7);
  h_lim->SetBinContent(1, obs_limit);
  h_lim->SetBinContent(2, med_limit);
  h_lim->SetBinContent(3, mu_up_p2);
  h_lim->SetBinContent(4, mu_up_p1);
  h_lim->SetBinContent(5, mu_up_n1);
  h_lim->SetBinContent(6, mu_up_n2);
  h_lim->SetBinContent(7, global_status);

  h_lim->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim->GetXaxis()->SetBinLabel(6, "-2sigma");
  h_lim->GetXaxis()->SetBinLabel(7, "Global status"); // do something with this later

  TH1D* h_lim_old = new TH1D("limit_old","limit_old",7,0,7); // include also old approximation of bands
  h_lim_old->SetBinContent(1, obs_limit);
  h_lim_old->SetBinContent(2, med_limit);
  h_lim_old->SetBinContent(3, mu_up_p2_approx);
  h_lim_old->SetBinContent(4, mu_up_p1_approx);
  h_lim_old->SetBinContent(5, mu_up_n1_approx);
  h_lim_old->SetBinContent(6, mu_up_n2_approx);
  h_lim_old->SetBinContent(7, global_status);

  h_lim_old->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim_old->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim_old->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim_old->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim_old->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim_old->GetXaxis()->SetBinLabel(6, "-2sigma");
  h_lim_old->GetXaxis()->SetBinLabel(7, "Global status"); 

  fout.Write();
  fout.Close();

  cout << "Finished with " << nrMinimize << " calls to minimize(nll)" << endl;
  timer.Print();
}

double getLimit(RooNLLVar* nll, double initial_guess)
{
  cout << "------------------------" << endl;
  cout << "Getting limit for nll: " << nll->GetName() << endl;
  //get initial guess based on muhat and sigma(muhat)
  firstPOI->setConstant(0);

  if (nll == asimov_0_nll) 
  {
    setMu(0);
    firstPOI->setConstant(1);
  }

  double muhat;
  if (map_nll_muhat.find(nll) == map_nll_muhat.end())
  {
    double nll_val = getNLL(nll);
    muhat = firstPOI->getVal();
    saveSnapshot(nll, muhat);
    map_muhat[nll] = muhat;
    if (muhat < 0 && doTilde)
    {
      setMu(0);
      firstPOI->setConstant(1);
      nll_val = getNLL(nll);
    }

    map_nll_muhat[nll] = nll_val;
  }
  else
  {
    muhat = map_muhat[nll];
  }

  if (muhat < 0.1 || initial_guess != 0) setMu(initial_guess);
  double qmu,qmuA;
  double sigma_guess = getSigma(asimov_0_nll, firstPOI->getVal(), 0, qmu);
  double sigma_b = sigma_guess;
  double mu_guess = findCrossing(sigma_guess, sigma_b, muhat);
  double pmu = calcPmu(qmu, sigma_b, mu_guess);
  double pb = calcPb(qmu, sigma_b, mu_guess);
  double CLs = calcCLs(qmu, sigma_b, mu_guess);
  double qmu95 = getQmu95(sigma_b, mu_guess);
  setMu(mu_guess);

  cout << "Initial guess:  " << mu_guess << endl;
  cout << "Sigma(obs):     " << sigma_guess << endl;
  cout << "Sigma(mu,0):    " << sigma_b << endl;
  cout << "muhat:          " << muhat << endl;
  cout << "qmu95:          " << qmu95 << endl;
  cout << "qmu:            " << qmu << endl;
  cout << "pmu:            " << pmu << endl;
  cout << "pb:             " << pb << endl;
  cout << "CLs:            " << CLs << endl;
  cout << endl;

  int nrDamping = 1;
  map<double, double> guess_to_corr;
  double damping_factor = 1.0;
  double damping_factor_pre = damping_factor;
  int nrItr = 0;
  double mu_pre = mu_guess-10*precision*mu_guess;
  while (fabs(mu_pre-mu_guess) > precision*mu_guess*direction)
  {
    cout << "----------------------" << endl;
    cout << "Starting iteration " << nrItr << " of " << nll->GetName() << endl;
    damping_factor_pre = damping_factor;
    if (nrItr != 0) loadSnapshot(nll, mu_pre); // do this to avoid comparing multiple minima in the conditional and unconditional fits
    sigma_guess=getSigma(nll, mu_guess, muhat, qmu);
    saveSnapshot(nll, mu_guess);


    if (nll != asimov_0_nll)
    {
      if (nrItr != 0) loadSnapshot(asimov_0_nll, mu_pre);
      sigma_b=getSigma(asimov_0_nll, mu_guess, 0, qmuA);
      saveSnapshot(asimov_0_nll, mu_guess);
    }
    else
    {
      sigma_b=sigma_guess;
      qmuA=qmu;
    }

    double corr = damping_factor*(mu_guess - findCrossing(sigma_guess, sigma_b, muhat));
    for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++)
    {
      if (fabs(itr->first - (mu_guess-corr)) < direction*mu_guess*0.02 && fabs(corr) > direction*mu_guess*precision) 
      {
	damping_factor *= 0.8;
	cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	if (nrDamping++ > 10)
	{
	  nrDamping = 1;
	  damping_factor = 1.0;
	}
	corr *= damping_factor;
	break;
      }
    }

    //subtract off the difference in the new and damped correction
    guess_to_corr[mu_guess] = corr;
    mu_pre = mu_guess;
    mu_guess -= corr;


    pmu = calcPmu(qmu, sigma_b, mu_pre);
    pb = calcPb(qmu, sigma_b, mu_pre);
    CLs = calcCLs(qmu, sigma_b, mu_pre);
    qmu95 = getQmu95(sigma_b, mu_pre);


    cout << "NLL:            " << nll->GetName() << endl;
    cout << "Previous guess: " << mu_pre << endl;
    cout << "Sigma(obs):     " << sigma_guess << endl;
    cout << "Sigma(mu,0):    " << sigma_b << endl;
    cout << "muhat:          " << muhat << endl;
    cout << "pmu:            " << pmu << endl;
    cout << "pb:             " << pb << endl;
    cout << "CLs:            " << CLs << endl;
    cout << "qmu95:          " << qmu95 << endl;
    cout << "qmu:            " << qmu << endl;
    cout << "qmuA:           " << qmuA << endl;
    cout << "Precision:      " << direction*mu_guess*precision << endl;
    cout << "Correction:    "  << (-corr<0?" ":"") << -corr << endl;
    cout << "New guess:      " << mu_guess << endl;
    cout << endl;

      

    nrItr++;
    if (nrItr > 25)
    {
      cout << "Infinite loop detected in getLimit(). Please intervene." << endl;
      break;
    }
  }


  cout << "Found limit for nll " << nll->GetName() << ": " << mu_guess << endl;
  cout << "Finished in " << nrItr << " iterations." << endl;
  cout << endl;
  return mu_guess;
}


double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu)
{
  qmu = getQmu(nll, mu);
  if (verbose) cout << "qmu = " << qmu << endl;
  if (mu*direction < muhat) return fabs(mu-muhat)/sqrt(qmu);
  else if (muhat < 0 && doTilde) return sqrt(mu*mu-2*mu*muhat*direction)/sqrt(qmu);
  else return (mu-muhat)*direction/sqrt(qmu);
}

double getQmu(RooNLLVar* nll, double mu)
{
  double nll_muhat = map_nll_muhat[nll];
  bool isConst = firstPOI->isConstant();
  firstPOI->setConstant(1);
  setMu(mu);
  double nll_val = getNLL(nll);
  firstPOI->setConstant(isConst);
  //cout << "qmu = 2 * (" << nll_val << " - " << nll_muhat << ")" << endl;
  return 2*(nll_val-nll_muhat);
}

void saveSnapshot(RooNLLVar* nll, double mu)
{
  stringstream snapshotName;
  snapshotName << nll->GetName() << "_" << mu;
  w->saveSnapshot(snapshotName.str().c_str(), *mc->GetNuisanceParameters());
}

void loadSnapshot(RooNLLVar* nll, double mu)
{
  stringstream snapshotName;
  snapshotName << nll->GetName() << "_" << mu;
  w->loadSnapshot(snapshotName.str().c_str());
}

double getNLL(RooNLLVar* nll)
{
  string snapshotName = map_snapshots[nll];
  if (snapshotName != "") w->loadSnapshot(snapshotName.c_str());
  minimize(nll);
  double val = nll->getVal();
  w->loadSnapshot("nominalGlobs");
  return val;
}


double findCrossing(double sigma_obs, double sigma, double muhat)
{
  double mu_guess = muhat + ROOT::Math::gaussian_quantile(1-target_CLs,1)*sigma_obs*direction;
  int nrItr = 0;
  int nrDamping = 1;

  map<double, double> guess_to_corr;
  double damping_factor = 1.0;
  double mu_pre = mu_guess - 10*mu_guess*precision;
  while (fabs(mu_guess-mu_pre) > direction*mu_guess*precision)
  {
    mu_pre = mu_guess;
    double qmu95 = getQmu95(sigma, mu_guess);
    double qmu = 1./sigma_obs/sigma_obs*(mu_guess-muhat)*(mu_guess-muhat);
    if (muhat < 0 && doTilde) qmu = 1./sigma_obs/sigma_obs*(mu_guess*mu_guess-2*mu_guess*muhat);

    double corr = damping_factor*(qmu-qmu95)/(2./sigma_obs/sigma_obs*(mu_guess-muhat));
    for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++)
    {
      if (fabs(itr->first - mu_guess) < direction*mu_guess*precision) 
      {
	damping_factor *= 0.8;
	if (verbose) cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	if (nrDamping++ > 10)
	{
	  nrDamping = 1;
	  damping_factor = 1.0;
	}
	corr *= damping_factor;
	break;
      }
    }
    guess_to_corr[mu_guess] = corr;

    mu_guess = mu_guess - corr;
    nrItr++;
    if (nrItr > 100)
    {
      cout << "Infinite loop detected in findCrossing. Please intervene." << endl;
      exit(1);
    }
    if (verbose) cout << "mu_guess = " << mu_guess << ", mu_pre = " << mu_pre << ", qmu = " << qmu << ", qmu95 = " << qmu95 << ", sigma_obs = " << sigma_obs << ", sigma = " << sigma << ", direction*mu*prec = " << direction*mu_guess*precision << endl;
  }

  return mu_guess;
}

void setMu(double mu)
{
  if (mu > 0 && firstPOI->getMax() < mu) firstPOI->setMax(2*mu);
  if (mu < 0 && firstPOI->getMin() > mu) firstPOI->setMin(2*mu);
  firstPOI->setVal(mu);
}


double getQmu95_brute(double sigma, double mu)
{
  double step_size = 0.001;
  double start = step_size;
  if (mu/sigma > 0.2) start = 0;
  for (double qmu=start;qmu<20;qmu+=step_size)
  {
    double CLs = calcCLs(qmu, sigma, mu);

    if (CLs < target_CLs) return qmu;
  }

  return 20;
}

double getQmu95(double sigma, double mu)
{
  double qmu95 = 0;
  //no sane man would venture this far down into |mu/sigma|
  double target_N = ROOT::Math::gaussian_cdf(1-target_CLs,1);
  if (fabs(mu/sigma) < 0.25*target_N)
  {
    qmu95 = 5.83/target_N;
  }
  else
  {
    map<double, double> guess_to_corr;
    double qmu95_guess = pow(ROOT::Math::gaussian_quantile(1-target_CLs,1),2);
    int nrItr = 0;
    int nrDamping = 1;
    double damping_factor = 1.0;
    double qmu95_pre = qmu95_guess - 10*2*qmu95_guess*precision;
    while (fabs(qmu95_guess-qmu95_pre) > 2*qmu95_guess*precision)
    {
      qmu95_pre = qmu95_guess;
      if (verbose)
      {
	cout << "qmu95_guess = " << qmu95_guess << endl;
	cout << "CLs = " << calcCLs(qmu95_guess, sigma, mu) << endl;
	cout << "Derivative = " << calcDerCLs(qmu95_guess, sigma, mu) << endl;
      }

      double corr = damping_factor*(calcCLs(qmu95_guess, sigma, mu)-target_CLs)/calcDerCLs(qmu95_guess, sigma, mu);
      for (map<double, double>::iterator itr=guess_to_corr.begin();itr!=guess_to_corr.end();itr++)
      {
	if (fabs(itr->first - qmu95_guess) < 2*qmu95_guess*precision) 
	{
	  damping_factor *= 0.8;
	  if (verbose) cout << "Changing damping factor to " << damping_factor << ", nrDamping = " << nrDamping << endl;
	  if (nrDamping++ > 10)
	  {
	    nrDamping = 1;
	    damping_factor = 1.0;
	  }
	  corr *= damping_factor;
	}
      }

      guess_to_corr[qmu95_guess] = corr;
      qmu95_guess = qmu95_guess - corr;

      if (verbose)
      {
	cout << "next guess = " << qmu95_guess << endl; 
	cout << "precision = " << 2*qmu95_guess*precision << endl;
	cout << endl;
      }
      nrItr++;
      if (nrItr > 200)
      {
	cout << "Infinite loop detected in getQmu95. Please intervene." << endl;
	exit(1);
      }
    }
    qmu95 = qmu95_guess;
  }

  if (qmu95 != qmu95) 
  {
    qmu95 = getQmu95_brute(sigma, mu);
  }
  if (verbose) cout << "Returning qmu95 = " << qmu95 << endl;

  return qmu95;
}

double calcCLs(double qmu_tilde, double sigma, double mu)
{
  double pmu = calcPmu(qmu_tilde, sigma, mu);
  double pb = calcPb(qmu_tilde, sigma, mu);
  if (verbose)
  {
    cout << "pmu = " << pmu << endl;
    cout << "pb = " << pb << endl;
  }
  if (pb == 1) return 0.5;
  return pmu/(1-pb);
}

double calcPmu(double qmu, double sigma, double mu)
{
  double pmu;
  if (qmu < mu*mu/(sigma*sigma) || !doTilde)
  {
    pmu = 1-ROOT::Math::gaussian_cdf(sqrt(qmu));
  }
  else
  {
    pmu = 1-ROOT::Math::gaussian_cdf((qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma)));
  }
  if (verbose) cout << "for pmu, qmu = " << qmu << ", sigma = " << sigma<< ", mu = " << mu << ", pmu = " << pmu << endl;
  return pmu;
}

double calcPb(double qmu, double sigma, double mu)
{
  if (qmu < mu*mu/(sigma*sigma) || !doTilde)
  {
    return 1-ROOT::Math::gaussian_cdf(fabs(mu/sigma) - sqrt(qmu));
  }
  else
  {
    return 1-ROOT::Math::gaussian_cdf((mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma)));
  }
}

double calcDerCLs(double qmu, double sigma, double mu)
{
  double dpmu_dq = 0;
  double d1mpb_dq = 0;

  if (qmu < mu*mu/(sigma*sigma))
  {
    double zmu = sqrt(qmu);
    dpmu_dq = -1./(2*sqrt(qmu*2*TMath::Pi()))*exp(-zmu*zmu/2);
  }
  else 
  {
    double zmu = (qmu+mu*mu/(sigma*sigma))/(2*fabs(mu/sigma));
    dpmu_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zmu*zmu/2);
  }

  if (qmu < mu*mu/(sigma*sigma))
  {
    double zb = fabs(mu/sigma)-sqrt(qmu);
    d1mpb_dq = -1./sqrt(qmu*2*TMath::Pi())*exp(-zb*zb/2);
  }
  else
  {
    double zb = (mu*mu/(sigma*sigma) - qmu)/(2*fabs(mu/sigma));
    d1mpb_dq = -1./(2*fabs(mu/sigma))*1./(sqrt(2*TMath::Pi()))*exp(-zb*zb/2);
  }

  double pb = calcPb(qmu, sigma, mu);
  return dpmu_dq/(1-pb)-calcCLs(qmu, sigma, mu)/(1-pb)*d1mpb_dq;
}

int minimize(RooNLLVar* nll)
{
  nrMinimize++;
  RooAbsReal* fcn = (RooAbsReal*)nll;
  return minimize(fcn);
}

int minimize(RooAbsReal* fcn)
{

//    cout << "Starting minimization. Using these global observables" << endl;
//    mc->GetGlobalObservables()->Print("v");


  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
  int save_strat = strat;
  RooMinimizer minim(*fcn);
  minim.setStrategy(strat);
  minim.setPrintLevel(printLevel);


  int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


//up the strategy
  if (status != 0 && status != 1 && strat < 2)
  {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }

  if (status != 0 && status != 1 && strat < 2)
  {
    strat++;
    cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
    minim.setStrategy(strat);
    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
  }


// //switch minuit version and try again
  if (status != 0 && status != 1)
  {
    string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
    string newMinType;
    if (minType == "Minuit2") newMinType = "Minuit";
    else newMinType = "Minuit2";
  
    cout << "Switching minuit type from " << minType << " to " << newMinType << endl;
  
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
    strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
    minim.setStrategy(strat);

    status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


    if (status != 0 && status != 1 && strat < 2)
    {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1 && strat < 2)
    {
      strat++;
      cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
      minim.setStrategy(strat);
      status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
  }

  if (status != 0 && status != 1)
  {
    cout << "Fit failed with status " << status << endl;
  }

  if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(msglevel);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(save_strat);

  return status;
}



RooDataSet* makeAsimovData2(RooDataSet* conditioningData, double mu_true, double mu_prof, string* mu_str, string* mu_prof_str)
{
  RooNLLVar* conditioningNLL = NULL;
  if (conditioningData)
  {
    conditioningNLL = (RooNLLVar*)mc->GetPdf()->createNLL(*conditioningData);
  }
  return makeAsimovData2(conditioningNLL, mu_true, mu_prof, mu_str, mu_prof_str);
}


RooDataSet* makeAsimovData2(RooNLLVar* conditioningNLL, double mu_true, double mu_prof, string* mu_str, string* mu_prof_str)
{
  if (mu_prof == -999) mu_prof = mu_true;
  bool doTest = 0;

  cout << "Creating asimov data at mu = " << mu_true << ", profiling at mu = " << mu_prof << endl;
  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();

  int test = 0;
  if (doTest) cout << "test = " << test++ << endl;

  int _printLevel = 0;

  stringstream muStr;
  muStr << setprecision(5);
  muStr << "_" << mu_true;
  if (mu_str) *mu_str = muStr.str();

  stringstream muStrProf;
  muStrProf << setprecision(5);
  muStrProf << "_" << mu_prof;
  if (mu_prof_str) *mu_prof_str = muStrProf.str();


  if (doTest) cout << "test = " << test++ << endl;
  const RooArgSet* globs = mc->GetGlobalObservables();
  const RooArgSet* nuis = mc->GetNuisanceParameters();
  const RooArgSet* obs = mc->GetObservables();
  const RooArgSet* pois = mc->GetParametersOfInterest();

  TIterator* gItr = globs->createIterator();
  TIterator* nItr = nuis->createIterator();
  RooRealVar* var;

  //cout << "test = " << test++ << endl;
  RooArgSet emptySet;
  RooArgSet params(*nuis);
  params.add(*globs);
  params.add(*pois);
  w->saveSnapshot("initial_params", params);
  
  if (doTest) cout << "test = " << test++ << endl;

//condition the MLEs
  if (conditioningNLL)
  {
    //get the conditional MLEs
    firstPOI->setVal(mu_prof);
    firstPOI->setConstant(1);
    minimize(conditioningNLL);
  }

  if (doTest) cout << "test = " << test++ << endl;
  w->saveSnapshot(("conditionalNuis" +muStrProf.str()).c_str(),*nuis);


//to find the conditional globs, do a fit to the constraint only pdf with the globs floating and the MLEs constant    
  RooArgSet obsCopy = *obs;
  RooArgSet nuisCopy = *nuis;
    
  RooArgSet constraints(*mc->GetPdf()->getAllConstraints(obsCopy, nuisCopy));
  RooRealVar minusOne("minusOne","minusOne",-1);
  constraints.add(minusOne);
  RooProduct constrFunc("constrFunc","constrFunc",constraints);
    
  if (doTest) cout << "test = " << test++ << endl;
  while ((var = (RooRealVar*)gItr->Next()))
  {
    var->setConstant(false);
  }
  gItr->Reset();

  while ((var = (RooRealVar*)nItr->Next()))
  {
    var->setConstant(true);
  }
  nItr->Reset();

  minimize(&constrFunc);

  while ((var = (RooRealVar*)gItr->Next()))
  {
    var->setConstant(true);
  }
  gItr->Reset();

  while ((var = (RooRealVar*)nItr->Next()))
  {
    var->setConstant(false);
  }
  nItr->Reset();

  w->saveSnapshot(("conditionalGlobs"+muStrProf.str()).c_str(),*globs);
  

  if (doTest) cout << "test = " << test++ << endl;
  


//make the asimov data
  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mc->GetObservables());

  RooRealVar* weightVar = NULL;
  if (!(weightVar = w->var(weightName)))
  {
    w->import(*(new RooRealVar(weightName, weightName, 1,0,10000000)));
    weightVar = w->var(weightName);
  }
  obsAndWeight.add(*w->var(weightName));


  if (doTest) cout << "test = " << test++ << endl;

  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());
  map<string, RooDataSet*> asimovDataMap;

  //try fix for sim pdf
  RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
  TIterator* iter = channelCat->typeIterator() ;
  RooCatType* tt = NULL;
  int nrIndices = 0;
  int iFrame=0;
  while((tt=(RooCatType*) iter->Next())) {
    nrIndices++;
  }
  for (int i=0;i<nrIndices;i++){
    channelCat->setIndex(i);
    iFrame++;
    // Get pdf associated with state from simpdf
    RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;
	
    // Generate observables defined by the pdf associated with this state
    RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;

    if (_printLevel >= 1)
    {
      obstmp->Print();
      cout << "on type " << channelCat->getLabel() << " " << iFrame << endl;
    }

    RooDataSet* obsDataUnbinned = new RooDataSet(Form("combAsimovData%d",iFrame),Form("combAsimovData%d",iFrame),RooArgSet(obsAndWeight,*channelCat),WeightVar(*weightVar));
    RooRealVar* thisObs = ((RooRealVar*)obstmp->first());
    double expectedEvents = pdftmp->expectedEvents(*obstmp);
    double thisNorm = 0;
    for(int jj=0; jj<thisObs->numBins(); ++jj){
      thisObs->setBin(jj);

      thisNorm=pdftmp->getVal(obstmp)*thisObs->getBinWidth(jj);
      if (thisNorm*expectedEvents > 0 && thisNorm*expectedEvents < pow(10.0, 18)) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
    }
    
    if (_printLevel >= 1)
    {
      obsDataUnbinned->Print();
      cout <<"sum entries "<<obsDataUnbinned->sumEntries()<<endl;
    }
    if(obsDataUnbinned->sumEntries()!=obsDataUnbinned->sumEntries()){
      cout << "sum entries is nan"<<endl;
      exit(1);
    }


    asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;

    if (_printLevel >= 1)
    {
      cout << "channel: " << channelCat->getLabel() << ", data: ";
      obsDataUnbinned->Print();
      cout << endl;
    }
  }

  if (doTest) cout << "test = " << test++ << endl;
  RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
  if (w->data(("asimovData"+muStr.str()).c_str()))
  {
    w->import(*asimovData, true);
  }
  else
  {
    w->import(*asimovData);
  }




  w->loadSnapshot("initial_params");
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

  if (doTest) cout << "test = " << test++ << endl;
  return asimovData;
  
}

