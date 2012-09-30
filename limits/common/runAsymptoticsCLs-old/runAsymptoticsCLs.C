/*
Author: Aaron Armbruster
Date:   2011-07-27 (Updated 2011-09-01)
Email:  armbrusa@umich.edu
Description: Script to run asymptotic CLs. 'conditionalSnapshot' and 'nominalSnapshot' are only required when running over workspaces 
             where the Asimov data has been generated at the conditional MLEs. If running over workspaces made by me, the default
             arguments are correct to use. Running over an asimov data from histfactory, these should be set to empty strings ("").


             The script uses an iterative method to find the upper limit by using the sequence:

       mu_N+1 = mu_N - (p_mu - (1-CL)) / [dp_mu / dmu]
              ~ mu_N - (p_mu - (1-CL)) / [(p_(mu*(1+epsilon)) - p_mu) / (mu*(1+epsilon) - mu)]

             epsilon is by default 0.05, but is adaptive.

             Most of the conditional statements in the main iterative loop are to avoid convergence problems, infs and nans. 99% of these
             have been worked out, but occasionally it can still fail to converge under extreme datasets. If this happens please adjust the print level
	     using ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1) and send me the output and the workspace.

             By default this uses \tilde{qmu}, but changing the range of mu in the iteration function to go negative will allow it to do qmu.
             Of course, with asymptotics both give equivalent results, but tilde avoids negative poisson means and is safer.

	     After running, the results will be printed as well as stored in a root file in the folder 'root-files/<folder>', where <folder>
	     is specified by you (default 'test')

	     The root file has a 6-bin TH1D, where each bin is filled with the upper limit values in this order:

	     1: Observed
	     2: Median
	     3: +2 sigma
	     4: +1 sigma
	     5: -1 sigma
	     6: -2 sigma

             
NOTE: The script runs significantly faster when compiled
*/


#include "TStopwatch.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "Math/MinimizerOptions.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooMinimizerFcn.h"
#include "RooMinimizer.h"

#include "TFile.h"
#include "TH1D.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace RooFit;
using namespace RooStats;

//global test stat objects
ProfileLikelihoodTestStat* asimov_testStat;
ProfileLikelihoodTestStat* obs_testStat;

//call this function to run!
double runAsymptoticsCLs(const char* infile,
			 const char* workspaceName = "combWS",
			 const char* modelConfigName = "ModelConfig",
			 const char* dataName = "combData",
			 const char* asimovDataName = "asimovData_0",
			 const char* conditionalSnapshot = "conditionalGlobs_0",
			 const char* nominalSnapshot = "nominalGlobs",
			 string folder = "test",
			 int mass = 130,
			 double CL = 0.95);


//common function to evaluate test stat
double evaluate(RooDataSet* data, RooDataSet* asimovData, RooArgSet& poi, RooWorkspace* ws, const char* conditionalSnapshot, const char* nominalSnapshot);

//make asimov data if needed
void makeAsimovData(ModelConfig* mcInWs, bool doConditional, RooWorkspace* combWS, RooAbsPdf* combPdf, RooDataSet* combData, bool b_only);

//helper for makeAsimovData
void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);

//calculate CLs
double calcCLs(double& mu_val, RooRealVar* mu, RooDataSet* data, RooDataSet* asimovData, const char* conditionalSnapshot, const char* nominalSnapshot, RooWorkspace* ws, double& testStat_val, double ncpScale);

//compute the limit once
double getLimit(RooWorkspace* ws,
		ModelConfig* mc,
		RooDataSet* data,
		RooDataSet* asimovData,
		double& testStat_val,
		const char* conditionalSnapshot,
		const char* nominalSnapshot,
		double CL,
		double initialGuess,
		double ncpScale);

//common function for setting up the test stat object
void setupTestStat(ProfileLikelihoodTestStat*& testStat, RooAbsPdf* pdf);

//get the +- N sigma quantiles according to frequentist recommendation paper
double getBandVal(int N, double alpha, double sigma);

//run the tool
double runAsymptoticsCLs(const char* infile,
			 const char* workspaceName,
			 const char* modelConfigName,
			 const char* dataName,
			 const char* asimovDataName,
			 const char* conditionalSnapshot,
			 const char* nominalSnapshot,
			 string folder,
			 int mass,
			 double CL)
{

//scale factor for non-centrality parameter (should be 1 in general, unless you've calibrated it)
  double ncpScale = 1.0;//1.35;

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

//check inputs
  TFile f(infile);
  RooWorkspace* ws = (RooWorkspace*)f.Get(workspaceName);
  if (!ws)
  {
    cout << "ERROR::Workspace: " << workspaceName << " doesn't exist!" << endl;
    return 0;
  }

  ModelConfig* mc = (ModelConfig*)ws->obj(modelConfigName);
  if (!mc)
  {
    cout << "ERROR::ModelConfig: " << modelConfigName << " doesn't exist!" << endl;
    return 0;
  }

  RooDataSet* data = (RooDataSet*)ws->data(dataName);
  if (!data)
  {
    cout << "ERROR::Dataset: " << dataName << " doesn't exist!" << endl;
    return 0;
  }

  RooDataSet* asimovData = (RooDataSet*)ws->data(asimovDataName);
  if (!asimovData)
  {
    bool doConditional = 1;
    cout << "Asimov data doesn't exist! Please, allow me to build one for you..." << endl;
    makeAsimovData(mc, doConditional, ws, mc->GetPdf(), data, 1);
    asimovData = (RooDataSet*)ws->data("asimovData_0");
    if (!asimovData)
    {
      cout << "Error in making asimov data." << endl;
      return 0;
    }
  }


//setup the test stats
  setupTestStat(asimov_testStat, mc->GetPdf());
  setupTestStat(obs_testStat, mc->GetPdf());


//get median limit, then use this as an initial guess for the observed
  double testStat_val = 1.0;
  double testStat_val_obs;
  double med_limit = getLimit(ws, mc, asimovData, asimovData, testStat_val,     conditionalSnapshot, nominalSnapshot, CL, 1,        ncpScale);
  double limit =     getLimit(ws, mc, data,       asimovData, testStat_val_obs, conditionalSnapshot, nominalSnapshot, CL, med_limit, ncpScale);

  if (testStat_val <= 0)
  {
    cout << "ERROR::test stat <= 0!" << endl;
    return 0.;
  }
  double sigma = med_limit/sqrt(testStat_val);

//print results
  cout << endl;
  cout << "Observed: " << limit << endl;
  cout << "Median: " << med_limit << endl;

  cout << "+2Sigma: " << getBandVal( 2, 1 - CL, sigma) << endl;
  cout << "+1Sigma: " << getBandVal( 1, 1 - CL, sigma) << endl;
  cout << "-1Sigma: " << getBandVal(-1, 1 - CL, sigma) << endl;
  cout << "-2Sigma: " << getBandVal(-2, 1 - CL, sigma) << endl;

//write out
  system(("mkdir -vp root-files/" + folder).c_str());
  stringstream fileName;
  fileName << "root-files/" << folder << "/" << mass << ".root";
  TFile fout(fileName.str().c_str(),"recreate");

  TH1D* h_lim = new TH1D("limit","limit",6,0,6);
  h_lim->SetBinContent(1, limit);
  h_lim->SetBinContent(2, med_limit);
  h_lim->SetBinContent(3, getBandVal( 2, 1 - CL, sigma));
  h_lim->SetBinContent(4, getBandVal( 1, 1 - CL, sigma));
  h_lim->SetBinContent(5, getBandVal(-1, 1 - CL, sigma));
  h_lim->SetBinContent(6, getBandVal(-2, 1 - CL, sigma));

  h_lim->GetXaxis()->SetBinLabel(1, "Observed");
  h_lim->GetXaxis()->SetBinLabel(2, "Expected");
  h_lim->GetXaxis()->SetBinLabel(3, "+2sigma");
  h_lim->GetXaxis()->SetBinLabel(4, "+1sigma");
  h_lim->GetXaxis()->SetBinLabel(5, "-1sigma");
  h_lim->GetXaxis()->SetBinLabel(6, "-2sigma");


  fout.Write();
  fout.Close();

  return limit;
}

//formula from frequentist recommendation paper
double getBandVal(int N, double alpha, double sigma)
{
  return sigma*(ROOT::Math::gaussian_quantile(1 - alpha*ROOT::Math::gaussian_cdf(N), 1) + N);
}

//function to do iterative loop to find limit. if data == asimovData, it will return median limit.
double getLimit(RooWorkspace* ws,
		ModelConfig* mc,
		RooDataSet* data,
		RooDataSet* asimovData,
		double& testStat_val,
		const char* conditionalSnapshot,
		const char* nominalSnapshot,
		double CL,
		double initialGuess,
		double ncpScale)
{
//setup
  double target_p = 1-CL; //target CLs value

  int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

  RooAbsPdf* pdf = mc->GetPdf();
  RooRealVar* mu = (RooRealVar*)mc->GetParametersOfInterest()->first();
  RooArgSet nuis = *mc->GetNuisanceParameters();



//get initial guess
  mu->setConstant(0);
  mu->setRange(0, 20); // set lower boundary of 0 for \tilde{qmu}. upper bound can change dynamically.
  double sigma = 0;
  double muhat = 0;
  
  if (data != asimovData && initialGuess == 0) 
  {
    RooAbsReal* nll = pdf->createNLL(*data, RooFit::Constrain(nuis));
    RooMinimizer minim(*nll);
    minim.setStrategy(2);
    minim.setPrintLevel(1);
    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    if (status != 0)
    {
      cout << "Fit failed for mu = " << mu->getVal() << " with status " << status << endl;
    }
    //pdf->fitTo(*data,Hesse(0),Minos(0),PrintLevel(0), Constrain(nuis), Strategy(ROOT::Math::MinimizerOptions::DefaultStrategy()));
    sigma = mu->getError();
    muhat = mu->getVal();
  }



  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(printLevel);

  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  cout << "Starting to run asymptoticCLs on dataset: " << data->GetName() << endl;

//try this first
  if (mu->getMax() <= initialGuess) mu->setMax(2*initialGuess);
  double mu_val;
  if (initialGuess == 0) 
  {
    if (muhat > 0.1) mu_val = muhat + sigma*1.64;
    else mu_val = 0.5;
  }
  else
  {
    mu_val = initialGuess;
  }

  RooArgSet poi(*mu);
  double ts = 0;

  mu->setVal(mu_val);

//make sure our mu val is low enough to avoid inf/nan in the test stat due to numerical errors
  int nrLower = 0;
  while (((ts = evaluate(data, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot)) > 5) || (ts != ts) /*nan*/)
  {
    double max_val = mu->getMax();
    if (data == asimovData) mu->setMax(min(1.1*mu_val, max_val));

    cout << "Lowering mu from " << mu_val << " to " << mu_val*0.5 << " (test stat = " << ts << ")" << endl;
    mu_val *= 0.5;
    mu->setVal(mu_val);

    nrLower++;
    if (nrLower > 10)
    {
      cout << "ERROR::Problem detected in data. Please intervene." << endl;
      exit(1);
    }

    if (data == asimovData) mu->setMax(max_val);
  }
  cout << "Test stat: " << ts << endl;


//Maybe muhat > initial guess, in which case tildeqmu = 0.
//If so, fit to find muhat and set initial guess
  if (data != asimovData)
  {
    if (ts < 1)
    {
      cout << "Initial test stat gave zero. Fitting to find muhat." << endl;

      int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
      RooNLLVar* fNll = (RooNLLVar*) pdf->createNLL(*data, RooFit::CloneData(kFALSE), RooFit::Constrain(nuis));
      RooMinimizer minim(*fNll);
      minim.setStrategy(strat);
      minim.setPrintLevel(1);
      int statusN = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

      if (statusN != 0)
      {
	cout << "WARNING::Minimizer returned status " << statusN << " when finding muhat." << endl;
      }

      delete fNll;

      sigma = mu->getError();
      muhat = mu->getVal();

      mu_val = muhat + sigma*1.64;
      cout << "sigma: " << sigma << ", muhat: " << muhat << ", mu_guess: " << mu_val << endl;
      mu->setVal(mu_val);
    }
  }
  else if (ts < 1) // if getting expected, raise mu until ts > 1 to avoid very small test stats
  {
    int nrRaise = 0;

    while (((ts = evaluate(asimovData, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot)) < 1) || (ts != ts) /*nan*/)
    {
      double max_val = mu->getMax();
      if (data == asimovData) mu->setMax(min(1.1*mu_val, max_val));

      nrRaise++;
      if (nrRaise > 15)
      {
	cout << "ERROR::Problem detected in asimov. Please intervene." << endl;
	exit(1);
      }

      double mu_save = mu_val;
      mu_val *= (ts < 0.1 ? 3 : ts < 0.3 ? 3 : ts < 0.5 ? 2 : 1.5);
      cout << "Raising mu from " << mu_save << " to " << mu_val << " (test stat = " << ts << ")" << endl;

      if (data == asimovData) mu->setMax(max(2*mu_val, max_val));
      mu->setVal(mu_val);

    }
    cout << "ts: " << ts << endl;
  }

//some variables used in avoiding infinite loop bouncing back and forth between two mu vals
  double p_nm2 = 0;
  bool even = true;
  double factor = 1.;

  int maxItr = 50; // maximum iterations
  int itr = 0; //current iteration
  double pmu=1; // CLs value
  double precision = 0.005; // precision requested in limit; defines loop cutoff
  static double epsilon = 0.05; // small value used to compute numerical derivative
  static int nrInf = 0; // inf counter

  double testStat_val_dummy;
  cout << "Epsilon:   " << epsilon << endl;
  cout << "Precision: " << precision << endl;
  cout << "ncpScale: " << ncpScale << endl;
  do
  {
//control if it's an even or odd iteration
    if (even) even = false;
    else even = true;
    if (even) p_nm2 = pmu;

    cout << endl << endl;
    cout << "---------------------------" << endl;
    cout << "Starting iteration: " << itr++ << endl;

//compute pmu/(1-pb) at mu
    pmu=calcCLs(mu_val, mu, data, asimovData, conditionalSnapshot, nominalSnapshot, ws, testStat_val, ncpScale);
    cout << "CLs(mu=" << mu_val << "): " << pmu << endl;

//try again if there's some numerical error
    if (testStat_val < 0)
    {
      mu_val *=1.1;
      itr++;
      continue;
    }
    cout << endl;
    cout << "Computing at mu*(1+epsilon)" << endl;

    double mu_val2 = mu_val*(1 + epsilon);

//compute p-val ratio at mu*(1+epsilon) and find difference with previous p-val ratio
    double delta_pmu=calcCLs(mu_val2, mu, data, asimovData, conditionalSnapshot, nominalSnapshot, ws, testStat_val_dummy, ncpScale) - pmu;
    cout << "CLs(mu=" << mu_val2 << "):" << delta_pmu+pmu << endl;
    cout << "delta(CLs): " << delta_pmu << endl;

//again, try again if there's a numerical error
    if (testStat_val_dummy < 0)
    {
      mu_val *= 1 + epsilon*2;
      itr++;
      continue;
    }

//if the N-2'th iteration == N'th iteration, it will end up in an infinite loop.
//try putting an adaptive scale factor in front of the next order correction factor to aid in convergence.
    if (p_nm2 > 0 && fabs(p_nm2 - pmu)/p_nm2 < precision) 
    {
      cout << "Changing factor: " << factor*0.8 << endl;
      factor *= 0.8;

      if (factor < pow(0.8, 5.0)) factor = 1;
    }

//get the next order correction and apply it to mu
    double diff = factor*(pmu-target_p)/(delta_pmu/(mu_val2 - mu_val));
    if (diff > mu_val) diff = 0.5*mu_val;
    cout << "diff: " << diff << endl;
    mu_val -= diff;

//check for inf. retry if it is, maximum of 3 times. adjust epsilon a little higher,
//sometimes it's too small
    if (mu_val > 10e9 || mu_val < -10e9)
    {
      nrInf++;
      if (nrInf > 3) return mu_val;
      epsilon *= 2;
      return getLimit(ws, mc, data, asimovData, testStat_val, conditionalSnapshot, nominalSnapshot, CL, initialGuess, ncpScale);
    }
    if ((testStat_val_dummy > 10e9 || testStat_val_dummy < -10e9) || 
	(testStat_val > 10e9 || testStat_val < -10e9))
    {
      nrInf++;
      if (nrInf > 3) return mu_val;
      epsilon *= 2;
      return getLimit(ws, mc, data, asimovData, testStat_val, conditionalSnapshot, nominalSnapshot, CL, initialGuess, ncpScale);
    }

//precision cutoff
    if (fabs(diff/mu_val)<precision) break;
    if (mu_val != mu_val || mu_val > 10e9 || mu_val < -10e9)
    {
      cout << "Error in computing limit. mu = " << mu_val << endl;
      break;
    }
  } while (true && itr < maxItr);

  return mu_val;
}


double calcCLs(double& mu_val, RooRealVar* mu, RooDataSet* data, RooDataSet* asimovData, const char* conditionalSnapshot, const char* nominalSnapshot, RooWorkspace* ws, double& testStat_val, double ncpScale)
{
  bool doQmu_tilde = !(mu->getMin() < 0); // check if we're doing \tilde{qmu}

//dynamically adjust the range to compensate for larger than expected limits
  cout << "mu_val: " << mu_val << endl;
  if (mu->getMax() < mu_val)
  {
    mu->setMax(2*mu_val);
  }
  mu->setVal(mu_val);

  RooArgSet poi(*mu);

  double CLsb = 1;
  double CLb = 0;
  double qmu_tilde;
  double sigma;

  if (data == asimovData) // this must be median expectation
  {
//help out the denominator fit by setting a tighter range.
    double max_val = mu->getMax();
    mu->setMax(min(1.1*mu_val, max_val));

    qmu_tilde=evaluate(asimovData, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot);
    cout << "qmu_A: " << qmu_tilde << endl;

    mu->setMax(max_val);

    sigma = ncpScale*mu_val/sqrt(qmu_tilde);
    testStat_val = qmu_tilde;
  }
  else
  {
    qmu_tilde=evaluate(data, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot);
    cout << "qmu_tilde: " << qmu_tilde << endl;

    if (qmu_tilde < -0.99)
    {
      cout << "qmu tilde gave " << qmu_tilde << ". Retrying with mu*1.01" << endl;

      mu_val *= 1.01;
      mu->setVal(mu_val);
      qmu_tilde=evaluate(data, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot);
    }

    double max_val = mu->getMax();
    mu->setMax(min(1.1*mu_val, max_val));
    double qmu_A = evaluate(asimovData, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot);
    if (qmu_A < 0) // just try again to see if it will work
    {
      qmu_A = evaluate(asimovData, asimovData, poi, ws, conditionalSnapshot, nominalSnapshot);
    }
    mu->setMax(max_val);
    cout << "qmu_A: " << qmu_A << endl;
    
    sigma = ncpScale*mu_val/sqrt(qmu_A);
  }

//compute p-values.
//eq 65/66 in Asimov paper (arxiv 1007.1727)
    CLsb=ROOT::Math::gaussian_cdf(sqrt(qmu_tilde));
  if (!doQmu_tilde || qmu_tilde <= mu_val*mu_val/(sigma*sigma))
  {
    CLsb=ROOT::Math::gaussian_cdf(sqrt(qmu_tilde));
    CLb =ROOT::Math::gaussian_cdf(mu_val / sigma - sqrt(qmu_tilde));  
    cout << "qmu <= mu^2/sigma^2" << endl;
  }
  else
  {
    CLsb=ROOT::Math::gaussian_cdf(0.5*(qmu_tilde+mu_val*mu_val/(sigma*sigma))/(mu_val/sigma));
    CLb =ROOT::Math::gaussian_cdf(0.5*(mu_val*mu_val/(sigma*sigma) - qmu_tilde)/(mu_val/sigma));  
    cout << "qmu > mu^2/sigma^2" << endl;
  }


  if (CLb == 0)
  {
    cout << "ERROR::CLb == 0. Returning 1" << endl;
    return 1;
  }

  cout << "p_s: " << 1-CLsb << endl;
  cout << "p_b: " << CLb << endl;

//return CLs
  return (1-CLsb)/CLb;
}




void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter)
{
  if (counter > 50)
  {
    cout << "ERROR::Couldn't unfold constraints!" << endl;
    cout << "Initial: " << endl;
    initial.Print("v");
    cout << endl;
    cout << "Final: " << endl;
    final.Print("v");
    exit(1);
  }
  TIterator* itr = initial.createIterator();
  RooAbsPdf* pdf;
  while ((pdf = (RooAbsPdf*)itr->Next()))
  {
    RooArgSet nuis_tmp = nuis;
    RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
    //if (constraint_set.getSize() > 1)
    //{
    string className(pdf->ClassName());
    if (className != "RooGaussian" && className != "RooLognormal" && className != "RooGamma" && className != "RooPoisson" && className != "RooBifurGauss")
    {
      counter++;
      unfoldConstraints(constraint_set, final, obs, nuis, counter);
    }
    else
    {
      final.add(*pdf);
    }
  }
  delete itr;
}


void makeAsimovData(ModelConfig* mcInWs, bool doConditional, RooWorkspace* combWS, RooAbsPdf* combPdf, RooDataSet* combData, bool b_only)
{
////////////////////
//make asimov data//
////////////////////

  int _printLevel = 0;
  stringstream muStr;
  muStr << "_" << !b_only;

  RooRealVar* mu = (RooRealVar*)mcInWs->GetParametersOfInterest()->first();//combWS->var("mu");
  mu->setVal(!b_only);

  RooArgSet mc_obs = *mcInWs->GetObservables();
  RooArgSet mc_globs = *mcInWs->GetGlobalObservables();
  RooArgSet mc_nuis = *mcInWs->GetNuisanceParameters();

//pair the nuisance parameter to the global observable
  RooArgSet mc_nuis_tmp = mc_nuis;
  RooArgList nui_list("ordered_nuis");
  RooArgList glob_list("ordered_globs");
  RooArgSet constraint_set_tmp(*combPdf->getAllConstraints(mc_obs, mc_nuis_tmp, false));
  RooArgSet constraint_set;
  int counter_tmp = 0;
  unfoldConstraints(constraint_set_tmp, constraint_set, mc_obs, mc_nuis_tmp, counter_tmp);

  TIterator* cIter = constraint_set.createIterator();
  RooAbsArg* arg;
  while ((arg = (RooAbsArg*)cIter->Next()))
  {
    RooAbsPdf* pdf = (RooAbsPdf*)arg;
    if (!pdf) continue;

    //pdf->Print();

    TIterator* nIter = mc_nuis.createIterator();
    RooRealVar* thisNui = NULL;
    RooAbsArg* nui_arg;
    while ((nui_arg = (RooAbsArg*)nIter->Next()))
    {
      if (pdf->dependsOn(*nui_arg))
      {
	thisNui = (RooRealVar*)nui_arg;
	break;
      }
    }
    delete nIter;


//need this incase the observable isn't fundamental. 
//in this case, see which variable is dependent on the nuisance parameter and use that.
    RooArgSet* components = pdf->getComponents();
    //components->Print();
    components->remove(*pdf);
    if (components->getSize())
    {
      TIterator* itr1 = components->createIterator();
      RooAbsArg* arg1;
      while ((arg1 = (RooAbsArg*)itr1->Next()))
      {
	TIterator* itr2 = components->createIterator();
	RooAbsArg* arg2;
	while ((arg2 = (RooAbsArg*)itr2->Next()))
	{
	  if (arg1 == arg2) continue;
	  if (arg2->dependsOn(*arg1))
	  {
	    components->remove(*arg1);
	  }
	}
	delete itr2;
      }
      delete itr1;
    }
    if (components->getSize() > 1)
    {
      cout << "ERROR::Couldn't isolate proper nuisance parameter" << endl;
      return;
    }
    else if (components->getSize() == 1)
    {
      thisNui = (RooRealVar*)components->first();
    }



    TIterator* gIter = mc_globs.createIterator();
    RooRealVar* thisGlob = NULL;
    RooAbsArg* glob_arg;
    while ((glob_arg = (RooAbsArg*)gIter->Next()))
    {
      if (pdf->dependsOn(*glob_arg))
      {
	thisGlob = (RooRealVar*)glob_arg;
	break;
      }
    }
    delete gIter;

    if (!thisNui || !thisGlob)
    {
      cout << "WARNING::Couldn't find nui or glob for constraint: " << pdf->GetName() << endl;
      //return;
      continue;
    }

    if (_printLevel >= 1) cout << "Pairing nui: " << thisNui->GetName() << ", with glob: " << thisGlob->GetName() << ", from constraint: " << pdf->GetName() << endl;

    nui_list.add(*thisNui);
    glob_list.add(*thisGlob);
    //thisNui->Print();
    //thisGlob->Print();
  }
  delete cIter;




//save the snapshots of nominal parameters
  combWS->saveSnapshot("nominalGlobs",glob_list);
  combWS->saveSnapshot("nominalNuis", nui_list);

  RooArgSet nuiSet_tmp(nui_list);

  mu->setVal(!b_only);
  mu->setConstant(1);

  if (doConditional)
  {
    cout << "Starting minimization.." << endl;
    RooAbsReal* nll = combPdf->createNLL(*combData, RooFit::Constrain(nuiSet_tmp));
    RooMinimizer minim(*nll);
    minim.setStrategy(2);
    minim.setPrintLevel(1);
    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    if (status != 0)
    {
      cout << "Fit failed for mu = " << mu->getVal() << " with status " << status << endl;
    }
    cout << "Done" << endl;
    //combPdf->fitTo(*combData,Hesse(false),Minos(false),PrintLevel(0),Extended(), Constrain(nuiSet_tmp));
  }
  mu->setConstant(0);



//loop over the nui/glob list, grab the corresponding variable from the tmp ws, and set the glob to the value of the nui
  int nrNuis = nui_list.getSize();
  if (nrNuis != glob_list.getSize())
  {
    cout << "ERROR::nui_list.getSize() != glob_list.getSize()!" << endl;
    return;
  }

  for (int i=0;i<nrNuis;i++)
  {
    RooRealVar* nui = (RooRealVar*)nui_list.at(i);
    RooRealVar* glob = (RooRealVar*)glob_list.at(i);

    if (_printLevel >= 1) cout << "nui: " << nui << ", glob: " << glob << endl;
    if (_printLevel >= 1) cout << "Setting glob: " << glob->GetName() << ", which had previous val: " << glob->getVal() << ", to conditional val: " << nui->getVal() << endl;

    glob->setVal(nui->getVal());
  }

//save the snapshots of conditional parameters
  //cout << "Saving conditional snapshots" << endl;
  combWS->saveSnapshot(("conditionalGlobs"+muStr.str()).c_str(),glob_list);
  combWS->saveSnapshot(("conditionalNuis" +muStr.str()).c_str(), nui_list);

  if (!doConditional)
  {
    combWS->loadSnapshot("nominalGlobs");
    combWS->loadSnapshot("nominalNuis");
  }

  //cout << "Making asimov" << endl;
//make the asimov data (snipped from Kyle)
  mu->setVal(!b_only);
  ModelConfig* mc = mcInWs;

  int iFrame=0;

  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  //cout << "adding obs" << endl;
  obsAndWeight.add(*mc->GetObservables());
  //cout << "adding weight" << endl;

  RooRealVar* weightVar = NULL;
  if (!(weightVar = combWS->var(weightName)))
  {
    combWS->import(*(new RooRealVar(weightName, weightName, 1,0,1000)));
    weightVar = combWS->var(weightName);
  }
  if (_printLevel >= 1) cout << "weightVar: " << weightVar << endl;
  obsAndWeight.add(*combWS->var(weightName));



  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////
  // MAKE ASIMOV DATA FOR OBSERVABLES

  // dummy var can just have one bin since it's a dummy
  if(combWS->var("ATLAS_dummyX"))  combWS->var("ATLAS_dummyX")->setBins(1);

  if (_printLevel >= 1) cout <<" check expectedData by category"<<endl;
  //RooDataSet* simData=NULL;
  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc->GetPdf());

  map<string, RooDataSet*> asimovDataMap;

    
  //try fix for sim pdf
  RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();//(RooCategory*)combWS->cat("master_channel");//(RooCategory*) (&simPdf->indexCat());
  //    TIterator* iter = simPdf->indexCat().typeIterator() ;
  TIterator* iter = channelCat->typeIterator() ;
  RooCatType* tt = NULL;
  int nrIndices = 0;
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
      if (thisNorm*expectedEvents <= 0)
      {
	cout << "WARNING::Detected bin with zero expected events! Please check your inputs." << endl;
      }
      if (thisNorm*expectedEvents > pow(10.0, -2) && thisNorm*expectedEvents < pow(10.0, 9)) obsDataUnbinned->add(*mc->GetObservables(), thisNorm*expectedEvents);
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

    //((RooRealVar*)obstmp->first())->Print();
    //cout << "expected events " << pdftmp->expectedEvents(*obstmp) << endl;
     

    asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;

    if (_printLevel >= 1)
    {
      cout << "channel: " << channelCat->getLabel() << ", data: ";
      obsDataUnbinned->Print();
      cout << endl;
    }
  }

  RooDataSet* asimovData = new RooDataSet(("asimovData"+muStr.str()).c_str(),("asimovData"+muStr.str()).c_str(),RooArgSet(obsAndWeight,*channelCat),Index(*channelCat),Import(asimovDataMap),WeightVar(*weightVar));
  combWS->import(*asimovData);


//bring us back to nominal for exporting
  combWS->loadSnapshot("nominalNuis");
  combWS->loadSnapshot("nominalGlobs");
}

void setupTestStat(ProfileLikelihoodTestStat*& testStat, RooAbsPdf* pdf)
{
  testStat = new ProfileLikelihoodTestStat(*pdf);
  testStat->SetOneSided(true);
}

double evaluate(RooDataSet* data, RooDataSet* asimovData, RooArgSet& poi, RooWorkspace* ws, const char* conditionalSnapshot, const char* nominalSnapshot)
{
  ProfileLikelihoodTestStat* testStat = data == asimovData ? asimov_testStat : obs_testStat;
  if (data == asimovData && string(conditionalSnapshot) != "") ws->loadSnapshot(conditionalSnapshot);

  RooRealVar* firstPOI = (RooRealVar*)poi.first();
  double mu_val = firstPOI->getVal();
  double ret = 2*testStat->Evaluate(*data, poi);
  firstPOI->setVal(mu_val); // reset to initial val. unconditional minimization is done last, so it will change the value of mu

  if (data == asimovData && string(nominalSnapshot) != "") ws->loadSnapshot(nominalSnapshot);
  return ret;
}
