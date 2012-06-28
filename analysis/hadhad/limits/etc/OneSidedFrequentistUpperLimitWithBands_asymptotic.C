/*
ProfileBands

Author: Kyle Cranmer
date: Jan. 2010

This is a standard demo that can be used with any ROOT file 
prepared in the standard way.  It runs the profile likelihood calculator
and generates background only pseudo-experiments to make bands
for the upper bound.  

Important! Currently the profile likelihood calculator is configured 
to give two-sided intervals, not one-sided upper-limits.  

You specify:
 - name for input ROOT file
 - name of workspace inside ROOT file that holds model and data
 - name of ModelConfig that specifies details for calculator tools
 - name of dataset 

With default parameters the macro will attempt to run the 
standard hist2workspace example and read the ROOT file
that it produces.

The actual heart of the demo is only about 10 lines long.

The ProfileLikelihoodCalculator is based on Wilks's theorem
and the asymptotic properties of the profile likeihood ratio
(eg. that it is chi-square distributed for the true value).
*/

#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooNLLVar.h"
#include "RooProfileLL.h"

#include "RooStats/RooStatsUtils.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ToyMCSampler.h"

using namespace RooFit;
using namespace RooStats;

double getEquivalentCL(int choice, RooAbsData* data, RooAbsPdf* pdf, 
		       RooRealVar* firstPOI, double confidenceLevel = 0.95){

  // This standard setup of ProfileLikleihoodCalculator 
  // will give a two-sided interval without regard to boundaries
  // on the POI.  This corresponds to a threshold on log(lambda) of
  // root [1] ROOT::Math::chisquared_quantile(.95,1)/2
  // (double)1.92072941034706179e+00

  // One may want to have one-sided upper limits instead, and
  // the test statistic will be slightly modified.
  // There are corresponding asymptotic distributions for both.
  // USE:
  //  - t_mu for 2-sided 'unified' interval
  //  - q_mu for 1-sided upper limit
  //
  //  If you have a boundary POI>=0, you can either:
  //   - artificially let POI<0 and formulae will work
  //   - or force POI>=0 and use tilde formulae

  // These different setups should all be done automatically
  // by the next version of ProfileLikleihoodCalculator
  // but here we do it manually.

  // This is the switch for below
  //  choice = 1; // t_mu:       2-sided upper limit, let POI<0
  //  choice = 2; // tilde t_mu: 2-sided upper limit, force POI>=0
  //  choice = 3; // q_mu:       1-sided upper limit, let POI<0
  //  choice = 4; // tilde q_mu: 1-sided upper limit, force POI>=0

  // We will calculate an equivalent CL for each case.
  // By this we mean, what artificial CL is needed 
  // to force the threshold to be what we want it to be.  

  // the current default is t_mu, so it is equivalent
  double tmuEquivalentCL = confidenceLevel;

  // This standard setup will give 2-sided limits
  // corresponds to a threshold at 
  // root [1] ROOT::Math::chisquared_quantile(.95,1)/2
  // (double)1.92072941034706179e+00
  // however if POI is constrained to be >0
  // this will not be completely accurate, would 
  // need to go to \tilde{t}_mu, but that is not
  // yet done.
  double tildetmuEquivalentCL = -1; //NOT YET DONE

  // To modify configuration of calculator for q_mu test statistic
  // and one-sided upper limit, need to let POI<0    
  // and artificially set confidence level to 90% to have
  // equivalent threshold for 1-sided 95% confidence level
  // With this setup, the -2-sigma band should be 0.
  // This one we can do in our heads, but let's do it by
  // inverting eqtn 59 to solve for qmu
  // qmu = (\Phi^{-1}(1-p))^2
  double qmuthreshold = pow(ROOT::Math::gaussian_quantile(confidenceLevel,1),2);
  // gives theshold = 2.7
  double qmuEquivalentCL = ROOT::Math::chisquared_cdf(qmuthreshold,1);
  // corresponding CL for the t_mu test is 90%


  if(choice==1) return tmuEquivalentCL;
  if(choice==2) return tildetmuEquivalentCL;
  if(choice==3) return qmuEquivalentCL;

  // To modify configuration of calculator for \tilede{q}_mu test statistic
  // and one-sided upper limit with constraint POI>=0    
  // We artificially set confidence level again, but the value
  // requires some more work (it's not 90%)
  // Now we need to invert eqtn 66 for tildeqmu.
  /* This is summary of equation 66
    if(0<tildeqmu<mu^2/sigma^2) // first gaussian
    F(tildeqmu|mu) = gaussian_cdf(sqrt(qmutilde));
    if(qmutilde>mu^2/sigma^2) // second gaussian
    F(tildeqmu|mu) = gaussian_cdf((qmutilde + mu^2/sigma^2)/(2*mu/sigma));
    //    invert to find  threshold;
    plTemp.SetConfidenceLevel(ROOT::Math::chisquared_cdf(threshold,1)); 
  */  
  // Here we need sigma.  In general, it depends on mu,
  // and we don't know what mu corresponds to mu95.
  // the problem requires a loop in general.
  // We can approximate by using sigma at muhat for this dataset
  // a better approach would be set mu= first estimate of mu95
  // then generate asimov data, and fit it
    
  // set POI to best fit
  firstPOI->setMin(-1.*firstPOI->getMax()); // let it go negative just for the error
  pdf->fitTo(*data);
  double muHat = firstPOI->getVal(); 
  if (muHat<0) cout << "NEGATIVE"<<endl;
  double sigmaAtMuHat=firstPOI->getError();
  // now for the mu
  double muOverSigma = muHat/sigmaAtMuHat;
  double muOverSigmaSq = pow(muOverSigma,2);
  //  cout << "sigma(mu=muhat) = " << sigmaAtMuHat << endl;
  //  cout << "muhat/sigma(mu=muhat) = " << muOverSigma << endl;


  double tildeqmuThreshTemp=-1;
  double tildeqmuEquivalentCL = -1;
  if(muHat < 0){
    cout << "muHat<0 without boundary "<<endl;
    

    RooNLLVar* nll = (RooNLLVar*) pdf->createNLL(*data);
    RooProfileLL* profile = (RooProfileLL*) nll->createProfile(*firstPOI);
    firstPOI->setVal(0);
    double deltaLR = profile->getVal();

    //    cout << "deltaLR = "<< deltaLR << " muHatOverSigmaSq = " << muOverSigmaSq << endl;
    tildeqmuThreshTemp=pow(ROOT::Math::gaussian_quantile(confidenceLevel,1),2);
    tildeqmuThreshTemp -= 2*deltaLR;
    // I guess we need profileLR(0), not approx muOverSigmaSq 
    tildeqmuEquivalentCL = ROOT::Math::chisquared_cdf(tildeqmuThreshTemp,1);
    delete profile;
    delete nll;
    if(tildeqmuThreshTemp<0)
      cout << "without boundary, upper limit<0 "<<endl;
    
  } else {
    tildeqmuThreshTemp =pow(ROOT::Math::gaussian_quantile(confidenceLevel,1),2);
    // gives theshold = 2.7
    tildeqmuEquivalentCL = ROOT::Math::chisquared_cdf(qmuthreshold,1);
    // corresponding CL for the t_mu test is 90%
  }

  firstPOI->setMin(0); // set it back to being bounded
  return tildeqmuEquivalentCL;


  // Notes to self, in more direct attempt to translate
  // Asimov paper to RooStast
  /*
  double tempThresh = ROOT::Math::gaussian_quantile(0.95,1);
  double tildeqmuThreshTemp = tempThresh*muOverSigma*2-muOverSigmaSq;

  //  double tildeqmuEquivalentCL = -1;
  if(tildeqmuThreshTemp> muOverSigmaSq){
    cout <<"tilde q_mu threshold in steeply falling tail due to boundary"<<endl;
    tildeqmuEquivalentCL = ROOT::Math::chisquared_cdf(tildeqmuThreshTemp,1);    
    // corresponding CL depends on sigma
    // sometimes CL~0 b/c mu near upper boundary    

  } else{
    cout <<"tilde q_mu threshold in main gaussian"<<endl;
    tildeqmuThreshTemp =pow(ROOT::Math::gaussian_quantile(confidenceLevel,1),2);
    // gives theshold = 2.7
    tildeqmuEquivalentCL = ROOT::Math::chisquared_cdf(qmuthreshold,1);
    // corresponding CL for the t_mu test is 90%

    
    // This turns out to be equivalent b/c prob in two tails have to be =
    // Here's the algorithm in words:
    // get tildeqmu quantile from second gaussian
    // check if it is < mu^2/sigma^2
    // if so get CDF up to boundary
    // then find cdf from second gaussian corresponding to boundary
    // then find cdf from first gaussian corresponding to boundary
    // find how much more prob we need (eg. CL-prob in second gaussian)
    // and add to that probability in first gaussian that is chopped by boundary
    // then find quantile of the first gaussian for that prob
    // and finally convert that threshold into the tmu CL

    
    // at boundary tildeqmu=muOverSigmaSq
    // and the argument to the second gaussian is
    // (tildeqmu+muOverSigmaSq)/2/muOverSigma = muOverSigma
    //tempThresh = muOverSigma; 
    //    double pInSecondTail = ROOT::Math::gaussian_cdf_c(tempThresh,1);
    //    double pInFirstTail = ROOT::Math::gaussian_cdf_c(muOverSigma,1);
    //    cout << " pInSecondTail " << pInSecondTail << " pInFirstTail " 
    //    << pInFirstTail << " additional prob needed = " 
    //    << (1-confidenceLevel-pInSecondTail) << endl;
    //    pInFirstTail += (1-confidenceLevel-pInSecondTail);
    //    tildeqmuThreshTemp = pow(ROOT::Math::gaussian_quantile_c(pInFirstTail,1),2);
    //    tildeqmuEquivalentCL = ROOT::Math::chisquared_cdf(tildeqmuThreshTemp,1);    
    //    cout << " boundary " << muOverSigma << " new tail " 
    //    << tildeqmuThreshTemp << " and equiv CL " << tildeqmuEquivalentCL << endl;
  }    
  if(muHat < 0){
    cout << "strange case.  part of CL is just going frm negative to 0, rest is penalty"<<endl;
    

    RooNLLVar* nll = (RooNLLVar*) pdf->createNLL(*data);
    RooProfileLL* profile = (RooProfileLL*) nll->createProfile(*firstPOI);
    firstPOI->setVal(0);
    double deltaLR = profile->getVal();

    cout << "deltaLR = "<< deltaLR << " muHatOverSigmaSq = " << muOverSigmaSq << endl;
    tildeqmuThreshTemp=pow(ROOT::Math::gaussian_quantile(confidenceLevel,1),2);
    tildeqmuThreshTemp -= 2*deltaLR;
    // I guess we need profileLR(0), not approx muOverSigmaSq 
    tildeqmuEquivalentCL = ROOT::Math::chisquared_cdf(tildeqmuThreshTemp,1);
    delete profile;
    delete nll;
    
  }
  
  firstPOI->setMin(0); // set it back to being bounded
  
  double ret = 0;
  if(choice==4) ret = tildeqmuEquivalentCL;
  return ret;
    */
}		 

void OneSidedFrequentistUpperLimitWithBands_asymptotic(const char* infile = "",
		      const char* workspaceName = "combined",
		      const char* modelConfigName = "ModelConfig",
		      const char* dataName = "obsData"){

  /////////////////////////////////////////////////////////////
  // First part is just to access a user-defined file 
  // or create the standard example file if it doesn't exist
  ////////////////////////////////////////////////////////////
  const char* filename = "";
  if (!strcmp(infile,""))
    filename = "results/example_combined_GaussExample_model.root";
  else
    filename = infile;
  // Check if example input file exists
  TFile *file = TFile::Open(filename);

  // if input file was specified byt not found, quit
  if(!file && strcmp(infile,"")){
    cout <<"file not found" << endl;
    return;
  } 

  // if default file not found, try to create it
  if(!file ){
    // Normally this would be run on the command line
    cout <<"will run standard hist2workspace example"<<endl;
    gROOT->ProcessLine(".! prepareHistFactory .");
    gROOT->ProcessLine(".! hist2workspace config/example.xml");
    cout <<"\n\n---------------------"<<endl;
    cout <<"Done creating example input"<<endl;
    cout <<"---------------------\n\n"<<endl;
  }

  // now try to access the file again
  file = TFile::Open(filename);
  if(!file){
    // if it is still not there, then we can't continue
    cout << "Not able to run hist2workspace to create example input" <<endl;
    return;
  }

  
  /////////////////////////////////////////////////////////////
  // Tutorial starts here
  ////////////////////////////////////////////////////////////

  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(workspaceName);
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }

  // get the modelConfig out of the file
  ModelConfig* mc = (ModelConfig*) w->obj(modelConfigName);

  // get the modelConfig out of the file
  RooAbsData* data = w->data(dataName);

  // make sure ingredients are found
  if(!data || !mc){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }

  RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();

  ////////////////////////////////////////////////////
  double confidenceLevel = 0.95;
  int choice = 4;
  //  choice = 1; // t_mu:       2-sided upper limit, let POI<0
  //  choice = 2; // tilde t_mu: 2-sided upper limit, force POI>=0
  //  choice = 3; // q_mu:       1-sided upper limit, let POI<0
  //  choice = 4; // tilde q_mu: 1-sided upper limit, force POI>=0
  double equivalentCL = getEquivalentCL(choice,data,mc->GetPdf(),firstPOI,confidenceLevel);
  cout << "getEquivalentCL = " << equivalentCL << endl;

  //  return;
  /////////////////////////////////////////////
  // create and use the ProfileLikelihoodCalculator
  // to find and plot the 95% confidence interval
  // on the parameter of interest as specified
  // in the model config
  ProfileLikelihoodCalculator pl(*data,*mc);

  if(choice == 1 ){
    // for t_mu
    firstPOI->setMin(-1*firstPOI->getMax());
    pl.SetConfidenceLevel(equivalentCL); 
  } else if( choice==2){
    // for tilde{t}_mu
    cout << "NOT YET IMPLEMENTED, will do tmu for now"<<endl;
    firstPOI->setMin(0);
    pl.SetConfidenceLevel(equivalentCL); 
  } else if( choice==3){
    // for q_mu
    firstPOI->setMin(-1*firstPOI->getMax());
    pl.SetConfidenceLevel(equivalentCL); 
  } else if( choice==4){
    // for tilde q_mu
    firstPOI->setMin(0);
    pl.SetConfidenceLevel(equivalentCL); 
  } else {
    cout << "make choice for test statistic" << endl;
  }

  LikelihoodInterval* interval = pl.GetInterval();

  // make a plot
  LikelihoodIntervalPlot plot(interval);
  plot.Draw();
  
  // print out the iterval on the first Parameter of Interest
  if(choice==1 || choice==2){
    cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
      interval->LowerLimit(*firstPOI) << ", "<<
      interval->UpperLimit(*firstPOI) <<"] "<<endl;
  } else if(choice==3 || choice==4){
    cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
      0 << ", "<<
      interval->UpperLimit(*firstPOI) <<"] "<<endl;
  } else {
    cout << "make choice for test statistic" << endl;
  }

  ////////////////////////////////////////////
  // Now find parameter point for sigma=0, 
  // with conditional MLEs for nuisance parameters
  RooAbsReal* nll = mc->GetPdf()->createNLL(*data);
  RooAbsReal* profile = nll->createProfile(*mc->GetParametersOfInterest());
  firstPOI->setVal(0.);
  profile->getVal(); // this will do fit and set nuisance parameters to profiled values
  RooArgSet* poiAndNuisance = new RooArgSet();
  poiAndNuisance->add(*mc->GetNuisanceParameters());
  poiAndNuisance->add(*mc->GetParametersOfInterest());
  w->saveSnapshot("paramsToGenerateData",*poiAndNuisance);
  RooArgSet* paramsToGenerateData = (RooArgSet*) poiAndNuisance->snapshot();
  cout << "\nWill use these parameter points to generate pseudo data for bkg only" << endl;
  paramsToGenerateData->Print("v");
  delete profile;
  delete nll;

  //  return;
  //////////////////////////////////////////////////////
  // To get sigma for qmutilde we need asimov
  

  //////////////////////////////////////////////////
  // Now do loop over background only pseudo-experiments
  int nToyMC=500;
  TH1F* histOfUL = new TH1F("histOfUL","",100,firstPOI->getMin(),firstPOI->getMax()+.2);
  for(int imc=0; imc<nToyMC; ++imc){

    // set parameters back to values for generating pseudo data
    //    cout << "\n get current nuis, set vals, print again" << endl;
    //poiAndNuisance->Print("v");
    w->loadSnapshot("paramsToGenerateData");
    // if you want to check the parameters used for generation
    //    poiAndNuisance->Print("v");

    // now generate a toy dataset
    //    RooDataSet* toyData = mc->GetPdf()->generate(*mc->GetObservables(),1);
    // if you want to check the toy data


    // in 5.30 there is a nicer way to generate toy data  & randomize global obs
    ToyMCSampler* toymcsampler = new ToyMCSampler();
    toymcsampler->SetPdf(*mc->GetPdf());
    toymcsampler->SetObservables(*mc->GetObservables());
    toymcsampler->SetGlobalObservables(*mc->GetGlobalObservables());
    RooAbsData* toyData = toymcsampler->GenerateToyData(*paramsToGenerateData);

    //    cout << "toy data = " << endl;
    //    toyData->Print("v");
    //    mc->GetObservables()->Print();
    //    mc->GetGlobalObservables()->Print("v");

    double thisUL = 0;

    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING); // lower message level

    ProfileLikelihoodCalculator plTemp(*toyData,*mc);

    equivalentCL = getEquivalentCL(choice,toyData,mc->GetPdf(),firstPOI,confidenceLevel);
    cout << "getEquivalentCL = " << equivalentCL << endl;
    
    // for testing
    //    if(equivalentCL < 0.9)
    //      continue;

    if(choice == 1 ){
      // for t_mu
      firstPOI->setMin(-1*firstPOI->getMax());
      plTemp.SetConfidenceLevel(equivalentCL); 
    } else if( choice==2){
      // for tilde{t}_mu
      cout << "NOT YET IMPLEMENTED, will do tmu for now"<<endl;
      firstPOI->setMin(0);
      plTemp.SetConfidenceLevel(equivalentCL); 
    } else if( choice==3){
      // for q_mu
      firstPOI->setMin(-1*firstPOI->getMax());
      plTemp.SetConfidenceLevel(equivalentCL); 
    } else if( choice==4){
      // for tilde q_mu
      firstPOI->setMin(0);
      plTemp.SetConfidenceLevel(equivalentCL); 
    } else {
      cout << "make choice for test statistic" << endl;
    }
    

    LikelihoodInterval* intervalTemp = plTemp.GetInterval();
    
    thisUL = intervalTemp->UpperLimit(*firstPOI);

    // special case
    if(equivalentCL ==0)
      thisUL = 0;

    cout << " thisUL = " << thisUL << endl;
    histOfUL->Fill(thisUL);

    /*
      // to illustraite what happens if muhat=0 (or <0 w/o constraint)
    if(equivalentCL <0.9){
      LikelihoodIntervalPlot* plotTemp = new LikelihoodIntervalPlot(intervalTemp);
      plotTemp->Draw("same");
    }
    */
    
    delete intervalTemp;
    delete toyData;
    //    delete plTemp;
  }
  histOfUL->Draw();
  Double_t* bins = histOfUL->GetIntegral();
  TH1F* cumulative = (TH1F*) histOfUL->Clone("cumulative");
  cumulative->SetContent(bins);
  double band2sigDown, band1sigDown, bandMedian, band1sigUp,band2sigUp;
  for(int i=1; i<=cumulative->GetNbinsX(); ++i){
    if(bins[i]<RooStats::SignificanceToPValue(2))
      band2sigDown=cumulative->GetBinCenter(i);
    if(bins[i]<RooStats::SignificanceToPValue(1))
      band1sigDown=cumulative->GetBinCenter(i);
    if(bins[i]<0.5)
      bandMedian=cumulative->GetBinCenter(i);
    if(bins[i]<RooStats::SignificanceToPValue(-1))
      band1sigUp=cumulative->GetBinCenter(i);
    if(bins[i]<RooStats::SignificanceToPValue(-2))
      band2sigUp=cumulative->GetBinCenter(i);
  }
  cout << "-2 sigma  band " << band2sigDown << endl;
  cout << "-1 sigma  band " << band1sigDown << endl;
  cout << "median of band " << bandMedian << endl;
  cout << "+1 sigma  band " << band1sigUp << endl;
  cout << "+2 sigma  band " << band2sigUp << endl;

    
}
