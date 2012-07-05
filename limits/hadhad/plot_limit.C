#include "TString.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TFile.h"
#include "TFile.h"
#include "TLine.h"
#include "TH1.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>


void plot_limit(float Lumi = 2052, bool showData = false, bool N95 = true , bool log = true , bool GEV = false, float ymin = 3 , float ymax = 250, int nTheoryPoints = 4 ){

//Enter number of mass points to evaluate and the their acceptance. 
int nMassPoints =31;//0.5 + n*50
//int nMassPoints =25;//0.5 + n*50

  bool doKKgluon=false;
//  const int numTheories = 4; 
  const int numTheories = 1; 
  
  //acceptance for KKgluon
  accEE_KK= new TF1("accKK","pol6",0.0,3000.0);
// non eta corrected acceptance fit
//  accEE_KK->SetParameters(0.136941,-0.00045355,1.75875e-06,-2.33786e-09,1.4147e-12,-4.04901e-16,4.47585e-20);
  // fit parameters with eta correction
  accEE_KK->SetParameters(0.0985968,-0.000383272,1.64764e-06,-2.23113e-09,1.36116e-12,-3.92057e-16,4.36006e-20);

  //acceptance for Zprime
  accEE_ZP= new TF1("accZP","pol6",0.0,3000.0);
//  accEE_ZP->SetParameters(0.18208,-0.000609399,1.6292e-06,-1.84565e-09,1.00645e-12,-2.66438e-16,2.76776e-20);
  accEE_ZP->SetParameters(-0.017067,0.000479189,-6.84504e-07,6.47602e-10,-4.37612e-13,1.62023e-16,-2.34613e-20);
  
  TGraph* graph = new TGraph();  
  TGraph* datagraph = new TGraph(); 

  TF1 *theory_uncertainty = new TF1("theory_uncertainty","pol1",0.07,3.07);
  theory_uncertainty->SetParameters(0.0,0.067);

  //Vector for each mass to store the PEs
  std::vector<double> vecLimit[nMassPoints];
  //vectors with info from PEs
  std::vector<double> Mass;   
  std::vector<double> vec68pos, vec95pos;
  std::vector<double> vec68neg, vec95neg;

  double pelimit = 0;
  double mass = 0;
  double datalimit = 0;
  double datamass = 0;
  double datainvmass = 0;

  TFile* zpsig = new TFile("LimitInputs/ZP/signalTemplate.root","READ");
  TH1D* zpacc = (TH1D*) zpsig->Get("acceptance");
  
  for (int j=0;j<nMassPoints;j++) {
    char buffer [50];
    sprintf (buffer, "%d",j);
    string iMassTxt = buffer;
    if ( j != 0 && j!= 4 && j!= 10) continue;
    if (doKKgluon)
        string inFile = ("zprime_ensembles_mass"+iMassTxt+"_run10_KK_eeSys_pe500.root").c_str();
    else
//        string inFile = ("zprime_ensembles_mass"+iMassTxt+"_run10_ZP_eeSys_pe500.root").c_str();
        string inFile = ("zprime_ensembles_mass"+iMassTxt+"_run99_ZP_eeSys.root").c_str();
    TFile* myfile = new TFile(inFile.c_str(),"READ");
    TTree* tree = (TTree*)myfile->Get("ensemble_test");     

    
    for (int i=0;i<tree->GetEntries();i++) 
    {
      tree->GetEntry(i);
      //double pelimit = tree->GetBranch("par_marg_quantile95_par_2")->GetL.eaf("par_marg_Quantile95_par_2")->GetValue();
      double pelimit = tree->GetBranch("95quantile_marginalized_1")->GetLeaf("marginalized 951antile of par. 23")->GetValue();
      vecLimit[j].push_back(pelimit);  
    }
      //templates come in 50GeV bins (starting at 500 GeV)
      double mass = 0.5+j*0.05;
      if(GEV)mass *= 1000; 
      Mass.push_back(mass);  
  }

  if (doKKgluon)
//    ifstream datafile ("./data_limit_run24_kkgluon_strip.txt");
    ifstream datafile ("./data_limit_run10_kk_strip.txt");
  else
//    ifstream datafile ("./data_limit_nosys_zprime_strip.txt");
    ifstream datafile ("./testcombdatarun_zp_strip.txt");
//    ifstream datafile ("./data_limit_run10_zprime_strip.txt");
    
  
  double xsecratiofactor = 1.0;
 
//get limit from data
//    double kk_eta_acc[39] = {0.855616,0.85939,0.877401,0.884387,0.877306,0.893388,0.902993,0.903293,0.905695,0.909199,0.915584,0.903082,0.919813,0.910765,0.918317,0.917066,0.921695,0.914809,0.924285,0.918122,0.923752,0.92022,0.929707,0.924938,0.926141,0.925097,0.920116,0.921875,0.92233,0.92497,0.926633,0.923063,0.923566,0.922518,0.91852,0.91479,0.920499,0.919744,0.917039};

for (int j=0;j<nMassPoints;j++) { 

    if ( j != 0 && j!= 4 && j!= 10) continue;
    datafile >> datamass >> datalimit;  
    int imass = int(datamass);
    double mass = 0.5+j*0.05;
    if (!N95)
    {
      double acceptance = 1.0;
      if (doKKgluon)
        acceptance = accEE_KK->Eval(mass*1000)*.104;
      else
        acceptance = zpacc->GetBinContent(j+1)*.104;
//        acceptance = accEE_ZP->Eval(mass*1000)*.104;

      cout<<" Mass "<<mass<<" accept "<<acceptance<<endl; 
      xsecratiofactor = 1./(Lumi*acceptance); 
    }
    else 
      xsecratiofactor = 1.; 

    //Display data limit on graph
    if(GEV)mass *= 1000; 
    datagraph->SetPoint(j, mass, xsecratiofactor*datalimit);    
    std::sort(vecLimit[j].begin(), vecLimit[j].end());  
    unsigned neg95 = unsigned(vecLimit[j].size() * (1 - .9545) / 2);
    unsigned neg68 = unsigned(vecLimit[j].size() * (1 - .6827) / 2);      
    unsigned central = unsigned(vecLimit[j].size() * .5);
    unsigned pos68 = unsigned(vecLimit[j].size() * (1 + .6827) / 2);
    unsigned pos95 = unsigned(vecLimit[j].size() * (1 + .9545) / 2);    

    double medianLimit = vecLimit[j][central];
    std::cout << " number of events @95% " << j << " bin : " << central << " : " << medianLimit << " : " << datalimit << std::endl;

    cout<<"Test mass "<<mass<<" Expected limit: "<<xsecratiofactor*medianLimit<<" Observed limit: "<<xsecratiofactor*datalimit<<" [counts] "<<endl;
    cout<<"Test mass "<<mass<<" Expected limit: "<<xsecratiofactor*medianLimit<<" Observed limit: "<<xsecratiofactor*datalimit<<" [pb] "<<endl;
  
    graph->SetPoint(j, Mass[j], xsecratiofactor*vecLimit[j][central]);    
    //put in points for expected limits
    vec68pos.push_back(xsecratiofactor*vecLimit[j][pos68]);
    vec68neg.push_back(xsecratiofactor*vecLimit[j][neg68]);
    vec95pos.push_back(xsecratiofactor*vecLimit[j][pos95]);
    vec95neg.push_back(xsecratiofactor*vecLimit[j][neg95]);  
    
 }

  datafile.close();


//  if (doKKgluon) numTheories = 4;
    

  if (doKKgluon)
//    ifstream smTheoryfile ("./kkgluon_theories.txt");
    ifstream smTheoryfile ("./kkgluon_theories_LO.txt");
  else
    ifstream smTheoryfile ("./zprimexsecLHC_cteq6l1_172.5_100M_20110210+20110505.dat");
  
  TGraph* smTheorygraph[numTheories];
  TGraph* smTheorygraphHigh[numTheories];
  TGraph* smTheorygraphLow[numTheories];
  TGraphErrors* smTheorygrapherrors[numTheories];
  TGraph* logsmTheorygraph[numTheories];
  TGraph* logsmTheorygraphLow[numTheories];

  for (int j=0;j<numTheories;j++) { 
    smTheorygraph[j] = new TGraph();
    smTheorygraphNLO[j] = new TGraph();
    smTheorygraphHigh[j] = new TGraph();
    smTheorygraphLow[j] = new TGraph();
    smTheorygrapherrors[j] = new TGraphErrors();
    logsmTheorygraph[j] = new TGraph();
    logsmTheorygraphNLO[j] = new TGraph();
    logsmTheorygraphLow[j] = new TGraph();
  }
    TGraph* smTheorygraphFill = new TGraph();
    smTheorygraphFill->SetLineWidth(0);
    smTheorygraphFill->SetLineStyle(0);
    smTheorygraphFill->SetLineColor(kWhite);
    smTheorygraphFill->SetFillColor(12);
    double theory[numTheories];

  if (doKKgluon)
  {
    nTheoryPoints = 16;
  }
  else
  {
    nTheoryPoints = 189;
  }
    
  for (int j=0;j<nTheoryPoints;j++) { 
    smTheoryfile >> datainvmass;  
    for (int k=0;k<numTheories;k++) {
      smTheoryfile >> theory[k];  
      if (doKKgluon) 
        theory[k] *= 1000000000; // input mb, plot pb
       
    //  if(k == 1 )theory[k] *= .95; // BR 
    }
    datamass = datainvmass;
    if(!GEV) datamass /= 1000.;
//    if (!doKKgluon) datamass /=1000;
    //datalogmass = log(1000.0/datainvmass)/log(10);
    for (int k=0;k<numTheories;k++) {
       smTheorygraph[k]->SetPoint(j, datamass, theory[k]);
       if (!doKKgluon) {
           smTheorygraphNLO[k]->SetPoint(j, datamass, theory[k]*1.3);
       }
       cout<<"Mass "<<datamass<<" xsec "<< theory[k]<<endl;
       smTheorygraphHigh[k]->SetPoint(j, datamass, (1.+theory_uncertainty->Eval(datamass))*theory[k]);
       smTheorygraphLow[k]->SetPoint(j, datamass, (1.-theory_uncertainty->Eval(datamass))*theory[k]);
       smTheorygrapherrors[k]->SetPoint(j, datamass, theory[k]);
       smTheorygrapherrors[k]->SetPointError(j, 0., 0.07*theory[k]);
       if(k==0){ 
         smTheorygraphFill->SetPoint(j, datamass, (1.+theory_uncertainty->Eval(datamass))*theory[k]);
         smTheorygraphFill->SetPoint(2*nTheoryPoints-j-1, datamass, (1.-theory_uncertainty->Eval(datamass))*theory[k]);
       }   
       logsmTheorygraph[k]->SetPoint(j, datamass, log(theory[k]));
       if (!doKKgluon) {
        logsmTheorygraph[k]->SetPoint(j, datamass, log(theory[k]*1.3));
       }
       logsmTheorygraphLow[k]->SetPoint(j, datamass, log((1.-theory_uncertainty->Eval(datamass))*theory[k]));

    }
 } 
  smTheoryfile.close();
  
 //-----Make the pretty plot-----//

   gDirectory->GetList()->Delete();
//   gROOT->ProcessLine(".x rootlogon_atlas.C");
   gROOT->ProcessLine(".L AtlasStyle.C");
   gROOT->SetStyle("Plain");
   SetAtlasStyle();
//   gROOT->SetStyle("ATLAS");
   gROOT->ForceStyle();

   TCanvas *c1 = new TCanvas("c1","C1",800,600);
   TGraph* hist95 = new TGraph();
   hist95->SetLineWidth(0);
   hist95->SetLineStyle(0);
   hist95->SetLineColor(kWhite);
   hist95->SetFillColor(kYellow);
   gPad->SetBottomMargin(.15);
   for (unsigned ibin = 0; ibin < vec95pos.size(); ++ibin)
   {
/*//-----Artificial 68% and 95% boundaries-----//
   std::cout <<"hey"<< ibin << " " << Mass[ibin] << " " << vec95pos[ibin] << std::endl;
   vec68pos[ibin] *= 1.3;
   vec68neg[ibin] /= 1.3.;
   vec95pos[ibin] *= 1.5;
   vec95neg[ibin] /= 1.5.;
*///---End: Artificial 68% and 95% boundaries---//

     hist95->SetPoint(ibin, Mass[ibin], vec95pos[ibin]);
   }
   unsigned counter95 = hist95->GetN();
   for (int bin = vec95neg.size() - 1; bin >= 0; --bin)
   {
     hist95->SetPoint(counter95++, Mass[bin], vec95neg[bin]);
   }
   hist95->GetYaxis()->SetTitleSize(0.045);
   hist95->GetYaxis()->SetLabelSize(0.045);
   hist95->GetXaxis()->SetTitleSize(0.05);
   hist95->GetXaxis()->SetLabelSize(0.045);
   hist95->GetXaxis()->SetTitleOffset(1.1);
   hist95->GetYaxis()->SetTitleOffset(1.1);
   if (doKKgluon)
     hist95->GetXaxis()->SetRangeUser(0.4,2.);
   else
     hist95->GetXaxis()->SetRangeUser(0.4,1.05);
  
   hist95->Draw("AF");  
   hist95->SetMaximum(ymax);
   hist95->SetMinimum(ymin);  
   hist95->SetTitle("");
   cout << "min and max of plots?" << hist95->GetMaximum() << " : " << hist95->GetMinimum() << endl;
   if (doKKgluon)
     hist95->GetXaxis()->SetTitle("m_{g_{KK}} [TeV]");
   else
     hist95->GetXaxis()->SetTitle("m_{Z'} [TeV]");
   hist95->GetYaxis()->SetTitle("#sigma B [pb]");

   //hist95->GetYaxis()->SetTitle("95% C.L. Limits on #sigma B(g_{KK}#rightarrow t #bar{t}) (pb)");
   if(N95)hist95->GetYaxis()->SetTitle("95% C.L. Limits on N(g_{KK})");
  
   TGraph* hist68 = new TGraph();
   hist68->SetLineWidth(0);
   hist68->SetLineStyle(0);
   hist68->SetLineColor(kWhite);
   hist68->SetFillColor(kGreen);
   for (unsigned gbin = 0; gbin < vec68pos.size(); ++gbin)
   {
      hist68->SetPoint(gbin, Mass[gbin], vec68pos[gbin]);
   }
   unsigned counter68 = hist68->GetN();
   for (int fbin = vec68neg.size() - 1; fbin >= 0; --fbin)
   {
      hist68->SetPoint(counter68++, Mass[fbin], vec68neg[fbin]);
   }
   hist68->Draw("F");

   graph->SetLineStyle(2);
   graph->SetMarkerStyle(28);
   graph->SetMarkerSize(.7);
   graph->Draw("PL"); 
   datagraph->SetMarkerStyle(20);
   datagraph->SetMarkerSize(.7);
   datagraph->SetMarkerColor(2);
   datagraph->SetLineColor(2);
   datagraph->SetLineWidth(2);
   datagraph->SetFillColor(2);
   if(showData)
     datagraph->Draw("samepl");   // UNCOMMENT this line to include datagraph !!!
   
   //int theoryColor[8] = {12,6,11,4,65,40,8};
   //int theoryMarker[8] = {24,25,22,26,23,27,24};
   int theoryColor[8] = {12,6,11,4,12,6,11,4};
   int theoryMarker[8] = {24,25,22,26,24,25,22,26};
   cout << theoryMarker << endl;
   if(!N95)
   for (int k=0;k<numTheories;k++) {
     double xIntersectionExp = findIntersection (graph, smTheorygraph[k], .500, 2.0000);
     double xIntersection = findIntersection (datagraph, smTheorygraph[k], .500, 2.000);
     std::cout << "Expected Mass Limit (TeV) = " << xIntersectionExp << " Observed Mass Limit (TeV) = " << xIntersection << std::endl;
     std::cout << "Expected Cross Section (pb) = " <<  smTheorygraph[k]->Eval(xIntersectionExp) << " Observed cross section (pb) = " << smTheorygraph[k]->Eval(xIntersection) << std::endl;    
     smTheorygraph[k]->SetMarkerStyle(theoryMarker[k]);
     //smTheorygraph[k]->SetMarkerSize(.7);
     smTheorygraph[k]->SetMarkerSize(0);
     smTheorygraph[k]->SetMarkerColor(theoryColor[k]);
     smTheorygraph[k]->SetLineColor(theoryColor[k]);
     smTheorygraph[k]->SetFillColor(theoryColor[k]);
     smTheorygraph[k]->Draw("samepl");
     smTheorygraph[k]->SetLineWidth(2);
    
     smTheorygrapherrors[k]->SetFillColor(kYellow);

     if (!doKKgluon) { 
         double xIntersectionExpNLO = findIntersection (graph, smTheorygraphNLO[k], .500, 2.0000);
         double xIntersectionNLO = findIntersection (datagraph, smTheorygraphNLO[k], .500, 2.000);
         std::cout << "NLO Expected Mass Limit (TeV) = " << xIntersectionExpNLO << " Observed Mass Limit (TeV) = " << xIntersectionNLO << std::endl;
         std::cout << "NLO Expected Cross Section (pb) = " <<  smTheorygraphNLO[k]->Eval(xIntersectionExpNLO) << " Observed cross section (pb) = " << smTheorygraphNLO[k]->Eval(xIntersectionNLO) << std::endl;    
         smTheorygraphNLO[k]->SetMarkerStyle(theoryMarker[k]);
         //smTheorygraph[k]->SetMarkerSize(.7);
         smTheorygraphNLO[k]->SetMarkerSize(0);
         smTheorygraphNLO[k]->SetMarkerColor(theoryColor[k]);
         smTheorygraphNLO[k]->SetLineColor(theoryColor[k]);
         smTheorygraph[k]->SetLineStyle(2);
         smTheorygraphNLO[k]->SetLineStyle(1);
         smTheorygraphNLO[k]->SetFillColor(theoryColor[k]);
         smTheorygraphNLO[k]->Draw("samepl");
         smTheorygraphNLO[k]->SetLineWidth(2);
     }   

   }   

   float legendXmin, legendXmax, legendYmin, legendYmax;
   legendXmin = 0.70; legendXmax = 0.90; legendYmin = 0.735; legendYmax = 0.915;
   TLegend* l = new TLegend(legendXmin, legendYmin, legendXmax, legendYmax);
   l->SetFillColor(10);
   l->SetLineColor(0);
   l->SetTextFont(42);
   l->AddEntry(graph,"Expected limit","lp");
   l->AddEntry(hist68,"Expected #pm 1#sigma","f");
   l->AddEntry(hist95,"Expected #pm 2#sigma","f");
   if(showData)
	l->AddEntry(datagraph,"Observed limit","lp");  // Change to "Observed limit" or "Pseudo-Data"
   l->SetLineColor(0);
   l->SetFillColor(10);
   l->SetBorderSize(0);
   l->SetTextSize(0.04);  

   legendXmin = 0.44; legendXmax = 0.66; legendYmin = 0.74; legendYmax = 0.90;
   TLegend* l1= new TLegend(legendXmin, legendYmin, legendXmax, legendYmax);
   l1->SetFillColor(10);
   l1->SetLineColor(0);
   l1->SetTextFont(42);
   l1->SetLineColor(0);
   l1->SetFillColor(10);
   l1->SetBorderSize(0);
   l1->SetTextSize(0.04);  

   if(!N95)
   for (int k=0;k<numTheories;k++) {
               // if (k==0) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.20","lp");
               // if (k==1) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.25","lp");
               // if (k==2) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.30","lp");
               // if (k==3) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.35","lp");
                if (k==0) 
		   if (doKKgluon) 
		     l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.20","lp");
                   else {
		     l1->AddEntry(smTheorygraph[k],"#sigma_{Z'} B","lp");
		     l1->AddEntry(smTheorygraphNLO[k],"NLO #sigma_{Z'} B","lp");
             }
		if (k==1) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.25","lp");
                if (k==2) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.30","lp");
                if (k==3) l1->AddEntry(smTheorygraph[k],"g_{qqg_{KK}}/g_{s} = -0.35","lp");
    }
     
   l->SetFillColor(10);
   l->SetLineColor(0);
   l->Draw();
   if(!N95){ 
   l1->SetFillColor(10);
   l1->SetLineColor(0);
   l1->Draw();   
   }
   char writetext1[100];
   char writetext2[100];
   char writetext3[100];
   char writetext4[100];
   TLatex *t = new TLatex();
   t->SetNDC(1);
   t->SetTextAlign(13);
   t->SetTextColor(kBlack);
//   sprintf(writetext1,"#font[72]{ATLAS} for Approval");
   sprintf(writetext1,"#font[72]{ATLAS} Internal");
   //     sprintf(writetext1,"#font[72]{ATLAS} Preliminary");
   //     sprintf(writetext1,"");
   if (doKKgluon)
     sprintf(writetext2,"g_{KK}#rightarrow t #bar{t}");
   else
     sprintf(writetext2,"Z'#rightarrow t #bar{t}");
   //     sprintf(writetext2,"Z'/Z* #rightarrow #mu#mu");
   //     sprintf(writetext2,"");
   sprintf(writetext3,"#sqrt{s} = 7 TeV");
   sprintf(writetext4,"#int L dt = 2.05 fb^{-1}");
   double xtext = 0.64;
   double ytext = 0.62;
   t->SetTextSize(0.05);
   t->DrawLatex(0.16      ,0.25      ,writetext1);
   t->SetTextSize(0.045);
   t->DrawLatex(xtext,ytext+0.08,writetext2);
   t->DrawLatex(xtext,ytext,writetext3);
   //t->DrawLatex(xtext-0.05,ytext     ,writetext4);
   t->DrawLatex(0.17,0.36     ,writetext4);



   TLatex *t = new TLatex();
   t->SetTextAlign(23);
   t->SetTextSize(0.035);
   int blackcolor=1;
   t->SetTextColor(blackcolor);
//   sprintf(writetext,"ATLAS preliminary                                          #int #font[52]{L dt} #approx ??? pb^{-1}");
//   sprintf(writetext,"                                                           #int #font[52]{L dt} #approx ??? pb^{-1}");

//   c1->SetGridx();
//   c1->SetGridy();
   c1->SetLogy();
 //  hist95->GetXaxis()->Seter(hist95->GetXaxis()->GetXmin(),hist95->GetXaxis()->GetXmax()*1.1);
   if (doKKgluon)
     hist95->GetXaxis()->SetRangeUser(0.41,1.985);
   else
     hist95->GetXaxis()->SetRangeUser(0.5,1.0);

   c1->RedrawAxis();
   c1->Update(); 
   c1->WaitPrimitive(); 
   if((!N95)&&log)c1->Print("Logmasslimit_zprimexsec_zoom.png");        
   if(N95&&log)c1->Print("Logmasslimit_N95_zoom.png");        
   if((!N95)&&(!log))c1->Print("masslimit_zprimexsec_zoom.png");       
   if(N95&&!log)c1->Print("masslimit_N95_zoom.png");        
   if((!N95)&&log)c1->Print("Logmasslimit_zprimexsec_zoom.eps");        
   if(N95&&log)c1->Print("Logmasslimit_N95_zoom.eps");        
   if((!N95)&&(!log))c1->Print("masslimit_zprimexsec_zoom.eps");       
   if(N95&&!log)c1->Print("masslimit_N95_zoom.eps");        

   gROOT->ProcessLine(".q") ;
}

Double_t findIntersection ( TGraph* graphOfData,  TGraph* graphOfTheory, Double_t xMin, Double_t xMax) {
  gDirectory->GetList()->ls();
  gDirectory->GetList()->Delete();
  double x = xMin;
  double intersection = -1;
  const int nSteps = 15000;
  double deltaX = (xMax-xMin)/nSteps;
  double deltaOld = 10;

  foundint = false;

  for (int i=0;i<nSteps;i++) {
    x += deltaX;
    double data = graphOfData->Eval(x);
    double theory = graphOfTheory->Eval(x);
    double delta = theory-data;
    if (delta<0. && !foundint ) {
      intersection = x;
      deltaOld = delta;
      foundint = true;
//    double delta = fabs(data-theory)/fabs(data);
//    if (delta<deltaOld) {
//      intersection = x;
//      deltaOld = delta;
    }
  }

  return intersection;
}



