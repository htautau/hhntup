
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TMath.h"

#include "RooWorkspace.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"

#include "SmoothHistogram.C"

using namespace RooFit;


// Append a weight var (with fixed value) to a RooDataSet
void AddWeight(RooRealVar* WeightVar, RooDataSet* Data ) {

  int NumRows = Data->numEntries();
  RooDataSet tmpWeightData("tmpWeightData", "", RooArgSet(*WeightVar) );
  
  for(int i = 0; i < NumRows; ++i) {
    RooArgSet Args(*WeightVar);
    tmpWeightData.add( Args );
  }
  
  Data->merge( &tmpWeightData );

}

void HistSmoothExample() {


  // Make the Example Variables
  RooWorkspace w("w");
  RooRealVar* VarA = (RooRealVar*) w.factory("VarA[5,0,10]");
  RooRealVar* VarB = (RooRealVar*) w.factory("VarB[5,-10,10]");
  RooRealVar* VarC = (RooRealVar*) w.factory("VarC[5,-10,10]");
  RooRealVar* VarD = (RooRealVar*) w.factory("VarD[5,0,60]");
  
  RooRealVar* Weight = (RooRealVar*) w.factory("Weight[1,0,10]");
  Weight->setConstant(true);

  // Make the "input" histograms

  TH1F* HistA = new TH1F("HistA", "HistA", 10,  0, 10);
  TH1F* HistB = new TH1F("HistB", "HistB", 50,-10, 10);
  TH1F* HistC = new TH1F("HistC", "HistC", 50,-10, 10);
  TH1F* HistD = new TH1F("HistD", "HistD", 50,  0, 60);


  // Fill the Histograms
  // Create Datasets in the process

  std::cout << "Creating PDF A" << std::endl;

  // Hist A is random:
  int NumEntriesA = 50;
  RooPolynomial PdfA("PdfA", "flat", *VarA, RooArgList() );
  w.import(PdfA, RecycleConflictNodes() );
  Weight->setVal(1.0);

  RooDataSet* DataA = w.pdf("PdfA")->generate( RooArgSet( *VarA), NumEntriesA, Name("DataA") );
  for(int i = 0; i < NumEntriesA; ++i) {
    const RooArgSet* Row = DataA->get(i);
    double RandVal = Row->getRealValue("VarA");
    HistA->Fill(RandVal, Weight->getVal() );
  }

  std::cout << "Creating PDF B" << std::endl;

  // Hist B is gaussian:
  int NumEntriesB = 200;

  w.factory("RooGaussian::PdfB(VarB,0,1)");
  Weight->setVal(1.0);
  RooDataSet* DataB = w.pdf("PdfB")->generate( RooArgSet( *VarB), NumEntriesB, Name("DataB") );
  for(int i = 0; i < NumEntriesB; ++i) {
    const RooArgSet* Row = DataB->get(i);
    double RandVal = Row->getRealValue("VarB");
    HistB->Fill(RandVal, Weight->getVal() );
  }

  std::cout << "Creating PDF C" << std::endl;

  // Hist C is two gaussians with different weights
  int NumEntriesC_First = 100;
  int NumEntriesC_Second = 100;
  double WeightC_First = 1;
  double WeightC_Second = 2;
  // Require that   


  std::cout << "Creating individual gaussians" << std::endl;

  w.factory("RooGaussian::GaussC_0(VarC, 3, 1)");
  w.factory("RooGaussian::GaussC_1(VarC, -3, 1)");
  w.factory("SUM::PdfC(1*GaussC_0,2*GaussC_1)");

  // Make the individual Datasets:

  Weight->setVal( WeightC_First );
  RooDataSet* DataC_0 = (RooDataSet*) w.pdf("GaussC_0")->generate(RooArgSet( *VarC ), NumEntriesC_First, Name("DataC_0"), Verbose(true) );
  AddWeight( Weight, DataC_0 );
  for(int i = 0; i < NumEntriesC_First; ++i) {
    const RooArgSet* Row = DataC_0->get(i);
    //Row->Print("V");
    double RandVal = Row->getRealValue("VarC");
    HistC->Fill(RandVal, Weight->getVal() );
  }
  Weight->setVal( WeightC_Second );
  RooDataSet* DataC_1 = (RooDataSet*) w.pdf("GaussC_1")->generate(RooArgSet( *VarC ), NumEntriesC_Second, Name("DataC_1"), Verbose(true)  );
  AddWeight( Weight, DataC_1 );
  for(int i = 0; i < NumEntriesC_Second; ++i) {
    const RooArgSet* Row = DataC_1->get(i);
    //Row->Print("V");
    double RandVal = Row->getRealValue("VarC");
    HistC->Fill(RandVal, Weight->getVal() );
  }

  std::cout << "Combining datasets" << std::endl;

  RooDataSet* DataC = new RooDataSet("DataC", "DataC", RooArgSet(*VarC, *Weight), "Weight" );

  Weight->setVal( WeightC_First );
  for(int i = 0; i < DataC_0->numEntries(); ++i) {
    const RooArgSet* Row = DataC_0->get(i);
    DataC->add(*Row, WeightC_First);
  }
  Weight->setVal( WeightC_Second );
  for(int i = 0; i < DataC_1->numEntries(); ++i) {
    const RooArgSet* Row = DataC_1->get(i);
    DataC->add(*Row, WeightC_Second);
  }

  // DataC->append( *DataC_0 );
  // DataC->append( *DataC_1 );
  //DataC->Print("V");

  std::cout << "Creating PDF D" << std::endl;

  // Hist D is two gaussians with different weights
  int NumEntriesD_First  = 100;
  int NumEntriesD_Second = 100;
  double WeightD_First = 3;
  double WeightD_Second = 1;

  w.factory("RooExponential::ExpD(VarD, -.1)");
  w.factory("RooGaussian::GaussD(VarD, 20, 3)");
  w.factory("SUM::PdfD(3*ExpD, 1*GaussD)");

  Weight->setVal( WeightD_First );
  RooDataSet* DataD_0 = (RooDataSet*) w.pdf("ExpD")->generate(RooArgSet( *VarD ), NumEntriesD_First, Name("DataD_0"), Verbose(true) );
  AddWeight( Weight, DataD_0 );
  for(int i = 0; i < NumEntriesD_First; ++i) {
    const RooArgSet* Row = DataD_0->get(i);
    //Row->Print("V");
    double RandVal = Row->getRealValue("VarD");
    HistD->Fill(RandVal, Weight->getVal() );
  }
  Weight->setVal( WeightD_Second );
  RooDataSet* DataD_1 = (RooDataSet*) w.pdf("GaussD")->generate(RooArgSet( *VarD ), NumEntriesD_Second, Name("DataD_1"), Verbose(true)  );
  AddWeight( Weight, DataD_1 );
  for(int i = 0; i < NumEntriesD_Second; ++i) {
    const RooArgSet* Row = DataD_1->get(i);
    //Row->Print("V");
    double RandVal = Row->getRealValue("VarD");
    HistD->Fill(RandVal, Weight->getVal() );
  }


  std::cout << "Combining datasets" << std::endl;

  RooDataSet* DataD = new RooDataSet("DataD", "DataD", RooArgSet(*VarD, *Weight), "Weight" );

  Weight->setVal( WeightD_First );
  for(int i = 0; i < DataD_0->numEntries(); ++i) {
    const RooArgSet* Row = DataD_0->get(i);
    DataD->add(*Row, WeightD_First);
  }
  Weight->setVal( WeightD_Second );
  for(int i = 0; i < DataD_1->numEntries(); ++i) {
    const RooArgSet* Row = DataD_1->get(i);
    DataD->add(*Row, WeightD_Second);
  }


  // Now make the smoothed histograms
  // Get the KEYS pdfs used in making
  // those histograms (via argument)

  RooKeysPdf* KeysA_Rho = NULL;
  RooKeysPdf* KeysB_Rho = NULL;
  RooKeysPdf* KeysC_Rho = NULL;
  RooKeysPdf* KeysD_Rho = NULL;

  RooKeysPdf* KeysA_no = NULL;
  RooKeysPdf* KeysB_no = NULL;
  RooKeysPdf* KeysC_no = NULL;
  RooKeysPdf* KeysD_no = NULL;

  RooKeysPdf::Mirror MirrorStyle = RooKeysPdf::MirrorLeft;

  bool VERBOSE = true;

  TH1F* HistA_sm = (TH1F*) SmoothHistogram( HistA, "HistA_sw", VarA, KeysA_Rho, true, MirrorStyle, VERBOSE );
  TH1F* HistB_sm = (TH1F*) SmoothHistogram( HistB, "HistB_sw", VarB, KeysB_Rho, true, MirrorStyle, VERBOSE );
  TH1F* HistC_sm = (TH1F*) SmoothHistogram( HistC, "HistC_sw", VarC, KeysC_Rho, true, MirrorStyle, VERBOSE );
  TH1F* HistD_sm = (TH1F*) SmoothHistogram( HistD, "HistD_sw", VarD, KeysD_Rho, true, MirrorStyle, VERBOSE );

  TH1F* HistA_sm_no = (TH1F*) SmoothHistogram( HistA, "HistA_sw_no", VarA, KeysA_no, false, MirrorStyle, VERBOSE );
  TH1F* HistB_sm_no = (TH1F*) SmoothHistogram( HistB, "HistB_sw_no", VarB, KeysB_no, false, MirrorStyle, VERBOSE );
  TH1F* HistC_sm_no = (TH1F*) SmoothHistogram( HistC, "HistC_sw_no", VarC, KeysC_no, false, MirrorStyle, VERBOSE );
  TH1F* HistD_sm_no = (TH1F*) SmoothHistogram( HistD, "HistD_sw_no", VarD, KeysD_no, false, MirrorStyle, VERBOSE );


  int smColor = 9; //38;
  int sm_noColor = 46;

  HistA_sm->SetLineColor(smColor);
  HistB_sm->SetLineColor(smColor);
  HistC_sm->SetLineColor(smColor);
  HistD_sm->SetLineColor(smColor);

  HistA_sm_no->SetLineColor(sm_noColor);
  HistB_sm_no->SetLineColor(sm_noColor);
  HistC_sm_no->SetLineColor(sm_noColor);
  HistD_sm_no->SetLineColor(sm_noColor);

  // Generate KEYS functions on the unbinned data for comparison
  // (Doesn't use rho)

  RooKeysPdf* KeysA_unbinned = new RooKeysPdf("KeysA_unbinned","", *VarA, *DataA, RooKeysPdf::MirrorLeft);
  RooKeysPdf* KeysB_unbinned = new RooKeysPdf("KeysB_unbinned","", *VarB, *DataB, RooKeysPdf::MirrorLeft);
  RooKeysPdf* KeysC_unbinned = new RooKeysPdf("KeysC_unbinned","", *VarC, *DataC, RooKeysPdf::MirrorLeft);
  RooKeysPdf* KeysD_unbinned = new RooKeysPdf("KeysD_unbinned","", *VarD, *DataD, RooKeysPdf::MirrorLeft);


  // Now draw the histograms

  TCanvas* Canvas = new TCanvas("Canvas","");
  Canvas->Divide(4,5);
  int CurrentPad = 0;


  RooPlot* PlotA = VarA->frame();
  RooPlot* PlotB = VarB->frame();
  RooPlot* PlotC = VarC->frame();
  RooPlot* PlotD = VarD->frame();

  RooPlot* PlotAKeys = VarA->frame();
  RooPlot* PlotBKeys = VarB->frame();
  RooPlot* PlotCKeys = VarC->frame();
  RooPlot* PlotDKeys = VarD->frame();

  // Make some Legends:
  int PdfColor = 44;
  int RawKeysColor = 8;
  int SmoothKeysColor = 46;
  int SmoothKeysRhoColor = 9; //38;

  // Dummy histograms
  TH1F* hPdf  = new TH1F("hPdf","", 1,0,1);      hPdf->SetLineColor( PdfColor );
  TH1F* hMark = new TH1F("hMark","", 1,0,1);     //hPdf->SetLineColor( PdfColor );

  TLegend* Leg1 = new TLegend(.1, .1, .9, .5);
  Leg1->SetFillColor(0);
  Leg1->AddEntry( hPdf,  "Pdf used to generate data", "l" );
  Leg1->AddEntry( hMark, "Raw Data", "P" );
  

  TH1F* hRawKeys = new TH1F("hRawKeys","", 1,0,1);             hRawKeys->SetLineColor( RawKeysColor );
  TH1F* hSmoothKeys = new TH1F("hSmoothKeys","", 1,0,1);       hSmoothKeys->SetLineColor( SmoothKeysColor );
  TH1F* hSmoothKeysRho = new TH1F("hSmoothKeysRho","", 1,0,1); hSmoothKeysRho->SetLineColor( SmoothKeysRhoColor );
 
  TLegend* Leg3 = new TLegend(.1, .1, .9, .5);
  Leg3->SetFillColor(0);
  Leg3->AddEntry( hRawKeys,       "Keys on Unbinned", "l" );
  Leg3->AddEntry( hSmoothKeys,    "Keys from histogram", "l" );
  Leg3->AddEntry( hSmoothKeysRho, "Keys with weight smoothing", "l" );


  TLegend* Leg4 = new TLegend(.1, .1, .9, .5);
  Leg4->SetFillColor(0);
  Leg4->AddEntry( (TObject*)0, "", "");
  Leg4->AddEntry( HistA_sm_no,    "Smoothed histogram from Keys", "l" );
  Leg4->AddEntry( HistA_sm,       "Histogram Smoothed with weight Keys", "l" );



  // Add the Titles

  Canvas->cd(++CurrentPad);
  TPaveText *pt_C = new TPaveText(.1,.6, .9,.9);
  pt_C->SetFillColor(0);
  pt_C->AddText("Raw, unbinned data");
  pt_C->AddText("and the underlying pdf");
  pt_C->Draw();
  Leg1->Draw();

  Canvas->cd(++CurrentPad);
  TPaveText *pt_A = new TPaveText(.1,.2, .9,.9);
  pt_A->SetFillColor(0);
  pt_A->AddText("Binned data");
  pt_A->AddText("This is the input histogram");
  pt_A->Draw();

  Canvas->cd(++CurrentPad);
  TPaveText *pt_D = new TPaveText(.1,.6, .9,.9);
  pt_D->SetFillColor(0);
  pt_D->AddText("Several Keys plots made");
  pt_D->AddText("on unbinned data, binned data,");
  pt_D->AddText("and binned data used weight parameter");
  pt_D->Draw();
  Leg3->Draw();

  Canvas->cd(++CurrentPad);
  TPaveText *pt_B = new TPaveText(.1,.6, .9,.9);
  pt_B->SetFillColor(0);
  pt_B->AddText("Smoothed histogram made from keys");
  pt_B->AddText("with and without weight parameter");
  pt_B->Draw();
  Leg4->Draw();

  ///////////// MAKE THE MAIN PLOTS HERE ///////////////


  // A
  Canvas->cd(++CurrentPad);
  DataA->plotOn( PlotA );
  w.pdf("PdfA")->plotOn(PlotA, LineColor(PdfColor) );
  PlotA->Draw();

  Canvas->cd(++CurrentPad);
  HistA->Draw();

  Canvas->cd(++CurrentPad);
  KeysA_unbinned->plotOn( PlotAKeys, LineColor(RawKeysColor) );
  KeysA_no->plotOn( PlotAKeys, LineColor(SmoothKeysColor) );
  KeysA_Rho->plotOn( PlotAKeys, LineColor(SmoothKeysRhoColor) );
  PlotAKeys->Draw();

  Canvas->cd(++CurrentPad);
  HistA_sm->Draw();
  HistA_sm_no->Draw("same");
  HistA_sm->Draw("same");



  // B
  Canvas->cd(++CurrentPad);
  DataB->plotOn( PlotB );
  w.pdf("PdfB")->plotOn(PlotB, LineColor(PdfColor) );
  PlotB->Draw();

  Canvas->cd(++CurrentPad);
  HistB->Draw();

  Canvas->cd(++CurrentPad);
  KeysB_unbinned->plotOn( PlotBKeys, LineColor(RawKeysColor) );
  KeysB_no->plotOn( PlotBKeys, LineColor(SmoothKeysColor) );
  KeysB_Rho->plotOn( PlotBKeys, LineColor(SmoothKeysRhoColor) );
  PlotBKeys->Draw();

  Canvas->cd(++CurrentPad);
  HistB_sm->Draw();
  HistB_sm_no->Draw("same");
  HistB_sm->Draw("same");

  // C
  Canvas->cd(++CurrentPad);
  DataC->plotOn( PlotC );
  w.pdf("PdfC")->plotOn(PlotC, LineColor(PdfColor) );
  PlotC->Draw();

  Canvas->cd(++CurrentPad);
  HistC->Draw();

  Canvas->cd(++CurrentPad);
  KeysC_unbinned->plotOn( PlotCKeys, LineColor(RawKeysColor) );
  KeysC_no->plotOn( PlotCKeys, LineColor(SmoothKeysColor) );
  KeysC_Rho->plotOn( PlotCKeys, LineColor(SmoothKeysRhoColor) );
  PlotCKeys->Draw();

  Canvas->cd(++CurrentPad);
  HistC_sm->Draw();
  HistC_sm_no->Draw("same");
  HistC_sm->Draw("same");


  // D
  Canvas->cd(++CurrentPad);
  DataD->plotOn( PlotD );
  w.pdf("PdfD")->plotOn(PlotD, LineColor(PdfColor) );
  PlotD->Draw();

  Canvas->cd(++CurrentPad);
  HistD->Draw();

  Canvas->cd(++CurrentPad);
  KeysD_unbinned->plotOn( PlotDKeys, LineColor(RawKeysColor) );
  KeysD_no->plotOn( PlotDKeys, LineColor(SmoothKeysColor) );
  KeysD_Rho->plotOn( PlotDKeys, LineColor(SmoothKeysRhoColor) );
  PlotDKeys->Draw();

  Canvas->cd(++CurrentPad);
  HistD_sm->Draw();
  HistD_sm_no->Draw("same");
  HistD_sm->Draw("same");

  Canvas->Print("KeysTest.eps");
  Canvas->Print("KeysTest.pdf");


}

