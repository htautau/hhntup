
/*SmoothHistogram

  A simple example to smooth a histogram using Kernel Esimation (KEYS).

  The smoothing algorithm it is based on:
  Kyle S. Cranmer, Comput.Phys.Commun. 136 (2001) 198-207
  e-Print: hep-ex/0011057


  Example Usage:

  // Smooth Histogram:
  TH1F* MyUglyHist = new TH1F(....); // Fill MyUglyHist
  TH1F* MySmoothHist = SmoothHistogram( MyUglyHist, "MySmoothHist");


  // Smooth all histograms in a file:
  SmoothHistogramsInFile( "Path/To/MyInputFile.root", "My/Output/FileName.root" );

  
  The functions do the following:

  TH1* SmoothHistogram( TH1* hist, const std::string& Name = "Hist", bool UseNEventsWeight=true, RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false );
  
     --> Return a TH1* representing the smoothed histogram from the input "hist".  The output TH1 has the name "Name", and the function has an optional argument
         to determine the amplitude parameter based on the total number of events used (and not on the weights of those events).  It can optionally use a mirror
         when constructing the KEYS pdf, and can print verbose information to the standard output.

  TH1* SmoothHistogram( TH1* hist, const std::string& Name, RooRealVar* Var, RooKeysPdf*& Keys_out, bool UseNEventsWeight=true, 
                        RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false );

      --> A similar method to the above, but it takes a RooRealVar as an input and constructes the Keys PDF as a function of that VAR.  The RooKeysPdf is then returned
          by setting the pointer argument "Keys_out" equal to the newly created RooKeysPdf (the user takes ownership of the pdf). 

  void SmoothHistogramsInFile( const std::string& InputFileName, const std::string& OutputFileName, bool CopyTree=false, bool UseNEventsWeight=true, 
                               RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false );

      --> Open the input file, loop over all histograms in the file (in all directories), smooth them, and save their smoothed versions in an output file.  The input
          file remains unchanged.  One can optionally copy over all trees that may exist in the input file (they are copied over identically).  One can optionally
	  use weights based on the number of evnets and can use mirrors; these are applied to all histograms that are smoothed.

*/


#include <string.h>

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TLegend.h"
#include "TKey.h"
#include "TRandom3.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TH1D.h"

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"

using namespace RooFit;

// Functions:

void SmoothHistogram(); 
TH1* SmoothHistogram( TH1* hist, const std::string& Name = "Hist", bool UseNEventsWeight=true, RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false );
TH1* SmoothHistogram( TH1* hist, const std::string& Name, RooRealVar* Var, RooKeysPdf*& Keys_out, bool UseNEventsWeight=true, RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false );
void SmoothHistogramsInFile(const std::string& InputFileName, const std::string& OutputFileName, bool CopyTree=false, bool UseNEventsWeight=true, RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false);


// "Internal" Functions
void RunExample();
void LoopAndSmoothHists( TDirectory *target, TList *sourcelist, bool CopyTree=false, bool UseNEventsWeight=true, RooKeysPdf::Mirror Mirror=RooKeysPdf::NoMirror, bool VERBOSE=false );

/////////////////////
// Implementations //
/////////////////////


void SmoothHistogram(){

// Run an example to create
// a simple gaussian histogram
// and to smooth it

  RunExample();

}

void SmoothHistogramsInFile(const std::string& InputFileName, const std::string& OutputFileName, bool CopyTree, bool UseNEventsWeight, RooKeysPdf::Mirror Mirror, bool VERBOSE) {

// Open Input file, collect all histograms, smooth them,
// and write them to the Output file:
// Optionally Copy any trees that may exist in the 
// original file


  TFile* Target = TFile::Open( OutputFileName.c_str(), "RECREATE" );
  TList* FileList = new TList();
  FileList->Add( TFile::Open(InputFileName.c_str()) );

  LoopAndSmoothHists( Target, FileList, CopyTree, UseNEventsWeight, Mirror, VERBOSE );

}


TH1* SmoothHistogram(TH1* hist, const std::string& Name, bool UseNEventsWeight, RooKeysPdf::Mirror Mirror,  bool VERBOSE) {

// Make a version that doesn't return the keys
// but instead deletes it (to prevent mem leaks)

  RooRealVar tmpVar("tmpVar","",hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
  RooKeysPdf* Keys_out = NULL;

  TH1* HistReturn = SmoothHistogram(hist, Name, &tmpVar, Keys_out, UseNEventsWeight, Mirror, VERBOSE);
  
  delete Keys_out;
  return HistReturn;
      
} 


// The main function:

TH1* SmoothHistogram(TH1* hist, const std::string& Name, RooRealVar* Var, RooKeysPdf*& Keys_out, bool UseNEventsWeight, RooKeysPdf::Mirror Mirror, bool VERBOSE) {

  // Create a histogram based on the input TH1
  // using KEYS to smooth the shape.

  // The input Var will be the variable
  // used to define the keys PDF
  // Return this keys via function argument

  
  if(VERBOSE)  std::cout << std::endl << "Begin SmoothHistogram() " << std::endl << std::endl;

  if( hist==NULL ) {
    std::cout << "Error: input histogram is NULL" << std::endl;
    if(VERBOSE) std::cout << "End SmoothHistogram() " << std::endl << std::endl;
    return NULL;
  }
  
  RooRealVar weight("weight","",0,100000);
  RooArgSet vars(*Var, weight);
  RooDataSet ds("ds","",vars,"weight");

  int NumBins = hist->GetNbinsX();
  double AxisMin = hist->GetXaxis()->GetXmin();
  double AxisMax = hist->GetXaxis()->GetXmax();

  int NumEntries = hist->GetEntries(); // Assumes it was filled properly...

  if( NumEntries == 0 ) {
    std::cout << "Histogram has no entries" << std::endl;
    if(VERBOSE) std::cout << "End SmoothHistogram() " << std::endl << std::endl;
    return hist;
  }


  int NumNonZeroBins = 0;

  for(int i=1; i<hist->GetNbinsX() + 1; ++i){
    if(hist->GetBinContent(i)<=0)
      continue;

    NumNonZeroBins++;
    Var->setVal(hist->GetBinCenter(i));
    //    weight->setVal(hist->GetBinContent(i));
    if(VERBOSE) std::cout << "adding bin at " << Var->getVal() <<" with weight " << hist->GetBinContent(i) << std::endl;
    ds.add(vars,hist->GetBinContent(i));
  }

  if( NumNonZeroBins < 3 ) {
    std::cout << "Warning: Making KEYS from histogram with < 3 non-zero entry" << std::endl;
    std::cout << "Currently this is not supported.  Returning original histogram." << std::endl;
    if(VERBOSE) std::cout << "End SmoothHistogram() " << std::endl << std::endl;
    return hist;
  }

  
  //  ds.Print();

  RooKeysPdf* keys;

  if( UseNEventsWeight ) {
    double Rho = TMath::Power( (double) NumNonZeroBins/NumEntries, .2); // ^1/5
    keys = new RooKeysPdf(Form("%s_keys_mirror",hist->GetName()),"",*Var,ds, Mirror, Rho);
    //keys = new RooKeysPdf(Form("%s_keys_nomirror",hist->GetName()),"",m,*ds, RooKeysPdf::NoMirror, Rho);
  }
  else {
    keys = new RooKeysPdf(Form("%s_keys_mirror",hist->GetName()),"",*Var,ds, Mirror, 1.0);
    //keys = new RooKeysPdf(Form("%s_keys_mirror",hist->GetName()),"",*Var,ds,RooKeysPdf::NoMirror);
  }

  if( keys == NULL ) {
    std::cout << "Error making Keys!  Returning NULL" << std::endl;
    if(VERBOSE) std::cout << "End SmoothHistogram() " << std::endl << std::endl;
    return NULL;
  }
  
  TH1* newHist = keys->createHistogram(Name.c_str(), *Var, Binning(NumBins, AxisMin, AxisMax));
  newHist->Scale(hist->Integral()/newHist->Integral());
  
  Keys_out = keys;

  if(VERBOSE) std::cout << "End SmoothHistogram() " << std::endl << std::endl;
  return newHist;

}



void LoopAndSmoothHists( TDirectory *target, TList *sourcelist, bool CopyTree, bool UseNEventsWeight, RooKeysPdf::Mirror Mirror, bool VERBOSE ) {

// Loop over all files in the source list, find their histograms,
// smooth them, and add them to the newly created target root file

  //  cout << "Target path: " << target->GetPath() << endl;
  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd( path );
  TDirectory *current_sourcedir = gDirectory;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;

  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    // read object from first source file
    first_source->cd( path );
    TObject *obj = key->ReadObj();

    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      // descendant of TH1 -> merge it
      
      TH1* h1 = (TH1*) obj;
      
      std::string Name = h1->GetName();
      // Make the new object:
      //      bool UseNEventsWeight = true;
      h1 = (TH1*) SmoothHistogram( h1, Name.c_str(), UseNEventsWeight, Mirror, VERBOSE );
      obj = h1;

    }
    else if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {

      // if it's a tree, continue
      if( CopyTree == false ) continue;

    }
    else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
      // it's a subdirectory

      if(VERBOSE) std::cout << "Found subdirectory " << obj->GetName() << endl;
      
      // create a new subdir of same name and title in the target file
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      // newdir is now the starting point of another round of merging
      // newdir still knows its depth within the target file via
      // GetPath(), so we can still figure out where we are in the recursion
      LoopAndSmoothHists( newdir, sourcelist, CopyTree, UseNEventsWeight, Mirror, VERBOSE );

    }
    if ( obj ) {
      target->cd();
      obj->Write( key->GetName() );
    }
    
  } // while ( ( TKey *key = (TKey*)nextkey() ) )
  
  // save modifications to target file
  target->SaveSelf(kTRUE);
  TH1::AddDirectory(status);

}


void RunExample() {

  // Make a toy histogram,
  // Smooth the histogram,
  // and plot them


  TH1F* HistRough = new TH1F("HistRough", "HistRough", 100,-10, 10);

  // Hist B is gaussian:
  int NumEntries = 200;
  TRandom3 RandomGauss(235325);
  for(int i = 0; i < NumEntries; ++i) {

    // Make a random gaussian value and fill the hist
    double RandVal = RandomGauss.Gaus( 0, 1 );
    HistRough->Fill(RandVal);

  }
  

  TH1F* HistSmooth = (TH1F*) SmoothHistogram( HistRough, "HistSmooth", true, RooKeysPdf::NoMirror, true  );
  
  TCanvas* Canvas = new TCanvas("Canvas","");
  Canvas->Divide(2,1);

  TLegend* RoughLeg  = new TLegend(.7, .8, .9, .9);
  RoughLeg->SetFillColor(0);
  TLegend* SmoothLeg = new TLegend(.7, .8, .9, .9);
  SmoothLeg->SetFillColor(0);


  Canvas->cd(1);
  HistRough->Draw();
  RoughLeg->AddEntry(HistRough, "Rough", "" );
  RoughLeg->Draw();

  Canvas->cd(2);
  HistSmooth->Draw();
  SmoothLeg->AddEntry(HistSmooth, "Smooth", "" );
  SmoothLeg->Draw();

}

