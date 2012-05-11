/******************************************************************************
Name:        TPileupReweighting

Author:      Will Buttinger
Created:     October 2011

Description: Tool to get the calculated MC pileup weight.
******************************************************************************/

// Preprocessor magic for debugging
#define XXX std::cout << " I am here: " << __FILE__ << ":" << __LINE__ << std::endl;

// This class' header
#include "TPileupReweighting.h"



// include math
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TLeaf.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TAxis.h>
#include <TString.h>
#include <TRandom3.h>


namespace AnalysisFramework
{
namespace External
{

//ClassImp(Root::TPileupReweighting)


//=============================================================================
// Constructor
//=============================================================================
Root::TPileupReweighting::TPileupReweighting(const char* name) :
  TNamed(name,"notitle"),
   m_SetWarnings(true),m_debugging(false),
   m_countingMode(true),m_defaultChannel(0),m_unrepresentedDataAction(0),m_isInitialized(false),m_lumiVectorIsLoaded(false),
   m_dataScaleFactorX(1.),m_dataScaleFactorY(1.),m_dataScaleFactorZ(1.),
   m_mcScaleFactorX(1.),m_mcScaleFactorY(1.),m_mcScaleFactorZ(1.),
   m_nextPeriodNumber(1),m_ignoreFilePeriods(false),m_metadatatree(0)
{
   m_random3 = new TRandom3(0);
   m_random3->SetSeed(1);
   //load the default pileup histogram binning
   m_emptyHistograms["pileup"] = new TH1D("pileup_default","pileup_default",100,0,50);
   m_emptyHistograms["pileup"]->SetDirectory(0);
   //set the lumi vector to size=11
   m_integratedLumiVector.resize(11,0.);
}

Double_t Root::TPileupReweighting::GetIntegratedLumiFraction(Int_t mcRunNumber, UInt_t start, UInt_t end) {
   if(!m_isInitialized) {
      Error("GetIntegratedLumiFraction", "Please initialize the tool before retrieving lumi fractions");
      throw std::runtime_error("Throwing 1");
   }
   if(!m_lumiVectorIsLoaded) {
      Error("GetIntegratedLumiFraction","No Lumicalc file loaded, so no lumi fraction possible, returning 0");
      return 0;
   }
   //check if the runNumber has been reassigned
   if(m_mcRemappings.find(mcRunNumber)!=m_mcRemappings.end()) mcRunNumber=m_mcRemappings[mcRunNumber];
   const std::map<Int_t, Double_t>::const_iterator it = periodTotals["pileup"][-1].find(mcRunNumber);
   if(it==periodTotals["pileup"][-1].end() || it->second == 0) {
      if(m_SetWarnings) Warning("GetIntegratedLumiFraction","%d has no associated lumi. Returning 0",mcRunNumber);
      return 0; /*throw 33;*/
   }
   //loop over assigned runs, if in range then include in numerator
   double total = 0;
   std::map<UInt_t, Double_t>::const_iterator itend = dataPeriodRunTotals["pileup"][mcRunNumber].end();
   for(std::map<UInt_t, Double_t>::const_iterator it2 = dataPeriodRunTotals["pileup"][mcRunNumber].begin(); it2 != itend; ++it2) {
      if(it2->first >= start && it2->first <= end) total += it2->second;
   }

   return total / it->second;

}

Double_t Root::TPileupReweighting::GetIntegratedLumi(UInt_t start, UInt_t end) {
   //look through dataPeriodRunTotals["pileup"][-1] for runs inside the given period
   double total = 0;
   std::map<UInt_t, Double_t>::const_iterator itend = dataPeriodRunTotals["pileup"][-1].end();
   for(std::map<UInt_t, Double_t>::const_iterator it2 = dataPeriodRunTotals["pileup"][-1].begin(); it2 != itend; ++it2) {
      if(it2->first >= start && it2->first <= end) total += it2->second;
   }
   return total/1E6;
}

Double_t Root::TPileupReweighting::GetIntegratedLumi(Int_t periodNumber, UInt_t start, UInt_t end) {
   //look through dataPeriodRunTotals["pileup"][periodNumber] for runs inside the given period
   double total = 0;
   std::map<UInt_t, Double_t>::const_iterator itend = dataPeriodRunTotals["pileup"][periodNumber].end();
   for(std::map<UInt_t, Double_t>::const_iterator it2 = dataPeriodRunTotals["pileup"][periodNumber].begin(); it2 != itend; ++it2) {
      if(it2->first >= start && it2->first <= end) total += it2->second;
   }
   return total/1E6;
}

//=============================================================================
// Destructor
//=============================================================================
Root::TPileupReweighting::~TPileupReweighting() {
   delete m_random3;
   //delete all histograms
   for(std::map<TString, std::map<Int_t,std::map<Int_t, TH1*> > >::iterator it0=m_inputHistograms.begin();it0!=m_inputHistograms.end();++it0) {
      for(std::map<Int_t,std::map<Int_t, TH1*> >::iterator it1=it0->second.begin();it1!=it0->second.end();++it1) {
         for(std::map<Int_t, TH1*>::iterator it2=it1->second.begin();it2!=it1->second.end();++it2) {
            if(it2->second) delete it2->second;
         }
         it1->second.clear();
      }
      it0->second.clear();
   }
   m_inputHistograms.clear();
   if(m_debugging) Info("destructor","Deleting primary distributions");
   for(std::map<TString, std::map<Int_t,std::map<Int_t, TH1D*> > >::iterator it0=primaryDistributions.begin();it0!=primaryDistributions.end();++it0) {
      for(std::map<Int_t,std::map<Int_t, TH1D*> >::iterator it1=it0->second.begin();it1!=it0->second.end();++it1) {
         for(std::map<Int_t, TH1D*>::iterator it2=it1->second.begin();it2!=it1->second.end();++it2) {
            if(it2->second) {
               if(m_debugging) { Info("destructor","weight=%s chan=%d run=%d",it0->first.Data(),it1->first,it2->first);
                  Info("destructor","Deleting %s",it2->second->GetName());
               }
               delete it2->second;
            }
         }
         it1->second.clear();
      }
      it0->second.clear();
   }
   primaryDistributions.clear();
   if(m_debugging) Info("destructor","Deleting secondary distributions");
   for(std::map<TString, std::map<Int_t,std::map<Int_t, TH2D*> > >::iterator it0=secondaryDistributions.begin();it0!=secondaryDistributions.end();++it0) {
      for(std::map<Int_t,std::map<Int_t, TH2D*> >::iterator it1=it0->second.begin();it1!=it0->second.end();++it1) {
         for(std::map<Int_t, TH2D*>::iterator it2=it1->second.begin();it2!=it1->second.end();++it2) {
            if(it2->second) {
               if(m_debugging) { Info("destructor","weight=%s chan=%d run=%d",it0->first.Data(),it1->first,it2->first);
                  Info("destructor","Deleting %s",it2->second->GetName());
               }
               delete it2->second;
            }
         }
         it1->second.clear();
      }
      it0->second.clear();
   }
   secondaryDistributions.clear();

}

Int_t Root::TPileupReweighting::UsePeriodConfig(const TString configName) {
   if(configName=="MC11a") {
      AddPeriod(180164, 177986,180481); //associates mc runnumber 180164 with data period 177986 to 180481 (period B-D)
      AddPeriod(183003, 180614,184169); //period E-H
      AddPeriod(185649, 185353,186934); //period I-K1. For I-K you would change the last number to 187815
      AddPeriod(185761, 186935,191933); //everything else. Thanks Ellie!
      Info("UsePeriodConfig","Using MC11a Period configuration");
      return 0;
   } else if(configName=="MC11b" || configName=="MC11c") {
      AddPeriod(180164, 177986, 180481);
      AddPeriod(183003, 180614, 184169);
      AddPeriod(186169, 185353, 187815);
      AddPeriod(189751, 188902, 191933);
      Info("UsePeriodConfig","Using MC11b/c Period configuration");
      return 0;
   }
   Error("UsePeriodConfig","Unrecognized period config");
   return -1;
}

std::vector<Double_t> Root::TPileupReweighting::getIntegratedLumiVector() {
   if(!m_lumiVectorIsLoaded) {
      Error("getIntegratedLumiVector","The vector has not been filled. You must AddDataDistribution, using the iLumiCalc file, to fill it");
      throw std::runtime_error("Throwing 99");
   }
   if(!m_isInitialized) {
      Error("getIntegratedLumiVector","Please initialize the tool before retrieving the lumi vector");
      throw std::runtime_error("Throwing 1");
   }
   return m_integratedLumiVector;
}

///returns a PeriodID. These count up from 1
Int_t Root::TPileupReweighting::AddPeriod(Int_t mcRunNumber, UInt_t start, UInt_t end) {
   if(m_isInitialized) {
      Error("AddPeriod","You cannot AddPeriod after initializing the tool. Reorder your code!");
      throw std::runtime_error("Throwing 1");
   }
   //check if this period is already defined
   for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods=m_periods[mcRunNumber].begin();periods!=m_periods[mcRunNumber].end();++periods) {
      std::pair<UInt_t,UInt_t> period = periods->second;
      if(period.first == start && period.second==end) {
         return 0;
      } else if((start <= period.second && start >= period.first) || (end <= period.second && end >= period.first) ) {
         //overlaps with existing period. reject this Add
         Error("AddPeriod","(%d,%d,%d) overlaps with existing period setup",mcRunNumber,start,end);
         throw std::runtime_error("Throwing 2");
      }
   }
   std::pair<UInt_t,UInt_t> thePair; thePair.first=start;thePair.second=end;
   m_periods[mcRunNumber][m_nextPeriodNumber] = thePair;
   //also put this in a handy backmapping for getting this period's periodWeight and mc distributions (uses same ones)
   m_periodToMCRun[m_nextPeriodNumber]=mcRunNumber;
   m_nextPeriodNumber++;
   return m_nextPeriodNumber-1;
}

Int_t Root::TPileupReweighting::GetPeriodNumber(Int_t mcRunNumber, UInt_t start, UInt_t end) {
   //try to locate this period number in the map
   for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods=m_periods[mcRunNumber].begin();periods!=m_periods[mcRunNumber].end();++periods) {
      std::pair<UInt_t,UInt_t> period = periods->second;
      if(period.first == start && period.second==end) {
         return periods->first;
      }
   }
   return -1;
}

Int_t Root::TPileupReweighting::AddBinning(const TString weightName,Int_t nbinsx, Double_t* xbins) {
   if(m_emptyHistograms.find(weightName)!=m_emptyHistograms.end()) {
        delete m_emptyHistograms[weightName];
   }
   m_emptyHistograms[weightName] = new TH1D(weightName+"_default","default",nbinsx,xbins);
   m_emptyHistograms[weightName]->SetDirectory(0);
   return 0;
}
Int_t Root::TPileupReweighting::AddBinning(const TString weightName,Int_t nbinsx, Double_t xlow, Double_t xup) {
   if(m_emptyHistograms.find(weightName)!=m_emptyHistograms.end()) {
        delete m_emptyHistograms[weightName];
   }
   m_emptyHistograms[weightName] = new TH1D(weightName+"_default","default",nbinsx,xlow,xup);
   m_emptyHistograms[weightName]->SetDirectory(0);
   return 0;
}
Int_t Root::TPileupReweighting::AddBinning(const TString weightName,Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) {
   if(m_emptyHistograms.find(weightName)!=m_emptyHistograms.end()) {
        delete m_emptyHistograms[weightName];
   }
   m_emptyHistograms[weightName] = new TH2D(weightName+"_default","default",nbinsx,xlow,xup,nbinsy,ylow,yup);
   m_emptyHistograms[weightName]->SetDirectory(0);
   return 0;
}

Int_t Root::TPileupReweighting::AddBinning(const TString weightName,TH1* hist) {
   if(!hist) return -1;
   if(m_emptyHistograms.find(weightName)!=m_emptyHistograms.end()) {
        delete m_emptyHistograms[weightName];
   }
   m_emptyHistograms[weightName] = dynamic_cast<TH1*>(hist->Clone(weightName+"_default"));
   m_emptyHistograms[weightName]->SetDirectory(0);
   return 0;
}


//combine run number's histograms, periods, everything, then remove trace of second run number
void Root::TPileupReweighting::MergeMCRunNumbers(Int_t mcRunNumber1, Int_t mcRunNumber2, Int_t newMCRunNumber) {
   if(m_isInitialized) {
      Error("MergeMCRunNumbers","You cannot MergeMCRunNumbers after initializing the tool. Reorder your code!");
      throw std::runtime_error("Throwing 4");
   }
   if(newMCRunNumber==0) newMCRunNumber=mcRunNumber1;
   m_mcMerges[newMCRunNumber][mcRunNumber1]=mcRunNumber2;

}

TH1* Root::TPileupReweighting::CloneEmptyHistogram(const TString weightName, Int_t runNumber, Int_t channelNumber) {

   TString s(weightName);
   if(channelNumber>=0) {
      s += "_chan"; s += channelNumber;
   } else {
      s += "_data";
   }
   s+="_run"; s+= runNumber;
   //get the empty histogram
   if(m_emptyHistograms.find(weightName)==m_emptyHistograms.end()) {
      Error("CloneEmptyHistogram","There is no binning info for weight %s",weightName.Data());
      throw std::runtime_error("Throwing 47");
   }

   TH1* out = dynamic_cast<TH1*>(m_emptyHistograms[weightName]->Clone(s));
   out->SetTitle(s);
   out->SetDirectory(0); //take ownership
   out->SetEntries(0);
   return out;

}

Int_t Root::TPileupReweighting::GenerateMetaDataFile(const TString fileName,const TString channelBranchName) {

   TTree inTree("in","in");
   inTree.ReadFile(fileName);
   TTree outTree("ChannelMetaData","ChannelMetaData");
   Int_t chanNum;
      std::map<TString, Double_t> data;
      TObjArray *leaves = inTree.GetListOfLeaves();
      if(leaves==0) {Error("GenerateMetaDataFile","No leaves"); return -1; }
      for(Int_t i=0;i<leaves->GetEntries();++i) {
         TLeaf *leaf = (TLeaf *)leaves->At(i);
         if(leaf) {
            TBranch *branch = leaf->GetBranch();
            if(strcmp(branch->GetName(),channelBranchName)==0) {
               //this is the channel branch
               if(strcmp(leaf->GetTypeName(),"Int_t")!=0) {
                  Error("GenerateMetaDataFile","Channel Branch must be type Int_t"); return -1;
               }
               branch->SetAddress(&chanNum);
               outTree.Branch(channelBranchName,&chanNum);
            } else if(strcmp(leaf->GetTypeName(),"Double_t")!=0) {
                  Warning("GenerateMetaDataFile","Cannot read non-double branch: %s",branch->GetName());
            } else {
               branch->SetAddress(&(data[branch->GetName()]));
               outTree.Branch(branch->GetName(),&(data[branch->GetName()]));
            }
         }
      }
      //loop over tree entries and read in the values according to property type
      for(Int_t i=0;i<inTree.GetEntries();++i) {
         inTree.GetEntry(i);
         outTree.Fill();
      }

      //remove the file extension and replace with .root
      TString outName = fileName(0,fileName.Last('.'));
      outName += ".prw.root";
      TFile f1(outName,"RECREATE");
      outTree.Write();
      f1.Close();
      Info("GenerateMetaDataFile","Succesfully Generated File %s",outName.Data());
      return 0;
}

Int_t Root::TPileupReweighting::AddMetaDataFile(const TString fileName,const TString channelBranchName) {
   TDirectory* origDir = gDirectory;
   TTree* tmp = 0;
   TFile* rootFile = 0;
   if(fileName.EndsWith(".root")) {
      //looks for ChannelMetaData ttree
      // Load the data ROOT file
      rootFile = new TFile( fileName, "READ" );
      if ( !rootFile ) {
         Error("AddMetaDataFile","Could not find file: %s",fileName.Data());
         throw std::runtime_error("Throwing 6");
      } else {
         //try to get the the known TTrees
         tmp = (TTree*)rootFile->Get( "ChannelMetaData" );
         if(!tmp) {
            Error("AddMetaDataFile","%s is not a valid metadata file. Should have a ChannelMetaData TTree",fileName.Data());
            throw std::runtime_error("Throwing 7");
         }
      }
   } else {
      //try to make the TTree by reading the file
      tmp = new TTree("ChannelMetaData","ChannelMetaData");
      tmp->ReadFile(fileName);
   }
   Int_t chanNum;
   std::map<TString, Double_t> data;
   TObjArray *leaves = tmp->GetListOfLeaves();
   if(leaves==0) {Error("AddMetaDataFile","No leaves"); return -1; }
   for(Int_t i=0;i<leaves->GetEntries();++i) {
      TLeaf *leaf = (TLeaf *)leaves->At(i);
      if(leaf) {
         TBranch *branch = leaf->GetBranch();
         if(strcmp(branch->GetName(),channelBranchName)==0) {
            //this is the channel branch
            if(strcmp(leaf->GetTypeName(),"Int_t")!=0) {
               Error("AddMetaDataFile","Channel Branch must be type Int_t");
               throw std::runtime_error("Throwing 7");
            }
            branch->SetAddress(&chanNum);
         } else if(strcmp(leaf->GetTypeName(),"Double_t")!=0) {
               Warning("AddMetaDataFile","Cannot read non-double branch: %s",branch->GetName());
         } else {
            branch->SetAddress(&(data[branch->GetName()]));
         }
      }
   }
   //loop over tree entries and read in the values
   for(Int_t i=0;i<tmp->GetEntries();++i) {
      tmp->GetEntry(i);
      for(std::map<TString, Double_t>::iterator it = data.begin();it!=data.end();++it) {
         if(m_metadata.find(it->first)!=m_metadata.end()&&m_metadata[it->first].find(chanNum)!=m_metadata[it->first].end() && m_metadata[it->first][chanNum]!=it->second) {
            Warning("AddMetaDataFile", "Overriding metadata [%s,%d]. %f becomes %f",(it->first).Data(),chanNum,m_metadata[it->first][chanNum],it->second);
         }
         m_metadata[it->first][chanNum] = it->second;
      }
   }

   if(tmp) delete tmp;
   if(rootFile) { rootFile->Close();delete rootFile;}

   // Return to the original ROOT directory
   gDirectory = origDir;

   return 0;
}

TTree* Root::TPileupReweighting::GetMetaDataTree() {
   if(m_metadatatree) delete m_metadatatree;
   m_metadatatree = new TTree(TString(this->GetName())+"MetaData",TString(this->GetName())+"MetaData");
   m_metadatatree->SetDirectory(0);
   Int_t channel=0;
   m_metadatatree->Branch("mc_channel_number",&channel);
   std::map<TString,Double_t> data;
   std::map<Int_t,bool> channels;
   //create a branch for each metadata item
   for(std::map<TString,std::map<Int_t,Double_t> >::iterator it=m_metadata.begin();it!=m_metadata.end();++it) {
      for(std::map<Int_t,Double_t>::iterator it2=(it->second).begin();it2!=(it->second).end();++it2) {
         channels[it2->first]=true; //records which channels we have metadata for
      }
      if(data.find(it->first)==data.end()) {
         //new branch
         data[it->first]=0.;
         m_metadatatree->Branch(it->first,&(data[it->first]));
      }
   }
   //also create branches for the NumberOfEvents and SumOfEventWeights
   data["NumberOfEvents"]=0.;
   m_metadatatree->Branch("NumberOfEvents",&(data["NumberOfEvents"]));
   data["SumOfEventWeights"]=0.;
   m_metadatatree->Branch("SumOfEventWeights",&(data["SumOfEventWeights"]));
   //and add to channels list any channels that have events
   for(std::map<Int_t,Double_t>::iterator it=globalNumberOfEntries["pileup"].begin();it!=globalNumberOfEntries["pileup"].end();++it) {
      if(it->first>=0 && it->second>0.) channels[it->first]=true;
   }

   //now loop over the channels and fill the TTree
   for(std::map<Int_t,bool>::iterator it=channels.begin();it!=channels.end();++it) {
      channel=it->first;
      for(std::map<TString,Double_t>::iterator it2=data.begin();it2!=data.end();++it2) {
         if(it2->first=="NumberOfEvents") {
            //check for this in globalNumberOfEntries
            if(globalNumberOfEntries["pileup"].find(channel)==globalNumberOfEntries["pileup"].end()) {
               data[it2->first]=0.;
               Warning("GetMetaDataTree","Channel %d does not have MetaData %s",it->first,(it2->first).Data());
            } else {
               data[it2->first]=globalNumberOfEntries["pileup"][channel];
            }
         } else if(it2->first=="SumOfEventWeights") {
            //check for this in globalTotals
            if(globalTotals["pileup"].find(channel)==globalTotals["pileup"].end()) {
               data[it2->first]=0.;
               Warning("GetMetaDataTree","Channel %d does not have MetaData %s",it->first,(it2->first).Data());
            } else {
               data[it2->first]=globalTotals["pileup"][channel];
            }
         } else {
            if(m_metadata[it2->first].find(channel)==m_metadata[it2->first].end()) {
               //this channel doesn't have this property
               data[it2->first]=0.;
               Warning("GetMetaDataTree","Channel %d does not have MetaData %s",it->first,(it2->first).Data());
            } else {
               data[it2->first]=m_metadata[it2->first][channel];
            }
         }
      }
      m_metadatatree->Fill();
   }

   m_metadatatree->BuildIndex("mc_channel_number");
   m_metadatatree->ResetBranchAddresses();

   return m_metadatatree;
}

Int_t Root::TPileupReweighting::AddDistribution(TH1* hist, Int_t runNumber, Int_t channelNumber) {
   return AddDistribution(hist,"pileup",runNumber,channelNumber);
}

void Root::TPileupReweighting::AddDistributionTree(TTree *tree, TFile *file) {

   bool isMC=true;
   //should have the standard structure
   Int_t channel = 0; UInt_t runNbr = 0;
   std::vector<UInt_t>* pStarts = 0;std::vector<UInt_t>* pEnds = 0;
   Char_t histName[150];
   Char_t customName[150];
   if(strcmp(tree->GetName(),"MCPileupReweighting")==0) {strcpy(customName,"pileup");isMC=true;}
   else {
      if(tree->SetBranchAddress("CustomName",&customName)!=0) {
         Error("AddDistributionTree","Could not find CustomName branch in TTree");throw std::runtime_error("Throwing 18");
      }
   }
   if(strcmp(tree->GetName(),"DataCustomReweighting")==0) {channel=-1;isMC=false;}
   else {
      if(tree->SetBranchAddress("Channel",&channel)!=0) {
         Error("AddDistributionTree","Could not find Channel branch in TTree");throw std::runtime_error("Throwing 18");
      }
   }
   if(tree->SetBranchAddress("RunNumber",&runNbr)!=0) {
      Error("AddDistributionTree","Could not find RunNumber branch in TTree");throw std::runtime_error("Throwing 18");
   }
   if(isMC) {
      if(tree->SetBranchAddress("PeriodStarts",&pStarts)!=0) {
         Error("AddDistributionTree","Could not find PeriodStarts branch in TTree");throw std::runtime_error("Throwing 18");
      }
      if(tree->SetBranchAddress("PeriodEnds",&pEnds)!=0) {
         Error("AddDistributionTree","Could not find PeriodEnds branch in TTree");throw std::runtime_error("Throwing 18");
      }
   }

   if(tree->SetBranchAddress("HistName",&histName)!=0) {
      Error("AddDistributionTree","Could not find HistName branch in TTree");throw std::runtime_error("Throwing 18");
   }
   long n = tree->GetEntries();
   std::map<TString,bool> loadedHistos; //avoid double-loading from this file
   for(long i=0;i<n;i++) {
      tree->GetEntry(i);
      TString sHistName(histName);
      TString weightName(customName);
      if(loadedHistos.find(sHistName)==loadedHistos.end()) {
         loadedHistos[sHistName]=true;
         if(!m_ignoreFilePeriods && isMC) {
            for(unsigned int j=0;j<pStarts->size();j++) {
               unsigned int start = pStarts->at(j);unsigned int end = pEnds->at(j);
               AddPeriod(runNbr,start,end);
            }
         }
         //load the histogram
         TH1 *histo = (TH1*)file->Get( sHistName );
         if(!histo) {
            Error("AddDistributionTree","Unable to find the histogram %s in the File %s",sHistName.Data(),file->GetName());
            throw std::runtime_error("Throwing 21");
         }
         //add it to the histograms
         AddDistribution(histo, weightName, runNbr, channel);
      }
   }

}

//data is signalled by a negative channel number
Int_t Root::TPileupReweighting::AddDistribution(TH1* hist,const TString weightName,  Int_t runNumber, Int_t channelNumber) {

   if(m_isInitialized) {
      Error("AddDistribution","You cannot AddDistribution after initializing the tool. Reorder your code!");
      throw std::runtime_error("Throwing 5");
   }


   //std::map<Int_t,TH1*>& inputHistograms = (channelNumber>=0) ? m_inputMCHistograms[weightName][channelNumber] : m_inputDataHistograms[weightName];
   std::map<Int_t,TH1*>& inputHistograms = m_inputHistograms[weightName][channelNumber];

   //check if we need to create an empty histogram
   if(inputHistograms.find(runNumber) == inputHistograms.end()) {
      //if no emptyHistogram specified, we will use this histogram as the structure for the empty;
      if(m_emptyHistograms.find(weightName)==m_emptyHistograms.end()) {
         TString s = weightName; s+= "_default";
         m_emptyHistograms[weightName] = dynamic_cast<TH1*>(hist->Clone(s));
         m_emptyHistograms[weightName]->Reset();//clear the histogram
         m_emptyHistograms[weightName]->SetEntries(0);
         m_emptyHistograms[weightName]->SetDirectory(0);
      }
      inputHistograms[runNumber] = CloneEmptyHistogram(weightName,runNumber,channelNumber);
   }
   //iterator over bins of histogram, filling the TH1 stored in the data map
   Int_t numEntries = (int)(inputHistograms[runNumber]->GetEntries());
   Int_t bin,binx,biny,binz;
   for(binz=1; binz<=hist->GetNbinsZ(); binz++) {
      for(biny=1; biny<=hist->GetNbinsY(); biny++) {
         for(binx=1; binx<=hist->GetNbinsX(); binx++) {
            bin = hist->GetBin(binx,biny,binz);
            Double_t value = hist->GetBinContent(bin);
            Double_t x = hist->GetXaxis()->GetBinCenter(binx);
            Double_t y = hist->GetYaxis()->GetBinCenter(biny);
            Double_t z = hist->GetZaxis()->GetBinCenter(binz);
            //shift x,y,z by the MCScaleFactors as appropriate
            if(channelNumber>=0) {x *= m_mcScaleFactorX; y *= m_mcScaleFactorY; z *= m_mcScaleFactorZ;}
            else { x *= m_dataScaleFactorX; y *= m_dataScaleFactorY; z *= m_dataScaleFactorZ; }
            Int_t inBin = inputHistograms[runNumber]->FindFixBin(x,y,z);
            Double_t inValue = inputHistograms[runNumber]->GetBinContent(inBin);
            inputHistograms[runNumber]->SetBinContent(inBin,inValue+value);
         }
      }
   }
   //also keep track of the number of entries
   //SetBinContent screws with the entry count, so had to record it before the loops above
   inputHistograms[runNumber]->SetEntries(numEntries+hist->GetEntries());
   m_countingMode=false;
   return 0;

}

Int_t Root::TPileupReweighting::AddDistributionFromFile(const TString fileName, const TString histName, const TString weightName, Int_t runNumber, Int_t channelNumber) {
   if(m_isInitialized) {
      Error("AddDistributionFromFile","You cannot AddDistributionFromFile after initializing the tool. Reorder your code!");
      throw std::runtime_error("Throwing 5");
   }

    // Cache the current directory in the root file
   TDirectory* origDir = gDirectory;
   // Load the data ROOT file
   TFile* rootFile = new TFile( fileName, "READ" );
   if ( !rootFile ) {
      Error("AddDistributionFromFile","Could not find file: %s",fileName.Data());
      throw std::runtime_error("Throwing 6");
   } else {
      //decide if this is a histogram or a ttree
      if(!rootFile->Get(histName)) {
         Error("AddDistributionFromFile","Could not find object %s in file %s",histName.Data(),fileName.Data());
         throw std::runtime_error("Throwing 7");
      }
      if(rootFile->Get(histName)->InheritsFrom("TH1")) {
         TH1 *tmp = (TH1*)rootFile->Get( histName );
         //tmp->SetDirectory(0);//take ownership
         AddDistribution(tmp, weightName, runNumber, channelNumber);
      } else {
         Error("AddDistributionFromFile","%s is not a Histogram",histName.Data());
         throw std::runtime_error("Throwing 8");
      }
      delete rootFile;
   }
   m_countingMode=false;
   // Return to the original ROOT directory
   gDirectory = origDir;

   return 0;

}

Int_t Root::TPileupReweighting::AddDistributionFromFile(const TString fileName, const TString histName, Int_t runNumber, Int_t channelNumber) {
   return AddDistributionFromFile(fileName,histName,"pileup",runNumber,channelNumber);
}

Int_t Root::TPileupReweighting::AddLumiCalcFile(const TString fileName) {

   if(m_isInitialized) {
      Error("AddLumiCalcFile","You cannot AddLumiCalcFile after initializing the tool. Reorder your code!");
      throw std::runtime_error("Throwing 5");
   }
   TDirectory* origDir = gDirectory;
   // Load the data ROOT file
   TFile* rootFile = new TFile( fileName, "READ" );
   if ( !rootFile ) {
      Error("AddConfigFile","Could not find file: %s",fileName.Data());
      throw std::runtime_error("Throwing 6");
   } else {
      //try to get the the known TTrees
      TTree *tmp = (TTree*)rootFile->Get( "LumiMetaData" );
      if(tmp) {
         Info("AddLumiCalcFile","Adding LumiMetaData...");
         //structure expected is as given by iLumiCalc:
         //   RunNbr, AvergeInteractionPerXing, IntLumi
         UInt_t runNbr=0;Float_t intLumi=0;TBranch *b_runNbr;TBranch *b_intLumi;
         Float_t mu=0.; TBranch *b_mu;
         if(tmp->SetBranchAddress("RunNbr",&runNbr,&b_runNbr)!=0) {
            Error("AddLumiCalcFile","Could not find RunNbr branch in Data TTree");throw std::runtime_error("Throwing 8");
         }
         if(tmp->SetBranchAddress("AvergeInteractionPerXing",&mu,&b_mu)!=0) {
            Error("AddLumiCalcFile","Could not find AvergeInteractionPerXing branch in Data TTree");throw std::runtime_error("Throwing 9");
         }
         if(tmp->SetBranchAddress("IntLumi",&intLumi,&b_intLumi)!=0) {
            Error("AddLumiCalcFile","Could not find IntLumi branch in Data TTree");throw std::runtime_error("Throwing 10");
         }
         long nEntries = tmp->GetEntries();
         //std::map<Int_t, TH1*>& histograms = m_inputDataHistograms["pileup"];
         std::map<Int_t, TH1*>& histograms = m_inputHistograms["pileup"][-1];
         for(long i=0;i<nEntries;i++) {
            b_runNbr->GetEntry(i);b_intLumi->GetEntry(i);b_mu->GetEntry(i);
            //rescale the mu value
            mu *= m_dataScaleFactorX;
            //fill into input data histograms
            //check if we need to create an empty histogram
            if(histograms.find(runNbr) == histograms.end()) {
               histograms[runNbr] = CloneEmptyHistogram("pileup",runNbr,-1);
            }
            histograms[runNbr]->Fill(mu,intLumi);
         }
         m_countingMode=false;
         m_lumiVectorIsLoaded=true;
      }
   }

   // Return to the original ROOT directory
   gDirectory = origDir;

   return 0;
}

Int_t Root::TPileupReweighting::AddConfigFile(const TString fileName) {

   //open the file and look for config TTrees
   //recognised TTrees are: ChannelMetaData, MCPileupReweighting, MCCustomReweighting, DataCustomReweighting

   if(m_isInitialized) {
      Error("AddConfigFile","You cannot AddLumiCalcFile after initializing the tool. Reorder your code!");
      throw std::runtime_error("Throwing 5");
   }

   TDirectory* origDir = gDirectory;
   // Load the data ROOT file
   TFile* rootFile = new TFile( fileName, "READ" );
   if ( !rootFile ) {
      Error("AddConfigFile","Could not find file: %s",fileName.Data());
      throw std::runtime_error("Throwing 6");
   } else {
      //try to get the the known TTrees
      TTree *tmp = (TTree*)rootFile->Get( "MCPileupReweighting" );
      if(tmp) {
         Info("AddConfigFile","Adding MCPileupReweighting...");
         AddDistributionTree(tmp,rootFile);
         m_countingMode=false;
      }
      tmp = 0;tmp = (TTree*)rootFile->Get( "MCCustomReweighting" );
      if(tmp) {
         Info("AddConfigFile","Adding MCCustomReweighting...");
         AddDistributionTree(tmp,rootFile);
         m_countingMode=false;
      }
      tmp = 0;tmp = (TTree*)rootFile->Get( "DataCustomReweighting" );
      if(tmp) {
         Info("AddConfigFile","Adding DataCustomReweighting...");
         AddDistributionTree(tmp,rootFile);
         m_countingMode=false;
      }
   }

   // Return to the original ROOT directory
   gDirectory = origDir;

   return 0;
}

void Root::TPileupReweighting::AddDataDistribution(const TString& dataRootFileName,
                    const TString& dataRootHistName,Int_t runNumber) {
   Warning("AddDataDistribution","Method Depricated. Please use either AddConfigFile or AddDistributionFromFile");
   if(dataRootHistName=="LumiMetaData") AddLumiCalcFile(dataRootFileName);
   else AddDistributionFromFile(dataRootFileName,dataRootHistName,"pileup",runNumber,-1);
}

void Root::TPileupReweighting::AddMCDistribution(const TString& mcRootFileName,
                    const TString& mcRootHistName, int runNumber, int channelNumber) {
   Warning("AddMCDistribution","Method Depricated. Please use either AddConfigFile or AddDistributionFromFile");
   if(mcRootHistName=="MCPileupReweighting") AddConfigFile(mcRootFileName);
   else AddDistributionFromFile(mcRootFileName,mcRootHistName,"pileup",runNumber,channelNumber);
}

void Root::TPileupReweighting::AddCustomDataDistribution(const TString& customRootFileName,
                    const TString& customRootHistName, TString customName, int runNumber) {
   Warning("AddCustomDataDistribution","Method Depricated. Please use either AddConfigFile or AddDistributionFromFile");
   if(customRootHistName=="DataCustomReweighting") AddConfigFile(customRootFileName);
   else AddDistributionFromFile(customRootFileName,customRootHistName,customName,runNumber,-1);
}


void Root::TPileupReweighting::AddCustomMCDistribution(const TString& customRootFileName,
                    const TString& customRootHistName, TString customName, int runNumber, int channelNumber) {
   Warning("AddCustomMCDistribution","Method Depricated. Please use either AddConfigFile or AddDistributionFromFile");
   if(customRootHistName=="MCCustomReweighting") AddConfigFile(customRootFileName);
   else AddDistributionFromFile(customRootFileName,customRootHistName,customName,runNumber,channelNumber);
}


//channelNumber=-1 is data
Int_t Root::TPileupReweighting::FactorizeDistribution(TH1* hist, const TString weightName,
                           Int_t channelNumber, Int_t periodNumber,bool includeInMCRun,bool includeInGlobal)   {

   //determine which mc runNumber this periodNumber is assigned to
   Int_t mcRunNumber = m_periodToMCRun[periodNumber];

   if(includeInGlobal) {
      globalTotals[weightName][channelNumber] += hist->GetSumOfWeights();
      globalNumberOfEntries[weightName][channelNumber] += hist->GetEntries();
   }
   if(includeInMCRun) {
      periodTotals[weightName][channelNumber][mcRunNumber] += hist->GetSumOfWeights();
      //initialize the corresponding data value if necessary
      if(periodTotals[weightName][-1].find(mcRunNumber)==periodTotals[weightName][-1].end()) {
         periodTotals[weightName][-1][mcRunNumber]=0;
      }
   }
   periodTotals[weightName][channelNumber][periodNumber] += hist->GetSumOfWeights();

   //initialize the corresponding data value if necessary
   if(periodTotals[weightName][-1].find(periodNumber)==periodTotals[weightName][-1].end()) {
      periodTotals[weightName][-1][periodNumber]=0;
   }

   //add to the primary/secondary/tertiary histograms
   if(hist->InheritsFrom("TH3")) {
      //not supported yet
      Error("initialize","3d not supported yet");throw 40;
   } else if(hist->InheritsFrom("TH2")) {
      TH2* histo2d = dynamic_cast<TH2*>(hist);
      //project out and add that to primary
      TH1D* proj = histo2d->ProjectionX();
      if(!primaryDistributions[weightName][channelNumber][periodNumber]) {
         //need to make a copy
         TH1D* projCopy = dynamic_cast<TH1D*>(proj->Clone(proj->GetName()));
         projCopy->SetDirectory(0);
         primaryDistributions[weightName][channelNumber][periodNumber] = projCopy;
         //initialize the corresponding data histogram if necessary
         if(!primaryDistributions[weightName][-1][periodNumber]) {
            TH1D* projEmpty = dynamic_cast<TH1D*>(proj->Clone(proj->GetName()));
            projEmpty->Reset();
            projEmpty->SetDirectory(0);
            primaryDistributions[weightName][-1][periodNumber]=projEmpty;
         }
      } else {
         primaryDistributions[weightName][channelNumber][periodNumber]->Add(proj);
      }
      if(!secondaryDistributions[weightName][channelNumber][periodNumber]) {
         secondaryDistributions[weightName][channelNumber][periodNumber] = dynamic_cast<TH2D*>(CloneEmptyHistogram(weightName,periodNumber,channelNumber));
         //initialize the corresponding data histogram if necessary
         if(!secondaryDistributions[weightName][-1][periodNumber]) {
            secondaryDistributions[weightName][-1][periodNumber]=dynamic_cast<TH2D*>(CloneEmptyHistogram(weightName,periodNumber,-1));
         }
      }
      secondaryDistributions[weightName][channelNumber][periodNumber]->Add(hist);

      if(includeInMCRun) {
         if(!primaryDistributions[weightName][channelNumber][mcRunNumber]) {
            //need to make a copy for the run number, to avoid double deleting
            TH1D* projCopy = dynamic_cast<TH1D*>(proj->Clone(proj->GetName()));
            projCopy->SetDirectory(0);
            primaryDistributions[weightName][channelNumber][mcRunNumber] = projCopy;
            //initialize the corresponding data histogram if necessary
            if(!primaryDistributions[weightName][-1][mcRunNumber]) {
               TH1D* projEmpty = dynamic_cast<TH1D*>(proj->Clone(proj->GetName()));
               projEmpty->Reset();
               projEmpty->SetDirectory(0);
               primaryDistributions[weightName][-1][mcRunNumber]=projEmpty;
            }
         } else {
            primaryDistributions[weightName][channelNumber][mcRunNumber]->Add(proj);
         }
         if(!secondaryDistributions[weightName][channelNumber][mcRunNumber]) {
            secondaryDistributions[weightName][channelNumber][mcRunNumber] = dynamic_cast<TH2D*>(CloneEmptyHistogram(weightName,mcRunNumber,channelNumber));
            //initialize the corresponding data histogram if necessary
            if(!secondaryDistributions[weightName][-1][mcRunNumber]) {
               secondaryDistributions[weightName][-1][mcRunNumber]=dynamic_cast<TH2D*>(CloneEmptyHistogram(weightName,mcRunNumber,-1));
            }
         }
         secondaryDistributions[weightName][channelNumber][mcRunNumber]->Add(hist);
      }
      delete proj;

   } else {
      //regular 1D distribution
      if(!primaryDistributions[weightName][channelNumber][periodNumber]) {
         primaryDistributions[weightName][channelNumber][periodNumber] = dynamic_cast<TH1D*>(CloneEmptyHistogram(weightName,periodNumber,channelNumber));
         //initialize the corresponding data histogram if necessary
         if(!primaryDistributions[weightName][-1][periodNumber]) {
            primaryDistributions[weightName][-1][periodNumber]=dynamic_cast<TH1D*>(CloneEmptyHistogram(weightName,periodNumber,-1));
         }
      }
      primaryDistributions[weightName][channelNumber][periodNumber]->Add(hist);
      if(includeInMCRun) {
         if(!primaryDistributions[weightName][channelNumber][mcRunNumber]) {
            primaryDistributions[weightName][channelNumber][mcRunNumber] = dynamic_cast<TH1D*>(CloneEmptyHistogram(weightName,mcRunNumber,channelNumber));
            //initialize the corresponding data histogram if necessary
            if(!primaryDistributions[weightName][-1][mcRunNumber]) {
               primaryDistributions[weightName][-1][mcRunNumber]=dynamic_cast<TH1D*>(CloneEmptyHistogram(weightName,mcRunNumber,-1));
            }
         }
         primaryDistributions[weightName][channelNumber][mcRunNumber]->Add(hist);
      }
   }

   return 0;
}

//=============================================================================
// Initialization
//=============================================================================
Int_t Root::TPileupReweighting::initialize( const TString dataRootFileName,
                                          const TString dataRootHistName,
                                          const TString mcRootFileName,
                                          const TString mcRootHistName, int runNumber, int channel )
{
   if(dataRootFileName != "") {
      AddMCDistribution(mcRootFileName,mcRootHistName,runNumber,channel);
      AddDataDistribution(dataRootFileName,dataRootHistName,runNumber);
      if(m_periods.size()==0) {
         Info("initialize","No Periods have been defined. I'm guessing you are doing old-style reweighting");
         Info("initialize","Creating dummy period and ignoring unrepresented data issues");
         AddPeriod(0,0,0); //assign mcRunNumber 0 to data runNumber 0-0.
         SetUnrepresentedDataAction(2);
      }
   }
   return Initialize();
}

Int_t Root::TPileupReweighting::Initialize() {

   //Need to calculate a period weight for each period
   //Need a reweighting histogram for each period
   //for merged mc run numbers, we must generate or modify the existing mc distributions


   //1. Deal with the registered mc mergers
   std::map<Int_t, bool> newMCRuns; //records if new mcRunNumbers have been created by merging
   for(std::map<Int_t, std::map<Int_t, Int_t> >::iterator mergers = m_mcMerges.begin();mergers != m_mcMerges.end();++mergers) {
      for(std::map<Int_t,Int_t>::iterator mergePair = (mergers->second).begin();mergePair!=(mergers->second).end();++mergePair) {
         Int_t run1 = mergePair->first;
         Int_t run2 = mergePair->second;
         Int_t runOut = mergers->first;
         if(run1==runOut) {
            if(m_debugging) Info("Initialize","Merging %d into %d",run2,run1);
            //merging run2 directly in to run1.
            //just combine run2 and run1 mc histos, and clear empty the run2 histograms
            for(std::map<TString, std::map<Int_t, std::map<Int_t,TH1*> > >::iterator mcHistos = m_inputHistograms.begin();mcHistos!=m_inputHistograms.end();++mcHistos) {
               TString weightName = mcHistos->first;
               for(std::map<Int_t, std::map<Int_t,TH1*> >::iterator channels = mcHistos->second.begin();channels!=mcHistos->second.end();++channels) {
                  if((channels->first)<0) continue; //skip data channels
                  if(channels->second.find(run2)!=channels->second.end()) {
                     if(channels->second.find(run1)!=channels->second.end()) {
                        channels->second[run1]->Add(channels->second[run2]);
                     } else {
                        channels->second[run1] = CloneEmptyHistogram(weightName,run1,channels->first);
                        channels->second[run1]->Add(channels->second[run2]);
                     }
                     channels->second[run2]->Reset();
                  }
               }
            }
            //merge the periods that were assigned
            //combine the periods
            if(m_periods.find(run2)!=m_periods.end()) {
               for(std::map<Int_t,std::pair<UInt_t,UInt_t> >::iterator periods2 = m_periods[run2].begin();periods2!=m_periods[run2].end();++periods2) {
                  AddPeriod(run1,periods2->second.first,periods2->second.second);
               }
               m_periods.erase(run2);
            }
            //and remap run2 in to run1
            m_mcRemappings[run2]=run1;
         } else {
            //merging run1 and run2 in to a new run (run3).
            newMCRuns[runOut]=true;
            for(std::map<TString, std::map<Int_t, std::map<Int_t,TH1*> > >::iterator mcHistos = m_inputHistograms.begin();mcHistos!=m_inputHistograms.end();++mcHistos) {
               TString weightName = mcHistos->first;
               for(std::map<Int_t, std::map<Int_t,TH1*> >::iterator channels = mcHistos->second.begin();channels!=mcHistos->second.end();++channels) {
                  if((channels->first)<0) continue;
                  if(channels->second.find(run2)!=channels->second.end()) {
                     if(channels->second.find(runOut)!=channels->second.end()) {
                        channels->second[runOut]->Add(channels->second[run2]);
                     } else {
                        channels->second[runOut] = CloneEmptyHistogram(weightName,runOut,channels->first);
                        channels->second[runOut]->Add(channels->second[run2]);
                     }
                  }
                  if(channels->second.find(run1)!=channels->second.end()) {
                     if(channels->second.find(runOut)!=channels->second.end()) {
                        channels->second[runOut]->Add(channels->second[run1]);
                     } else {
                        channels->second[runOut] = CloneEmptyHistogram(weightName,runOut,channels->first);
                        channels->second[runOut]->Add(channels->second[run1]);
                     }
                  }
               }
            }
            //merge the periods that were assigned, but do not clear the period assignment
            if(m_periods.find(run2)!=m_periods.end()) {
               for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods2 = m_periods[run2].begin();periods2!=m_periods[run2].end();++periods2) {
                  AddPeriod(runOut,periods2->second.first,periods2->second.second);
               }
            }
            if(m_periods.find(run1)!=m_periods.end()) {
               for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods2 = m_periods[run1].begin();periods2!=m_periods[run1].end();++periods2) {
                  AddPeriod(runOut,periods2->second.first,periods2->second.second);
               }
            }
         }
      }
   }

   //print out the period assignments if in debug mode
   if(m_debugging) {
      for(std::map<Int_t,std::map<Int_t, std::pair<UInt_t,UInt_t> > >::iterator mcruns = m_periods.begin();mcruns!=m_periods.end();++mcruns) {
         Info("Initialize","MCRunNumber %d has following periods:",mcruns->first);
         for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods = mcruns->second.begin();periods!=mcruns->second.end();++periods) {
            Info("Initialize","  PeriodNumber %d: %d - %d",periods->first,periods->second.first,periods->second.second);
         }
      }
   }

   //Factorize the MC distributions to build up the totals and primary/seconday/tertiary distributions
   for(std::map<TString, std::map<Int_t, std::map<Int_t,TH1*> > >::iterator mcHistos = m_inputHistograms.begin();mcHistos!=m_inputHistograms.end();++mcHistos) {
      TString weightName = mcHistos->first;
      for(std::map<Int_t, std::map<Int_t,TH1*> >::iterator channels = mcHistos->second.begin();channels!=mcHistos->second.end();++channels) {
         Int_t channel = channels->first;
         if(channel<0) continue;
         for(std::map<Int_t,TH1*>::iterator mcruns = channels->second.begin();mcruns!=channels->second.end();++mcruns) {
            Int_t mcRunNum = mcruns->first;
            //all the periods associated to the mcrun will get the same distribution
            //only include the first subperiod in the totals and mcrunnum's distributions
            //and only include genuine mcrunnums (not ones generated in mergers) in global totals
            bool isFirst = true;
            for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator period = m_periods[mcRunNum].begin(); period!=m_periods[mcRunNum].end();++period) {
               FactorizeDistribution(mcruns->second,weightName,channel,period->first,isFirst,isFirst&&(!newMCRuns[mcRunNum]));
               isFirst=false;
            }
         }
      }
   }

   std::map<TString,Double_t> unrepData;

   for(std::map<Int_t, std::map<Int_t,std::pair<UInt_t,UInt_t> > >::iterator runPeriods = m_periods.begin();runPeriods != m_periods.end();++runPeriods) {
      Int_t thisMCRunNumber = runPeriods->first;
      for(std::map<Int_t,std::pair<UInt_t,UInt_t> >::iterator period = runPeriods->second.begin();period!=runPeriods->second.end();++period) {
         Int_t periodNumber = period->first;
         UInt_t pStart = period->second.first; UInt_t pEnd = period->second.second;
         //for(std::map<TString,std::map<Int_t,TH1*> >::iterator dataHists = m_inputDataHistograms.begin();dataHists!=m_inputDataHistograms.end();++dataHists) {
         for(std::map<TString,std::map<Int_t,std::map<Int_t,TH1*> > >::iterator dataHists = m_inputHistograms.begin();dataHists!=m_inputHistograms.end();++dataHists) {
            TString weightName = dataHists->first;
            //get list of loaded runs associated to this period
            std::vector<UInt_t> associatedRuns;
            //std::map<Int_t,TH1*> hists = m_inputDataHistograms[weightName];
            std::map<Int_t,TH1*> hists = (dataHists->second)[-1];
            //look for loaded runs
            for(UInt_t run=pStart;run<=pEnd;run++) {
               if(hists.find(run)!=hists.end() && hists[run]) {
                  associatedRuns.push_back(run);
               }
            }
            //loop over data runs in this period and build up the data distributions
            for(std::vector<UInt_t>::iterator run = associatedRuns.begin();run!=associatedRuns.end();++run) {
               TH1* thisHist = hists[*run];
               //look for unrepresented data in any of the bins of any of the corresponding channels
               for(std::map<Int_t,std::map<Int_t,TH1*> >::iterator channels=m_inputHistograms[weightName].begin();channels!=m_inputHistograms[weightName].end();++channels) {
                  if((channels->first)<0) continue; //skip data channels
                  if((channels->second)[thisMCRunNumber]) {
                     Int_t bin,binx,biny,binz;
                     for(binz=1; binz<=thisHist->GetNbinsZ(); binz++) {
                        for(biny=1; biny<=thisHist->GetNbinsY(); biny++) {
                           for(binx=1; binx<=thisHist->GetNbinsX(); binx++) {
                              bin = thisHist->GetBin(binx,biny,binz);
                              Double_t value = thisHist->GetBinContent(bin);
                              if(value==0.) continue;
                              Double_t mcValue = (channels->second)[thisMCRunNumber]->GetBinContent(bin);
                              if(mcValue==0. && !m_badbins[weightName][*run][bin]) { //we don't want to double-count as we loop over channels
                                 //unrepresented data bin
                                 m_badbins[weightName][*run][bin]=true;
                                 unrepData[weightName] += value;
                                 if(m_unrepresentedDataAction==0) {
                                    if(m_debugging) Error("Initialize","Unrepresented data at coords [%f,%f,%f] caused by mcrun %d in channel %d",thisHist->GetXaxis()->GetBinCenter(binx),thisHist->GetYaxis()->GetBinCenter(biny),thisHist->GetZaxis()->GetBinCenter(binz),thisMCRunNumber,channels->first);
                                 }
                                 else if(m_unrepresentedDataAction==1) thisHist->SetBinContent(bin,0);
                              }
                           }
                        }
                     }
                  }
               }
                  dataPeriodRunTotals[weightName][periodNumber][*run] += thisHist->GetSumOfWeights(); //for random runnumber
                  dataPeriodRunTotals[weightName][thisMCRunNumber][*run] += thisHist->GetSumOfWeights(); //for random runnumber
                  //don't double-count the runs. we track the run totals in dataPeriodRunTotals["pileup"][-1]
                  if(weightName=="pileup" && dataPeriodRunTotals["pileup"][-1].find(*run)==dataPeriodRunTotals["pileup"][-1].end()) {
                     dataPeriodRunTotals["pileup"][-1][*run]=thisHist->GetSumOfWeights();
                     int lumiVectorPosition = -1;
                     int runNbr = *run;
                     if(runNbr>=177986&&runNbr<=178109) lumiVectorPosition=0;//period B
                     else if(runNbr>=179710&&runNbr<=180481) lumiVectorPosition=1;//period D
                     else if(runNbr>=180614&&runNbr<=180776) lumiVectorPosition=2;//period E
                     else if(runNbr>=182013&&runNbr<=182519) lumiVectorPosition=3;//period F
                     else if(runNbr>=182726&&runNbr<=183462) lumiVectorPosition=4;//period G
                     else if(runNbr>=183544&&runNbr<=184169) lumiVectorPosition=5;//period H
                     else if(runNbr>=185353&&runNbr<=186493) lumiVectorPosition=6;//period I
                     else if(runNbr>=186516&&runNbr<=186755) lumiVectorPosition=7;//period J
                     else if(runNbr>=186873&&runNbr<=187815) lumiVectorPosition=8;//period K
                     else if(runNbr>=188902&&runNbr<=190343) lumiVectorPosition=9; //period L
                     else if(runNbr>=190503&&runNbr<=191933) lumiVectorPosition=10; //period M
                     if(lumiVectorPosition>=0) { m_integratedLumiVector[lumiVectorPosition]+=(thisHist->GetSumOfWeights())/1E6; }
                  }
                  //factorize out the distribution and add it to existing distributions
                  //include this distribution in the mcRunNumber distribution and in the global total (only if mcRunNumber is a genuine run num. Don't want to double-count)
                  FactorizeDistribution(thisHist,weightName,-1,periodNumber,true,!newMCRuns[thisMCRunNumber]);
            } //end loop over data in period
         } //loop over weights
      } //loop over periods in mcrun
   } //loop over mcruns

   //we should also check if there are any data runs that are not covered by the period assignments
   for(std::map<TString, std::map<Int_t,std::map<Int_t, TH1*> > >::iterator weights=m_inputHistograms.begin();weights!=m_inputHistograms.end();++weights) {
      Double_t ignoredData = 0;
      if(weights->second.find(-1)==weights->second.end()) continue; //no data
      for(std::map<Int_t,TH1*>::iterator runs = weights->second[-1].begin();runs!=weights->second[-1].end();++runs) {
         //look for period assignment for run
         Int_t periodNum=-1;
         for(std::map<Int_t, std::map<Int_t, std::pair<UInt_t,UInt_t> > >::iterator mcruns = m_periods.begin();mcruns!=m_periods.end();++mcruns) {
            for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods=mcruns->second.begin();periods!=mcruns->second.end();++periods) {
               if((UInt_t)runs->first >= periods->second.first && (UInt_t)runs->first <= periods->second.second) {
                  periodNum = periods->first;break;
               }
            }
            if(periodNum!=-1) break;
         }
         if(periodNum==-1) {
            //this run isn't assigned to a period. Give warning
            Warning("Initialize","%s weight contains loaded data in run %d that is not covered by period assignments",weights->first.Data(),runs->first);
            ignoredData += (runs->second)->GetSumOfWeights();
         }
      }
      if(ignoredData>0.) Warning("Initialize", "Period Assignments missed %f data in %s weight",ignoredData,weights->first.Data());
   }

   bool displayOptions = false;
   for(std::map<TString,Double_t>::iterator it=unrepData.begin();it!=unrepData.end();++it) {
      if(it->second>0.) {
           if(m_unrepresentedDataAction==0) {
               Error("Initialize","%s weight has %f%% unrepresented data",it->first.Data(),100.*it->second/globalTotals[it->first][-1]);
               //sweep back over the bad bins and add up the unrepresented data grouped by channel
               std::map<Int_t,Double_t> unrepChannelTotals;
               for(std::map<Int_t, std::map<Int_t, Bool_t> >::iterator badRuns=m_badbins[it->first].begin();badRuns!=m_badbins[it->first].end();++badRuns) {
                  //find out which mcrunNum this run belongs to
                  Int_t mcrunNum=-1;
                  for(std::map<Int_t, std::map<Int_t, std::pair<UInt_t,UInt_t> > >::iterator mcruns = m_periods.begin();mcruns!=m_periods.end();++mcruns) {
                     for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::iterator periods=mcruns->second.begin();periods!=mcruns->second.end();++periods) {
                        if((UInt_t)badRuns->first >= periods->second.first && (UInt_t)badRuns->first <= periods->second.second) {
                           mcrunNum = mcruns->first;break;
                        }
                     }
                     if(mcrunNum!=-1) break;
                  }
                  //now sweep over mc channels and check this periodNum for the bad bins
                  for(std::map<Int_t, std::map<Int_t, TH1*> >::iterator channels=m_inputHistograms[it->first].begin();channels!=m_inputHistograms[it->first].end();++channels) {
                     if(channels->first < 0) continue;
                     for(std::map<Int_t, Bool_t>::iterator bins=badRuns->second.begin();bins!=badRuns->second.end();++bins) {
                        if(bins->second && (channels->second)[mcrunNum]->GetBinContent(bins->first)==0.) {
                           //found one of the bad bins.
                           unrepChannelTotals[channels->first] += m_inputHistograms[it->first][-1][badRuns->first]->GetBinContent(bins->first);
                        }
                     }
                  }
               }
               //display the results of the scan
               for(std::map<Int_t,Double_t>::iterator badChans = unrepChannelTotals.begin();badChans!=unrepChannelTotals.end();++badChans) {
                  Error("Initialize","   Channel %d has %f%% of the unrepresented",badChans->first,badChans->second*100./it->second);
               }
               displayOptions=true;
           } else if(m_unrepresentedDataAction==1) {
               Warning("Initialize","%s weight has %f%% unrepresented data. This was removed (UnrepresentedDataAction=1)",it->first.Data(),100.*it->second/(globalTotals[it->first][-1]+it->second));
           } else if(m_unrepresentedDataAction==2) {
               Warning("Initialize","%s weight has %f%% unrepresented data. This was kept in (UnrepresentedDataAction=2)",it->first.Data(),100.*it->second/(globalTotals[it->first][-1]));
           }
      }
   }

   if(displayOptions) {
      Error("Initialize","Exiting. You must decide how to proceed...");
      Error("Initialize","1) use AddPeriod or MergeMCRunNumbers to redefine the mc periods to include this data. You should not need to regenerate your mc config file");
      Error("Initialize","2) use SetUnrepresentedDataAction(1) to remove this data from the weight calculations. You should also veto such data events (using IsUnrepresentedData(..,..) method)");
      Error("Initialize","3) use SetUnrepresentedDataAction(2) to leave this data in the calculation. I hope you know what you're doing!!");
      throw std::runtime_error("Throwing exception 22");
   }

   if(m_debugging) Info("Initialize","Normalizing 1D histograms...");
   //now that all the distributions are built. Normalize them all
   for(std::map<TString, std::map<Int_t, std::map<Int_t, TH1D*> > >::iterator weights=primaryDistributions.begin();weights!=primaryDistributions.end();++weights) {
      for(std::map<Int_t, std::map<Int_t, TH1D*> >::iterator channels=weights->second.begin();channels!=weights->second.end();++channels) {
         for(std::map<Int_t,TH1D*>::iterator mcrunsAndPeriods=channels->second.begin();mcrunsAndPeriods!=channels->second.end();++mcrunsAndPeriods) {
            normalizeHistogram(mcrunsAndPeriods->second);
         }
      }
   }
   if(m_debugging) Info("Initialize","Normalizing 2D histograms...");
   for(std::map<TString, std::map<Int_t, std::map<Int_t, TH2D*> > >::iterator weights=secondaryDistributions.begin();weights!=secondaryDistributions.end();++weights) {
      for(std::map<Int_t, std::map<Int_t, TH2D*> >::iterator channels=weights->second.begin();channels!=weights->second.end();++channels) {
         for(std::map<Int_t,TH2D*>::iterator mcrunsAndPeriods=channels->second.begin();mcrunsAndPeriods!=channels->second.end();++mcrunsAndPeriods) {
            normalizeHistogram(mcrunsAndPeriods->second);
         }
      }
   }

   if(m_debugging) Info("Initialize","Deleting input histograms...");
   //to save space, delete all the input histograms
   for(std::map<TString, std::map<Int_t,std::map<Int_t, TH1*> > >::iterator it0=m_inputHistograms.begin();it0!=m_inputHistograms.end();++it0) {
      for(std::map<Int_t,std::map<Int_t, TH1*> >::iterator it1=it0->second.begin();it1!=it0->second.end();++it1) {
         for(std::map<Int_t, TH1*>::iterator it2=it1->second.begin();it2!=it1->second.end();++it2) {
            delete it2->second;
         }
         it1->second.clear();
      }
      it0->second.clear();
   }
   m_inputHistograms.clear();
   if(m_debugging) Info("Initialize","...Done");
   //if no input histograms were added, we are in counting mode
   if(m_countingMode) {
      Info("Initialize","In Config File Generating mode. Remember to call WriteToFile!");
     //histograms are instantiated in the event loop as the channels occur
      //but since we are here, check that the user has at least defined some periods
      if(m_periods.size()==0) {
         Error("Initialize", "In Config File Generating mode, but no periods have been defined. This makes no sense!? Define some periods please (with AddPeriod method)!");
         throw std::runtime_error("Throwing 43");
      }
   }

   m_isInitialized=true;

  return 0;
}

Bool_t Root::TPileupReweighting::IsUnrepresentedData(const TString weightName, Int_t runNumber, Float_t x, Float_t y, Float_t z) {
   if(m_emptyHistograms.find(weightName)==m_emptyHistograms.end()) {
      Error("IsUnrepresentedData", "Unknown weight %s", weightName.Data());
      throw std::runtime_error("Throwing 23");
   }
   TH1* hist = m_emptyHistograms[weightName];
   Int_t bin = hist->GetXaxis()->FindFixBin(x);
   if(hist->GetDimension()==2) {
      Int_t biny = hist->GetYaxis()->FindFixBin(y);
      Int_t nx = hist->GetXaxis()->GetNbins()+2;
      bin += nx*biny;
   }
   if(hist->GetDimension()==3) {
      Int_t biny = hist->GetYaxis()->FindFixBin(y);
      Int_t nx = hist->GetXaxis()->GetNbins()+2;
      Int_t ny = hist->GetYaxis()->GetNbins()+2;
      Int_t binz = hist->GetZaxis()->FindFixBin(z);
      bin += nx*(biny+ny*binz);
   }
   return m_badbins[weightName][runNumber][bin];
}

Bool_t Root::TPileupReweighting::IsUnrepresentedData(Int_t runNumber, Float_t x, Float_t y, Float_t z) {
   return IsUnrepresentedData("pileup",runNumber,x,y,z);
}

Int_t Root::TPileupReweighting::WriteToFile(const TString fname) {

   if(!m_countingMode) {Warning("WriteToFile","Not in counting mode, so no file will be written");return 0;}

   TDirectory* origDir = gDirectory;

   //build a TTree with the correct mc structure, and dump the mc histogram info in to it
   //also build aggregate histograms across all channels - shouldn't be used as input histograms unless for a single channel
   TString filename = fname;
   filename += (filename=="") ? TString(this->GetName()) + ".prw.root" : "";

   //loop over the weights. "pileup" gets it own MCPileupReweighting ttree.
   //data goes in to DataCustomReweighting. other goes in to MCCustomReweighting

   TFile outFile(filename,"RECREATE");
   TTree *outTreeMC=0;
   TTree *outTreeCustomMC=0;
   TTree *outTreeData=0;
   TTree *outTreeCustomData=0;
   TTree *outTree = 0;//points to tree currently being written to
   Int_t channel = 0;UInt_t runNumber = 0;
   std::vector<UInt_t>* pStarts = 0;std::vector<UInt_t>* pEnds = 0;
   Char_t histName[150];Char_t customName[150];

   bool lastChannelWasData = true;

   for(std::map<TString, std::map<Int_t, std::map<Int_t, TH1*> > >::iterator weights=m_inputHistograms.begin();weights!=m_inputHistograms.end();++weights) {
      TString weightName = weights->first;
      strcpy(customName,weightName);
      if(outTree) {
         outTree=0;
      }
      for(std::map<Int_t, std::map<Int_t, TH1*> >::iterator channels=weights->second.begin();channels!=weights->second.end();++channels) {
         channel = channels->first;
         if(lastChannelWasData&&channel>=0) {
            if(outTree) {outTree=0;}
            lastChannelWasData=false;
         }
         if(channel<0) lastChannelWasData=true;
         if(weightName=="pileup") {
            if(channel>=0 && !outTree) {
               if(!outTreeMC) {
                  outTreeMC = new TTree("MCPileupReweighting","MCPileupReweighting");
                  outTreeMC->Branch("Channel",&channel);
                  outTreeMC->Branch("RunNumber",&runNumber);
                  outTreeMC->Branch("PeriodStarts",&pStarts);
                  outTreeMC->Branch("PeriodEnds",&pEnds);
                  outTreeMC->Branch("HistName",&histName,"HistName/C");
               }
               outTree = outTreeMC;
            } else if(channel<0 && !outTree) {
               if(!outTreeData) {
                  outTreeData = new TTree("DataPileupReweighting","DataPileupReweighting");
                  outTreeData->Branch("RunNumber",&runNumber);
                  outTreeData->Branch("HistName",&histName,"HistName/C");
               }
               outTree = outTreeData;
            }
         } else {
            if(channel>=0 && !outTree) {
               if(!outTreeCustomMC) {
                  outTreeCustomMC = new TTree("MCCustomReweighting","MCCustomReweighting");
                  outTreeCustomMC->Branch("CustomName",&customName,"CustomName/C");
                  outTreeCustomMC->Branch("Channel",&channel);
                  outTreeCustomMC->Branch("RunNumber",&runNumber);
                  outTreeCustomMC->Branch("PeriodStarts",&pStarts);
                  outTreeCustomMC->Branch("PeriodEnds",&pEnds);
                  outTreeCustomMC->Branch("HistName",&histName,"HistName/C");
               }
               outTree = outTreeCustomMC;
            } else if(channel<0 && !outTree) {
               if(!outTreeCustomData) {
                  outTreeCustomData = new TTree("DataCustomReweighting","DataCustomReweighting");
                  outTreeCustomData->Branch("CustomName",&customName,"CustomName/C");
                  outTreeCustomData->Branch("RunNumber",&runNumber);
                  outTreeCustomData->Branch("HistName",&histName,"HistName/C");
               }
               outTree = outTreeCustomData;
            }
         }
         for(std::map<Int_t, TH1*>::iterator runs=channels->second.begin();runs!=channels->second.end();++runs) {
            runNumber = runs->first;
            if(!lastChannelWasData) { //load periods only for mc
               pStarts = new std::vector<UInt_t>;pEnds = new std::vector<UInt_t>;
               for(std::map<Int_t,std::pair<UInt_t,UInt_t> >::iterator periods=m_periods[runNumber].begin();periods!=m_periods[runNumber].end();++periods) {
                  pStarts->push_back(periods->second.first);pEnds->push_back(periods->second.second);
               }
            }
            strcpy(histName,runs->second->GetName());
            outTree->Fill();
            runs->second->Write();
            if(!lastChannelWasData) { //delete periods only for mc
               delete pStarts;delete pEnds;
            }
         }
      }
   }
   //write the non-zero ttrees
   if(outTreeMC) {outTreeMC->Write();delete outTreeMC;outTreeMC=0;}
   if(outTreeData) {outTreeData->Write();delete outTreeData;outTreeData=0;}
   if(outTreeCustomMC) {outTreeCustomMC->Write();delete outTreeCustomMC;outTreeCustomMC=0;}
   if(outTreeCustomData) {outTreeCustomData->Write();delete outTreeCustomData;outTreeCustomData=0;}

   outFile.Close();

   Info("WriteToFile", "Successfully generated config file: %s",filename.Data());
   Info("WriteToFile", "Happy Reweighting :-)");

   gDirectory = origDir;

   return 0;
}

UInt_t Root::TPileupReweighting::GetRandomRunNumber(Int_t mcRunNumber) {
   if(m_countingMode) { return 0; } //do nothing when in counting mode
   //check if the runNumber has been reassigned
   if(m_mcRemappings.find(mcRunNumber)!=m_mcRemappings.end()) mcRunNumber=m_mcRemappings[mcRunNumber];
   const std::map<Int_t, Double_t>::const_iterator it = periodTotals["pileup"][-1].find(mcRunNumber);
   if(it==periodTotals["pileup"][-1].end()|| it->second == 0) {
      if(m_SetWarnings) Warning("GetRandomRunNumber","%d has no associated lumi. Returning 0",mcRunNumber);
      return 0; /*throw 33;*/
   }
   double lumi = it->second * m_random3->Rndm();
    double lumisum = 0.;
    std::map<UInt_t, Double_t>::const_iterator itend = dataPeriodRunTotals["pileup"][mcRunNumber].end();
    for(std::map<UInt_t, Double_t>::const_iterator it2 = dataPeriodRunTotals["pileup"][mcRunNumber].begin(); it2 != itend; ++it2)
      {
	lumisum += it2->second;
	if (lumisum >= lumi)
	  return it2->first;
      }
    Error("GetRandomRunNumber","overran integrated luminosity for mcRunNumber=%d",mcRunNumber);
     throw std::runtime_error("Throwing 46");
    return 0;
}

Int_t Root::TPileupReweighting::GetRandomPeriodNumber(Int_t mcRunNumber) {
   if(m_countingMode) { return 0; } //do nothing when in counting mode
   //check if the runNumber has been reassigned
   if(m_mcRemappings.find(mcRunNumber)!=m_mcRemappings.end()) mcRunNumber=m_mcRemappings[mcRunNumber];

   if(m_periods.find(mcRunNumber)==m_periods.end()) {
      Error("GetRandomPeriodNumber","Unknown MC RunNumber: %d",mcRunNumber);
      throw std::runtime_error("Throwing 46");
   }

   //if this runnumber has only one period, return itself
   if(m_periods[mcRunNumber].size()==1) {
      return mcRunNumber;
   }

   const std::map<Int_t, Double_t>::const_iterator it = periodTotals["pileup"][-1].find(mcRunNumber);
   if(it==periodTotals["pileup"][-1].end()) {
      if(m_SetWarnings) Warning("GetRandomPeriodNumber","%d has no associated lumi. Returning 0",mcRunNumber);
      return 0; /*throw 33;*/
   }
   Double_t lumi = it->second * m_random3->Rndm();
   Double_t lumisum = 0.;
   for(std::map<Int_t, std::pair<UInt_t,UInt_t> >::const_iterator it2 = m_periods[mcRunNumber].begin(); it2 != m_periods[mcRunNumber].end(); ++it2) {
      lumisum += periodTotals["pileup"][-1][it2->first];
      if (lumisum >= lumi)
	  return it2->first;
      }
    Error("GetRandomPeriodNumber","overran integrated luminosity for mcRunNumber=%d",mcRunNumber);
     throw std::runtime_error("Throwing 46");
    return 0;
}

Float_t Root::TPileupReweighting::getPileupWeight(Float_t mu, Int_t runNumber, Int_t channel ) {
   if(m_SetWarnings) Warning("getPileupWeight","Depricated. Please use GetCombinedWeight(runNumber,channelNumber,mu)");
   return GetCombinedWeight(runNumber,channel,mu);
}

Float_t Root::TPileupReweighting::GetCombinedWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber, Float_t x, Float_t y, Float_t /*z*/) {
   if(m_countingMode) return 0.;
   if(m_countingMode) {
      Error("GetCombinedWeight","You cannot get a weight when in config file generating mode. Please use FillWeight method");
      throw std::runtime_error("Throwing 2");
   }
   //decide how many dimensions this weight has - use the emptyHistogram to tell...
   TH1* hist = m_emptyHistograms[weightName];
   if(!hist) {
      Error("GetCombinedWeight","Weight not recognised: %s",weightName.Data());
   }
   Float_t out = GetPeriodWeight(weightName,periodNumber,channelNumber)*GetPrimaryWeight(weightName,periodNumber,channelNumber,x);
   if(hist->GetDimension()>1) out *= GetSecondaryWeight(weightName,periodNumber,channelNumber,x,y);
   //if(hist->GetDimension()>2) out *= GetTertiaryWeight(weightName,periodNumber,channelNumber,x,y,z);

   return out;
}

Float_t Root::TPileupReweighting::GetPeriodWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber) {
   if(m_countingMode) return 0.;
   //check if the runNumber has been reassigned
   //only do this for mc, never for data
   if(channelNumber>=0) {
      if(m_mcRemappings.find(periodNumber)!=m_mcRemappings.end()) periodNumber=m_mcRemappings[periodNumber];
   }

   if(periodTotals[weightName][-1].find(periodNumber)==periodTotals[weightName][-1].end()) {
      Error("GetPeriodWeight","Unrecognised periodNumber/mcRunNumber: %d",periodNumber);
      throw std::runtime_error("Throwing 1");
   }

   Double_t data = periodTotals[weightName][-1][periodNumber];
   Double_t totalData = globalTotals[weightName][-1];

   if(totalData==0.) {
      Error("GetPeriodWeight","No data for weight %s",weightName.Data());
      throw std::runtime_error("Throwing 2");
   }

   if(periodTotals[weightName].find(channelNumber)==periodTotals[weightName].end()) {
      //use the default
      if(periodTotals[weightName].find(m_defaultChannel)==periodTotals[weightName].end()) {
         Error("GetPeriodWeight","Unrecognised channel %d and no default found: %d",channelNumber,m_defaultChannel);
         throw std::runtime_error("Throwing 3");
      }
      channelNumber=m_defaultChannel;
   }
   if(periodTotals[weightName][channelNumber].find(periodNumber)==periodTotals[weightName][channelNumber].end()) {
      Error("GetPeriodWeight","Unrecognised periodNumber/mcRunNumber %d in channelNumber %d",periodNumber,channelNumber);
      throw std::runtime_error("Throwing 1");
   }

   Double_t mc = periodTotals[weightName][channelNumber][periodNumber];
   Double_t totalMC = globalTotals[weightName][channelNumber];

   if(totalMC==0.) {
      if(m_SetWarnings) Warning("GetPeriodWeight","No mc for weight %s in channelNumber %d",weightName.Data(),channelNumber);
      return 0.;
   }
   if(mc==0.) {
      if(m_SetWarnings) Warning("GetPeriodWeight","No mc for weight %s in channelNumber %d and periodNumber/runNumber %d",weightName.Data(),channelNumber,periodNumber);
      return 0.;
   }

   return (data/totalData) / (mc/totalMC) ;

}

Float_t Root::TPileupReweighting::GetPrimaryWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber, Float_t x) {
   if(m_countingMode) return 0.;
   //check if the runNumber has been reassigned
   //only do this for mc, never for data
   if(channelNumber>=0) {
      if(m_mcRemappings.find(periodNumber)!=m_mcRemappings.end()) periodNumber=m_mcRemappings[periodNumber];
   }

   //determine bin number
   TH1D* data = primaryDistributions[weightName][-1][periodNumber];
   if(!data) {
      Error("GetPrimaryWeight","Unrecognised periodNumber/mcRunNumber: %d",periodNumber);
      throw std::runtime_error("Throwing 4");
   }
   if(primaryDistributions[weightName].find(channelNumber)==primaryDistributions[weightName].end()) {
      //try switching to default channel
      if(primaryDistributions[weightName].find(m_defaultChannel)==primaryDistributions[weightName].end()) {
         Error("GetPrimaryWeight","Unrecognised channelNumber %d, and no distribution found for defaultChannel %d",channelNumber,m_defaultChannel);
         throw std::runtime_error("Throwing 5");
      }
      channelNumber=m_defaultChannel;
   }
   TH1D* mc = primaryDistributions[weightName][channelNumber][periodNumber];
   if(!mc) {
      Error("GetPrimaryWeight","Unable to find distribution for [%s,%d,%d]",weightName.Data(),channelNumber,periodNumber);
      throw std::runtime_error("Throwing 6");
   }

   //MC10b correction
   if ( fabs(x - 100.0) < 0.00001 && mc->GetXaxis()->GetXmax() < 99.0 ) x = 0.0;

   Int_t bin = data->GetXaxis()->FindFixBin(x);

   if(mc->GetBinContent(bin)==0.) {
      if(m_SetWarnings) Warning("GetPrimaryWeight","No mc for weight %s in channelNumber %d and periodNumber/runNumber %d",weightName.Data(),channelNumber,periodNumber);
      return 0.;
   }

   return (data->GetBinContent(bin)/mc->GetBinContent(bin));
}

Float_t Root::TPileupReweighting::GetSecondaryWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber,Float_t x, Float_t y) {
   if(m_countingMode) return 0.;
   //check if the runNumber has been reassigned
   //only do this for mc, never for data
   if(channelNumber>=0) {
      if(m_mcRemappings.find(periodNumber)!=m_mcRemappings.end()) periodNumber=m_mcRemappings[periodNumber];
   }

   //determine bin number
   TH2D* data = secondaryDistributions[weightName][-1][periodNumber];
   if(!data) {
      Error("GetSecondaryWeight","Unrecognised periodNumber/mcRunNumber: %d",periodNumber);
      throw std::runtime_error("Throwing 4");
   }
   if(secondaryDistributions[weightName].find(channelNumber)==secondaryDistributions[weightName].end()) {
      //try switching to default channel
      if(secondaryDistributions[weightName].find(m_defaultChannel)==secondaryDistributions[weightName].end()) {
         Error("GetSecondaryWeight","Unrecognised channelNumber %d, and no distribution found for defaultChannel %d",channelNumber,m_defaultChannel);
         throw std::runtime_error("Throwing 5");
      }
      channelNumber=m_defaultChannel;
   }
   TH2D* mc = secondaryDistributions[weightName][channelNumber][periodNumber];
   if(!mc) {
      Error("GetSecondaryWeight","Unable to find distribution for [%s,%d,%d]",weightName.Data(),channelNumber,periodNumber);
      throw std::runtime_error("Throwing 6");
   }
   Int_t nx = data->GetXaxis()->GetNbins()+2;
   Int_t binx = data->GetXaxis()->FindFixBin(x);
   Int_t biny = data->GetYaxis()->FindFixBin(y);
   Int_t bin = binx + nx*biny;

   if(mc->GetBinContent(bin)==0.) {
      if(m_SetWarnings) Warning("GetSecondaryWeight","No mc for weight %s in channelNumber %d and periodNumber/runNumber %d",weightName.Data(),channelNumber,periodNumber);
      return 0.;
   }

   return (data->GetBinContent(bin)/mc->GetBinContent(bin));

}

Int_t Root::TPileupReweighting::Fill(Int_t runNumber,Int_t channelNumber,Float_t w, Float_t x, Float_t y, Float_t z) {
   return Fill("pileup",runNumber,channelNumber,w,x,y,z);
}

//fills the appropriate inputHistograms
Int_t Root::TPileupReweighting::Fill(const TString weightName,Int_t runNumber,Int_t channelNumber,Float_t w,Float_t x, Float_t y, Float_t z) {

   //should only be given genuine mcRunNumbers if mc (channel>=0). We don't fill periodNumber distributions
   if(channelNumber>=0) {
      if(m_mcRemappings.find(runNumber)!=m_mcRemappings.end()) runNumber=m_mcRemappings[runNumber];

      //also check that channel histograms exist for all defined periods. If not, define them
      if(m_inputHistograms[weightName].find(channelNumber)==m_inputHistograms[weightName].end()) {
         //loop over mcruns, get empty histogram for each
         for(std::map<Int_t,std::map<Int_t, std::pair<UInt_t,UInt_t> > >::iterator it=m_periods.begin();it!=m_periods.end();++it) {
            m_inputHistograms[weightName][channelNumber][it->first] = CloneEmptyHistogram(weightName,it->first,channelNumber);
         }
      }
   } else {
      //create histograms for data too, only for this particular run
      if(m_inputHistograms[weightName][channelNumber].find(runNumber)==m_inputHistograms[weightName][channelNumber].end()) {
         m_inputHistograms[weightName][channelNumber][runNumber] = CloneEmptyHistogram(weightName,runNumber,channelNumber);
      }
   }

   TH1* hist = m_inputHistograms[weightName][channelNumber][runNumber];

   if(!hist) {
      Error("Fill","Unknown [weight,run,channel] = [%s,%d,%d]",weightName.Data(),runNumber,channelNumber);
      throw std::runtime_error("Throwing 45");
   }

   if(hist->GetDimension()==1) {
      return hist->Fill(x,w);
   } else if(hist->GetDimension()==2) {
      return (dynamic_cast<TH2*>(hist))->Fill(x,y,w);
   } else if(hist->GetDimension()==3) {
      return (dynamic_cast<TH3*>(hist))->Fill(x,y,z,w);
   }
   return -1;
}

//=============================================================================
// Method to calculate the ratio histogram
//=============================================================================
void Root::TPileupReweighting::normalizeHistogram(TH1* hist){
   // normalize the data histogram based on which sort of histogram it is
   if(hist){
      if(hist->InheritsFrom("TH3")) {
         Error("normalizeHistogram","3D reweighting not supported yet");
         throw std::runtime_error("Throwing 3");
      }
      else if(hist->InheritsFrom("TH2")) {
         bool skipNorm=false;
         //normalize each bin according to the projection in x
         TH1D* proj = ((TH2*)hist)->ProjectionX();
         Int_t bin,binx,biny,binz;
         for(binz=1; binz<=hist->GetNbinsZ(); binz++) {
            for(biny=1; biny<=hist->GetNbinsY(); biny++) {
               for(binx=1; binx<=hist->GetNbinsX(); binx++) {
                  bin = hist->GetBin(binx,biny,binz);
                  Double_t value = hist->GetBinContent(bin);
                  Double_t normalizer = proj->GetBinContent(binx);
                  if(fabs(normalizer)>0.00001) {
                     hist->SetBinContent(bin,value/normalizer);
                  } else {
                     skipNorm=true;
                  }
               }
            }
         }
         delete proj;
         if(skipNorm && m_debugging) Warning("normalizeHistogram","Skipped normalization in hist %s",hist->GetName());
      } else {
         //normalize to the sum of weights
         if(fabs(hist->GetSumOfWeights())>0.00001) {
            hist->Scale(1.0/hist->GetSumOfWeights());
         } else {
            if (m_debugging) Warning("normalizeHistogram","Skipping Normalizing histogram %s to ~zero: %f",hist->GetName(),hist->GetSumOfWeights());
         }
      }
   } else {
      Error("normalizeHistogram","Non existent histogram for normalizing");throw std::runtime_error("Throwing 56");
   }
}


//This allows PROOF to merge the generated histograms before the WriteToFile or WriteCustomToFile calls
Int_t Root::TPileupReweighting::Merge(TCollection *coll) {
   if(!coll) return 0;
   if(coll->IsEmpty()) return 0;

   // Iterate over the elements of the collection:
   TIter next( coll );
   TObject* obj = 0;
   while( ( obj = next() ) ) {

      // Check that the element is of the right type:
      Root::TPileupReweighting* vobj = dynamic_cast< Root::TPileupReweighting* >( obj );
      if( ! vobj ) {
         Error( "Merge", "Unknown object type encountered: %s",obj->ClassName() );
         return 0;
      }

      //merge the inputHistograms
      for(std::map<TString, std::map<Int_t,std::map<Int_t, TH1*> > >::iterator weights=vobj->GetInputHistograms().begin();weights!=vobj->GetInputHistograms().end();++weights) {
         TString weightName = weights->first;
         for(std::map<Int_t,std::map<Int_t, TH1*> >::iterator channels=weights->second.begin();channels!=weights->second.end();++channels) {
            Int_t channel = channels->first;
            for(std::map<Int_t, TH1*>::iterator runs=channels->second.begin();runs!=channels->second.end();++runs) {
               Int_t run=runs->first;
               if(runs->second) {
                  if(m_inputHistograms[weightName][channel][run]) {
                     m_inputHistograms[weightName][channel][run]->Add(runs->second);
                  } else {
                     m_inputHistograms[weightName][channel][run] = dynamic_cast<TH1*>(runs->second->Clone(runs->second->GetName()));
                     m_inputHistograms[weightName][channel][run]->SetDirectory(0);
                  }
               }
            }
         }
      }
   }

   return 1;
}

}
}

