// Dear emacs, this is -*-c++-*-

#ifndef __MYTPILEUPREWEIGHTING__
#define __MYTPILEUPREWEIGHTING__

/**
   @class TPileupReweighting
   @brief Tool to get the calculated MC pileup weight. Also does custom weights and other useful stuff.

   @author  Will Buttinger
            Based on the original tool by Karsten Koeneke
            Thanks to Nick Edwards for testing and useful feedback

*/

/* Developer notes:
* Could have incorporated the pileup histos in to the generic custom maps
* but decided to keep them seperate for probable slight performance improvement
* Allow the user to change the histo binning through AddBinning with name="pileup"
* */

#include "TNamed.h"
#include <TFile.h>
#include <TString.h>
#include "TVectorD.h"
#include <map>
#include <vector>
#include <TRandom3.h>

// STL includes
#include <iostream>
#include <stdexcept>

class TH1;
class TH1D;
class TH2D;
class TTree;
class TFile;

namespace AnalysisFramework
{
namespace External
{

namespace Root {
  class TPileupReweighting : public TNamed {

  public: 
      /** Standard constructor */
      TPileupReweighting(const char* name="TPileupReweighting");

      /** Standard destructor */
      ~TPileupReweighting();

  public:
      /** Initialize this class once before the event loop starts 
          If distribution information is provided, it is assumed to be 
          for the standard pileup reweighting */
      Int_t initialize(const TString dataRootFileName="",
                    const TString dataRootHistName="",
                    const TString mcRootFileName="",
                    const TString mcRootHistName="", int runNumber=0, int channel=0);

      Int_t Initialize();

      //-------------------------------------------------------------
      //Bonus method for retrieving the int lumi in each data period
      //the MuonEfficiencyCorrections package wants this
      //-------------------------------------------------------------
      std::vector<Double_t> getIntegratedLumiVector();

      /** total luminosity loaded and accepted by the tool */
      inline Double_t GetIntegratedLumi() { 
         if(!m_isInitialized) {
            Error("GetIntegratedLumi", "Please initialize the tool before retrieving the total lumi");
            throw std::runtime_error("Throwing 1");
         }
         if(!m_lumiVectorIsLoaded) {
            Error("GetIntegratedLumi","No Lumicalc file loaded, so cannot get integrated lumi, returning 0");
            return 0;
         }
         return globalTotals["pileup"][-1]/1E6; 
      }
      /** totalMC maps should also hold the total number of entries for each channel */
      inline Double_t GetNumberOfEvents(Int_t channel) {
         std::map<Int_t, Double_t>& tMap = globalNumberOfEntries["pileup"];
         if(tMap.find(channel)==tMap.end()) {
            Error("GetNumberOfEvents", "Unknown channel: %d",channel);
            return 0;
         }
         return tMap[channel];
      }
      inline Double_t GetSumOfEventWeights(Int_t channel) {
         std::map<Int_t, Double_t>& tMap = globalTotals["pileup"];
         if(tMap.find(channel)==tMap.end()) {
            Error("GetSumOfEventWeights", "Unknown channel: %d",channel);
            return 0;
         }
         return tMap[channel];
      }
      /** return fraction of lumi assigned to periodNumber (or runNumber) that is between start and end data run numbers*/
      Double_t GetIntegratedLumiFraction(Int_t periodNumber, UInt_t start, UInt_t end);

      /** get the amount of integrated lumi between the two given runs */
      Double_t GetIntegratedLumi(UInt_t start, UInt_t end);
      /** similar to above, but for only the given mcRunNumber/periodNumber */
      Double_t GetIntegratedLumi(Int_t periodNumber, UInt_t start, UInt_t end);

      //-----------------------------------------------------
      //General Tool Configuring methods
      //-----------------------------------------------------
      /** Indicate if warnings should be suppressed */
      inline void DisableWarnings(Bool_t in) { m_SetWarnings = !in;}
      /** Indicate if additional debugging information should be output */
      inline void EnableDebugging(Bool_t in) { m_debugging = in;}
      /** Set which channel should be used as a default when specific mc channel distributions cannot be found */
      inline void SetDefaultChannel(Int_t channel) { m_defaultChannel=channel;}
      /** Set how to handle configurations where some of your data has no corresponding mc to describe it 
          0=Default (Throw exception), 1=Subtract lumi from normalizations (i.e. discard data), 2=keep lumi and continue*/
      inline void SetUnrepresentedDataAction(Int_t action) {
         if(action<0 || action>2) {
            Fatal("SetUnrepresentedDataAction","Set action=%d - I'm afraid I can't let you do that, Dave",action);
            throw std::runtime_error("Throwing 2001");
         } 
         m_unrepresentedDataAction=action;
      }
      /** Should the tool ignore period assignments encoded in config file */
      inline void IgnoreConfigFilePeriods(Bool_t in) { m_ignoreFilePeriods=in; }

      /** Assign an mc RunNumber to a data period */
      Int_t AddPeriod(Int_t runNumber, UInt_t start, UInt_t end);
      /** Get the PeriodNumber of the given runNumber,start and end. If not found, returns -1 */
      Int_t GetPeriodNumber(Int_t runNumber, UInt_t start, UInt_t end);
      /** use a hardcoded period configuration */
      Int_t UsePeriodConfig(const TString configName);
      /** Add a histogram binning config. To modify the pileup histo binning, use "pileup" as name */
      Int_t AddBinning(const TString weightName,Int_t nbinsx, Double_t* xbins);
      Int_t AddBinning(const TString weightName,Int_t nbinsx, Double_t xlow, Double_t xup);
      Int_t AddBinning(const TString weightName,Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);
      Int_t AddBinning(const TString weightName,TH1* hist);

      /** Combine two mc RunNumbers. Histos are merged and the second number will be redirected to the first */
      void MergeMCRunNumbers(Int_t runNumber1, Int_t runNumber2, Int_t newRunNumber=0);
      /** Scale the LumiMetaData mu values by this factor */
      void SetDataScaleFactors(Float_t x,Float_t y=1., Float_t z=1.) { m_dataScaleFactorX=x;m_dataScaleFactorY=y;m_dataScaleFactorZ=z; }
      void SetMCScaleFactors(Float_t x,Float_t y=1.,Float_t z=1.) { m_mcScaleFactorX=x;m_mcScaleFactorY=y;m_mcScaleFactorZ=z;}

      //-----------------------------------------------------
      //Methods to load config files or histograms
      //-----------------------------------------------------
      Int_t AddConfigFile(const TString fileName);
      Int_t AddLumiCalcFile(const TString fileName);
      Int_t AddMetaDataFile(const TString fileName,const TString channelBranchName="mc_channel_number");

      Int_t AddDistribution(TH1* hist, const TString weightName, Int_t runNumber, Int_t channelNumber);
      Int_t AddDistribution(TH1* hist, Int_t runNumber, Int_t channelNumber);
      Int_t AddDistributionFromFile(const TString fileName, const TString histName, const TString weightName, Int_t runNumber=0, Int_t channelNumber=0);
      Int_t AddDistributionFromFile(const TString fileName, const TString histName, Int_t runNumber=0, Int_t channelNumber=0);

      //-----------------------------------------------------
      //Methods to get the various weights
      //-----------------------------------------------------
      Float_t GetCombinedWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber,Float_t x,Float_t y=0., Float_t z=0.);
      Float_t GetPeriodWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber);
      Float_t GetPrimaryWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber,Float_t x);
      Float_t GetSecondaryWeight(const TString weightName, Int_t periodNumber, Int_t channelNumber,Float_t x, Float_t y);


      //inline methods for above methods with the default "pileup" weight
      inline Float_t GetCombinedWeight(Int_t periodNumber, Int_t channelNumber,Float_t x, Float_t y=0., Float_t z=0.) {
         return GetCombinedWeight("pileup",periodNumber,channelNumber,x,y,z);
      }
      inline Float_t GetPeriodWeight(Int_t periodNumber, Int_t channelNumber) {
         return GetPeriodWeight("pileup",periodNumber,channelNumber);
      }
      inline Float_t GetPrimaryWeight(Int_t periodNumber, Int_t channelNumber,Float_t x) {
         return GetPrimaryWeight("pileup",periodNumber,channelNumber,x);
      }
      inline Float_t GetSecondaryWeight(Int_t periodNumber, Int_t channelNumber,Float_t x, Float_t y) {
         return GetSecondaryWeight("pileup",periodNumber,channelNumber,x,y);
      }

      //-----------------------------------------------------
      //Methods to work with the metadata
      //-----------------------------------------------------
      Double_t GetMetaData(const TString metadataName,Int_t channelNumber) {
         if(m_metadata.find(metadataName)==m_metadata.end()) {
            Error("GetMetaData","Metadata %s not known",metadataName.Data());
            return 0;
         }
         if(m_metadata[metadataName].find(channelNumber)==m_metadata[metadataName].end()) {
            Error("GetMetaData","Metadata %s not known in channel %d",metadataName.Data(), channelNumber);
            return 0;
         }
         return m_metadata[metadataName][channelNumber];
      }
      /** combines loaded metadata with channel sumsofweights and entry counts */
      TTree* GetMetaDataTree(); 
      Int_t GenerateMetaDataFile(const TString fileName,const TString channelBranchName="mc_channel_number");


      //-----------------------------------------------------
      //Methods to generate config files
      //-----------------------------------------------------
      Int_t Fill(const TString weightName,Int_t runNumber,Int_t channelNumber,Float_t w,Float_t x, Float_t y=0., Float_t z=0.);
      Int_t Fill(Int_t runNumber,Int_t channelNumber,Float_t w,Float_t x, Float_t y=0., Float_t z=0.);
      Int_t WriteToFile(const TString filename=""); //if no name given, will use tool name


      //-----------------------------------------------------
      //Methods to veto data in Action=1 mode
      //-----------------------------------------------------
      Bool_t IsUnrepresentedData(const TString weightName, Int_t runNumber, Float_t x, Float_t y=0., Float_t z=0.);
      Bool_t IsUnrepresentedData(Int_t runNumber, Float_t x, Float_t y=0., Float_t z=0.);


      //-----------------------------------------------------
      //RandomDataPeriod functionality stuff
      //numbers generated seperately for each mc period
      //-----------------------------------------------------
      void SetRandomSeed(int seed) {m_random3->SetSeed(seed);}
      int GetRandomSeed() {return m_random3->GetSeed();}
      /** Gets a random data run number according to the integrated lumi distribution associated to this mcRunNumber */
      UInt_t GetRandomRunNumber(Int_t mcRunNumber);
      /** Get random period number from the sub-periods assigned to this run number */
      Int_t GetRandomPeriodNumber(Int_t mcRunNumber);


      //-----------------------------------------------------
      //Methods for PROOF cluster merging of generated histos
      //-----------------------------------------------------
      Int_t Merge(TCollection *coll);
      std::map<TString, std::map<Int_t,std::map<Int_t, TH1*> > >& GetInputHistograms() { return m_inputHistograms;}


      //-----------------------------------------------------
      //Methods to inspect the input and weighting histograms
      //-----------------------------------------------------
      TH1* GetInputHistogram(const TString weightName, Int_t channelNumber, Int_t runNumber) {
         if(m_inputHistograms.find(weightName)==m_inputHistograms.end()||
            m_inputHistograms[weightName].find(channelNumber)==m_inputHistograms[weightName].end()||
            m_inputHistograms[weightName][channelNumber].find(runNumber)==m_inputHistograms[weightName][channelNumber].end()) {
            return 0;
         }
         return m_inputHistograms[weightName][channelNumber][runNumber];
      }

      TH1D* GetPrimaryDistribution(const TString weightName, Int_t channelNumber, Int_t periodNumber) {
         if(primaryDistributions.find(weightName)==primaryDistributions.end()||
            primaryDistributions[weightName].find(channelNumber)==primaryDistributions[weightName].end()||
            primaryDistributions[weightName][channelNumber].find(periodNumber)==primaryDistributions[weightName][channelNumber].end()) {
            return 0;
         }
         return primaryDistributions[weightName][channelNumber][periodNumber];
      }

      TH2D* GetSecondaryDistribution(const TString weightName, Int_t channelNumber, Int_t periodNumber) {
         if(secondaryDistributions.find(weightName)==secondaryDistributions.end()||
            secondaryDistributions[weightName].find(channelNumber)==secondaryDistributions[weightName].end()||
            secondaryDistributions[weightName][channelNumber].find(periodNumber)==secondaryDistributions[weightName][channelNumber].end()) {
            return 0;
         }
         return secondaryDistributions[weightName][channelNumber][periodNumber];
      }


      //-----------------------------------------------------
      //Depricated methods
      //-----------------------------------------------------
      /** Add a data distribution histogram or ttree (do this after mc) */
      void AddDataDistribution(const TString& dataRootFileName,const TString& dataRootHistName, int runNumber=0);
      /** Add a mc distribution histogram or ttree */
      void AddMCDistribution(const TString& mcRootFileName,const TString& mcRootHistName, int runNumber=0, int channel=0);
      /** Add a custom data distribution histogram or ttree (do this after mc) */
      void AddCustomDataDistribution(const TString& customRootFileName,const TString& customRootHistName, TString customName="", int runNumber=0);
      /** Add a custom mc distribution histogram or ttree */
      void AddCustomMCDistribution(const TString& customRootFileName,const TString& customRootHistName, TString customName="", int runNumber=0, int channel=0);
      Float_t getPileupWeight( float mu, int runNumber=0, int channel=0 );

  private:

      TH1* CloneEmptyHistogram(const TString weightName, Int_t runNumber, Int_t channelNumber);
      /** Normalize histograms */
      void normalizeHistogram(TH1* histo);
      void AddDistributionTree(TTree *tree, TFile *file);
      Int_t FactorizeDistribution(TH1* hist, const TString weightName, Int_t channelNumber, Int_t periodNumber,bool includeInMCRun,bool includeInGlobal);


      //********** Private members*************************
      Bool_t m_SetWarnings;
      Bool_t m_debugging;
      Bool_t m_countingMode;
      Int_t m_defaultChannel;
      Int_t m_unrepresentedDataAction;
      Bool_t m_isInitialized;
      Bool_t m_lumiVectorIsLoaded;
      Float_t m_dataScaleFactorX;Float_t m_dataScaleFactorY;Float_t m_dataScaleFactorZ;
      Float_t m_mcScaleFactorX;Float_t m_mcScaleFactorY;Float_t m_mcScaleFactorZ;
      Int_t m_nextPeriodNumber;
      Bool_t m_ignoreFilePeriods;
      TTree *m_metadatatree;

      //-----------------------------------------------------
      //Bonus method members
      //-----------------------------------------------------
      std::vector<Double_t> m_integratedLumiVector;

      //-----------------------------------------------------
      //Shared private data members
      //-----------------------------------------------------
      /** Map from mc RunNumber to (periodNumber,[start,end]) period definitions. One mc runNumber can have many disconnected periods */
      std::map<Int_t,std::map<Int_t, std::pair<UInt_t,UInt_t> > > m_periods;
      /** Quick map back from periodNumber to mcRunNumber */
      std::map<Int_t,Int_t> m_periodToMCRun;
      /** Map of empty histograms to use for custom weights. Also holds the standard pileup histogram in "pileup" */
      std::map<TString, TH1*> m_emptyHistograms;
      /** Remappings of mc run numbers. Used through the MergeMCRunNumbers method */
      std::map<Int_t,Int_t> m_mcRemappings;
      std::map<Int_t,std::map<Int_t,Int_t> > m_mcMerges;

      //-----------------------------------------------------
      //The standard pileup private members
      //-----------------------------------------------------
      /** [weightName,channelNumber] -> total (N,L,D) */
      std::map<TString, std::map<Int_t, Double_t> > globalTotals;
      /** [weightName,channelNumber,periodNumber] -> total (N_A,L_A,D_A) */
      std::map<TString, std::map<Int_t, std::map<Int_t, Double_t> > > periodTotals;
      /** [weightName,channelNumber,periodNumber] -> 1D distribution */
      std::map<TString, std::map<Int_t, std::map<Int_t, TH1D*> > > primaryDistributions;
      /** [weightName,channelNumber,periodNumber] -> 2D distribution */
      std::map<TString, std::map<Int_t, std::map<Int_t, TH2D*> > > secondaryDistributions;

      /** [weightName,channelNumber] -> NumberOfEntries from histograms */
      std::map<TString, std::map<Int_t, Double_t> > globalNumberOfEntries;

      //data is stored in channel=-1
      std::map<TString, std::map<Int_t,std::map<Int_t, TH1*> > > m_inputHistograms;

      /** channel metadata map */
      std::map<TString, std::map<Int_t, Double_t> > m_metadata; 

      //-----------------------------------------------------
      //RandomDataPeriod functionality stuff
      //numbers generated seperately for each mc period
      //-----------------------------------------------------
      TRandom3 *m_random3;
      /** Integrated lumi in each run, split in to the different mc periods */
      std::map<TString, std::map<Int_t, std::map<UInt_t, Double_t> > > dataPeriodRunTotals;
      /** [weightName,datarunnum,binnum] -> badbin flag */
      std::map<TString, std::map<Int_t, std::map<Int_t, Bool_t> > > m_badbins;

      ClassDef(TPileupReweighting,1)


  }; // End: class definition

} // End: namespace Root

}
}

#endif
