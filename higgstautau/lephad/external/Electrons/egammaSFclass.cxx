// Efficiency scale factors macro
// Date: 11/01/2010
// Author: Olivier Arnaez <olivier.arnaez@cern.ch>
//         Jan Kretzschmar <jan.kretzschmar@cern.ch>
//
// Usage: 
// std::pair<float,float> sf_error = objsf->scaleFactor(eta(cluster), ET(MeV), set, range, rel, etcorrection)
//
// Please note: alternative accessors scaleFactorLoose, scaleFactorMedium, ... 
//              have been disabled, as the amount of sets is expanding and this was not maintainable
//              Read the "set=X" documentation below and use the correct arguments for scaleFactor() function
//
// The first number (sf_error.first) returns the efficiency scale factor,
// the second number is its uncertainty (sf_errror.second)
// 
// The combined (W/Z/jpsi) scale factor and uncertainty vs eta and ET (in MeV) are given
// 
// different "sets" of numbers available (not all with all sets):
//    * Loose SF (set=0)
//    * Medium SF (set=1)
//    * Tight SF (set=2)
//    * e20_medium trigger SF (set=3) (use set >=8 for release 17 2011 data/MC11a)
//    * reco+trkqual SF (set=4)
//    * Loose++ SF (set=5)
//    * Medium++ SF (set=6)
//    * Tight++ SF (set=7)
//    * e20_medium trigger SF w.r.t Medium++ offline (set=8)
//    * e20_medium MC efficiency w.r.t Medium++ offline (set=9)
//    * e20_medium trigger SF w.r.t Tight++ offline (set=10)
//    * e20_medium MC efficiency w.r.t Tight++ offline (set=11)
//    * e22_medium trigger SF w.r.t Medium++ offline (set=12)
//    * e22_medium MC efficiency w.r.t Medium++ offline (set=13)
//    * e22_medium trigger SF w.r.t Tight++ offline (set=14)
//    * e22_medium MC efficiency w.r.t Tight++ offline (set=15)
//    * e22vh_medium1 trigger SF (using e22_medium1 on MC11a) w.r.t Medium++ offline (set=16)
//    * e22_medium1 MC efficiency w.r.t Medium++ offline (set=17)
//    * e22vh_medium1 trigger SF (using e22_medium1 on MC11a) w.r.t Tight++ offline (set=18)
//    * e22_medium1 MC efficiency w.r.t Tight++ offline (set=19)
//    * e20_medium MC efficiency w.r.t Loose++ offline (set=20)
//    * e22_medium1 MC efficiency w.r.t Loose++ offline (set=21)
//    * e22vh_medium1 MC efficiency w.r.t Loose++ offline (set=22)
// data and MC release selection:
//    * release 15 2010 data/MC09 (rel=0)
//    * release 16 2010 data/MC10 (rel=1)
//    * release 16.6 2011 data/MC10ab (estimated from 2010 data) (rel=2)
//    * release 16.6 estimated from 2011 data "EPS recommendations" and MC10b (rel=3)
//    * release 16.6 estimated from 2011 data "EPS recommendations" including Jpsi measurements (rel=4)
//    * release 17 estimated from 2011 data/MC11a "CERN council recommendations" (rel=5)
//    * release 17 estimated from 2011 data/MC11a/b/c "Moriond recommendations" G4 FullSim MC (rel=6)
//    * release 17 estimated from 2011 data/MC11a/b/c "Moriond recommendations" AFII MC (rel=7)
// measured with probes in the ranges:
//    * 20-50 GeV range (range=0)
//    * 30-50 GeV (range=1)
// and correcting (etcorrection=1) or not (etcorrection=0) for the ET-dependence
//
// Eta binning is changing from release to release
//
// Note that for rel>=4 range should be left at 0 and etcorrection=1 always
// 
// For now separete function for Forward Electrons (|eta|>=2.5)
// and ONLY for release 16.6, 2011 data in analogy to function for central electrons
// std::pair<float,float> sf_error = objsf->scaleFactorForward(eta, set)
// where eta = electron eta = cluster eta (as no track)
// cuts are ForwardLoose (set=0) or ForwardTight (set=2)
// 

#include "egammaSFclass.h"
#include <cmath>

namespace AnalysisFramework
{
namespace External
{
egammaSFclass::egammaSFclass()
{
  //Definition of the eta binning
  m_Etabins.push_back(-2.47);
  m_Etabins.push_back(-2.01); 
  m_Etabins.push_back(-1.52); 
  m_Etabins.push_back(-1.37); 
  m_Etabins.push_back(-0.8); 
  m_Etabins.push_back(0); 
  m_Etabins.push_back(0.8); 
  m_Etabins.push_back(1.37); 
  m_Etabins.push_back(1.52); 
  m_Etabins.push_back(2.01); 
  m_Etabins.push_back(2.47);
  //Definition of the fine eta binning
  m_FineEtabins.push_back(-2.47);
  m_FineEtabins.push_back(-2.37);
  m_FineEtabins.push_back(-2.01);
  m_FineEtabins.push_back(-1.81);
  m_FineEtabins.push_back(-1.52);
  m_FineEtabins.push_back(-1.37);
  m_FineEtabins.push_back(-1.15);
  m_FineEtabins.push_back(-0.8 );
  m_FineEtabins.push_back(-0.6 );
  m_FineEtabins.push_back(-0.1 );
  m_FineEtabins.push_back( 0.  );
  m_FineEtabins.push_back( 0.1 );
  m_FineEtabins.push_back( 0.6 );
  m_FineEtabins.push_back( 0.8 );
  m_FineEtabins.push_back( 1.15);
  m_FineEtabins.push_back( 1.37);
  m_FineEtabins.push_back( 1.52);
  m_FineEtabins.push_back( 1.81);
  m_FineEtabins.push_back( 2.01);
  m_FineEtabins.push_back( 2.37);
  m_FineEtabins.push_back( 2.47);
  //Definition of the eta binning with 11 bins
  m_11Etabins.push_back(-2.47);
  m_11Etabins.push_back(-2.01); 
  m_11Etabins.push_back(-1.52); 
  m_11Etabins.push_back(-1.37); 
  m_11Etabins.push_back(-0.8); 
  m_11Etabins.push_back(-0.1); 
  m_11Etabins.push_back(0.1); 
  m_11Etabins.push_back(0.8); 
  m_11Etabins.push_back(1.37); 
  m_11Etabins.push_back(1.52); 
  m_11Etabins.push_back(2.01); 
  m_11Etabins.push_back(2.47);
  //Definition of the ET binning
  m_ETbins.push_back(0.);
  m_ETbins.push_back(20000.); 
  m_ETbins.push_back(25000.); 
  m_ETbins.push_back(30000.); 
  m_ETbins.push_back(35000.); 
  m_ETbins.push_back(40000.); 
  m_ETbins.push_back(45000.); 
  m_ETbins.push_back(500000000.); 
  //Definition of the ET binning on the full range
  m_ETbinsFullRange.push_back(    0.);
  m_ETbinsFullRange.push_back( 7000.);
  m_ETbinsFullRange.push_back(10000.);
  m_ETbinsFullRange.push_back(15000.);
  m_ETbinsFullRange.push_back(20000.); 
  m_ETbinsFullRange.push_back(25000.); 
  m_ETbinsFullRange.push_back(30000.); 
  m_ETbinsFullRange.push_back(35000.); 
  m_ETbinsFullRange.push_back(40000.); 
  m_ETbinsFullRange.push_back(45000.); 
  m_ETbinsFullRange.push_back(500000000.); 
  //Definition of the ET binning for trigger
  m_ETbinsTrigger.push_back(21000.);
  m_ETbinsTrigger.push_back(23000.);
  m_ETbinsTrigger.push_back(25000.); 
  m_ETbinsTrigger.push_back(30000.); 
  m_ETbinsTrigger.push_back(35000.); 
  m_ETbinsTrigger.push_back(40000.); 
  m_ETbinsTrigger.push_back(500000000.); 


  //For the scale factors of the standard egamma cuts 

  //Release 15
  //Probes between 30 and 50 GeV (plateau region)
  //Loose
  efficienciesRel15Loose3050.push_back(98.1); 
  efficienciesRel15Loose3050.push_back(99.0); 
  efficienciesRel15Loose3050.push_back(0.); 
  efficienciesRel15Loose3050.push_back(98.6); 
  efficienciesRel15Loose3050.push_back(99.5); 
  efficienciesRel15Loose3050.push_back(99.1); 
  efficienciesRel15Loose3050.push_back(98.8); 
  efficienciesRel15Loose3050.push_back(0.); 
  efficienciesRel15Loose3050.push_back(99.9); 
  efficienciesRel15Loose3050.push_back(98.2);
  uncertaintiesRel15Loose3050.push_back(1.6); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back(0.); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back(0.); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.6);
  //Medium
  efficienciesRel15Medium3050.push_back(95.4); 
  efficienciesRel15Medium3050.push_back(98.7);
  efficienciesRel15Medium3050.push_back(0.); 
  efficienciesRel15Medium3050.push_back(97.9);
  efficienciesRel15Medium3050.push_back(98.1);
  efficienciesRel15Medium3050.push_back(97.7); 
  efficienciesRel15Medium3050.push_back(97.9); 
  efficienciesRel15Medium3050.push_back(0.); 
  efficienciesRel15Medium3050.push_back(99.9); 
  efficienciesRel15Medium3050.push_back(97.4);
  uncertaintiesRel15Medium3050.push_back(1.7);
  uncertaintiesRel15Medium3050.push_back( 1.6);
  uncertaintiesRel15Medium3050.push_back(0.); 
  uncertaintiesRel15Medium3050.push_back( 1.6);
  uncertaintiesRel15Medium3050.push_back( 1.5); 
  uncertaintiesRel15Medium3050.push_back( 1.5); 
  uncertaintiesRel15Medium3050.push_back( 1.5); 
  uncertaintiesRel15Medium3050.push_back(0.); 
  uncertaintiesRel15Medium3050.push_back( 1.6);
  uncertaintiesRel15Medium3050.push_back( 1.7);
  //Tight
  efficienciesRel15Tight3050.push_back(92.3); 
  efficienciesRel15Tight3050.push_back(99.2); 
  efficienciesRel15Tight3050.push_back(0.);
  efficienciesRel15Tight3050.push_back(101.5); 
  efficienciesRel15Tight3050.push_back(98.9); 
  efficienciesRel15Tight3050.push_back(99.9);
  efficienciesRel15Tight3050.push_back(104.2); 
  efficienciesRel15Tight3050.push_back(0.);
  efficienciesRel15Tight3050.push_back(102.6); 
  efficienciesRel15Tight3050.push_back(95.5);
  uncertaintiesRel15Tight3050.push_back(3.3);
  uncertaintiesRel15Tight3050.push_back( 2.3); 
  uncertaintiesRel15Tight3050.push_back(0.);
  uncertaintiesRel15Tight3050.push_back( 2.0); 
  uncertaintiesRel15Tight3050.push_back( 1.8); 
  uncertaintiesRel15Tight3050.push_back( 1.8);
  uncertaintiesRel15Tight3050.push_back( 2.5); 
  uncertaintiesRel15Tight3050.push_back(0.); 
  uncertaintiesRel15Tight3050.push_back( 5.0); 
  uncertaintiesRel15Tight3050.push_back( 3.2);

  //Probes between 20 and 50 GeV
  //Loose
  efficienciesRel15Loose2050.push_back(97.6); 
  efficienciesRel15Loose2050.push_back(99.0); 
  efficienciesRel15Loose2050.push_back(0.); 
  efficienciesRel15Loose2050.push_back(98.2); 
  efficienciesRel15Loose2050.push_back(99.1); 
  efficienciesRel15Loose2050.push_back(98.8); 
  efficienciesRel15Loose2050.push_back(98.2); 
  efficienciesRel15Loose2050.push_back(0.); 
  efficienciesRel15Loose2050.push_back(99.6); 
  efficienciesRel15Loose2050.push_back(97.4);
  uncertaintiesRel15Loose2050.push_back(1.6); 
  uncertaintiesRel15Loose2050.push_back(1.5); 
  uncertaintiesRel15Loose2050.push_back(0.); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back(0.); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.6);
  //Medium
  efficienciesRel15Medium2050.push_back(94.5); 
  efficienciesRel15Medium2050.push_back(98.8);
  efficienciesRel15Medium2050.push_back(0.); 
  efficienciesRel15Medium2050.push_back(97.2);
  efficienciesRel15Medium2050.push_back(97.4);
  efficienciesRel15Medium2050.push_back(97.2); 
  efficienciesRel15Medium2050.push_back(96.7); 
  efficienciesRel15Medium2050.push_back(0.); 
  efficienciesRel15Medium2050.push_back(99.5); 
  efficienciesRel15Medium2050.push_back(96.1);
  uncertaintiesRel15Medium2050.push_back(1.7);
  uncertaintiesRel15Medium2050.push_back( 1.6);
  uncertaintiesRel15Medium2050.push_back(0.); 
  uncertaintiesRel15Medium2050.push_back( 1.6);
  uncertaintiesRel15Medium2050.push_back( 1.5); 
  uncertaintiesRel15Medium2050.push_back( 1.5); 
  uncertaintiesRel15Medium2050.push_back( 1.5); 
  uncertaintiesRel15Medium2050.push_back(0.); 
  uncertaintiesRel15Medium2050.push_back( 2.9);
  uncertaintiesRel15Medium2050.push_back( 1.7);
  //Tight
  efficienciesRel15Tight2050.push_back(92.5); 
  efficienciesRel15Tight2050.push_back(99.5); 
  efficienciesRel15Tight2050.push_back(0.);
  efficienciesRel15Tight2050.push_back(100.6); 
  efficienciesRel15Tight2050.push_back(98.2); 
  efficienciesRel15Tight2050.push_back(98.7);
  efficienciesRel15Tight2050.push_back(103.3); 
  efficienciesRel15Tight2050.push_back(0.);
  efficienciesRel15Tight2050.push_back(102.8); 
  efficienciesRel15Tight2050.push_back(93.6);
  uncertaintiesRel15Tight2050.push_back(3.4);
  uncertaintiesRel15Tight2050.push_back( 2.4); 
  uncertaintiesRel15Tight2050.push_back(0.);
  uncertaintiesRel15Tight2050.push_back( 2.1); 
  uncertaintiesRel15Tight2050.push_back( 1.8); 
  uncertaintiesRel15Tight2050.push_back( 1.8);
  uncertaintiesRel15Tight2050.push_back( 2.5); 
  uncertaintiesRel15Tight2050.push_back(0.); 
  uncertaintiesRel15Tight2050.push_back( 4.5); 
  uncertaintiesRel15Tight2050.push_back( 3.4);


  //Release 16
  //Probes between 30 and 50 GeV (plateau region)
  //Medium
  efficienciesRel16Medium3050.push_back(98.8); 
  efficienciesRel16Medium3050.push_back(98.0);
  efficienciesRel16Medium3050.push_back(96.9); 
  efficienciesRel16Medium3050.push_back(98.0);
  efficienciesRel16Medium3050.push_back(97.4);
  efficienciesRel16Medium3050.push_back(98.1); 
  efficienciesRel16Medium3050.push_back(98.1); 
  efficienciesRel16Medium3050.push_back(98.3); 
  efficienciesRel16Medium3050.push_back(98.6); 
  efficienciesRel16Medium3050.push_back(97.5);
  uncertaintiesRel16Medium3050.push_back(0.8);
  uncertaintiesRel16Medium3050.push_back(0.9);
  uncertaintiesRel16Medium3050.push_back(2.5); 
  uncertaintiesRel16Medium3050.push_back(0.8);
  uncertaintiesRel16Medium3050.push_back(0.7); 
  uncertaintiesRel16Medium3050.push_back(0.7); 
  uncertaintiesRel16Medium3050.push_back(0.8); 
  uncertaintiesRel16Medium3050.push_back(2.6); 
  uncertaintiesRel16Medium3050.push_back(0.8);
  uncertaintiesRel16Medium3050.push_back(0.8);
  //Tight
  efficienciesRel16Tight3050.push_back(102.0); 
  efficienciesRel16Tight3050.push_back(102.7); 
  efficienciesRel16Tight3050.push_back(114.4);
  efficienciesRel16Tight3050.push_back(106.7); 
  efficienciesRel16Tight3050.push_back( 99.0); 
  efficienciesRel16Tight3050.push_back(100.1);
  efficienciesRel16Tight3050.push_back(105.7); 
  efficienciesRel16Tight3050.push_back(110.8);
  efficienciesRel16Tight3050.push_back(104.2); 
  efficienciesRel16Tight3050.push_back(102.7);
  uncertaintiesRel16Tight3050.push_back(3.0);
  uncertaintiesRel16Tight3050.push_back(1.1); 
  uncertaintiesRel16Tight3050.push_back(3.9);
  uncertaintiesRel16Tight3050.push_back(1.1); 
  uncertaintiesRel16Tight3050.push_back(0.8); 
  uncertaintiesRel16Tight3050.push_back(0.8);
  uncertaintiesRel16Tight3050.push_back(0.9); 
  uncertaintiesRel16Tight3050.push_back(4.6); 
  uncertaintiesRel16Tight3050.push_back(2.6); 
  uncertaintiesRel16Tight3050.push_back(1.2);

  //Probes between 20 and 50 GeV
  //Medium
  efficienciesRel16Medium2050.push_back(97.6); 
  efficienciesRel16Medium2050.push_back(96.8);
  efficienciesRel16Medium2050.push_back(97.7); 
  efficienciesRel16Medium2050.push_back(97.1);
  efficienciesRel16Medium2050.push_back(96.8);
  efficienciesRel16Medium2050.push_back(97.6); 
  efficienciesRel16Medium2050.push_back(97.2); 
  efficienciesRel16Medium2050.push_back(98.2); 
  efficienciesRel16Medium2050.push_back(97.9); 
  efficienciesRel16Medium2050.push_back(96.2);
  uncertaintiesRel16Medium2050.push_back(1.0);
  uncertaintiesRel16Medium2050.push_back(1.0);
  uncertaintiesRel16Medium2050.push_back(3.3); 
  uncertaintiesRel16Medium2050.push_back(1.1);
  uncertaintiesRel16Medium2050.push_back(0.8); 
  uncertaintiesRel16Medium2050.push_back(0.8); 
  uncertaintiesRel16Medium2050.push_back(0.9); 
  uncertaintiesRel16Medium2050.push_back(3.2); 
  uncertaintiesRel16Medium2050.push_back(1.0);
  uncertaintiesRel16Medium2050.push_back(2.8);
  //Tight
  efficienciesRel16Tight2050.push_back(100.2); 
  efficienciesRel16Tight2050.push_back(101.5); 
  efficienciesRel16Tight2050.push_back(117.9);
  efficienciesRel16Tight2050.push_back(105.7); 
  efficienciesRel16Tight2050.push_back( 98.1); 
  efficienciesRel16Tight2050.push_back( 99.1);
  efficienciesRel16Tight2050.push_back(105.2); 
  efficienciesRel16Tight2050.push_back(113.9);
  efficienciesRel16Tight2050.push_back(103.8); 
  efficienciesRel16Tight2050.push_back(101.2);
  uncertaintiesRel16Tight2050.push_back(1.1);
  uncertaintiesRel16Tight2050.push_back(1.2); 
  uncertaintiesRel16Tight2050.push_back(4.4);
  uncertaintiesRel16Tight2050.push_back(1.5); 
  uncertaintiesRel16Tight2050.push_back(0.9); 
  uncertaintiesRel16Tight2050.push_back(1.0);
  uncertaintiesRel16Tight2050.push_back(1.1); 
  uncertaintiesRel16Tight2050.push_back(5.2); 
  uncertaintiesRel16Tight2050.push_back(3.0); 
  uncertaintiesRel16Tight2050.push_back(1.3);


  //For the ET-corrections of the scale factors
  //Medium
  ETCorrectionsMediumRel16.push_back( 79.6);
  ETCorrectionsMediumRel16.push_back( 93.9);
  ETCorrectionsMediumRel16.push_back( 96.2);
  ETCorrectionsMediumRel16.push_back( 99.7);
  ETCorrectionsMediumRel16.push_back(100.6);
  ETCorrectionsMediumRel16.push_back(100.4);
  ETCorrectionsMediumRel16.push_back(101.00);
  uncertaintiesETCorrectionsMediumRel16.push_back( 9.4);
  uncertaintiesETCorrectionsMediumRel16.push_back( 3.6);
  uncertaintiesETCorrectionsMediumRel16.push_back( 1.4);
  uncertaintiesETCorrectionsMediumRel16.push_back( 0.7);
  uncertaintiesETCorrectionsMediumRel16.push_back( 0.5);
  uncertaintiesETCorrectionsMediumRel16.push_back( 0.7);
  uncertaintiesETCorrectionsMediumRel16.push_back( 1.7);
  //Medium
  ETCorrectionsTightRel16.push_back( 76.7);
  ETCorrectionsTightRel16.push_back( 93.6);
  ETCorrectionsTightRel16.push_back( 95.1);
  ETCorrectionsTightRel16.push_back( 99.9);
  ETCorrectionsTightRel16.push_back(100.4);
  ETCorrectionsTightRel16.push_back(100.0);
  ETCorrectionsTightRel16.push_back(100.7);
  uncertaintiesETCorrectionsTightRel16.push_back(10.0);
  uncertaintiesETCorrectionsTightRel16.push_back( 3.7);
  uncertaintiesETCorrectionsTightRel16.push_back( 1.6);
  uncertaintiesETCorrectionsTightRel16.push_back( 0.9);
  uncertaintiesETCorrectionsTightRel16.push_back( 0.7);
  uncertaintiesETCorrectionsTightRel16.push_back( 0.9);
  uncertaintiesETCorrectionsTightRel16.push_back( 1.8);



  //Release 16.6 Data 2010
  //Probes between 30 and 50 GeV (plateau region)
  //Medium
  efficienciesRel166Data2010Medium3050.push_back(98.44); 
  efficienciesRel166Data2010Medium3050.push_back(96.93);
  efficienciesRel166Data2010Medium3050.push_back(96.61); 
  efficienciesRel166Data2010Medium3050.push_back(96.87);
  efficienciesRel166Data2010Medium3050.push_back(97.06);
  efficienciesRel166Data2010Medium3050.push_back(97.49); 
  efficienciesRel166Data2010Medium3050.push_back(97.04); 
  efficienciesRel166Data2010Medium3050.push_back(97.17); 
  efficienciesRel166Data2010Medium3050.push_back(97.31); 
  efficienciesRel166Data2010Medium3050.push_back(97.51);
  uncertaintiesRel166Data2010Medium3050.push_back(2.14);
  uncertaintiesRel166Data2010Medium3050.push_back(2.20);
  uncertaintiesRel166Data2010Medium3050.push_back(2.84); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.13);
  uncertaintiesRel166Data2010Medium3050.push_back(2.18); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.10); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.13); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.89); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.13);
  uncertaintiesRel166Data2010Medium3050.push_back(2.21);
  //Tight
  efficienciesRel166Data2010Tight3050.push_back(101.47); 
  efficienciesRel166Data2010Tight3050.push_back(104.02); 
  efficienciesRel166Data2010Tight3050.push_back(112.70);
  efficienciesRel166Data2010Tight3050.push_back(106.82); 
  efficienciesRel166Data2010Tight3050.push_back( 99.35); 
  efficienciesRel166Data2010Tight3050.push_back(100.13);
  efficienciesRel166Data2010Tight3050.push_back(105.94); 
  efficienciesRel166Data2010Tight3050.push_back(113.57);
  efficienciesRel166Data2010Tight3050.push_back(105.48); 
  efficienciesRel166Data2010Tight3050.push_back(101.99);
  uncertaintiesRel166Data2010Tight3050.push_back(3.46);
  uncertaintiesRel166Data2010Tight3050.push_back(2.65); 
  uncertaintiesRel166Data2010Tight3050.push_back(3.65);
  uncertaintiesRel166Data2010Tight3050.push_back(2.49); 
  uncertaintiesRel166Data2010Tight3050.push_back(2.33); 
  uncertaintiesRel166Data2010Tight3050.push_back(2.28);
  uncertaintiesRel166Data2010Tight3050.push_back(2.45); 
  uncertaintiesRel166Data2010Tight3050.push_back(3.72); 
  uncertaintiesRel166Data2010Tight3050.push_back(3.38); 
  uncertaintiesRel166Data2010Tight3050.push_back(2.70);

  //Probes between 20 and 50 GeV
  //Medium
  efficienciesRel166Data2010Medium2050.push_back(97.35); 
  efficienciesRel166Data2010Medium2050.push_back(95.86);
  efficienciesRel166Data2010Medium2050.push_back(96.25); 
  efficienciesRel166Data2010Medium2050.push_back(95.80);
  efficienciesRel166Data2010Medium2050.push_back(96.01);
  efficienciesRel166Data2010Medium2050.push_back(96.84); 
  efficienciesRel166Data2010Medium2050.push_back(96.04); 
  efficienciesRel166Data2010Medium2050.push_back(96.54); 
  efficienciesRel166Data2010Medium2050.push_back(96.59); 
  efficienciesRel166Data2010Medium2050.push_back(96.33);
  uncertaintiesRel166Data2010Medium2050.push_back(2.21);
  uncertaintiesRel166Data2010Medium2050.push_back(2.25);
  uncertaintiesRel166Data2010Medium2050.push_back(3.22); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.27);
  uncertaintiesRel166Data2010Medium2050.push_back(2.23); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.13); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.17); 
  uncertaintiesRel166Data2010Medium2050.push_back(3.20); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.24);
  uncertaintiesRel166Data2010Medium2050.push_back(2.41);
  //Tight
  efficienciesRel166Data2010Tight2050.push_back(99.90); 
  efficienciesRel166Data2010Tight2050.push_back(103.11); 
  efficienciesRel166Data2010Tight2050.push_back(116.16);
  efficienciesRel166Data2010Tight2050.push_back(105.70); 
  efficienciesRel166Data2010Tight2050.push_back( 97.98); 
  efficienciesRel166Data2010Tight2050.push_back( 99.08);
  efficienciesRel166Data2010Tight2050.push_back(105.23); 
  efficienciesRel166Data2010Tight2050.push_back(115.12);
  efficienciesRel166Data2010Tight2050.push_back(104.91); 
  efficienciesRel166Data2010Tight2050.push_back(101.99);
  uncertaintiesRel166Data2010Tight2050.push_back(2.28);
  uncertaintiesRel166Data2010Tight2050.push_back(2.89); 
  uncertaintiesRel166Data2010Tight2050.push_back(4.35);
  uncertaintiesRel166Data2010Tight2050.push_back(2.72); 
  uncertaintiesRel166Data2010Tight2050.push_back(2.40); 
  uncertaintiesRel166Data2010Tight2050.push_back(2.24);
  uncertaintiesRel166Data2010Tight2050.push_back(2.48); 
  uncertaintiesRel166Data2010Tight2050.push_back(4.17); 
  uncertaintiesRel166Data2010Tight2050.push_back(2.45); 
  uncertaintiesRel166Data2010Tight2050.push_back(3.29);
  //For the ET-corrections of the scale factors
  //Medium
  ETCorrectionsMediumRel166Data2010.push_back(80.60);
  ETCorrectionsMediumRel166Data2010.push_back(92.07);
  ETCorrectionsMediumRel166Data2010.push_back(96.34);
  ETCorrectionsMediumRel166Data2010.push_back(100.19);
  ETCorrectionsMediumRel166Data2010.push_back(101.54);
  ETCorrectionsMediumRel166Data2010.push_back(101.25);
  ETCorrectionsMediumRel166Data2010.push_back(102.29);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 9.60);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 3.27);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 1.40);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 0.70);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 0.53);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 0.74);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 1.59);
  //Tight
  ETCorrectionsTightRel166Data2010.push_back(77.78);
  ETCorrectionsTightRel166Data2010.push_back(91.84);
  ETCorrectionsTightRel166Data2010.push_back(95.67);
  ETCorrectionsTightRel166Data2010.push_back(100.86);
  ETCorrectionsTightRel166Data2010.push_back(101.83);
  ETCorrectionsTightRel166Data2010.push_back(101.33);
  ETCorrectionsTightRel166Data2010.push_back(102.10);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back(10.29);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 3.47);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 1.52);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 1.04);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 0.66);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 0.92);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 1.90);


  //Release 16.6 Data 2011 EPS recommendations
  //Identification for probes between 20 and 50 GeV
  //Medium
  efficienciesRel166EPSMedium2050.push_back(95.7273);
  efficienciesRel166EPSMedium2050.push_back(95.5243);
  efficienciesRel166EPSMedium2050.push_back(96.403);
  efficienciesRel166EPSMedium2050.push_back(96.3494);
  efficienciesRel166EPSMedium2050.push_back(97.9518);
  efficienciesRel166EPSMedium2050.push_back(96.3292);
  efficienciesRel166EPSMedium2050.push_back(97.0952);
  efficienciesRel166EPSMedium2050.push_back(96.3317);
  efficienciesRel166EPSMedium2050.push_back(97.1977);
  efficienciesRel166EPSMedium2050.push_back(97.8678);
  efficienciesRel166EPSMedium2050.push_back(96.5697);
  efficienciesRel166EPSMedium2050.push_back(96.7783);
  efficienciesRel166EPSMedium2050.push_back(97.0532);
  efficienciesRel166EPSMedium2050.push_back(96.4621);
  efficienciesRel166EPSMedium2050.push_back(95.3501);
  efficienciesRel166EPSMedium2050.push_back(97.9656);
  efficienciesRel166EPSMedium2050.push_back(96.3031);
  efficienciesRel166EPSMedium2050.push_back(97.3978);
  efficienciesRel166EPSMedium2050.push_back(95.7546);
  efficienciesRel166EPSMedium2050.push_back(97.2443);
  uncertaintiesRel166EPSMedium2050.push_back(0.758538);
  uncertaintiesRel166EPSMedium2050.push_back(1.48083);
  uncertaintiesRel166EPSMedium2050.push_back(0.778086);
  uncertaintiesRel166EPSMedium2050.push_back(0.496963);
  uncertaintiesRel166EPSMedium2050.push_back(1.0011);
  uncertaintiesRel166EPSMedium2050.push_back(0.694056);
  uncertaintiesRel166EPSMedium2050.push_back(0.603261);
  uncertaintiesRel166EPSMedium2050.push_back(0.719089);
  uncertaintiesRel166EPSMedium2050.push_back(0.635625);
  uncertaintiesRel166EPSMedium2050.push_back(0.825545);
  uncertaintiesRel166EPSMedium2050.push_back(0.777055);
  uncertaintiesRel166EPSMedium2050.push_back(0.655198);
  uncertaintiesRel166EPSMedium2050.push_back(0.736623);
  uncertaintiesRel166EPSMedium2050.push_back(0.633197);
  uncertaintiesRel166EPSMedium2050.push_back(1.04172);
  uncertaintiesRel166EPSMedium2050.push_back(0.612204);
  uncertaintiesRel166EPSMedium2050.push_back(0.47725);
  uncertaintiesRel166EPSMedium2050.push_back(1.32532);
  uncertaintiesRel166EPSMedium2050.push_back(0.74313);
  uncertaintiesRel166EPSMedium2050.push_back(1.44683);
  //Tight
  efficienciesRel166EPSTight2050.push_back( 99.9569);
  efficienciesRel166EPSTight2050.push_back( 99.1664);
  efficienciesRel166EPSTight2050.push_back(103.421 );
  efficienciesRel166EPSTight2050.push_back(102.688 );
  efficienciesRel166EPSTight2050.push_back(113.028 );
  efficienciesRel166EPSTight2050.push_back(111.078 );
  efficienciesRel166EPSTight2050.push_back(103.481 );
  efficienciesRel166EPSTight2050.push_back( 99.5783);
  efficienciesRel166EPSTight2050.push_back( 98.4303);
  efficienciesRel166EPSTight2050.push_back(100.837 );
  efficienciesRel166EPSTight2050.push_back( 99.1868);
  efficienciesRel166EPSTight2050.push_back( 98.1188);
  efficienciesRel166EPSTight2050.push_back(100.492 );
  efficienciesRel166EPSTight2050.push_back(102.816 );
  efficienciesRel166EPSTight2050.push_back(109.09  );
  efficienciesRel166EPSTight2050.push_back(113.772 );
  efficienciesRel166EPSTight2050.push_back(103.355 );
  efficienciesRel166EPSTight2050.push_back(103.454 );
  efficienciesRel166EPSTight2050.push_back( 98.4376);
  efficienciesRel166EPSTight2050.push_back(102.174 );
  uncertaintiesRel166EPSTight2050.push_back(2.82899);
  uncertaintiesRel166EPSTight2050.push_back(1.47076);
  uncertaintiesRel166EPSTight2050.push_back(2.64305);
  uncertaintiesRel166EPSTight2050.push_back(0.692373);
  uncertaintiesRel166EPSTight2050.push_back(2.0146 );
  uncertaintiesRel166EPSTight2050.push_back(0.967662);
  uncertaintiesRel166EPSTight2050.push_back(0.714802);
  uncertaintiesRel166EPSTight2050.push_back(0.807023);
  uncertaintiesRel166EPSTight2050.push_back(0.686988);
  uncertaintiesRel166EPSTight2050.push_back(1.4562);
  uncertaintiesRel166EPSTight2050.push_back(0.984975);
  uncertaintiesRel166EPSTight2050.push_back(0.703155);
  uncertaintiesRel166EPSTight2050.push_back(0.80346);
  uncertaintiesRel166EPSTight2050.push_back(0.742777);
  uncertaintiesRel166EPSTight2050.push_back(1.78409);
  uncertaintiesRel166EPSTight2050.push_back(1.13598);
  uncertaintiesRel166EPSTight2050.push_back(0.716145);
  uncertaintiesRel166EPSTight2050.push_back(2.28302);
  uncertaintiesRel166EPSTight2050.push_back(1.13891);
  uncertaintiesRel166EPSTight2050.push_back(2.02877);
  //Identification for low ET probes
  //Medium
  efficienciesRel166EPSMediumLowET.push_back(91.16);
  efficienciesRel166EPSMediumLowET.push_back(99.84);
  efficienciesRel166EPSMediumLowET.push_back( 0.00);
  efficienciesRel166EPSMediumLowET.push_back(101.4);
  efficienciesRel166EPSMediumLowET.push_back(96.76);
  efficienciesRel166EPSMediumLowET.push_back(98.11);
  efficienciesRel166EPSMediumLowET.push_back(96.75);
  efficienciesRel166EPSMediumLowET.push_back( 0.00);
  efficienciesRel166EPSMediumLowET.push_back(86.38);
  efficienciesRel166EPSMediumLowET.push_back(84.37);
  uncertaintiesRel166EPSMediumLowET.push_back(11.0);
  uncertaintiesRel166EPSMediumLowET.push_back( 8.5);
  uncertaintiesRel166EPSMediumLowET.push_back( 0.0);
  uncertaintiesRel166EPSMediumLowET.push_back(10.8);
  uncertaintiesRel166EPSMediumLowET.push_back( 6.7);
  uncertaintiesRel166EPSMediumLowET.push_back( 7.0);
  uncertaintiesRel166EPSMediumLowET.push_back( 7.2);
  uncertaintiesRel166EPSMediumLowET.push_back( 0.0);
  uncertaintiesRel166EPSMediumLowET.push_back(10.1);
  uncertaintiesRel166EPSMediumLowET.push_back(10.2);
  //Tight
  efficienciesRel166EPSTightLowET.push_back(91.67);
  efficienciesRel166EPSTightLowET.push_back(100.6);
  efficienciesRel166EPSTightLowET.push_back( 0.00);
  efficienciesRel166EPSTightLowET.push_back(101.1);
  efficienciesRel166EPSTightLowET.push_back(96.88);
  efficienciesRel166EPSTightLowET.push_back(98.14);
  efficienciesRel166EPSTightLowET.push_back(98.23);
  efficienciesRel166EPSTightLowET.push_back( 0.00);
  efficienciesRel166EPSTightLowET.push_back(86.59);
  efficienciesRel166EPSTightLowET.push_back(84.39);
  uncertaintiesRel166EPSTightLowET.push_back(10.9);
  uncertaintiesRel166EPSTightLowET.push_back( 9.6);
  uncertaintiesRel166EPSTightLowET.push_back( 0.0);
  uncertaintiesRel166EPSTightLowET.push_back(10.5);
  uncertaintiesRel166EPSTightLowET.push_back( 6.1);
  uncertaintiesRel166EPSTightLowET.push_back( 6.1);
  uncertaintiesRel166EPSTightLowET.push_back( 9.5);
  uncertaintiesRel166EPSTightLowET.push_back( 0.0);
  uncertaintiesRel166EPSTightLowET.push_back(11.3);
  uncertaintiesRel166EPSTightLowET.push_back( 8.6);
  //For the ET-corrections of the identification scale factors
  //Medium
  ETCorrectionsMediumRel166EPS.push_back( 87.0781);
  ETCorrectionsMediumRel166EPS.push_back( 90.9091);
  ETCorrectionsMediumRel166EPS.push_back( 97.3568);
  ETCorrectionsMediumRel166EPS.push_back(100.453);
  ETCorrectionsMediumRel166EPS.push_back(101.55);
  ETCorrectionsMediumRel166EPS.push_back(101.365);
  ETCorrectionsMediumRel166EPS.push_back(102.087);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(6.00538);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(2.62057);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.93479);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.94788);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.43064);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.40351);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.53891);
  //Tight
  ETCorrectionsTightRel166EPS.push_back( 84.3469);
  ETCorrectionsTightRel166EPS.push_back( 89.3899);
  ETCorrectionsTightRel166EPS.push_back( 97.1825);
  ETCorrectionsTightRel166EPS.push_back(100.33);
  ETCorrectionsTightRel166EPS.push_back(101.319);
  ETCorrectionsTightRel166EPS.push_back(101.238);
  ETCorrectionsTightRel166EPS.push_back(101.552);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(6.52625);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(2.75939);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(1.6303);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(1.29104);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(0.420933);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(0.435997);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(1.05739);
  //For the low ET electrons
  //Medium
  ETCorrectionsMediumRel166EPSFullRange.push_back(0.000/0.9666);
  ETCorrectionsMediumRel166EPSFullRange.push_back(97.36/0.9666);
  ETCorrectionsMediumRel166EPSFullRange.push_back(93.55/0.9666);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(7.25/0.9666);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(7.41/0.9666);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(8.57/0.9666);
  ETCorrectionsMediumRel166EPSFullRange.push_back( 87.0781);
  ETCorrectionsMediumRel166EPSFullRange.push_back( 90.9091);
  ETCorrectionsMediumRel166EPSFullRange.push_back( 97.3568);
  ETCorrectionsMediumRel166EPSFullRange.push_back(100.453);
  ETCorrectionsMediumRel166EPSFullRange.push_back(101.55);
  ETCorrectionsMediumRel166EPSFullRange.push_back(101.365);
  ETCorrectionsMediumRel166EPSFullRange.push_back(102.087);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(9.18078);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(2.62057);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.93479);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.94788);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.43064);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.40351);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.53891);
  //Tight
  ETCorrectionsTightRel166EPSFullRange.push_back(0.000/0.9673);
  ETCorrectionsTightRel166EPSFullRange.push_back(105.8/0.9673);
  ETCorrectionsTightRel166EPSFullRange.push_back(98.8/0.9673);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.24/0.9673);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.43/0.9673);  
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.50/0.9673);
  ETCorrectionsTightRel166EPSFullRange.push_back( 84.3469);
  ETCorrectionsTightRel166EPSFullRange.push_back( 89.3899);
  ETCorrectionsTightRel166EPSFullRange.push_back( 97.1825);
  ETCorrectionsTightRel166EPSFullRange.push_back(100.33);
  ETCorrectionsTightRel166EPSFullRange.push_back(101.319);
  ETCorrectionsTightRel166EPSFullRange.push_back(101.238);
  ETCorrectionsTightRel166EPSFullRange.push_back(101.552);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.1599);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(2.75939);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(1.6303);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(1.29104);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(0.420933);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(0.435997);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(1.05739);
  //Trigger efficiency scale factors
  efficienciesRel166EPSTrigger.push_back(96.5517);
  efficienciesRel166EPSTrigger.push_back(97.3861);
  efficienciesRel166EPSTrigger.push_back(98.4245);
  efficienciesRel166EPSTrigger.push_back(98.6712);
  efficienciesRel166EPSTrigger.push_back(97.7936);
  efficienciesRel166EPSTrigger.push_back(99.7033);
  efficienciesRel166EPSTrigger.push_back(98.9571);
  efficienciesRel166EPSTrigger.push_back(98.4703);
  efficienciesRel166EPSTrigger.push_back(99.3016);
  efficienciesRel166EPSTrigger.push_back(99.1186);
  efficienciesRel166EPSTrigger.push_back(99.2838);
  efficienciesRel166EPSTrigger.push_back(99.2266);
  efficienciesRel166EPSTrigger.push_back(99.709);
  efficienciesRel166EPSTrigger.push_back(99.1478);
  efficienciesRel166EPSTrigger.push_back(99.5733);
  efficienciesRel166EPSTrigger.push_back(98.9866);
  efficienciesRel166EPSTrigger.push_back(99.8198);
  efficienciesRel166EPSTrigger.push_back(97.821);
  efficienciesRel166EPSTrigger.push_back(97.862);
  efficienciesRel166EPSTrigger.push_back(97.901);
  uncertaintiesRel166EPSTrigger.push_back(0.645476);
  uncertaintiesRel166EPSTrigger.push_back(0.588429);
  uncertaintiesRel166EPSTrigger.push_back(0.432384);
  uncertaintiesRel166EPSTrigger.push_back(0.43052);
  uncertaintiesRel166EPSTrigger.push_back(0.579508);
  uncertaintiesRel166EPSTrigger.push_back(0.410817);
  uncertaintiesRel166EPSTrigger.push_back(0.457);
  uncertaintiesRel166EPSTrigger.push_back(0.515013);
  uncertaintiesRel166EPSTrigger.push_back(0.402588);
  uncertaintiesRel166EPSTrigger.push_back(0.418344);
  uncertaintiesRel166EPSTrigger.push_back(0.415669);
  uncertaintiesRel166EPSTrigger.push_back(0.404291);
  uncertaintiesRel166EPSTrigger.push_back(0.407594);
  uncertaintiesRel166EPSTrigger.push_back(0.460203);
  uncertaintiesRel166EPSTrigger.push_back(0.410275);
  uncertaintiesRel166EPSTrigger.push_back(0.53542);
  uncertaintiesRel166EPSTrigger.push_back(0.425722);
  uncertaintiesRel166EPSTrigger.push_back(0.667037);
  uncertaintiesRel166EPSTrigger.push_back(0.426163);
  uncertaintiesRel166EPSTrigger.push_back(0.976323);
  //Reco+trackquality efficiencies
  efficienciesRel166EPSRecoTrkQual.push_back( 97.59);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back( 97.59);
  uncertaintiesRel166EPSRecoTrkQual.push_back(1.84);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(1.84);

  // 2011 data with rel. 17 and MC11a ("CERN Council SF")
  // for technical reasons the values are first stored in float[],
  // then converted to vector<float>
  // Raw Rel17CC Reco+TQ SF
const float Sf_RecoTrkQ_Eta[]={ 102.01, 100.67, 100.97, 100.17, 99.40, 99.16, 99.25, 100.13, 100.73, 100.57, 102.30};
const float Sf_RecoTrkQ_Eta_err[]={ 0.70, 0.57, 0.70, 0.57, 1.11, 1.16, 0.99, 0.55, 0.90, 0.60, 0.71};
  // Raw Rel17CC Identification SF
const float sfLoosePP_Combined_eta[] = {0.978162, 0.989691, 0.9892, 1.00281, 0.993113, 0.994409, 0.995224, 1.00113, 0.9927, 0.990337, 0.98053};
const float errsfLoosePP_Combined_eta[] = {0.0184629, 0.015968, 0.00871837, 0.00385742, 0.00430604, 0.00414063, 0.00707358, 0.003712, 0.00843564, 0.0164764, 0.0178917};
const float sfLoosePP_Jpsi_eta[] = {0.928469, 0.876753, 0.947689, 0.940677, 0.933882, 0.932504, 0.943054, 0.924861, 1.07193, 0.909942, 0.94223};
const float errsfLoosePP_Jpsi_eta[] = {0.0442547, 0.0651155, 0.100367, 0.0459643, 0.0318983, 0.0337912, 0.0316421, 0.0362685, 0.0843151, 0.0566668, 0.0470655};
const float sfLoosePP_Combined_pt[] = {0., 1.04564, 1.02127, 0.950536, 0.956266, 0.985196, 1.00014, 1.00734, 1.00668, 1.00266};
const float errsfLoosePP_Combined_pt[] = {1., 0.0577688, 0.0532959, 0.0192058, 0.0159554, 0.0120122, 0.00643931, 0.00608316, 0.00608894, 0.00670763};
const float sfMediumPP_Combined_eta[] = {0.956077, 0.984517, 0.9933, 0.998451, 0.998374, 1.01566, 0.999115, 0.995048, 0.9972, 0.98697, 0.957895};
const float errsfMediumPP_Combined_eta[] = {0.013147, 0.0124841, 0.00889719, 0.00400233, 0.00446367, 0.00438371, 0.00441865, 0.00390813, 0.0090824, 0.0131541, 0.0154712};
const float sfMediumPP_Jpsi_eta[] = {0.913436, 0.892599, 0.981171, 0.918171, 0.939638, 0.935174, 0.934618, 0.907705, 1.09734, 0.874291, 0.903363};
const float errsfMediumPP_Jpsi_eta[] = {0.0451658, 0.0664901, 0.10942, 0.0456451, 0.0314451, 0.0334607, 0.0309696, 0.0362455, 0.0959015, 0.0564854, 0.0509632};
const float sfMediumPP_Combined_pt[] = {0., 1.06787, 1.0114, 0.949246, 0.940358, 0.974558, 0.994974, 1.0084, 1.00916, 1.0066};
const float errsfMediumPP_Combined_pt[] = {1., 0.0569981, 0.0482483, 0.0216574, 0.0173227, 0.0114571, 0.00633696, 0.00606375, 0.00609331, 0.00677809};
const float sfTightPP_Combined_eta[] = {0.970385, 1.00039, 1.0294, 1.02121, 1.00159, 1.01284, 1.00105, 1.01674, 1.0349, 1.00659, 0.971479};
const float errsfTightPP_Combined_eta[] = {0.0144101, 0.0116894, 0.00947048, 0.00424625, 0.00453523, 0.00451146, 0.00448671, 0.00412469, 0.00987978, 0.0116082, 0.0147539};
const float sfTightPP_Jpsi_eta[] = {0.961754, 0.913472, 1.00017, 0.920565, 0.940924, 0.930151, 0.934168, 0.898207, 1.19533, 0.887737, 0.949335};
const float errsfTightPP_Jpsi_eta[] = {0.0504488, 0.0706803, 0.117501, 0.0490665, 0.0352094, 0.0382616, 0.035019, 0.0403119, 0.109184, 0.0611907, 0.055913};
const float sfTightPP_Combined_pt[] = {0., 1.067, 1.0142, 0.953088, 0.94455, 0.974825, 0.995567, 1.00683, 1.00781, 1.00327};
const float errsfTightPP_Combined_pt[] = {1., 0.0635091, 0.0501458, 0.0228872, 0.0181984, 0.0118053, 0.00635714, 0.00609709, 0.00613041, 0.00679589};
  // Raw Rel17CC Trigger SF and eff
//////////////////////////////////////
/////// trigger e20_medium
//////////////////////////////////////
///// MC Efficiencies  vs Et
const float mcEff_e20_loo1_Et[6]={0.864088, 0.89944, 0.939482, 0.97356, 0.988114, 1.0045};
const float mcEff_e20_med1_Et[6]={0.932642, 0.959893, 0.98097, 0.990982, 0.995659, 1.00106};
const float mcEff_e20_tig1_Et[6]={0.940908, 0.965677, 0.984083, 0.992407, 0.996765, 1.00114};
//////  SF vs Et
const float SF_e20_med1_Et[6]={1.00245, 0.999368, 0.998914, 0.999601, 1.00048, 0.999809};
const float  SF_e20_med1_Et_toterror[6]={0.00597661,0.00298427,0.00154714,0.00105366,0.000587446,0.000577095};
const float SF_e20_tig1_Et[6]={1.0025, 0.995979, 0.997492, 0.999587, 1.0005, 1.00004}; 
const float  SF_e20_tig1_Et_toterror[6]={0.00570995,0.00346066,0.00152597,0.000993986,0.000568664,0.000466317};
//// MC Efficiencies vs eta
const float  mcEff_e20_loo1_eta[20]={0.83473,0.949204,0.944577,0.944431,0.828458,0.964727,0.970421,0.965379,0.955024,0.950054,0.940576,0.955574,0.966962,0.971468,0.963756,0.840913,0.95043,0.944207,0.951791,0.840843};
const float  mcEff_e20_med1_eta[20]={0.872355,0.969229,0.977964,0.973087,0.849374,0.985378,0.986858,0.988503,
0.979529,0.980543,0.974007,0.980666,0.98927,0.988601,0.98482,0.860735,0.978141,0.977674,0.970493,0.875797};
const float  mcEff_e20_tig1_eta[20]={0.886963,0.975401,0.983324,0.979107,0.853658,0.988534,0.990043,
0.991176,0.981667,0.982593,0.975875,0.982813,0.992012,0.99155,0.987791,0.865031,0.983154,0.982475,
0.975659,0.888129};
//// SF vs eta
const float  SF_e20_med1_eta[20]={1.01132,0.988154,0.98865,0.987197,1.03248,1.00244,0.994201,0.981812,1.00942,
0.923586,0.97601,1.00496,1.00136,0.995827,1.00157,1.02941,0.988385,0.983446,0.99408,1.01889};
const float SF_e20_med1_eta_toterror[20]={0.0138381, 0.010173, 0.0101778, 0.0105304, 0.0112517, 0.0100553, 
0.0100303, 0.0100703, 0.0100135, 0.0108708, 0.0105664, 0.0100226, 0.0100271, 0.0100309, 
0.0100894, 0.0114349, 0.0101388, 0.010223, 0.010139, 0.0129749};
const float  SF_e20_tig1_eta[20]={1.00768,0.98708,0.988363,0.986254,1.03095,1.00141,0.993596,0.980864,
1.00878,0.926747,0.975567,1.00439,0.999926,0.995449,1.00165,1.02746,0.989032,0.984948,0.994878,1.01632};
const float SF_e20_tig1_eta_toterror[20]={0.0133341, 0.0102052, 0.0101942, 0.0104491, 0.0114524, 0.0100434, 
0.0100267, 0.010074, 0.0100116, 0.0109426, 0.0105505, 0.010022, 0.0100252, 0.0100245, 0.0100841, 
0.0115263, 0.0101251, 0.0101966, 0.0101346, 0.0131077};
 /////////////////////////////////////////
////////  Trigger e22medium
/////////////////////////////////////////
///// MC efficiencies vs Et

float mcEff_e22_loo1_Et[6]={0., 0.877805, 0.933197, 0.973957, 0.990786, 1.00939}; 
float mcEff_e22_med1_Et[6]={0., 0.938168, 0.97598, 0.990546, 0.996888, 1.00369}; 
float mcEff_e22_tig1_Et[6]={0., 0.945008, 0.978626, 0.991608, 0.997479, 1.00331};
///   SF vs Et
const float SF_e22_med1_Et[6]={0., 1.00106, 0.997813, 1.00152, 1.00105, 0.999557};
const float  SF_e22_med1_Et_toterror[6]={1.,0.00788436,0.00348746,0.00196079,0.00128428,0.000730529};
const float SF_e22_tig1_Et[6]={0., 0.997016, 0.996317, 1.00213, 1.00131, 0.999715};
const float  SF_e22_tig1_Et_toterror[6]={1.,0.00795901,0.00301526,0.00186497,0.00135485,0.000586026};
///  MC efficiencies vs eta
const float  mcEff_e22_loo1_eta[20]={0.80987,0.935756,0.930274,0.933041,0.774026,0.955098,0.962545,0.958426,0.946747,0.940417,0.932461,0.946313,0.957602,0.961308,0.951713,0.796754,0.937367,0.935115,0.938189,0.820834};
const float  mcEff_e22_med1_eta[20]={0.850781,0.957769,0.965943,0.962989,0.796618,0.976691,0.98031,0.981781,0.9729,0.973739,0.966018,0.972751,0.981814,0.98023,0.973673,0.818503,0.967127,0.968771,0.959493,0.855468};
const float  mcEff_e22_tig1_eta[20]={0.866554,0.963796,0.97137,0.968614,0.800967,0.98014,0.983543,
0.98483,0.975155,0.975365,0.967583,0.974961,0.984649,0.983594,0.977511,0.821548,0.972552,
0.973661,0.965113,0.867067};
/// SF vs eta
const float  SF_e22_med1_eta[20]={1.0429,0.993361,0.990606,0.983569,1.07278,1.00356,0.99341,0.983135,
1.00858,0.922439,0.975137,1.00435,1.00485,0.998841,1.00251,1.05195,0.988802,0.974716,0.998945,1.03681};
const float SF_e22_med1_eta_toterror[20]={0.0170972, 0.0104995, 0.010535, 0.0113039, 0.0164707, 0.0101648, 
0.0101347, 0.0102433, 0.0100746, 0.01125, 0.0113536, 0.0100524, 0.010093, 0.0100712, 0.0102185, 
0.0129822, 0.0106089, 0.0108853, 0.0103578, 0.0177669};
const float  SF_e22_tig1_eta[20]={1.04354,0.9923,0.991662,0.983054,1.07381,1.00123,0.993014,0.982268,
1.00798,0.923537,0.977255,1.00411,1.00372,0.997861,1.00232,1.05236,0.988432,0.977005,0.997971,1.03317};
const float SF_e22_tig1_eta_toterror[20]={0.0167165, 0.0104359, 0.0106771, 0.0114704, 0.0167785, 0.0101834, 
0.0101823, 0.0102809, 0.0100571, 0.0113901, 0.011439, 0.0100538, 0.0100842, 0.0100766, 0.0101922, 
0.0131861, 0.0105969, 0.010921, 0.0103099, 0.0168021};
 
//////////////////////////////////////////////
///////  trigger e22vh medium1
////////////////////////////////////////////
//////////  MC efficiencies vs Et

const float mcEff_e22vh_loo1_Et[6]={0., 0.867613, 0.925484, 0.971542, 0.996169, 1.02254};
const float mcEff_e22vh_med1_Et[6]={0., 0.935306, 0.976425, 0.989257, 0.998406, 1.00712};
const float mcEff_e22vh_tig1_Et[6]={0., 0.946316, 0.984708, 0.993602, 1.00022, 1.00552};
///  SF vs Et
const float SF_e22vh_med1_Et[6]={0., 0.976255, 0.990213, 1.00065, 0.999608, 1.00088};
const float  SF_e22vh_med1_Et_toterror[6]={1.,0.00839289,0.00617629,0.00476838,0.00253987,0.00111677};
const float SF_e22vh_tig1_Et[6]={0., 0.971424, 0.986322, 0.998, 0.999025, 1.00115};
const float  SF_e22vh_tig1_Et_toterror[6]={1.,0.00840815,0.00611464,0.00418324,0.00233599,0.00100704};

/// MC efficiencies vs eta
const float  mcEff_e22vh_loo1_eta[20]={0.70588,0.850925,0.885496,0.88813,0.708504,0.898123,0.918406,0.916392,0.898269,0.815695,0.81267,0.897603,0.920762,0.919525,0.902682,0.726406,0.872302,0.884492,0.855155,0.72748};

const float  mcEff_e22vh_med1_eta[20]={0.806431,0.934885,0.950239,0.936858,0.783633,0.948002,0.9631,0.960062,
0.94823,0.909058,0.905434,0.948098,0.967413,0.963304,0.946392,0.793201,0.933363,0.945796,0.94782,0.828012};
const float  mcEff_e22vh_tig1_eta[20]={0.840512,0.949399,0.959665,0.95542,0.801649,0.969443,0.972239,
0.970461,0.959291,0.956205,0.956761,0.960181,0.977778,0.973809,0.965969,0.804513,0.953652,0.956252,
0.957446,0.847437};
// SF vs eta
const float  SF_e22vh_med1_eta[20]={0.984147,0.980365,0.970567,0.984624,0.97203,1.01202,0.999753,0.991051,
1.02403,0.990278,1.02291,1.01998,1.00556,1.00275,1.00858,1.02538,0.993383,0.965577,0.973939,0.953943};
const float SF_e22vh_med1_eta_toterror[20]={0.0305649, 0.0112576, 0.0118168, 0.0127388, 0.0190914, 0.0112009, 
0.0103002, 0.0105304, 0.0102212, 0.0128408, 0.0125235, 0.0103448, 0.0104107, 0.010408, 0.0108098, 
0.016527, 0.0125478, 0.0111028, 0.0107405, 0.0223372};
const float  SF_e22vh_tig1_eta[20]={0.964279,0.977745,0.975187,0.978157,0.961275,0.999952,0.995353,0.984084,
1.01769,0.956689,0.985402,1.01267,0.998659,0.99684,0.997916,1.0241,0.985229,0.970021,0.976287,0.947805};
const float SF_e22vh_tig1_eta_toterror[20]={0.033432, 0.0204931, 0.0208272, 0.0212713, 0.026605, 0.0203521, 
0.0201282, 0.020186, 0.0200903, 0.0213355, 0.0208187, 0.020117, 0.0201333, 0.0201497, 0.0203914, 
0.0245803, 0.0209498, 0.0205718, 0.0205268, 0.0272371};


 copyToVector(Sf_RecoTrkQ_Eta, 11, efficienciesRel17CCRecoTrkQual, 1.);
 copyToVector(Sf_RecoTrkQ_Eta_err, 11, uncertaintiesRel17CCRecoTrkQual, 1.);
 
 // Identification eta for probes between 15 and 50 GeV
 copyToVector(sfLoosePP_Combined_eta, 11, efficienciesRel17CCLoosePP1550);
 copyToVector(errsfLoosePP_Combined_eta, 11, uncertaintiesRel17CCLoosePP1550);
 copyToVector(sfMediumPP_Combined_eta, 11, efficienciesRel17CCMediumPP1550);
 copyToVector(errsfMediumPP_Combined_eta, 11, uncertaintiesRel17CCMediumPP1550);
 copyToVector(sfTightPP_Combined_eta, 11, efficienciesRel17CCTightPP1550);
 copyToVector(errsfTightPP_Combined_eta, 11, uncertaintiesRel17CCTightPP1550);
 //Identification eta for low ET probes
 copyToVector(sfLoosePP_Jpsi_eta, 11, efficienciesRel17CCLoosePP415);
 copyToVector(errsfLoosePP_Jpsi_eta, 11, uncertaintiesRel17CCLoosePP415);
 copyToVector(sfMediumPP_Jpsi_eta, 11, efficienciesRel17CCMediumPP415);
 copyToVector(errsfMediumPP_Jpsi_eta, 11, uncertaintiesRel17CCMediumPP415);
 copyToVector(sfTightPP_Jpsi_eta, 11, efficienciesRel17CCTightPP415);
 copyToVector(errsfTightPP_Jpsi_eta, 11, uncertaintiesRel17CCTightPP415);
 // ET correction
 copyToVector(sfLoosePP_Combined_pt, 10, ETCorrectionsRel17CCLoosePP);
 copyToVector(errsfLoosePP_Combined_pt, 10, uncertaintiesETCorrectionsRel17CCLoosePP);
 copyToVector(sfMediumPP_Combined_pt, 10, ETCorrectionsRel17CCMediumPP);
 copyToVector(errsfMediumPP_Combined_pt, 10, uncertaintiesETCorrectionsRel17CCMediumPP);
 copyToVector(sfTightPP_Combined_pt, 10, ETCorrectionsRel17CCTightPP);
 copyToVector(errsfTightPP_Combined_pt, 10, uncertaintiesETCorrectionsRel17CCTightPP);
 
 // Trigger efficiencies
 // e20_medium B-J
  copyToVector(SF_e20_med1_eta, 20, efficienciesRel17CCe20_mediumMediumPP);
  copyToVector(SF_e20_med1_eta_toterror, 20, uncertaintiesRel17CCe20_mediumMediumPP);
  copyToVector(SF_e20_med1_Et, 6, efficienciesRel17CCe20_mediumMediumPPET);
  copyToVector(SF_e20_med1_Et_toterror, 6, uncertaintiesRel17CCe20_mediumMediumPPET);

  copyToVector(mcEff_e20_med1_eta, 20, MCefficienciesRel17CCe20_mediumMediumPP);
  copyToVector(mcEff_e20_med1_Et, 6, MCefficienciesRel17CCe20_mediumMediumPPET);

  copyToVector(SF_e20_tig1_eta, 20, efficienciesRel17CCe20_mediumTightPP);
  copyToVector(SF_e20_tig1_eta_toterror, 20, uncertaintiesRel17CCe20_mediumTightPP);
  copyToVector(SF_e20_tig1_Et, 6, efficienciesRel17CCe20_mediumTightPPET);
  copyToVector(SF_e20_tig1_Et_toterror, 6, uncertaintiesRel17CCe20_mediumTightPPET);

  copyToVector(mcEff_e20_tig1_eta, 20, MCefficienciesRel17CCe20_mediumTightPP);
  copyToVector(mcEff_e20_tig1_Et, 6, MCefficienciesRel17CCe20_mediumTightPPET);

  copyToVector(mcEff_e20_loo1_eta, 20, MCefficienciesRel17CCe20_mediumLoosePP);
  copyToVector(mcEff_e20_loo1_Et, 6, MCefficienciesRel17CCe20_mediumLoosePPET);


  // e22_medium K
  copyToVector(SF_e22_med1_eta, 20, efficienciesRel17CCe22_mediumMediumPP);
  copyToVector(SF_e22_med1_eta_toterror, 20, uncertaintiesRel17CCe22_mediumMediumPP);
  copyToVector(SF_e22_med1_Et, 6, efficienciesRel17CCe22_mediumMediumPPET);
  copyToVector(SF_e22_med1_Et_toterror, 6, uncertaintiesRel17CCe22_mediumMediumPPET);

  copyToVector(mcEff_e22_med1_eta, 20, MCefficienciesRel17CCe22_mediumMediumPP);
  copyToVector(mcEff_e22_med1_Et, 6, MCefficienciesRel17CCe22_mediumMediumPPET);

  copyToVector(SF_e22_tig1_eta, 20, efficienciesRel17CCe22_mediumTightPP);
  copyToVector(SF_e22_tig1_eta_toterror, 20, uncertaintiesRel17CCe22_mediumTightPP);
  copyToVector(SF_e22_tig1_Et, 6, efficienciesRel17CCe22_mediumTightPPET);
  copyToVector(SF_e22_tig1_Et_toterror, 6, uncertaintiesRel17CCe22_mediumTightPPET);

  copyToVector(mcEff_e22_tig1_eta, 20, MCefficienciesRel17CCe22_mediumTightPP);
  copyToVector(mcEff_e22_tig1_Et, 6, MCefficienciesRel17CCe22_mediumTightPPET);

  copyToVector(mcEff_e22_loo1_eta, 20, MCefficienciesRel17CCe22_mediumLoosePP);
  copyToVector(mcEff_e22_loo1_Et, 6, MCefficienciesRel17CCe22_mediumLoosePPET);


  // e22vh_medium1 L-M
  copyToVector(SF_e22vh_med1_eta, 20, efficienciesRel17CCe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_eta_toterror, 20, uncertaintiesRel17CCe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_Et, 6, efficienciesRel17CCe22vh_medium1MediumPPET);
  copyToVector(SF_e22vh_med1_Et_toterror, 6, uncertaintiesRel17CCe22vh_medium1MediumPPET);

  copyToVector(mcEff_e22vh_med1_eta, 20, MCefficienciesRel17CCe22vh_medium1MediumPP);
  copyToVector(mcEff_e22vh_med1_Et, 6, MCefficienciesRel17CCe22vh_medium1MediumPPET);

  copyToVector(SF_e22vh_tig1_eta, 20, efficienciesRel17CCe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_eta_toterror, 20, uncertaintiesRel17CCe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_Et, 6, efficienciesRel17CCe22vh_medium1TightPPET);
  copyToVector(SF_e22vh_tig1_Et_toterror, 6, uncertaintiesRel17CCe22vh_medium1TightPPET);

  copyToVector(mcEff_e22vh_tig1_eta, 20, MCefficienciesRel17CCe22vh_medium1TightPP);
  copyToVector(mcEff_e22vh_tig1_Et, 6, MCefficienciesRel17CCe22vh_medium1TightPPET);


  copyToVector(mcEff_e22vh_loo1_eta, 20, MCefficienciesRel17CCe22vh_medium1LoosePP);
  copyToVector(mcEff_e22vh_loo1_Et, 6, MCefficienciesRel17CCe22vh_medium1LoosePPET);



  // 2011 data with rel. 17 and MC11a/b/c G4 FullSim ("Moriond SF")
const float sfLoosePP_Combined_Moriond_eta[] = {0.983727, 0.992645, 0.9892, 1.00425, 0.994238, 0.99583, 0.995955, 1.00307, 0.9907, 0.994061, 0.984558};
const float errsfLoosePP_Combined_Moriond_eta[] = {0.0190335, 0.0158532, 0.00838153, 0.00394904, 0.00432761, 0.00425582, 0.00429913, 0.00357115, 0.0079806, 0.0164504, 0.0178845};
const float sfLoosePP_Jpsi_Moriond_eta[] = {0.928469, 0.876753, 0.947689, 0.940677, 0.933882, 0.932504, 0.943054, 0.924861, 1.07193, 0.909942, 0.94223};
const float errsfLoosePP_Jpsi_Moriond_eta[] = {0.0442547, 0.0651155, 0.100367, 0.0459643, 0.0318983, 0.0337912, 0.0316421, 0.0362685, 0.0843151, 0.0566668, 0.0470655};

const float sfLoosePP_Combined_Moriond_pt[] = {0., 1.04564, 1.02127, 0.952484, 0.958019, 0.986368, 1.00122, 1.00889, 1.00725, 1.00164};
const float errsfLoosePP_Combined_Moriond_pt[] = {1., 0.0584143, 0.0539949, 0.0220578, 0.0191925, 0.011236, 0.00582933, 0.00661811, 0.00532031, 0.00642074};


const float sfMedium_Combined_Moriond_eta[] = {0.980389, 0.984739, 0.987, 0.995007, 0.992199, 0.999625, 0.992335, 0.991465, 0.9896, 0.990364, 0.982351};
const float errsfMedium_Combined_Moriond_eta[] = {0.0187827, 0.01235, 0.00807527, 0.00418403, 0.0044624, 0.00429145, 0.00446588, 0.00383729, 0.00793473, 0.0138236, 0.0179662};
const float sfMedium_Jpsi_Moriond_eta[] = {0.942597, 0.974524, 1.08641, 1.01779, 0.975308, 0.97664, 0.982064, 0.971485, 1.07755, 0.981647, 0.910087};
const float errsfMedium_Jpsi_Moriond_eta[] = {0.042, 0.085, 0.091, 0.049, 0.028, 0.03, 0.028, 0.049, 0.095, 0.06, 0.049};
const float sfMedium_Combined_Moriond_pt[] = {0., 1.02283, 0.980082, 0.954413, 0.950976, 0.982899, 0.998409, 1.00998, 1.00825, 1.00394};
const float errsfMedium_Combined_Moriond_pt[] = {1., 0.0562668, 0.0482163, 0.0225589, 0.0206991, 0.00776764, 0.00602853, 0.00685129, 0.00541447, 0.00645837};


const float sfMediumPP_Combined_Moriond_eta[] = {0.967144, 0.990423, 0.9926, 0.998398, 1.00044, 1.01609, 1.00082, 0.995716, 0.9922, 0.994387, 0.967035};
const float errsfMediumPP_Combined_Moriond_eta[] = {0.0152302, 0.0126984, 0.00866025, 0.00419149, 0.00447153, 0.00445324, 0.00447221, 0.00391371, 0.00843564, 0.0148758, 0.0164505};
const float sfMediumPP_Jpsi_Moriond_eta[] = {0.913436, 0.892599, 0.981171, 0.918171, 0.939638, 0.935174, 0.934618, 0.907705, 1.09734, 0.874291, 0.903363};
const float errsfMediumPP_Jpsi_Moriond_eta[] = {0.0451658, 0.0664901, 0.10942, 0.0456451, 0.0314451, 0.0334607, 0.0309696, 0.0362455, 0.0959015, 0.0564854, 0.0509632};
const float sfMediumPP_Combined_Moriond_pt[] = {0., 1.06787, 1.0114, 0.952377, 0.942511, 0.977914, 0.995868, 1.00973, 1.01033, 1.0053};
const float errsfMediumPP_Combined_Moriond_pt[] = {1., 0.0576523, 0.0490194, 0.0245875, 0.0217735, 0.00821047, 0.00626121, 0.00701316, 0.00551153, 0.00653204};


const float sfTightPP_Combined_Moriond_eta[] = {0.987679, 1.00704, 1.027, 1.02319, 1.00345, 1.01298, 1.00215, 1.0204, 1.0276, 1.01574, 0.985344};
const float errsfTightPP_Combined_Moriond_eta[] = {0.00770372, 0.0116938, 0.00920923, 0.00471006, 0.00464143, 0.0051928, 0.0045983, 0.00446568, 0.0090824, 0.0140841, 0.0168724};
const float sfTightPP_Jpsi_Moriond_eta[] = {0.961754, 0.913472, 1.00017, 0.920565, 0.940924, 0.930151, 0.934168, 0.898207, 1.19533, 0.887737, 0.949335};
const float errsfTightPP_Jpsi_Moriond_eta[] = {0.0504488, 0.0706803, 0.117501, 0.0490665, 0.0352094, 0.0382616, 0.035019, 0.0403119, 0.109184, 0.0611907, 0.055913};
const float sfTightPP_Combined_Moriond_pt[] = {0., 1.067, 1.0142, 0.957576, 0.949916, 0.979087, 0.996125, 1.00763, 1.00897, 1.0016};
const float errsfTightPP_Combined_Moriond_pt[] = {1., 0.0640969, 0.0508882, 0.0259462, 0.0232052, 0.0085569, 0.00700376, 0.00719543, 0.00578251, 0.00667381};



  // reco+TQ are by choice identical to the "Cern council" results!
  copyToVector(Sf_RecoTrkQ_Eta, 11, efficienciesRel17MoriondRecoTrkQual, 1.);
  copyToVector(Sf_RecoTrkQ_Eta_err, 11, uncertaintiesRel17MoriondRecoTrkQual, 1.);
 
  // Identification eta for probes between 15 and 50 GeV
  copyToVector(sfLoosePP_Combined_Moriond_eta, 11, efficienciesRel17MoriondLoosePP1550);
  copyToVector(errsfLoosePP_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondLoosePP1550);
  copyToVector(sfMedium_Combined_Moriond_eta, 11, efficienciesRel17MoriondMedium1550);
  copyToVector(errsfMedium_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondMedium1550);
  copyToVector(sfMediumPP_Combined_Moriond_eta, 11, efficienciesRel17MoriondMediumPP1550);
  copyToVector(errsfMediumPP_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondMediumPP1550);
  copyToVector(sfTightPP_Combined_Moriond_eta, 11, efficienciesRel17MoriondTightPP1550);
  copyToVector(errsfTightPP_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondTightPP1550);
  //Identification eta for low ET probes
  copyToVector(sfLoosePP_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondLoosePP415);
  copyToVector(errsfLoosePP_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondLoosePP415);
  copyToVector(sfMedium_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondMedium415);
  copyToVector(errsfMedium_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondMedium415);
  copyToVector(sfMediumPP_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondMediumPP415);
  copyToVector(errsfMediumPP_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondMediumPP415);
  copyToVector(sfTightPP_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondTightPP415);
  copyToVector(errsfTightPP_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondTightPP415);
  // ET correction
  copyToVector(sfLoosePP_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondLoosePP);
  copyToVector(errsfLoosePP_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondLoosePP);
  copyToVector(sfMedium_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondMedium);
  copyToVector(errsfMedium_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondMedium);
  copyToVector(sfMediumPP_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondMediumPP);
  copyToVector(errsfMediumPP_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondMediumPP);
  copyToVector(sfTightPP_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondTightPP);
  copyToVector(errsfTightPP_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondTightPP);
 
 // Trigger efficiencies
  // by choice identical to the "Cern council" results - to be updated!
 // e20_medium B-J
  copyToVector(SF_e20_med1_eta, 20, efficienciesRel17Morionde20_mediumMediumPP);
  copyToVector(SF_e20_med1_eta_toterror, 20, uncertaintiesRel17Morionde20_mediumMediumPP);
  copyToVector(SF_e20_med1_Et, 6, efficienciesRel17Morionde20_mediumMediumPPET);
  copyToVector(SF_e20_med1_Et_toterror, 6, uncertaintiesRel17Morionde20_mediumMediumPPET);

  copyToVector(mcEff_e20_med1_eta, 20, MCefficienciesRel17Morionde20_mediumMediumPP);
  copyToVector(mcEff_e20_med1_Et, 6, MCefficienciesRel17Morionde20_mediumMediumPPET);

  copyToVector(SF_e20_tig1_eta, 20, efficienciesRel17Morionde20_mediumTightPP);
  copyToVector(SF_e20_tig1_eta_toterror, 20, uncertaintiesRel17Morionde20_mediumTightPP);
  copyToVector(SF_e20_tig1_Et, 6, efficienciesRel17Morionde20_mediumTightPPET);
  copyToVector(SF_e20_tig1_Et_toterror, 6, uncertaintiesRel17Morionde20_mediumTightPPET);

  copyToVector(mcEff_e20_tig1_eta, 20, MCefficienciesRel17Morionde20_mediumTightPP);
  copyToVector(mcEff_e20_tig1_Et, 6, MCefficienciesRel17Morionde20_mediumTightPPET);

  copyToVector(mcEff_e20_loo1_eta, 20, MCefficienciesRel17Morionde20_mediumLoosePP);
  copyToVector(mcEff_e20_loo1_Et, 6, MCefficienciesRel17Morionde20_mediumLoosePPET);


  // e22_medium K
  copyToVector(SF_e22_med1_eta, 20, efficienciesRel17Morionde22_mediumMediumPP);
  copyToVector(SF_e22_med1_eta_toterror, 20, uncertaintiesRel17Morionde22_mediumMediumPP);
  copyToVector(SF_e22_med1_Et, 6, efficienciesRel17Morionde22_mediumMediumPPET);
  copyToVector(SF_e22_med1_Et_toterror, 6, uncertaintiesRel17Morionde22_mediumMediumPPET);

  copyToVector(mcEff_e22_med1_eta, 20, MCefficienciesRel17Morionde22_mediumMediumPP);
  copyToVector(mcEff_e22_med1_Et, 6, MCefficienciesRel17Morionde22_mediumMediumPPET);

  copyToVector(SF_e22_tig1_eta, 20, efficienciesRel17Morionde22_mediumTightPP);
  copyToVector(SF_e22_tig1_eta_toterror, 20, uncertaintiesRel17Morionde22_mediumTightPP);
  copyToVector(SF_e22_tig1_Et, 6, efficienciesRel17Morionde22_mediumTightPPET);
  copyToVector(SF_e22_tig1_Et_toterror, 6, uncertaintiesRel17Morionde22_mediumTightPPET);

  copyToVector(mcEff_e22_tig1_eta, 20, MCefficienciesRel17Morionde22_mediumTightPP);
  copyToVector(mcEff_e22_tig1_Et, 6, MCefficienciesRel17Morionde22_mediumTightPPET);

  copyToVector(mcEff_e22_loo1_eta, 20, MCefficienciesRel17Morionde22_mediumLoosePP);
  copyToVector(mcEff_e22_loo1_Et, 6, MCefficienciesRel17Morionde22_mediumLoosePPET);


  // e22vh_medium1 L-M
  copyToVector(SF_e22vh_med1_eta, 20, efficienciesRel17Morionde22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_eta_toterror, 20, uncertaintiesRel17Morionde22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_Et, 6, efficienciesRel17Morionde22vh_medium1MediumPPET);
  copyToVector(SF_e22vh_med1_Et_toterror, 6, uncertaintiesRel17Morionde22vh_medium1MediumPPET);

  copyToVector(mcEff_e22vh_med1_eta, 20, MCefficienciesRel17Morionde22vh_medium1MediumPP);
  copyToVector(mcEff_e22vh_med1_Et, 6, MCefficienciesRel17Morionde22vh_medium1MediumPPET);

  copyToVector(SF_e22vh_tig1_eta, 20, efficienciesRel17Morionde22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_eta_toterror, 20, uncertaintiesRel17Morionde22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_Et, 6, efficienciesRel17Morionde22vh_medium1TightPPET);
  copyToVector(SF_e22vh_tig1_Et_toterror, 6, uncertaintiesRel17Morionde22vh_medium1TightPPET);

  copyToVector(mcEff_e22vh_tig1_eta, 20, MCefficienciesRel17Morionde22vh_medium1TightPP);
  copyToVector(mcEff_e22vh_tig1_Et, 6, MCefficienciesRel17Morionde22vh_medium1TightPPET);


  copyToVector(mcEff_e22vh_loo1_eta, 20, MCefficienciesRel17Morionde22vh_medium1LoosePP);
  copyToVector(mcEff_e22vh_loo1_Et, 6, MCefficienciesRel17Morionde22vh_medium1LoosePPET);





  // 2011 data with rel. 17 and MC11a/b/c AFII ("Moriond SF")
const float sfLoosePP_Combined_Moriond_AFII_eta[] = {0.990782, 0.984431, 1.0473, 1.00123, 0.999666, 0.998526, 1.00158, 1.0011, 1.0422, 0.985262, 0.990121};
const float errsfLoosePP_Combined_Moriond_AFII_eta[] = {0.0206951, 0.0160733, 0.00920923, 0.00394199, 0.00433968, 0.00427718, 0.00434684, 0.00354063, 0.00832827, 0.0159877, 0.0185366};
const float sfLoosePP_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float errsfLoosePP_Jpsi_Moriond_AFII_eta[] = {1.00045, 1.00045, 1.00045, 1.00045, 1.00045, 1.00045, 1.00045, 1.00045, 1.00045, 1.00045, 1.00045};
const float sfLoosePP_Combined_Moriond_AFII_pt[] = {0., 1, 1, 0.94271, 0.952656, 0.985613, 1.00222, 1.01042, 1.00727, 1.00107}; 
const float errsfLoosePP_Combined_Moriond_AFII_pt[] = {1.00005, 1.00005, 1.00005, 0.0233024, 0.0191394, 0.0122292, 0.00584829, 0.00683113, 0.00533981, 0.0064188};

const float sfMedium_Combined_Moriond_AFII_eta[] = {0.993134, 0.973003, 1.0982, 0.997358, 1.0093, 1.01942, 1.01548, 0.998222, 1.0948, 0.97789, 0.992183};
const float errsfMedium_Combined_Moriond_AFII_eta[] = {0.0209136, 0.0138568, 0.00914549, 0.00418142, 0.00855772, 0.00459038, 0.0089923, 0.00384268, 0.00838153, 0.0146522, 0.0187527};
const float sfMedium_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float errsfMedium_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float sfMedium_Combined_Moriond_AFII_pt[] = {0., 1, 1, 0.94817, 0.938252, 0.976575, 0.997979, 1.01049, 1.01606, 1.01232};
const float errsfMedium_Combined_Moriond_AFII_pt[] = {1.00005, 1.00005, 1.00005, 0.0240558, 0.0205825, 0.00782366, 0.00611601, 0.00752741, 0.00545777, 0.00649618};

const float sfMediumPP_Combined_Moriond_AFII_eta[] = {0.985436, 0.979718, 1.1627, 1.00117, 1.01279, 1.02311, 1.01838, 1.00231, 1.154, 0.982628, 0.983485};
const float errsfMediumPP_Combined_Moriond_AFII_eta[] = {0.0182925, 0.0136327, 0.00914549, 0.00414907, 0.00754674, 0.00472048, 0.00843378, 0.00387936, 0.00866025, 0.0150049, 0.0178148};
const float sfMediumPP_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float errsfMediumPP_Jpsi_Moriond_AFII_eta[] = {1.00042, 1.00042, 1.00042, 1.00042, 1.00042, 1.00042, 1.00042, 1.00042, 1.00042, 1.00042, 1.00042};
const float sfMediumPP_Combined_Moriond_AFII_pt[] = {0., 1, 1, 0.933341, 0.932939, 0.973115, 0.996294, 1.00996, 1.01637, 1.0127};
const float errsfMediumPP_Combined_Moriond_AFII_pt[] = {1.00005, 1.00005, 1.00005, 0.0264644, 0.0217461, 0.00823226, 0.0063997, 0.00788586, 0.00555109, 0.00656721};


const float sfTightPP_Combined_Moriond_AFII_eta[] = {1.00032, 0.997738, 1.2036, 1.0263, 1.01493, 1.03012, 1.01828, 1.02795, 1.1994, 1.00119, 0.999282};
const float errsfTightPP_Combined_Moriond_AFII_eta[] = {0.0180631, 0.013335, 0.0120847, 0.00472454, 0.00465211, 0.00587229, 0.00919978, 0.00446069, 0.00974115, 0.0128199, 0.0176403};
const float sfTightPP_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float errsfTightPP_Jpsi_Moriond_AFII_eta[] = {1.00054, 1.00054, 1.00054, 1.00054, 1.00054, 1.00054, 1.00054, 1.00054, 1.00054, 1.00054, 1.00054};
const float sfTightPP_Combined_Moriond_AFII_pt[] = {0., 1, 1, 0.942954, 0.938838, 0.974907, 0.995292, 1.0075, 1.01237, 1.01166};
const float errsfTightPP_Combined_Moriond_AFII_pt[] = {1.00005, 1.00005, 1.00005, 0.0279654, 0.0231078, 0.00852419, 0.00717741, 0.00797035, 0.00573223, 0.00678183};

  // reco+TQ are by choice identical to the "Cern council" results!
  copyToVector(Sf_RecoTrkQ_Eta, 11, efficienciesRel17MoriondAFIIRecoTrkQual, 1.);
  copyToVector(Sf_RecoTrkQ_Eta_err, 11, uncertaintiesRel17MoriondAFIIRecoTrkQual, 1.);
 
 // Identification eta for probes between 15 and 50 GeV
  copyToVector(sfLoosePP_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIILoosePP1550);
  copyToVector(errsfLoosePP_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIILoosePP1550);
  copyToVector(sfMedium_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMedium1550);
  copyToVector(errsfMedium_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMedium1550);
  copyToVector(sfMediumPP_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMediumPP1550);
  copyToVector(errsfMediumPP_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMediumPP1550);
  copyToVector(sfTightPP_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIITightPP1550);
  copyToVector(errsfTightPP_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIITightPP1550);
  //Identification eta for low ET probes
  copyToVector(sfLoosePP_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIILoosePP415);
  copyToVector(errsfLoosePP_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIILoosePP415);
  copyToVector(sfMedium_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMedium415);
  copyToVector(errsfMedium_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMedium415);
  copyToVector(sfMediumPP_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMediumPP415);
  copyToVector(errsfMediumPP_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMediumPP415);
  copyToVector(sfTightPP_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIITightPP415);
  copyToVector(errsfTightPP_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIITightPP415);
  // ET correction
  copyToVector(sfLoosePP_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIILoosePP);
  copyToVector(errsfLoosePP_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIILoosePP);
  copyToVector(sfMedium_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIIMedium);
  copyToVector(errsfMedium_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIIMedium);
  copyToVector(sfMediumPP_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIIMediumPP);
  copyToVector(errsfMediumPP_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIIMediumPP);
  copyToVector(sfTightPP_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIITightPP);
  copyToVector(errsfTightPP_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIITightPP);
  
 // Trigger efficiencies
  // for now by choice identical to the "Cern council" results - to be updated!
 // e20_medium B-J
  copyToVector(SF_e20_med1_eta, 20, efficienciesRel17MoriondAFIIe20_mediumMediumPP);
  copyToVector(SF_e20_med1_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe20_mediumMediumPP);
  copyToVector(SF_e20_med1_Et, 6, efficienciesRel17MoriondAFIIe20_mediumMediumPPET);
  copyToVector(SF_e20_med1_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe20_mediumMediumPPET);

  copyToVector(mcEff_e20_med1_eta, 20, MCefficienciesRel17MoriondAFIIe20_mediumMediumPP);
  copyToVector(mcEff_e20_med1_Et, 6, MCefficienciesRel17MoriondAFIIe20_mediumMediumPPET);

  copyToVector(SF_e20_tig1_eta, 20, efficienciesRel17MoriondAFIIe20_mediumTightPP);
  copyToVector(SF_e20_tig1_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe20_mediumTightPP);
  copyToVector(SF_e20_tig1_Et, 6, efficienciesRel17MoriondAFIIe20_mediumTightPPET);
  copyToVector(SF_e20_tig1_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe20_mediumTightPPET);

  copyToVector(mcEff_e20_tig1_eta, 20, MCefficienciesRel17MoriondAFIIe20_mediumTightPP);
  copyToVector(mcEff_e20_tig1_Et, 6, MCefficienciesRel17MoriondAFIIe20_mediumTightPPET);

  copyToVector(mcEff_e20_loo1_eta, 20, MCefficienciesRel17MoriondAFIIe20_mediumLoosePP);
  copyToVector(mcEff_e20_loo1_Et, 6, MCefficienciesRel17MoriondAFIIe20_mediumLoosePPET);


  // e22_medium K
  copyToVector(SF_e22_med1_eta, 20, efficienciesRel17MoriondAFIIe22_mediumMediumPP);
  copyToVector(SF_e22_med1_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22_mediumMediumPP);
  copyToVector(SF_e22_med1_Et, 6, efficienciesRel17MoriondAFIIe22_mediumMediumPPET);
  copyToVector(SF_e22_med1_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22_mediumMediumPPET);

  copyToVector(mcEff_e22_med1_eta, 20, MCefficienciesRel17MoriondAFIIe22_mediumMediumPP);
  copyToVector(mcEff_e22_med1_Et, 6, MCefficienciesRel17MoriondAFIIe22_mediumMediumPPET);

  copyToVector(SF_e22_tig1_eta, 20, efficienciesRel17MoriondAFIIe22_mediumTightPP);
  copyToVector(SF_e22_tig1_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22_mediumTightPP);
  copyToVector(SF_e22_tig1_Et, 6, efficienciesRel17MoriondAFIIe22_mediumTightPPET);
  copyToVector(SF_e22_tig1_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22_mediumTightPPET);

  copyToVector(mcEff_e22_tig1_eta, 20, MCefficienciesRel17MoriondAFIIe22_mediumTightPP);
  copyToVector(mcEff_e22_tig1_Et, 6, MCefficienciesRel17MoriondAFIIe22_mediumTightPPET);

  copyToVector(mcEff_e22_loo1_eta, 20, MCefficienciesRel17MoriondAFIIe22_mediumLoosePP);
  copyToVector(mcEff_e22_loo1_Et, 6, MCefficienciesRel17MoriondAFIIe22_mediumLoosePPET);


  // e22vh_medium1 L-M
  copyToVector(SF_e22vh_med1_eta, 20, efficienciesRel17MoriondAFIIe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_Et, 6, efficienciesRel17MoriondAFIIe22vh_medium1MediumPPET);
  copyToVector(SF_e22vh_med1_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPPET);

  copyToVector(mcEff_e22vh_med1_eta, 20, MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPP);
  copyToVector(mcEff_e22vh_med1_Et, 6, MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPPET);

  copyToVector(SF_e22vh_tig1_eta, 20, efficienciesRel17MoriondAFIIe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_Et, 6, efficienciesRel17MoriondAFIIe22vh_medium1TightPPET);
  copyToVector(SF_e22vh_tig1_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22vh_medium1TightPPET);

  copyToVector(mcEff_e22vh_tig1_eta, 20, MCefficienciesRel17MoriondAFIIe22vh_medium1TightPP);
  copyToVector(mcEff_e22vh_tig1_Et, 6, MCefficienciesRel17MoriondAFIIe22vh_medium1TightPPET);


  copyToVector(mcEff_e22vh_loo1_eta, 20, MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePP);
  copyToVector(mcEff_e22vh_loo1_Et, 6, MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePPET);


}

std::pair<float,float> egammaSFclass::scaleFactor(float eta, float ET, int set, int range, int rel, bool etcorrection) {

   std::vector<float> * vectEff=0;
   std::vector<float> * vectUnc=0;
   std::vector<float> * vectEtaBinning=0;

   if (rel==7) { //release 17 for 2011 data and AFII MC11a/b/c, "Moriond recommendations"
     // range is ignored here
     vectEtaBinning = &m_11Etabins;
     if (set == 0 || set == 2 || set == 3 || set > 22) {
       std::cout << "egammaSFclass: ERROR : only Reco+TrackQuality, Medium, Loose++, Medium++, Tight++, and 3 single electron triggers exist" << std::endl;
       return make_pair(-1.,-1.);
     }
     else if (set==4) {//Reco + track quality requirements
       // this has implicit ET dependence, so don't confuse the user
       etcorrection = false;
       vectEff = &efficienciesRel17MoriondAFIIRecoTrkQual;
       vectUnc = &uncertaintiesRel17MoriondAFIIRecoTrkQual;
       if (ET<15000.) {
	 float eff=1.,unc=0.05;
	 if (fabs(eta)<1.37) {
	   eff=1.;unc=0.02;
	 }
	 return make_pair(eff,unc);
       }
     }
     else if (set==1) {//Medium
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIILoosePP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIILoosePP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIILoosePP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIILoosePP415;
       }
     }
     else if (set==5) {//Loose++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIILoosePP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIILoosePP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIILoosePP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIILoosePP415;
       }
     }
     else if (set==6) {//Medium++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIIMediumPP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIIMediumPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIIMediumPP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIIMediumPP415;
       }
     }
     else if (set==7) {//Tight++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIITightPP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIITightPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIITightPP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIITightPP415;
       }
     }
     else if (set==8) {//e20_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe20_mediumMediumPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe20_mediumMediumPP;
     }
     else if (set==9) {//e20_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe20_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==10) {//e20_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe20_mediumTightPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe20_mediumTightPP;
     }
     else if (set==11) {//e20_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe20_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==12) {//e22_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22_mediumMediumPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22_mediumMediumPP;
     }
     else if (set==13) {//e22_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==14) {//e22_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22_mediumTightPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22_mediumTightPP;
     }
     else if (set==15) {//e22_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==16) {//e22vh_medium1 trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22vh_medium1MediumPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPP;
     }
     else if (set==17) {//e22vh_medium1 MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPP;
       vectUnc = 0; // no error
     }
     else if (set==18) {//e22vh_medium1 trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22vh_medium1TightPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22vh_medium1TightPP;
     }
     else if (set==19) {//e22vh_medium1 MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22vh_medium1TightPP;
       vectUnc = 0; // no error
     }
     else if (set==20) {//e20_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe20_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==21) {//e22_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==22) {//e22vh_medium1 MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePP;
       vectUnc = 0; // no error
     }
   }
   else if (rel==6) { //release 17 for 2011 data and G4 FullSim MC11a/b/c, "Moriond recommendations"
     // range is ignored here
     vectEtaBinning = &m_11Etabins;
     if (set == 0 || set == 2 || set == 3 || set > 22) {
       std::cout << "egammaSFclass: ERROR : only Reco+TrackQuality, Medium, Loose++, Medium++, Tight++, and 3 single electron triggers exist" << std::endl;
       return make_pair(-1.,-1.);
     }
     else if (set==4) {//Reco + track quality requirements
       // this has implicit ET dependence, so don't confuse the user
       etcorrection = false;
       vectEff = &efficienciesRel17MoriondRecoTrkQual;
       vectUnc = &uncertaintiesRel17MoriondRecoTrkQual;
       if (ET<15000.) {
	 float eff=1.,unc=0.05;
	 if (fabs(eta)<1.37) {
	   eff=1.;unc=0.02;
	 }
	 return make_pair(eff,unc);
       }
     }
     else if (set==1) {//Medium
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondLoosePP1550;
	 vectUnc = &uncertaintiesRel17MoriondLoosePP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondLoosePP415;
	 vectUnc = &uncertaintiesRel17MoriondLoosePP415;
       }
     }
     else if (set==5) {//Loose++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondLoosePP1550;
	 vectUnc = &uncertaintiesRel17MoriondLoosePP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondLoosePP415;
	 vectUnc = &uncertaintiesRel17MoriondLoosePP415;
       }
     }
     else if (set==6) {//Medium++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondMediumPP1550;
	 vectUnc = &uncertaintiesRel17MoriondMediumPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondMediumPP415;
	 vectUnc = &uncertaintiesRel17MoriondMediumPP415;
       }
     }
     else if (set==7) {//Tight++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondTightPP1550;
	 vectUnc = &uncertaintiesRel17MoriondTightPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondTightPP415;
	 vectUnc = &uncertaintiesRel17MoriondTightPP415;
       }
     }
     else if (set==8) {//e20_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde20_mediumMediumPP;
       vectUnc = &uncertaintiesRel17Morionde20_mediumMediumPP;
     }
     else if (set==9) {//e20_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde20_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==10) {//e20_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde20_mediumTightPP;
       vectUnc = &uncertaintiesRel17Morionde20_mediumTightPP;
     }
     else if (set==11) {//e20_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde20_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==12) {//e22_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22_mediumMediumPP;
       vectUnc = &uncertaintiesRel17Morionde22_mediumMediumPP;
     }
     else if (set==13) {//e22_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==14) {//e22_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22_mediumTightPP;
       vectUnc = &uncertaintiesRel17Morionde22_mediumTightPP;
     }
     else if (set==15) {//e22_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==16) {//e22vh_medium1 trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22vh_medium1MediumPP;
       vectUnc = &uncertaintiesRel17Morionde22vh_medium1MediumPP;
     }
     else if (set==17) {//e22vh_medium1 MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22vh_medium1MediumPP;
       vectUnc = 0; // no error
     }
     else if (set==18) {//e22vh_medium1 trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22vh_medium1TightPP;
       vectUnc = &uncertaintiesRel17Morionde22vh_medium1TightPP;
     }
     else if (set==19) {//e22vh_medium1 MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22vh_medium1TightPP;
       vectUnc = 0; // no error
     }
     else if (set==20) {//e20_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde20_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==21) {//e22_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==22) {//e22vh_medium1 MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22vh_medium1LoosePP;
       vectUnc = 0; // no error
     }
   }
   else if (rel==5) { //release 17 for 2011 data and MC11a, "CERN council recommendations"
     // range is ignored here
     vectEtaBinning = &m_11Etabins;
     if (set < 4 || set > 22) {
       std::cout << "egammaSFclass: ERROR : only Reco+TrackQuality, IsEM++ menu, and 3 single electron triggers exist" << std::endl;
       return make_pair(-1.,-1.);
     }
     else if (set==4) {//Reco + track quality requirements
       // this has implicit ET dependence, so don't confuse the user
       etcorrection = false;
       vectEff = &efficienciesRel17CCRecoTrkQual;
       vectUnc = &uncertaintiesRel17CCRecoTrkQual;
       if (ET<15000.) {
	 float eff=1.,unc=0.05;
	 if (fabs(eta)<1.37) {
	   eff=1.;unc=0.02;
	 }
	 return make_pair(eff,unc);
       }
     }
     else if (set==5) {//Loose++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17CCLoosePP1550;
	 vectUnc = &uncertaintiesRel17CCLoosePP1550;
       } else {
	 vectEff = &efficienciesRel17CCLoosePP415;
	 vectUnc = &uncertaintiesRel17CCLoosePP415;
       }
     }
     else if (set==6) {//Medium++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17CCMediumPP1550;
	 vectUnc = &uncertaintiesRel17CCMediumPP1550;
       } else {
	 vectEff = &efficienciesRel17CCMediumPP415;
	 vectUnc = &uncertaintiesRel17CCMediumPP415;
       }
     }
     else if (set==7) {//Tight++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17CCTightPP1550;
	 vectUnc = &uncertaintiesRel17CCTightPP1550;
       } else {
	 vectEff = &efficienciesRel17CCTightPP415;
	 vectUnc = &uncertaintiesRel17CCTightPP415;
       }
     }
     else if (set==8) {//e20_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe20_mediumMediumPP;
       vectUnc = &uncertaintiesRel17CCe20_mediumMediumPP;
     }
     else if (set==9) {//e20_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe20_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==10) {//e20_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe20_mediumTightPP;
       vectUnc = &uncertaintiesRel17CCe20_mediumTightPP;
     }
     else if (set==11) {//e20_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe20_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==12) {//e22_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22_mediumMediumPP;
       vectUnc = &uncertaintiesRel17CCe22_mediumMediumPP;
     }
     else if (set==13) {//e22_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==14) {//e22_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22_mediumTightPP;
       vectUnc = &uncertaintiesRel17CCe22_mediumTightPP;
     }
     else if (set==15) {//e22_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==16) {//e22vh_medium1 trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22vh_medium1MediumPP;
       vectUnc = &uncertaintiesRel17CCe22vh_medium1MediumPP;
     }
     else if (set==17) {//e22vh_medium1 MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22vh_medium1MediumPP;
       vectUnc = 0; // no error
     }
     else if (set==18) {//e22vh_medium1 trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22vh_medium1TightPP;
       vectUnc = &uncertaintiesRel17CCe22vh_medium1TightPP;
     }
     else if (set==19) {//e22vh_medium1 MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22vh_medium1TightPP;
       vectUnc = 0; // no error
     }
     else if (set==20) {//e20_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe20_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==21) {//e22_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==22) {//e22vh_medium1 MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22vh_medium1LoosePP;
       vectUnc = 0; // no error
     }
   }
   else if (rel==4) { //release 16.6 numbers estimated from 2011 data, "EPS recommendations" including the low ET region
     vectEtaBinning = &m_FineEtabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>4) {
	 std::cout << "egammaSFclass: ERROR : only Medium, Tight and trigger scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 if (ET>=15000.) {
	   vectEff = &efficienciesRel166EPSMedium2050;
	   vectUnc = &uncertaintiesRel166EPSMedium2050;
	 }
	 else {
	   vectEtaBinning = &m_Etabins;
	   vectEff = &efficienciesRel166EPSMediumLowET;
	   vectUnc = &uncertaintiesRel166EPSMediumLowET;
	 }
       }
       else if (set==2) {//Tight
	 if (ET>=15000.) {
	   vectEff = &efficienciesRel166EPSTight2050;
	   vectUnc = &uncertaintiesRel166EPSTight2050;
	 }
	 else {
	   vectEtaBinning = &m_Etabins;
	   vectEff = &efficienciesRel166EPSTightLowET;
	   vectUnc = &uncertaintiesRel166EPSTightLowET;
	 }
       }
       else if (set==3) {//Trigger
	 vectEff = &efficienciesRel166EPSTrigger;
	 vectUnc = &uncertaintiesRel166EPSTrigger;
       }
       else if (set==4) {//Reco + track quality requirements
	 vectEff = &efficienciesRel166EPSRecoTrkQual;
	 vectUnc = &uncertaintiesRel166EPSRecoTrkQual;
	 if (ET<15000.) {
	   float eff=1.,unc=0.05;
	   if (fabs(eta)<1.37) {
	     eff=1.;unc=0.02;
	   }
	   return make_pair(eff,unc);
	 }
       }
     }//endif 20-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   } 
   else if (rel==3) { //release 16.6 numbers estimated from 2011 data, "EPS recommendations"
     vectEtaBinning = &m_FineEtabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>4) {
	 std::cout << "egammaSFclass: ERROR : only Medium, Tight and trigger scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel166EPSMedium2050;
	 vectUnc = &uncertaintiesRel166EPSMedium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel166EPSTight2050;
	 vectUnc = &uncertaintiesRel166EPSTight2050;
       }
       else if (set==3) {//Trigger
	 vectEff = &efficienciesRel166EPSTrigger;
	 vectUnc = &uncertaintiesRel166EPSTrigger;
       }
       else if (set==4) {//Reco + track quality requirements
	 vectEff = &efficienciesRel166EPSRecoTrkQual;
	 vectUnc = &uncertaintiesRel166EPSRecoTrkQual;
       }
     }//endif 20-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   } 
   else if (rel==2) { //release 16.6 numbers estimated from 2010 data
     vectEtaBinning = &m_Etabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel166Data2010Medium2050;
	 vectUnc = &uncertaintiesRel166Data2010Medium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel166Data2010Tight2050;
	 vectUnc = &uncertaintiesRel166Data2010Tight2050;
       }
     }//endif 20-50 GeV
     else if (range==1) { //30-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel166Data2010Medium3050;
	 vectUnc = &uncertaintiesRel166Data2010Medium3050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel166Data2010Tight3050;
	 vectUnc = &uncertaintiesRel166Data2010Tight3050;
       }
     }//endif 30-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   } 
   else if (rel==1) { //release 16 numbers
     vectEtaBinning = &m_Etabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel16Medium2050;
	 vectUnc = &uncertaintiesRel16Medium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel16Tight2050;
	 vectUnc = &uncertaintiesRel16Tight2050;
       }
     }//endif 20-50 GeV
     else if (range==1) { //30-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel16Medium3050;
	 vectUnc = &uncertaintiesRel16Medium3050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel16Tight3050;
	 vectUnc = &uncertaintiesRel16Tight3050;
       }
     }//endif 30-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   }
   else { //release 15 numbers
     vectEtaBinning = &m_Etabins;
     if (range==0) { //20-50 GeV region
       if (set==0) {//Loose
	 vectEff = &efficienciesRel15Loose2050;
	 vectUnc = &uncertaintiesRel15Loose2050;
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel15Medium2050;
	 vectUnc = &uncertaintiesRel15Medium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel15Tight2050;
	 vectUnc = &uncertaintiesRel15Tight2050;
       }
       else {
	 std::cout << "egammaSFclass: ERROR : invalid set of cuts" << std::endl;
	 return make_pair(-1.,-1.);
       }
     }//endif 20-50 GeV
     else if (range==1) { //30-50 GeV region
       if (set==0) {//Loose
	 vectEff = &efficienciesRel15Loose3050;
	 vectUnc = &uncertaintiesRel15Loose3050;
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel15Medium3050;
	 vectUnc = &uncertaintiesRel15Medium3050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel15Tight3050;
	 vectUnc = &uncertaintiesRel15Tight3050;
       }
       else {
	 std::cout << "egammaSFclass: ERROR : invalid set of cuts" << std::endl;
	 return make_pair(-1.,-1.);
       }
     }//endif 30-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   }//endif rel15

   //Choice of the eta bin
   int ietabin=-1;
   while (ietabin<((int)vectEtaBinning->size()-1) && eta>=vectEtaBinning->at(ietabin+1)) ietabin++;
   if (ietabin<0 || ietabin>=((int)vectEtaBinning->size()-1)) {
     std::cout << "egammaSFclass: ERROR : given eta value outside range of existing scale factors" << std::endl;
     return make_pair(-1.,-1.);
   }


   float effvseta = vectEff->at(ietabin)/100.;
   float uncvseta = 0.;
   if (vectUnc)
     uncvseta = vectUnc->at(ietabin)/100.;

   float eff = effvseta;
   float unc = uncvseta;

   if (etcorrection) {
     std::pair<float,float> corr = etCorrection(ET, set, rel);
     if (corr.first <= 0)
       unc = 1.;
     else
       unc = eff*corr.first * sqrt( unc*unc/(eff*eff) + corr.second*corr.second/(corr.first*corr.first) );
     eff *= corr.first;
   }

   return make_pair(eff,unc);
}


//Returns the ET-correction factor (and uncertainty) to the scale factor for the correspond ET bin 
std::pair<float,float> egammaSFclass::etCorrection(float ET, int set, int rel) {
  //for backport of rel16 SF-ET-dependence to rel15
  if (rel==0) rel=1;  
  
  std::vector<float> * vectCorr=0;
  std::vector<float> * vectUncCorr=0;
  std::vector<float> * vectETBinning=0;
  vectETBinning = &m_ETbinsFullRange;
  
  if (rel==1) {
    vectETBinning = &m_ETbins;
    if (set==1) {
      vectCorr = &ETCorrectionsMediumRel16;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel16;
    }
    else if (rel==2) { // tight
      vectCorr = &ETCorrectionsTightRel16;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel16;
    }
  }
  else if (rel==2) {
    vectETBinning = &m_ETbins;
    if (set==1) { //Medium
      vectCorr = &ETCorrectionsMediumRel166Data2010;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel166Data2010;
    }
    else if (set==2) {
      vectCorr = &ETCorrectionsTightRel166Data2010;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel166Data2010;
    }
  } 
  else if (rel==3) {
    vectETBinning = &m_ETbins;
    if (set==1) { //Medium
      vectCorr = &ETCorrectionsMediumRel166EPS;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel166EPS;
    }
    else if (set==2) {
      vectCorr = &ETCorrectionsTightRel166EPS;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel166EPS;
    }
  }
  else if (rel==4) {
    if (set==1) { //Medium
      vectCorr = &ETCorrectionsMediumRel166EPSFullRange;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel166EPSFullRange;
    }
    else if (set==2) {
      vectCorr = &ETCorrectionsTightRel166EPSFullRange;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel166EPSFullRange;
    }
  }
  else if (rel==5) {// Rel17CC
    if (set==5) { // Loose++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17CCLoosePP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17CCLoosePP;
    }
    else if (set==6) { // Medium++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17CCMediumPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17CCMediumPP;
    }
    else if (set==7) { // Tight++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17CCTightPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17CCTightPP;
    }
    else if (set==8) { // e20_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe20_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17CCe20_mediumMediumPPET;
    }
    else if (set==9) { // e20_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe20_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==10) { // e20_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe20_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17CCe20_mediumTightPPET;
    }
    else if (set==11) { // e20_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe20_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==12) { // e22_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17CCe22_mediumMediumPPET;
    }
    else if (set==13) { // e22_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==14) { // e22_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17CCe22_mediumTightPPET;
    }
    else if (set==15) { // e22_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==16) { // e22vh_medium1 trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22vh_medium1MediumPPET;
      vectUncCorr = &uncertaintiesRel17CCe22vh_medium1MediumPPET;
    }
    else if (set==17) { // e22vh_medium1 trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22vh_medium1MediumPPET;
      vectUncCorr = 0;
    }
    else if (set==18) { // e22vh_medium1 trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22vh_medium1TightPPET;
      vectUncCorr = &uncertaintiesRel17CCe22vh_medium1TightPPET;
    }
    else if (set==19) { // e22vh_medium1 trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22vh_medium1TightPPET;
      vectUncCorr = 0;
    }
    else if (set==20) { // e20_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe20_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==21) { // e22_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==22) { // e22vh_medium1 trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22vh_medium1LoosePPET;
      vectUncCorr = 0;
    }
  }
  else if (rel==6) {// Rel17 Moriond G4FullSim
    if (set==1) { // Medium
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondMedium;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondMedium;
    }
    else if (set==5) { // Loose++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondLoosePP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondLoosePP;
    }
    else if (set==6) { // Medium++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondMediumPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondMediumPP;
    }
    else if (set==7) { // Tight++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondTightPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondTightPP;
    }
    else if (set==8) { // e20_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde20_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17Morionde20_mediumMediumPPET;
    }
    else if (set==9) { // e20_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde20_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==10) { // e20_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde20_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17Morionde20_mediumTightPPET;
    }
    else if (set==11) { // e20_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde20_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==12) { // e22_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22_mediumMediumPPET;
    }
    else if (set==13) { // e22_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==14) { // e22_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22_mediumTightPPET;
    }
    else if (set==15) { // e22_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==16) { // e22vh_medium1 trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22vh_medium1MediumPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22vh_medium1MediumPPET;
    }
    else if (set==17) { // e22vh_medium1 trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22vh_medium1MediumPPET;
      vectUncCorr = 0;
    }
    else if (set==18) { // e22vh_medium1 trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22vh_medium1TightPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22vh_medium1TightPPET;
    }
    else if (set==19) { // e22vh_medium1 trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22vh_medium1TightPPET;
      vectUncCorr = 0;
    }
    else if (set==20) { // e20_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde20_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==21) { // e22_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==22) { // e22vh_medium1 trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22vh_medium1LoosePPET;
      vectUncCorr = 0;
    }
  }
  else if (rel==7) {// Rel17 Moriond AFII
    if (set==1) { // Medium
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIIMedium;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIIMedium;
    }
    else if (set==5) { // Loose++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIILoosePP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIILoosePP;
    }
    else if (set==6) { // Medium++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIIMediumPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIIMediumPP;
    }
    else if (set==7) { // Tight++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIITightPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIITightPP;
    }
    else if (set==8) { // e20_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe20_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe20_mediumMediumPPET;
    }
    else if (set==9) { // e20_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe20_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==10) { // e20_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe20_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe20_mediumTightPPET;
    }
    else if (set==11) { // e20_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe20_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==12) { // e22_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22_mediumMediumPPET;
    }
    else if (set==13) { // e22_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==14) { // e22_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22_mediumTightPPET;
    }
    else if (set==15) { // e22_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==16) { // e22vh_medium1 trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22vh_medium1MediumPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPPET;
    }
    else if (set==17) { // e22vh_medium1 trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPPET;
      vectUncCorr = 0;
    }
    else if (set==18) { // e22vh_medium1 trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22vh_medium1TightPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22vh_medium1TightPPET;
    }
    else if (set==19) { // e22vh_medium1 trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22vh_medium1TightPPET;
      vectUncCorr = 0;
    }
    else if (set==20) { // e20_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe20_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==21) { // e22_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==22) { // e22vh_medium1 trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePPET;
      vectUncCorr = 0;
    }
  }


  if (vectCorr == 0) { // catch all missing cases
    std::cout << "egammaSFclass: ERROR : ET-correction factors not implemented for given selection" << std::endl;
    return make_pair(-1.,-1.);
  }

  int iETbin=-1;
  while (iETbin < int(vectETBinning->size()-1)
	 && ET>=vectETBinning->at(iETbin+1))
    iETbin++;
  if (iETbin<0 || iETbin>= int(vectETBinning->size()-1)) {
    std::cout << "egammaSFclass: ERROR : given ET value (" 
	      << ET << ") outside range of existing ET-correction factors" << std::endl;
    return make_pair(-1.,-1.);
  }

  float eff=vectCorr->at(iETbin)/100.;
  float unc=0;
  if (vectUncCorr)
    unc=vectUncCorr->at(iETbin)/100.;
  return make_pair(eff, unc);
}

std::pair<float,float> egammaSFclass::scaleFactorForward(float eta, int set)
{
  if (set == 0) {
    // Forward loose
    float abseta = std::abs(eta);
    if (2.5 <= abseta && abseta <= 3.16)
      return make_pair(1.010, 0.039);
    else if (3.35 <= abseta && abseta <= 4.9)
      return make_pair(1.020, 0.059);
    else {
      std::cout << "egammaSFclass: ERROR : Out of eta range for forward electrons" << std::endl;
      return make_pair(-1.,-1.);
    }

  } else if (set == 2) {
    // Forward tight
    float abseta = std::abs(eta);
    if (2.5 <= abseta && abseta <= 3.16)
      return make_pair(0.893, 0.047);
    else if (3.35 <= abseta && abseta <= 4.9)
      return make_pair(0.940, 0.061);
    else {
      std::cout << "egammaSFclass: ERROR : Out of eta range for forward electrons" << std::endl;
      return make_pair(-1.,-1.);
    }

  } else {
    std::cout << "egammaSFclass: ERROR : Forward electrons only have Loose or Tight cuts" << std::endl;
    return make_pair(-1.,-1.);
  }

}

void egammaSFclass::copyToVector(const float *myarray, int n, std::vector<float> &dest, double renorm)
{
  for (int i = 0; i < n; i++)
    dest.push_back(myarray[i]*renorm);
}

}
}
