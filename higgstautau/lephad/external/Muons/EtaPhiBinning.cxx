#include <cmath>
using std::fabs;

#include "EtaPhiBinning.h"

namespace AnalysisFramework
{
namespace External
{

int EtaPhiBinning::getSector(double phi) const {

  int sector = -9;
  if (phi>2.905 || phi<=-2.905) sector = 9;
  else if (phi>2.59 && phi<=2.905) sector = 8;
  else if (phi>2.12 && phi<=2.59) sector = 7;
  else if (phi>1.805 && phi<=2.12) sector = 6;
  else if (phi>1.335 && phi<=1.805) sector = 5;
  else if (phi>1.02 && phi<=1.335) sector = 4;
  else if (phi>0.55 && phi<=1.02) sector = 3;
  else if (phi>0.235 && phi<=0.55) sector = 2;
  else if (phi>-0.235 && phi<=0.235) sector = 1;
  else if (phi>-0.55 && phi<=-0.235) sector = 16;
  else if (phi>-1.02 && phi<=-0.55) sector = 15;
  else if (phi>-1.335 && phi<=-1.02) sector = 14;
  else if (phi>-1.805 && phi<=-1.335) sector = 13;
  else if (phi>-2.12 && phi<=-1.805) sector = 12;
  else if (phi>-2.59 && phi<=-2.12) sector = 11;
  else if (phi>-2.905 && phi<=-2.59) sector = 10;
  return sector;

}

int EtaPhiBinning::getECSector(double phi) const {

  int sector = -9;
  if (phi>3.011 || phi<=-3.011) sector = 9;
  else if (phi>2.487 && phi<=3.011) sector = 8;
  else if (phi>2.225 && phi<=2.487) sector = 7;
  else if (phi>1.702 && phi<=2.225) sector = 6;
  else if (phi>1.440 && phi<=1.702) sector = 5;
  else if (phi>0.916 && phi<=1.440) sector = 4;
  else if (phi>0.655 && phi<=0.916) sector = 3;
  else if (phi>0.131 && phi<=0.655) sector = 2;
  else if (phi>-0.131 && phi<=0.131) sector = 1;
  else if (phi>-0.655 && phi<=-0.131) sector = 16;
  else if (phi>-0.916 && phi<=-0.655) sector = 15;
  else if (phi>-1.440 && phi<=-0.916) sector = 14;
  else if (phi>-1.702 && phi<=-1.440) sector = 13;
  else if (phi>-2.225 && phi<=-1.702) sector = 12;
  else if (phi>-2.487 && phi<=-2.225) sector = 11;
  else if (phi>-3.011 && phi<=-2.487) sector = 10;
  return sector;

}

int EtaPhiBinning::getCoarseNSector(double phi) const {

    if ( (fabs(phi)>=0.18 && fabs(phi)<0.285)
            || (fabs(phi)>=0.5 && fabs(phi)<0.605)
            || (fabs(phi)>=0.965 && fabs(phi)<1.07)
            || (fabs(phi)>=1.285 && fabs(phi)<1.39)
            || (fabs(phi)>=1.75 && fabs(phi)<1.855)
            || (fabs(phi)>=2.07 && fabs(phi)<2.175)
            || (fabs(phi)>=2.535 && fabs(phi)<2.64)
            || (fabs(phi)>=2.855 && fabs(phi)<2.96) )
        return 2;
    else return 1;

}

int EtaPhiBinning::symmetricBin(const TLorentzVector* mst) const
{

  //Region bin based on eta and phi
  double mu_phi = mst->Phi();
  double mu_eta = mst->Eta();
  int mu_sector = getSector(mu_phi);
  bool isSmSect = mu_sector%2==0 ? true : false;
  int mu_nsect = getCoarseNSector(mu_phi);
  int mu_ec_sector = getECSector(mu_phi);
  bool isSmECSect = mu_ec_sector%2==0 ? true : false;
  
  //Feet
  if ( fabs(mu_eta)<0.97
       && ((mu_phi<-1.01 && mu_phi>-1.36)
	   || (mu_phi<-1.79 && mu_phi>-2.14)) )
    return binFEET;
  if ( ( fabs(mu_eta)>0.51 && fabs(mu_eta)<0.81)
       && (mu_phi<-1.36 && mu_phi>-1.79) )
    return binFEET;
  //To put tiny regions of overlap into feet 
  //if ( mu_sector==13 && mst->barrelSectors()>1) return binFEET;

  //Non-Feet Barrel
  if ( fabs(mu_eta)<0.97 ) {
    if (mu_nsect==2) return bin2BARREL;
    else if (mu_nsect==1) {
      if (isSmSect) return bin1BARRELSM;
      else return bin1BARRELLG;
    }
  }
  
  //Transition and BIS78
  if ( fabs(mu_eta)>=0.97 && fabs(mu_eta)<1.11 )
    return binTRANSITION;
  
  else if ( fabs(mu_eta)>=1.11 && fabs(mu_eta)<1.19 ) {
    if ( isSmSect || mu_nsect==2 ) return binTRANSITION;
    else {
      if (isSmECSect) return binENDCAPSM;
      else return binENDCAPLG;
    }
    
  } else if ( fabs(mu_eta)>=1.19 && fabs(mu_eta)<1.25 ) {
    if ( isSmSect && mu_nsect==1 ) return binTRANSITION;
    else if ( isSmECSect ) return binENDCAPSM;
    else return binENDCAPLG;
    
  }
  
  //BEE
  if ( fabs(mu_eta)>=1.42 && fabs(mu_eta)<1.72
       && mu_nsect==1 && isSmSect )
    return binBEE;
  
  //Endcap
  if ( fabs(mu_eta)>=1.25 && fabs(mu_eta)<1.97 ) {
    if (isSmECSect) return binENDCAPSM;
    else return binENDCAPLG;
  }
  
  //Forward
  if ( fabs(mu_eta)>=1.97 && fabs(mu_eta)<2.01 ) {
    if ( mu_nsect==1 && !isSmSect ) {
      if (isSmECSect) return binENDCAPSM;
      else return binENDCAPLG;
    } else return binFORWARDSM;
  }
  
  if ( fabs(mu_eta)>=2.01 ) {
    if ( mu_nsect>1 || isSmSect ) return binFORWARDSM;
    else {
      if (isSmECSect) return binFORWARDSM;
      else return binFORWARDLG;
    }
  }

  //control flow should now be here!!
  std::cout<<"Control flow should not be here: MuSelectionToolsLine="<<__LINE__<<std::endl;
  return binUNKNOWN;

}

}
}
