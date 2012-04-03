#include "_cleaning.h"
#include <math.h>


bool is_bad_jet(BADLEVEL criteria,
	double quality, double NegE,
	double emf,     double hecf,
	double time,    double fmax,
	double eta,     double chf ,
	double HecQ,    double LArQmean )
  {

    if( criteria==LooseMinusBad || criteria==LooseBad || criteria==MediumBad || criteria==TightBad )
      {
	// HEC spike
	if( hecf>0.5 && fabs(HecQ)>0.5 && LArQmean>0.8)         return true;
	if( fabs(NegE)>60000./*MeV*/)                          return true;
	// EM coherent noise
	if( emf>0.95 && fabs(quality)>0.8 && LArQmean>0.8 && fabs(eta)<2.8 )   return true;
	// Cosmics and Beam background
	if( emf<0.05 && chf<0.05 && fabs(eta)<2. )             return true;
	if( emf<0.05 && fabs(eta)>2. )                         return true;
	if( fmax>0.99&&fabs(eta)<2)                            return true;
      }



    if( criteria==LooseBad || criteria==MediumBad || criteria==TightBad )
      {
	// HEC spike
	if( hecf>0.5 && fabs(HecQ)>0.5 )                       return true;
	// EM coherent noise
	if( emf>0.95 && fabs(quality)>0.8 && fabs(eta)<2.8 )   return true;
	// Cosmics and Beam background
	if( fabs(time)>25. )                                   return true;
      }


    if(criteria==MediumBad || criteria==TightBad)
      {
	// HEC spike
	if( hecf> 1-fabs(HecQ))                               return true;
	// EM coherent noise
	if( emf>0.9 && fabs(quality)>0.8 && fabs(eta)<2.8 )   return true;
	// Cosmics and Beam background
	if( fabs(time)>10. )                                  return true;
	if( emf < 0.05 && chf < 0.1  && fabs(eta)<2. )        return true;
	if( emf > 0.95 && chf < 0.05 && fabs(eta)<2. )        return true;
      }

    if( criteria==TightBad)
      {
	// EM coherent noise
	if ( fabs(quality)>0.95 )                   return true ;
	if (emf>0.98 && fabs(quality)>0.05)         return true ;
	// Cosmics and Beam background
	if (emf<0.1 && chf < 0.2 && fabs(eta)<2.5 ) return true;
	if (emf>0.9 && chf < 0.1 && fabs(eta)<2.5 ) return true;
	if (chf<0.01 && fabs(eta)<2.5 )             return true;
	if (emf<0.1 && fabs(eta)>2.5 )              return true;


      }
    return false;
  }
