// CLASS DERIVED FROM HARVARD GROUP's CLASS //

#ifndef MYETAPHI_BINNING_CLASS_H
#define MYETAPHI_BINNING_CLASS_H

#include "TLorentzVector.h"
#include <iostream>

namespace AnalysisFramework
{
namespace External
{

class EtaPhiBinning {

    public:
        EtaPhiBinning(){ };
        virtual ~EtaPhiBinning() {;};

        virtual int bin(const TLorentzVector *m) const { 
            if(m->Eta() > 0)
                return this->symmetricBin(m) + 11;

            return 11 - this->symmetricBin(m);
        }
        virtual int symmetricBin(const TLorentzVector *m) const;

        enum binregion{ binUNKNOWN=0, bin1BARRELLG=1, bin1BARRELSM=2, bin2BARREL=3, binFEET=4,
			binTRANSITION=5, binENDCAPLG=6, binENDCAPSM=7, binBEE=8, binFORWARDLG=9, 
			binFORWARDSM=10};

    private:
        int getCoarseNSector(double phi) const;
        int getSector(double phi) const;
        int getECSector(double phi) const;

};

}
}

#endif
