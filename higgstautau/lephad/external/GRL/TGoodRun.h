
/**********************************************************************************
 * Class  : TGoodRun                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Max Baak <mbaak@cern.ch> - CERN, Switzerland                              *
 **********************************************************************************/

#ifndef __MyTGoodRun__
#define __MyTGoodRun__

#include "TLumiBlockRange.h"
#include "TObject.h"
#include <vector>

#ifdef __CINT__
#pragma link C++ class vector<AnalysisFramework::External::TLumiBlockRange>+;
#else
template class std::vector<AnalysisFramework::External::TLumiBlockRange>;
#endif

namespace AnalysisFramework
{
namespace External
{
   class TGoodRun : public std::vector< TLumiBlockRange >, public TObject {

   public:

      TGoodRun() ;
      TGoodRun( const Int_t& runnr );
      virtual ~TGoodRun();

      TGoodRun(const TGoodRun& other) ;
      TGoodRun& operator=(const TGoodRun& other) ;

      const TGoodRun GetOverlapWith(const TGoodRun& other) const ;
      const TGoodRun GetSumWith(const TGoodRun& other) const ;
      const TGoodRun GetPartOnlyIn(const TGoodRun& other) const ;
      const TGoodRun GetPartNotIn(const TGoodRun& other) const ;

      Bool_t IsEmpty() const;
      Bool_t HasLB( const Int_t& lumiblocknr )  const;
      std::vector<TLumiBlockRange>::iterator Find( const Int_t& lumiblocknr );
      std::vector< TLumiBlockRange >::const_iterator Find( const Int_t& lumiblocknr ) const;
      inline Int_t GetRunNumber() const { return m_runnr; }
      inline void SetRunNumber( const Int_t& runnr ) { m_runnr=runnr; }

      void Summary() const ;

      void Sort();
      void Compress();
      void AddLB( const Int_t& lumiblocknr );

   private:

      Int_t m_runnr;

      // sorter function for lumiblock ranges
      struct SorterL2H {
        SorterL2H () {}
        bool operator() (const TLumiBlockRange& p1, const TLumiBlockRange& p2) {
          return (p1.Begin()<p2.Begin());
        }
      };

      ClassDef(TGoodRun,1)
   };
}
}

#endif

