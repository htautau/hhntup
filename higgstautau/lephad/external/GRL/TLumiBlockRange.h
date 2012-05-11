
/**********************************************************************************
 * Class  : TLumiBlockRange                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Max Baak <mbaak@cern.ch> - CERN, Switzerland                              *
 **********************************************************************************/

#ifndef __MyTLumiBlockRange__
#define __MyTLumiBlockRange__

#include "TObject.h"
#include <vector>

namespace AnalysisFramework
{
namespace External
{
   class TLumiBlockRange : public TObject {
      
   public:

      TLumiBlockRange();      
      TLumiBlockRange( const Int_t& start, const Int_t& end = 2147483647 );
      virtual ~TLumiBlockRange();

      TLumiBlockRange(const TLumiBlockRange& other) ;
      TLumiBlockRange& operator=(const TLumiBlockRange& other) ;

      const TLumiBlockRange GetOverlapWith(const TLumiBlockRange& other) const ; 
      const std::vector<AnalysisFramework::External::TLumiBlockRange> GetPartOnlyIn(const TLumiBlockRange& other) const ; 
      const std::vector<AnalysisFramework::External::TLumiBlockRange> GetPartNotIn(const TLumiBlockRange& other) const ;

      Bool_t Contains( const Int_t& number ) const;

      inline Int_t Begin() const { return m_begin; }
      inline Int_t End() const { return m_end; }
      inline Bool_t IsEmpty() const { return (Begin()>End()); }

      inline void SetBegin(const Int_t& begin) { m_begin=begin; }
      inline void SetEnd(const Int_t& end) { m_end=end; }

      void Summary() const ;
      
   private:

      Int_t m_begin;
      Int_t m_end;

      ClassDef(TLumiBlockRange,1)
   };
}
}

#endif

