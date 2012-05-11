
#include "TLumiBlockRange.h"
#include "TMsgLogger.h"

// ClassImp(TLumiBlockRange)

namespace AnalysisFramework
{
namespace External
{

TLumiBlockRange::TLumiBlockRange()
  : TObject()
  , m_begin(0)
  , m_end(2147483647)
{
}


TLumiBlockRange::TLumiBlockRange( const Int_t& begin, const Int_t& end )
  : TObject()
  , m_begin(begin)
  , m_end(end)
{
}


TLumiBlockRange::~TLumiBlockRange()
{
}


TLumiBlockRange::TLumiBlockRange(const TLumiBlockRange& other)
  : TObject(other)
  , m_begin(other.m_begin)
  , m_end(other.m_end)
{
}


TLumiBlockRange&
TLumiBlockRange::operator=(const TLumiBlockRange& other)
{
  if (&other==this) {
    return *this ;
  } 
  m_begin = other.m_begin ;
  m_end = other.m_end ;

  return *this ;
}


const TLumiBlockRange 
TLumiBlockRange::GetOverlapWith(const TLumiBlockRange& other) const
{
  Int_t begin = ( this->Begin()>=other.Begin() ? this->Begin() : other.Begin() );
  Int_t end   = ( this->End()<=other.End() ? this->End() : other.End() );
  return TLumiBlockRange(begin,end);
}


const std::vector<TLumiBlockRange> 
TLumiBlockRange::GetPartOnlyIn(const TLumiBlockRange& other) const
{
  return other.GetPartNotIn(*this);
}


const std::vector<TLumiBlockRange> 
TLumiBlockRange::GetPartNotIn(const TLumiBlockRange& other) const
{
  std::vector<TLumiBlockRange> notinother;

  Int_t begin(-1), end(-1);

  if ( other.Begin()>this->Begin() ) {
    begin = this->Begin() ;
    end = ( other.Begin()>this->End() ? this->End() : other.Begin()-1 ) ; 
    if (end>=begin) notinother.push_back( TLumiBlockRange(begin,end) ) ;
  }

  if ( other.End()<this->End() ) {
    begin = ( other.End()>=this->Begin() ? other.End()+1 : this->Begin() ) ;
    end = this->End();
    if (end>=begin) notinother.push_back( TLumiBlockRange(begin,end) ) ;
  }

  return notinother;
}


Bool_t 
TLumiBlockRange::Contains( const Int_t& number )  const
{
  if (number>=m_begin && number<=m_end) return kTRUE;
  else return kFALSE;
}


void
TLumiBlockRange::Summary() const 
{
  TMsgLogger m_logger("TLumiBlockRange");
  m_logger << kINFO << "Begin: " << m_begin << " End: " << m_end << GEndl;  
}

}
}
