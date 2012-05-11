
#include "TGoodRun.h"
#include "TMsgLogger.h"

#include <algorithm>

// ClassImp(TGoodRun)

namespace AnalysisFramework
{
namespace External
{

TGoodRun::TGoodRun( )
  : std::vector<TLumiBlockRange>()
  , TObject()
  , m_runnr(-1)
{
}


TGoodRun::TGoodRun( const Int_t& runnr )
  : std::vector<TLumiBlockRange>()
  , TObject()
  , m_runnr(runnr)
{
}


TGoodRun::~TGoodRun()
{
}


TGoodRun::TGoodRun(const TGoodRun& other)
  : std::vector<TLumiBlockRange>(other)
  , TObject(other)
  , m_runnr(other.m_runnr)
{
}


TGoodRun&
TGoodRun::operator=(const TGoodRun& other)
{
  if (&other==this) {
    return *this ;
  } 
  std::vector<TLumiBlockRange>::operator=(other);
  m_runnr = other.m_runnr;

  return *this ;
}


const TGoodRun
TGoodRun::GetOverlapWith(const TGoodRun& other) const
{
  Int_t runnr = ( this->GetRunNumber()==other.GetRunNumber() ? this->GetRunNumber() : -1 ) ;

  TGoodRun overlaprun(runnr);
  if (runnr==-1) return overlaprun; // no overlap -> empty overlap

  std::vector<TLumiBlockRange>::const_iterator itr1 = this->begin();
  std::vector<TLumiBlockRange>::const_iterator itr2 = other.begin();

  for (; itr1!=this->end(); ++itr1) {
    for (itr2=other.begin(); itr2!=other.end(); ++itr2) {
      const TLumiBlockRange lbr = itr1->GetOverlapWith(*itr2) ;
      if (!lbr.IsEmpty()) overlaprun.push_back(lbr);
    }
  } 

  return overlaprun;
}


const TGoodRun
TGoodRun::GetSumWith(const TGoodRun& other) const
{
  Int_t runnr = ( this->GetRunNumber()==other.GetRunNumber() ? this->GetRunNumber() : -1 ) ;

  TGoodRun sumrun(runnr);
  if (runnr==-1) return sumrun; // only add runs with same runnumber

  TGoodRun overlap     = this->GetOverlapWith(other);
  TGoodRun onlyinother = this->GetPartOnlyIn(other);
  TGoodRun notinother  = this->GetPartNotIn(other);

  std::vector<TLumiBlockRange>::const_iterator itr;
  for (itr=overlap.begin(); itr!=overlap.end(); ++itr) sumrun.push_back(*itr);
  for (itr=onlyinother.begin(); itr!=onlyinother.end(); ++itr) sumrun.push_back(*itr);
  for (itr=notinother.begin(); itr!=notinother.end(); ++itr) sumrun.push_back(*itr);
  sumrun.Sort();     // sort lb ranges

  return sumrun;
}


const TGoodRun 
TGoodRun::GetPartOnlyIn(const TGoodRun& other) const
{
  return other.GetPartNotIn(*this);
}


const TGoodRun 
TGoodRun::GetPartNotIn(const TGoodRun& other) const
{
  TGoodRun notinother(*this);

  if (this->GetRunNumber()==other.GetRunNumber()) {

    TGoodRun tmprun(this->GetRunNumber());
    std::vector<TLumiBlockRange>::const_iterator itro,itrt,lbvecitr;

    for (itro=other.begin(); itro!=other.end(); ++itro) {
      tmprun.clear();
      for (itrt=notinother.begin(); itrt!=notinother.end(); ++itrt) {
        std::vector<TLumiBlockRange> lbvec = itrt->GetPartNotIn(*itro);
        for (lbvecitr=lbvec.begin(); lbvecitr!=lbvec.end(); ++lbvecitr) 
          tmprun.push_back(*lbvecitr);
      }
      notinother=tmprun;
    }
  } 

  return notinother;
}


Bool_t 
TGoodRun::HasLB( const Int_t& lumiblocknr )  const
{
  std::vector<TLumiBlockRange>::const_iterator itr = this->begin();
  std::vector<TLumiBlockRange>::const_iterator end_ = this->end();

  Bool_t pass(kFALSE);
  for (; itr!=end_; ++itr) {
    pass = pass || itr->Contains(lumiblocknr) ;
    if (pass) break;
  }

  return pass;
}


std::vector< TLumiBlockRange >::iterator
TGoodRun::Find( const Int_t& lumiblocknr )
{
  std::vector<TLumiBlockRange>::iterator itr = this->begin();

  Bool_t pass(kFALSE);
  for (; itr!=this->end(); ++itr) {
    pass = itr->Contains(lumiblocknr) ;
    if (pass) break;
  }

  return itr;
}


std::vector< TLumiBlockRange >::const_iterator
TGoodRun::Find( const Int_t& lumiblocknr ) const
{
  std::vector<TLumiBlockRange>::const_iterator itr = this->begin();

  Bool_t pass(kFALSE);
  for (; itr!=this->end(); ++itr) {
    pass = itr->Contains(lumiblocknr) ;
    if (pass) break;
  }

  return itr;
}


void
TGoodRun::Summary () const
{
  TMsgLogger m_logger("TGoodRun");
  m_logger << kINFO << "Run: " << m_runnr << GEndl;

  std::vector<TLumiBlockRange>::const_iterator itr = this->begin();
  std::vector<TLumiBlockRange>::const_iterator end_ = this->end();

  for (; itr!=end_; ++itr)
    itr->Summary() ;
}


void
TGoodRun::Sort ()
{
  std::sort(this->begin(),this->end(),SorterL2H());
}


void 
TGoodRun::Compress()
{
  // this function throws away redundant lb ranges.

  // need to be sure goodrun is sorted in lumiblock ranges
  this->Sort();

  std::vector<TLumiBlockRange>::iterator itr, cbegin, cend;
  int fbegin(-1), fend(-1), tbegin(-1), tend(-1), lbegin(-1), lend(-1);

  for (itr=this->begin(); itr!=this->end();) {
    cbegin = itr;
    cend = itr;

    // first lumiblock of merge
    fbegin = cbegin->Begin();
    fend = cbegin->End();

    // search all connecting lumi ranges
    bool go_on(true);
    for (++itr; itr!=this->end() && go_on; ++itr) {
      tbegin = cend->Begin();
      tend   = cend->End();
      lbegin = itr->Begin();
      lend   = itr->End();
      if ( lbegin>=tbegin && lend>=tend ) { // should always be true
        if (tend+1<lbegin) { go_on=false; }
        else { cend = itr; }
      } 
    }

    // reset iterator
    itr=cend; ++itr;

    if (cbegin!=cend) {
      // first lumiblock of merge
      fbegin = cbegin->Begin();
      fend = cbegin->End();
      // last lumiblock of merge
      lbegin = cend->Begin();
      lend = cend->End();
      // remove [cbegin,itr)
      itr = this->erase(cbegin,itr); 
      // insert merge lbrange before itr
      this->insert(itr,TLumiBlockRange(fbegin,lend));
    }
  }
}


Bool_t 
TGoodRun::IsEmpty() const
{
  if (this->GetRunNumber()<0) return kTRUE;
  if (this->empty()) return kTRUE;

  Bool_t isEmpty(kTRUE);
  std::vector< TLumiBlockRange >::const_iterator litr = this->begin();
  for (; litr!=this->end() && isEmpty; ++litr)
    isEmpty = isEmpty && litr->IsEmpty();

  return isEmpty;
}


void 
TGoodRun::AddLB( const Int_t& lumiblocknr )
{
  if ( lumiblocknr<0 ) return;
  if ( this->HasLB(lumiblocknr) ) return;

  // add at end
  std::vector< TLumiBlockRange >::iterator eitr = this->Find(lumiblocknr-1);
  if (eitr!=this->end()) { 
    eitr->SetEnd(lumiblocknr);
    this->Compress();
    return;
  } 
  // add at beginning
  std::vector< TLumiBlockRange >::iterator bitr = this->Find(lumiblocknr+1);
  if (bitr!=this->end()) {
    bitr->SetBegin(lumiblocknr);
    this->Compress();
    return;
  } 
  // add as new range 
  this->push_back(TLumiBlockRange(lumiblocknr,lumiblocknr)); 
  this->Sort();
  return;
}

}
}
