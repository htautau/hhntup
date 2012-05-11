
#include "TGRLCollection.h"
#include "TMsgLogger.h"

#include <algorithm>
#include <iostream>

//ClassImp(TGRLCollection)
namespace AnalysisFramework
{
namespace External
{

TGRLCollection::TGRLCollection( Bool_t checkGRLInfo )
 : std::vector<TGoodRunsList>()
 , TObject()
 , m_checkGRLInfo( checkGRLInfo )
{
}


TGRLCollection::~TGRLCollection()
{
  this->Reset();
}


TGRLCollection::TGRLCollection(const TGRLCollection& other)
 : std::vector<TGoodRunsList>(other)
 , TObject(other)
 , m_checkGRLInfo(other.m_checkGRLInfo)
{
}


TGRLCollection&
TGRLCollection::operator=(const TGRLCollection& other)
{
  if (&other==this) {
    return *this ;
  }
  std::vector<TGoodRunsList>::operator=(other);
  m_checkGRLInfo=other.m_checkGRLInfo;

  return *this ;
}


void
TGRLCollection::Reset()
{
  this->clear();
  m_checkGRLInfo=kFALSE;
}


void
TGRLCollection::SetVersion( const TString& version )
{
  std::vector<TGoodRunsList>::iterator itr = this->begin();
  std::vector<TGoodRunsList>::iterator end_ = this->end();
  for (; itr!=end_; ++itr) itr->SetVersion(version);
}


void
TGRLCollection::SetMetaData( const std::map<TString,TString>& metadata )
{
  std::vector<TGoodRunsList>::iterator itr = this->begin();
  std::vector<TGoodRunsList>::iterator end_ = this->end();
  for (; itr!=end_; ++itr) itr->SetMetaData(metadata);
}

void
TGRLCollection::Summary (Bool_t verbose /*= kFALSE*/) const
{
  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  std::vector<TGoodRunsList>::const_iterator end_ = this->end();

  for (; itr!=end_; ++itr)
    itr->Summary(verbose) ;
}


Bool_t
TGRLCollection::HasRun( const Int_t& runnr )  const
{
  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  std::vector<TGoodRunsList>::const_iterator end_ = this->end();

  Bool_t pass(kFALSE);
  for (; itr!=end_ && !pass; ++itr)
    pass = itr->HasRun(runnr);

  return pass;
}


Bool_t
TGRLCollection::HasRunLumiBlock( const Int_t& runnr, const Int_t& lumiblocknr ) const
{
  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  std::vector<TGoodRunsList>::const_iterator end_ = this->end();

  Bool_t pass(kFALSE);
  for (; itr!=end_ && !pass; ++itr)
    pass = itr->HasRunLumiBlock(runnr,lumiblocknr);

  return pass;
}


Bool_t
TGRLCollection::IsEmpty() const
{
  if (this->empty()) return kTRUE;

  Bool_t isEmpty(kTRUE);
  std::vector< TGoodRunsList >::const_iterator litr = this->begin();
  for (; litr!=this->end() && isEmpty; ++litr)
    isEmpty = isEmpty && litr->IsEmpty();

  return isEmpty;
}


const TGoodRunsList
TGRLCollection::GetMergedGoodRunsList( const BoolOperation& operation ) const
{
  // nothing interpreted. Return empty grl.
  if (this->empty()) return TGoodRunsList();

  // set first goodrunslist
  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  TGoodRunsList grl(*itr);
  if (this->size()==1) {
    grl.Compress();
    return grl;
  }

  TMsgLogger mylogger( "TGRLCollection" );
  mylogger << kINFO << "Now merging GRLs." << GEndl;

  // check version and metadata when merging goodrunslists?
  grl.SetCheckGRLInfo(m_checkGRLInfo);

  if (!m_checkGRLInfo)
    mylogger << kINFO << "Metadata and other info not required to be identical between GRLs." << GEndl;

  // start AND-ing or OR-ring with following goodrunslists
  for (++itr; itr!=this->end(); ++itr) {
    switch (operation) {
      case OR :
        if ( grl.HasOverlapWith(*itr,false/*verbose*/) ) { // MB 22-june: LB splitting across files, turn off warning.
          //mylogger << kWARNING << "Merging GRLs with overlapping lumi-blocks! Overlapping LBs rejected." << GEndl;
          //mylogger << kWARNING << "IMPORTANT : Check your analysis for possible duplicate events!" << GEndl;
        }
        grl.AddGRL( *itr );
      break;
      case AND :
        grl = grl.GetOverlapWith(*itr);
      break;
    }
  }
  grl.Compress(); // cleanup, safe space
  return grl;
}


const TGoodRunsList
TGRLCollection::GetGoodRunsList( unsigned int idx ) const
{
  // invalid idx. Return empty grl.
  if (idx>=this->size()) return TGoodRunsList();

  return (*this)[idx];
}


const TGRLCollection
TGRLCollection::GetMergedGRLCollection( const BoolOperation& operation ) const
{
  if (this->empty() /*|| this->size()==1*/) return *this;  // nothing to merge, return this

  TMsgLogger mylogger( "TGRLCollection" );
  mylogger << kINFO << "Now merging GRLs where possible. Metadata required to be identical." << GEndl;

  TGRLCollection mergevec;

  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  std::vector<TGoodRunsList>::const_iterator end_ = this->end();
  std::vector<TGoodRunsList>::iterator mitr;

  for (; itr!=end_; ++itr) {
    bool matchFound(false);
    for (mitr=mergevec.begin(); mitr!=mergevec.end() && !matchFound ; ++mitr) {
      if (mitr->HasSameGRLInfo(*itr)) {
        matchFound = true;
        switch (operation) {
          case OR :
            if ( mitr->HasOverlapWith(*itr,false/*verbose*/) ) { // // MB 22-june: LB splitting across files, turn off warning.
              //mylogger << kWARNING << "Merging GRLs with overlapping lumi-blocks! Overlapping LBs rejected." << GEndl;
              //mylogger << kWARNING << "IMPORTANT : Check your analysis for possible duplicate events!" << GEndl;
            }
            mitr->AddGRL( *itr );
          break;
          case AND :
            *mitr = mitr->GetOverlapWith( *itr );
          break;
        }
        mitr->Compress(); // safe space
      }
    }
    if (!matchFound) {
      mergevec.push_back(*itr);
      mergevec.rbegin()->Compress(); // safe space
    }
  }

  return mergevec;
}


std::vector<TGoodRunsList>::iterator
TGRLCollection::find( const TString& name )
{
  Bool_t found(false);

  std::vector<TGoodRunsList>::iterator itr = this->begin();
  std::vector<TGoodRunsList>::iterator end_ = this->end();

  for (; itr!=end_; ++itr) {
    found = ( name==TString(itr->GetName()) ) ;
    if (found) break;
  }

  return itr;
}


std::vector<TGoodRunsList>::const_iterator
TGRLCollection::find( const TString& name ) const
{
  Bool_t found(false);

  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  std::vector<TGoodRunsList>::const_iterator end_ = this->end();

  for (; itr!=end_; ++itr) {
    found = ( name==TString(itr->GetName()) ) ;
    if (found) break;
  }

  return itr;
}


Bool_t
TGRLCollection::HasGoodRunsList( const TString& name ) const
{
  return (this->find(name)!=this->end());
}


const TGRLCollection
TGRLCollection::GetOverlapWith( const TGoodRunsList& other ) const
{
  TGRLCollection overlapvec;

  std::vector<TGoodRunsList>::const_iterator itr = this->begin();
  for (; itr!=this->end(); ++itr) {
    TGoodRunsList overlapgrl = itr->GetOverlapWith(other);
    overlapgrl.SetName(itr->GetName());
    overlapgrl.SetVersion(itr->GetVersion());
    overlapgrl.SetMetaData(itr->GetMetaData());
    overlapgrl.Compress();
    overlapvec.push_back(overlapgrl); // also push_back if empty!
  }

  return overlapvec;
}

}
}

