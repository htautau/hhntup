
/**********************************************************************************
 * Class  : TGoodRunsList                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Max Baak <mbaak@cern.ch> - CERN, Switzerland                              *
 **********************************************************************************/

#ifndef __MyTGoodRunsList__
#define __MyTGoodRunsList__

#include "TGoodRun.h"
#include "TNamed.h"
#include "TString.h"
#include <map>
#include <vector>
#include <string>

namespace AnalysisFramework
{
namespace External
{
   class TGoodRunsList : public std::map< Int_t, TGoodRun >, public TNamed {

   public:

      TGoodRunsList();
      TGoodRunsList(const char* name);
      virtual ~TGoodRunsList();

      TGoodRunsList(const TGoodRunsList& other) ;
      TGoodRunsList& operator=(const TGoodRunsList& other) ;

      void AddGRL(const TGoodRunsList& other);
      const TGoodRunsList GetOverlapWith(const TGoodRunsList& other) const ;
      const TGoodRunsList GetSumWith(const TGoodRunsList& other) const ;
      const TGoodRunsList GetPartOnlyIn(const TGoodRunsList& other) const ;
      const TGoodRunsList GetPartNotIn(const TGoodRunsList& other) const ;

      Bool_t HasTriggerInfo() const;
      Bool_t HasRun( const Int_t& runnr )  const;
      Bool_t HasRunLumiBlock( const Int_t& runnr, const Int_t& lumiblocknr ) const ;
      Bool_t HasSameGRLInfo( const TGoodRunsList& other ) const;
      Bool_t HasOverlapWith( const TGoodRunsList& other, bool verb=false ) const;

      void AddRunLumiBlock( const Int_t& runnr, const Int_t& lumiblocknr );
      inline void SetVersion(const TString& version) { m_version = version; }
      inline void AddMetaData(const TString& key, const TString& value) { m_metadata[key] = value; }
      inline void SetMetaData(const std::map<TString,TString>& metadata) { m_metadata = metadata; }
      inline void SetCheckGRLInfo(Bool_t check=kTRUE) { m_checkGRLInfo=check; }

      inline const Bool_t& GetCheckGRLInfo() const { return m_checkGRLInfo; }
      inline const TString& GetVersion() const { return m_version; }
      inline const std::map<TString,TString>& GetMetaData() const { return m_metadata; }
      inline unsigned int GetMetaDataSize() const { return m_metadata.size(); }

      void Summary(Bool_t verbose = kFALSE) const;
      Bool_t IsEmpty() const;

      const std::vector<int> GetRunlist() const;
      const std::vector<AnalysisFramework::External::TGoodRun> GetGoodRuns() const;
      const std::vector<std::string> GetTriggerList() const;
      const std::vector<std::string> GetStreamList() const;

      const TString GetSuggestedName() const;
      void Compress();

   private:

      TString m_version;
      std::map<TString,TString> m_metadata;
      Bool_t m_checkGRLInfo;

      mutable Bool_t m_hasRun;
      mutable Bool_t m_hasLB;
      mutable Int_t m_prevRun;
      mutable Int_t m_prevLB;

      ClassDef(TGoodRunsList,1)
   };
}
}

#endif

