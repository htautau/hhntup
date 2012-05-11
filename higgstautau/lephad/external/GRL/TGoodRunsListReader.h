/**********************************************************************************
 * Class  : TGoodRun                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Max Baak <mbaak@cern.ch> - CERN, Switzerland                              *
 **********************************************************************************/

#ifndef _TGoodRunsListReader
#define _TGoodRunsListReader

#include <vector>

#include "TROOT.h"
#include "TList.h"
#include "TObject.h"
#include "TString.h"

#include "TMsgLogger.h"
#include "TGRLCollection.h"

#include <vector>

class TXMLNode;

namespace AnalysisFramework
{
namespace External
{
   class TStringList;
   class TLumiBlockRange;
   class TGoodRun;
   class TGoodRunsList;

   class TGoodRunsListReader : public TObject {

   public:

      TGoodRunsListReader( Bool_t checkGRLInfo=kFALSE );
      TGoodRunsListReader( const TString& dataCardName, Bool_t checkGRLInfo=kFALSE );
      ~TGoodRunsListReader();

      Bool_t Interpret();

      // accessor
      inline const TString& GetXMLString()   const { return m_xmlstring; }
      inline const TString& GetXMLFilename() const { return m_dataCardName; }

      inline void AddXMLFile( const TString& xmlfile )          { m_dataCardList.push_back(xmlfile); }
      inline void AddXMLString( const TString& xmlstring )      { m_xmlstringList.push_back(xmlstring); }
      inline void SetXMLString( const TString& xmlstring )      { Reset(); m_xmlstringList.push_back(xmlstring); }
      inline void SetXMLFile( const TString& xmlfile )          { Reset(); m_dataCardList.push_back(xmlfile); }
      inline void SetCheckGRLInfo( Bool_t check=kTRUE )         { m_grlvec.SetCheckGRLInfo( check ); }

      const TGoodRunsList  GetMergedGoodRunsList( const BoolOperation& operation = OR ) const ;
      const TGoodRunsList  GetGoodRunsList( unsigned int idx ) const ;
      inline const TGRLCollection GetGRLCollection() const      { return m_grlvec; }
      const TGRLCollection GetMergedGRLCollection( const BoolOperation& operation = OR ) const ;

      void Reset();

   private:

      void ReadNamedLumiRange ( TXMLNode* );
      void ReadLumiBlockCollection ( TXMLNode* );
      void ReadAttribs ( TXMLNode* );

      TGoodRun GetLumiBlockCollection( TXMLNode* dataNode ) ;

      TString m_xmlstring;
      TString m_dataCardName; // current xmlfile processed
      std::vector<TString> m_dataCardList;
      std::vector<TString> m_xmlstringList;
      TMsgLogger m_logger;
      TGRLCollection m_grlvec;

      ClassDef(TGoodRunsListReader,0)
   };
}
}
#endif
