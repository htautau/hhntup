
/**********************************************************************************
 * Class  : TGRLCollection                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Max Baak <mbaak@cern.ch> - CERN, Switzerland                              *
 **********************************************************************************/

#ifndef __MYTGRLCollection__
#define __MYTGRLCollection__

#include "TGoodRunsList.h"
#include "TObject.h"
#include <vector>

namespace AnalysisFramework
{
namespace External
{
   enum BoolOperation { OR, AND };

   class TGRLCollection : public std::vector< TGoodRunsList >, public TObject {

   public:

      TGRLCollection( Bool_t checkGRLInfo=kFALSE ) ;
      virtual ~TGRLCollection();

      TGRLCollection(const TGRLCollection& other) ;
      TGRLCollection& operator=(const TGRLCollection& other) ;

      void SetVersion( const TString& version );
      void SetMetaData( const std::map<TString,TString>& metadata );
      inline void SetCheckGRLInfo( Bool_t check=kTRUE ) { m_checkGRLInfo=check; }

      Bool_t HasRun( const Int_t& runnr )  const;
      Bool_t HasRunLumiBlock( const Int_t& runnr, const Int_t& lumiblocknr ) const ;
      Bool_t IsEmpty() const;
      Bool_t HasGoodRunsList( const TString& name ) const;

      const TGoodRunsList GetMergedGoodRunsList( const BoolOperation& operation = OR ) const ;
      const TGoodRunsList GetGoodRunsList( unsigned int idx ) const ;
      std::vector<TGoodRunsList>::iterator find( const TString& name );
      std::vector<TGoodRunsList>::const_iterator find( const TString& name ) const ;
      const TGRLCollection GetMergedGRLCollection( const BoolOperation& operation = OR) const ;
      const TGRLCollection GetOverlapWith(const TGoodRunsList& other) const ;

      void Reset() ;
      void Summary(Bool_t verbose=kFALSE) const ;

   private:

      Bool_t m_checkGRLInfo;

      ClassDef(TGRLCollection,1)
   };
}
}

#endif

