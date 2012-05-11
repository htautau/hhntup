
/**********************************************************************************
 * Class  : TGoodRunsReader                                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Max Baak <mbaak@cern.ch> - CERN, Switzerland                              *
 **********************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <string>

#include "TFormula.h"
#include "Riostream.h"
#include "TObjString.h"
#include "TDOMParser.h"
#include "TXMLNode.h"
#include "TXMLDocument.h"
#include "TXMLAttr.h"

#include "TGoodRunsListReader.h"
#include "TGoodRunsList.h"
#include "TGoodRun.h"
#include "TLumiBlockRange.h"
#include "StrUtil.h"

//ClassImp(TGoodRunsListReader)

using namespace std;

namespace AnalysisFramework
{
namespace External
{

TGoodRunsListReader::TGoodRunsListReader( Bool_t checkGRLInfo )
 : TObject()
 , m_logger( "TGoodRunsListReader" )
{
  m_grlvec.SetCheckGRLInfo( checkGRLInfo );
}


TGoodRunsListReader::TGoodRunsListReader( const TString& dataCardName, Bool_t checkGRLInfo )
 : TObject()
 , m_logger( "TGoodRunsListReader" )
{
  m_dataCardList.push_back(dataCardName);
  m_grlvec.SetCheckGRLInfo( checkGRLInfo );
}


TGoodRunsListReader::~TGoodRunsListReader()
{
  this->Reset();
}


void
TGoodRunsListReader::Reset()
{
  m_grlvec.Reset();
  m_dataCardList.clear();
  m_xmlstringList.clear();
  m_xmlstring=""; //Clear() only works in root5.24
  m_dataCardName="";
}


Bool_t
TGoodRunsListReader::Interpret()
{
  Bool_t xmlInterpret(kTRUE);

  if (m_dataCardList.empty() && m_xmlstringList.empty()) {
    m_logger << kWARNING << "No xml data-card or string set. Return false." << GEndl;
    return kFALSE;
  }

  Int_t parseCode(1);
  TDOMParser* xmlparser = new TDOMParser();

  //////////////////////////////////////////////////////////////////////

  // --------------- xml file read
  for (unsigned int j=0; j<m_dataCardList.size() && xmlInterpret; ++j) {
    m_dataCardName = m_dataCardList[j];

    if (!m_dataCardName.IsNull()) {
      m_logger << kDEBUG << "Read xml data-card: \"" << m_dataCardName << "\"" << GEndl;
      xmlparser->SetValidate(kFALSE); // MB 14/4/'10 : don't validate structure of dtd file in case runquery down
      parseCode = xmlparser->ParseFile( m_dataCardName );
    } else {
      m_logger << kWARNING << "No xml data-card set. Skip." << GEndl;
      continue;
    }

    m_logger << kDEBUG << "XML parser returned code: " << parseCode << GEndl;
    if (parseCode != 0) {
      m_logger << kERROR << "loading of xml document failed" << GEndl;
      xmlInterpret = kFALSE;
    } else {
      // --------------- parse JobConfiguration
      TXMLDocument* xmldoc = xmlparser->GetXMLDocument();

      TXMLNode* jobConfig_node = xmldoc->GetRootNode();
      TXMLNode* jobConfig_elem = jobConfig_node->GetChildren();

      while (jobConfig_elem != 0) {
        if (jobConfig_elem->GetNodeName() == TString("NamedLumiRange")) {
          this->ReadNamedLumiRange ( jobConfig_elem );
        }
        // crawl on...
        jobConfig_elem = jobConfig_elem->GetNextNode();
      }
    }
  }
  m_dataCardList.clear(); // Can now add fresh xml files

  //////////////////////////////////////////////////////////////////////

  // --------------- xml string read
  for (unsigned int j=0; j<m_xmlstringList.size() && xmlInterpret; ++j) {
    m_xmlstring = m_xmlstringList[j];

    if (!m_xmlstring.IsNull()) {
      m_logger << kDEBUG << "Read xml string." << GEndl;
      xmlparser->SetValidate(kFALSE); // MB 14/4/'10 : don't validate structure of dtd file in case runquery down
      parseCode = xmlparser->ParseBuffer( m_xmlstring.Data(), m_xmlstring.Length() );
    } else {
      m_logger << kWARNING << "No xml string set. Skip." << GEndl;
      continue;
    }

    m_logger << kDEBUG << "XML parser returned code: " << parseCode << GEndl;
    if (parseCode != 0) {
      m_logger << kERROR << "loading of xml document failed" << GEndl;
      xmlInterpret = kFALSE;
    } else {
      // --------------- parse JobConfiguration
      TXMLDocument* xmldoc = xmlparser->GetXMLDocument();

      TXMLNode* jobConfig_node = xmldoc->GetRootNode();
      TXMLNode* jobConfig_elem = jobConfig_node->GetChildren();

      while (jobConfig_elem != 0) {
        if (jobConfig_elem->GetNodeName() == TString("NamedLumiRange")) {
          this->ReadNamedLumiRange ( jobConfig_elem );
        }
        // crawl on...
        jobConfig_elem = jobConfig_elem->GetNextNode();
      }
    }
  }
  m_xmlstringList.clear(); // Can now add fresh xml strings

  //////////////////////////////////////////////////////////////////////

  delete xmlparser;

  return xmlInterpret;
}


void
TGoodRunsListReader::ReadAttribs( TXMLNode* node )
{
   if (!node->HasAttributes()) return;

   TListIter attribIt( node->GetAttributes() );
   TXMLAttr* curAttr( 0 );
   while ((curAttr = (TXMLAttr*)attribIt()) != 0) {
     m_logger << kDEBUG << node->GetNodeName() << ": " << curAttr->GetName()
              << "  =  \"" << curAttr->GetValue() << "\"" << GEndl;
     // add variable
   }
}


void
TGoodRunsListReader::ReadLumiBlockCollection( TXMLNode* dataNode )
{
  // retrieve childrens of node
  TXMLNode* node = dataNode->GetChildren();

  for (; node != 0; node = node->GetNextNode()) {

    // read run
    if (TString("Run")==node->GetNodeName()) {
      m_logger << kDEBUG << "subchild node value: " << node->GetText() << GEndl;
    }

    // read lumi block
    if (TString("LBRange")==node->GetNodeName() && node->HasAttributes()) {
      TXMLAttr* curAttr( 0 );
      TListIter attribIt(node->GetAttributes());
      while ((curAttr = (TXMLAttr*)attribIt()) != 0) {
        if (TString("Start")==curAttr->GetName()) {
          m_logger << kDEBUG << node->GetNodeName() << ": " << curAttr->GetName()
                   << "  =  \"" << curAttr->GetValue() << "\"" << GEndl;
        } else if (TString("End")==curAttr->GetName()) {
          m_logger << kDEBUG << node->GetNodeName() << ": " << curAttr->GetName()
                   << "  =  \"" << curAttr->GetValue() << "\"" << GEndl;
        }
      }
    }
  }
}


TGoodRun
TGoodRunsListReader::GetLumiBlockCollection( TXMLNode* dataNode )
{
  TGoodRun goodrun;

  // sanity check
  if (!dataNode->HasChildren()) {
     m_logger << kWARNING << "<Data> ... </Data> Part does not contain any parameters" << GEndl;
     return goodrun;
  }

  // retrieve childrens of node
  TXMLNode* node = dataNode->GetChildren();

  for (; node != 0; node = node->GetNextNode()) {
    // read run
    if (TString("Run")==node->GetNodeName()) {
      m_logger << kDEBUG << "subchild node value: " << node->GetText() << GEndl;
      goodrun.SetRunNumber(atoi(node->GetText()));
    }
    // read lumi block
    if (TString("LBRange")==node->GetNodeName() && node->HasAttributes()) {
      TLumiBlockRange lbr;
      TXMLAttr* curAttr( 0 );
      TListIter attribIt(node->GetAttributes());
      while ((curAttr = (TXMLAttr*)attribIt()) != 0) {
        if (TString("Start")==curAttr->GetName()) {
          m_logger << kDEBUG << node->GetNodeName() << ": " << curAttr->GetName()
                   << "  =  \"" << curAttr->GetValue() << "\"" << GEndl;
          lbr.SetBegin(atoi(curAttr->GetValue()));
          lbr.SetEnd(2147483647); // set in case lb turns out to be open-ended
        } else if (TString("End")==curAttr->GetName()) {
          m_logger << kDEBUG << node->GetNodeName() << ": " << curAttr->GetName()
                   << "  =  \"" << curAttr->GetValue() << "\"" << GEndl;
          lbr.SetEnd(atoi(curAttr->GetValue()));
        }
      }
      if (!lbr.IsEmpty()) goodrun.push_back(lbr);
    }
  }

  goodrun.Sort();
  return goodrun;
}


void
TGoodRunsListReader::ReadNamedLumiRange( TXMLNode* dataNode )
{
  // sanity check
  if (!dataNode->HasChildren()) {
     m_logger << kWARNING << "<Data> ... </Data> Part does not contain any parameters" << GEndl;
     return;
  }

  // retrieve childrens of node
  TXMLNode* node = dataNode->GetChildren();
  TGoodRunsList grl(m_dataCardName.Data());
  std::string nameStr, valueStr;

  for (; node != 0; node = node->GetNextNode()) {
    /// set name
    if (TString("Name") == node->GetNodeName()) {
      if (node->GetText()!=0) {
        m_logger << kDEBUG << "child node value: " << node->GetText() << GEndl;
        nameStr=node->GetText();
      } else { nameStr=""; }
      GRLStrUtil::trim(nameStr);
      grl.SetName(nameStr.c_str());
    }
    /// set version
    else if (TString("Version") == node->GetNodeName()) {
      if (node->GetText()!=0) {
        m_logger << kDEBUG << "child node value: " << node->GetText() << GEndl;
        valueStr=node->GetText();
      } else { valueStr=""; }
      GRLStrUtil::trim(valueStr);
      grl.SetVersion(valueStr);
    }
    /// set metadata
    else if (TString("Metadata") == node->GetNodeName()) {
      if (node->GetText()!=0) m_logger << kDEBUG << node->GetNodeName() << " value: " << node->GetText() << GEndl;
      this->ReadAttribs(node);

      if (node->HasAttributes()) {
        TListIter attribIt( node->GetAttributes() );
        TXMLAttr* curAttr( 0 );
        while ((curAttr = (TXMLAttr*)attribIt()) != 0) {
          if (curAttr->GetValue()!=0) { nameStr=curAttr->GetValue(); } else { nameStr=""; }
          if (node->GetText()!=0) { valueStr=node->GetText(); } else { valueStr=""; }
          GRLStrUtil::trim(nameStr); GRLStrUtil::trim(valueStr);
          if (!nameStr.empty() && !valueStr.empty()) grl.AddMetaData(nameStr,valueStr);
        }
      }
    }
    /// set run and lb range
    else if (TString("LumiBlockCollection") == node->GetNodeName()) {
      TGoodRun goodrun = this->GetLumiBlockCollection(node);
      if (!goodrun.IsEmpty()) grl[goodrun.GetRunNumber()] = goodrun ;
    }
  }

  if (!grl.IsEmpty()) m_grlvec.push_back(grl);
}


const TGoodRunsList
TGoodRunsListReader::GetMergedGoodRunsList( const BoolOperation& operation ) const
{
  return m_grlvec.GetMergedGoodRunsList(operation);
}


const TGoodRunsList
TGoodRunsListReader::GetGoodRunsList( unsigned int idx ) const
{
  return m_grlvec.GetGoodRunsList(idx);
}


const TGRLCollection
TGoodRunsListReader::GetMergedGRLCollection( const BoolOperation& operation ) const
{
  return m_grlvec.GetMergedGRLCollection(operation);
}

}
}
