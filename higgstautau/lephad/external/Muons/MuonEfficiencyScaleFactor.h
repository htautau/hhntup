//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 07.12.2010, MCP working group
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef MyMuonEfficiencyScaleFactorH
#define MyMuonEfficiencyScaleFactorH

/////////////////////////////////////
// CLASS MuonEfficiencyScaleFactor //
/////////////////////////////////////

/// \class MuonEfficiencyScaleFactor
///
/// This class is the base class for muon efficiency scale factors.
///
/// \date 27.04.2011
///
/// \author Oliver.Kortner@CERN.CH

//////////////////
// HEADER FILES //
//////////////////

// ROOT //
#include "TLorentzVector.h"

namespace AnalysisFramework
{
namespace External
{
  class MuonEfficiencyScaleFactor {
  public:
    MuonEfficiencyScaleFactor() {}
    virtual ~MuonEfficiencyScaleFactor() {}
    // Methods //
    virtual double scaleFactor(const TLorentzVector & tlv) const = 0;
    ///< Get the efficiency scale factor for the given
    ///< fourmomentum.
    virtual double scaleFactorUncertainty(const TLorentzVector & tlv) const = 0;
    ///< Get the uncertainty of the efficiency scale
    ///< factor for the given fourmomentum.
    virtual double scaleFactorSystematicUncertainty(
                                        const TLorentzVector & tlv) const = 0;
  };
}
}
#endif
