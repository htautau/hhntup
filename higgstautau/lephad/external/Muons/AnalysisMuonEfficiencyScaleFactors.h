//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 31.10.2011, MCP working group
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef MyAnalysisMuonEfficiencyScaleFactorH
#define MyAnalysisMuonEfficiencyScaleFactorH

//////////////////////////////////
// CLASS AnalysisMuonEfficiencyScaleFactors //
//////////////////////////////////

/// \class AnalysisMuonEfficiencyScaleFactors
///
/// This class provides efficiency scale factors and their uncertainties
/// for physics analyses using muons. The muon type is specified in the
/// constructor as well as the integrated luminosity per period used in the
/// analysis.
///
/// For release 17 data analysis with MC11a and MC11b.
///
/// \date 31.10.2011
///
/// \author Oliver.Kortner@CERN.CH
/// \author Marco.Vanadia@CERN>CH

//////////////////
// HEADER FILES //
//////////////////

// ROOT //
#include "TLorentzVector.h"

// STL //
#include <vector>
#include <string>

// Base class //
#include "MuonEfficiencyScaleFactor.h"
#include "EtaPhiBinning.h"

namespace AnalysisFramework
{
namespace External
{

class AnalysisMuonEfficiencyScaleFactors : public MuonEfficiencyScaleFactor {

public:
    //! Constructor
    AnalysisMuonEfficiencyScaleFactors(void);
    ///< Default constructor.
    ///< Please do not use the default constructor!

    AnalysisMuonEfficiencyScaleFactors(const std::string & alg,
                                       const std::vector<double> int_lum,
                                       const std::string & unit,
                                       std::string dir="");
    ///< Constructor.
    ///< \param alg Name of the reconstruction algorithm. Possible values are:
    ///<            STACO_CB (STACO combined muons),
    ///<            STACO_CB_plus_ST (STACO combined plus segment tagged muons),
    ///<            Muid_CB (Muid combined muons),
    ///<            Muid_CB_plus_ST (Muid combined plus segment tagged muons)
    ///< \param int_lum Vector containing the integrated luminosity per period
    ///<                used in the physics analysis.
    ///<                int_lum[0] = luminosity for period B,
    ///<                int_lum[1] = luminosity for period D,
    ///<                int_lum[2] = luminosity for period E,
    ///<                int_lum[3] = luminosity for period F,
    ///<                int_lum[4] = luminosity for period G,
    ///<                int_lum[5] = luminosity for period H,
    ///<                int_lum[6] = luminosity for period I,
    ///<                int_lum[7] = luminosity for period J,
    ///<                int_lum[8] = luminosity for period K,
    ///<                int_lum[9] = luminosity for period L,
    ///<                int_lum[10] = luminosity for period M.
    ///< \param unit MeV if muon momenta are provided in MeV,
    ///<             GeV if muon momenta are provided in GeV.
    ///< \param dir Directory containing the scale factors.

    virtual ~AnalysisMuonEfficiencyScaleFactors() {}

    // Methods //
    double scaleFactor(const TLorentzVector & tlv) const;
    ///< Get the efficiency scale factor for the given
    ///< fourmomentum. Scale factors for period after B are provided.
    double scaleFactorUncertainty(const TLorentzVector & tlv) const;
    ///< Get the uncertainty of the efficiency scale
    ///< factor for the given fourmomentum.
    double scaleFactorSystematicUncertainty(const TLorentzVector & tlv) const;
    ///< Get the systematic uncertainty of the scale factor. The momentum
    ///< is assumed to be given in MeV or GeV (see above).
    void PrintValues() const; //debug methos to èrint the values of the vectors.

private:
// unit //
    double m_unit; // m_unit = 1 for GeV, 0.001 for MeV

// pt-eta-phi map of the efficiency scale factors //
    std::vector<int> m_region;
    std::vector<double> m_eta_min, m_eta_max; // bin boundaries in eta
    std::vector<double> m_phi_min, m_phi_max; // bin boundaries in phi
    std::vector<double> m_pt_min, m_pt_max; // bin boundaries in pt
    std::vector<double> m_sf; // scale factors
    std::vector<double> m_sf_stat_err; // statistical errors of the scale
                                       // factors
    std::vector<double> m_sf_syst_err; // systematic errors of the scale 
                                       // factors

    void init(const std::string & alg,
                                       const std::vector<double> int_lum,
                                       const std::string & unit,
                                       std::string dir="");

    void read_file(const std::string & name,
                   std::vector<int> & region,
                   std::vector<double> & eta_min,
                   std::vector<double> & eta_max,
                   std::vector<double> & phi_min,
                   std::vector<double> & phi_max,
                   std::vector<double> & pt_min,
                   std::vector<double> & pt_max,
                   std::vector<double> & sf,
                   std::vector<double> & sf_stat_err,
                   std::vector<double> & st_syst_err) const;

    int get_pt_eta_phi_bin_index(const TLorentzVector & tlv) const;
    EtaPhiBinning m_EPbin;
 };

}
}
#endif
