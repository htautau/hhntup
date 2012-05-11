//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// 31.10.2011, MCP working group
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//////////////////
// HEADER FILES //
//////////////////

#include <fstream>
#include <iostream>
#include "AnalysisMuonEfficiencyScaleFactors.h"

using namespace std;

//*****************************************************************************

/////////////////////////
// DEFAULT CONSTRUCTOR //
/////////////////////////

namespace AnalysisFramework
{
namespace External
{


AnalysisMuonEfficiencyScaleFactors::AnalysisMuonEfficiencyScaleFactors(void) {

  //OLD APPROACH. 


  //The values used here are the ones from the AllGood GRL
  /*
    vector<double> int_lum(11);
    int_lum[0] = 11.991;
    int_lum[1] = 166.467;
    int_lum[2] = 50.418;
    int_lum[3] = 136.759;
    int_lum[4] = 517.988;
    int_lum[5] = 264.752;
    int_lum[6] = 333.797;
    int_lum[7] = 233.152;
    int_lum[8] = 576.290;
    int_lum[9] = 1415.898;
    int_lum[10] = 1031.540;
  */    
    //this->init("STACO_CB", int_lum, "GeV", "");

    //initializing to standard values all the vectors

  //NOW USING INIT WITH STANDARD VALUES. STACO_CB values are used

  //using MeV by default
  m_unit = 0.001;

  //side A
    for(int i=1; i<11; i++)  
      {
        m_pt_min.push_back(3);
	m_pt_min.push_back(7);
	m_pt_min.push_back(10);
	m_pt_min.push_back(7000);

	m_pt_max.push_back(7);
	m_pt_max.push_back(10);
	m_pt_max.push_back(7000);
	m_pt_max.push_back(1e+09);

	//syst error
	m_sf_syst_err.push_back(0.02);
	m_sf_syst_err.push_back(0.01);
	m_sf_syst_err.push_back(0.002);
	m_sf_syst_err.push_back(0.002);

	for(int j=0; j<4; j++) 
	  {
	    m_region.push_back(i);
	    m_eta_min.push_back(100);
	    m_eta_max.push_back(100);
	    m_phi_min.push_back(100);
	    m_phi_max.push_back(100);
	  }
      }

    //side C
    for(int i=1; i<11; i++)  
      {
        m_pt_min.push_back(3);
	m_pt_min.push_back(7);
	m_pt_min.push_back(10);
	m_pt_min.push_back(7000);

	m_pt_max.push_back(7);
	m_pt_max.push_back(10);
	m_pt_max.push_back(7000);
	m_pt_max.push_back(1e+09);

	//syst error
	m_sf_syst_err.push_back(0.02);
	m_sf_syst_err.push_back(0.01);
	m_sf_syst_err.push_back(0.002);
	m_sf_syst_err.push_back(0.002);

	for(int j=0; j<4; j++) 
	  {
	    m_region.push_back(-i);
	    m_eta_min.push_back(100);
	    m_eta_max.push_back(100);
	    m_phi_min.push_back(100);
	    m_phi_max.push_back(100);
	  }
      }
	    //statical error
    for(int j=0; j<4; j++) 	    m_sf_stat_err.push_back(0.000779433);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000869424);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000893747);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00140777);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00135603);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000694567);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000702335);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000980779);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00120622);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000830213);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000824911);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000831151);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000930575);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00137499);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00141248);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000667718);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.000693062);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00104101);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00115552);
    for(int j=0; j<4; j++)	    m_sf_stat_err.push_back(0.00115552);
	    //sfs
    for(int j=0; j<4; j++)	    m_sf.push_back(0.991434);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.97241);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.984447);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.978189);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.92935);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.990144);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.990442);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.992722);
    for(int j=0; j<4; j++)	    m_sf.push_back(1.0029);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.996391);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.989659);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.988802);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.991695);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.98634);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.915818);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.993528);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.994529);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.993491);
    for(int j=0; j<4; j++)	    m_sf.push_back(1.00459);
    for(int j=0; j<4; j++)	    m_sf.push_back(0.996809);


    

}



//*****************************************************************************

/////////////////
// CONSTRUCTOR //
/////////////////


AnalysisMuonEfficiencyScaleFactors::AnalysisMuonEfficiencyScaleFactors(
                                    const std::string & alg,
                                    const std::vector<double> int_lum,
                                    const std::string & unit,
                                    std::string dir)
{
 
  this->init(alg, int_lum,unit,dir);
}

/////////////////
//Print method
////////////////

void AnalysisMuonEfficiencyScaleFactors::PrintValues() const
{
  cout<<"Print method"<<endl;
  cout<<"size region "<<m_region.size()<<endl;
  cout<<"size eta_min "<<m_eta_min.size()<<endl;
  cout<<"size eta_max "<<m_eta_max.size()<<endl;
  cout<<"size phi_min "<<m_phi_min.size()<<endl;
  cout<<"size phi_max "<<m_phi_max.size()<<endl;
  cout<<"size pt_min "<<m_pt_min.size()<<endl;
  cout<<"size pt_max "<<m_pt_max.size()<<endl;
  cout<<"m_sf_stat_err size "<<m_sf_stat_err.size()<<endl;
  cout<<"m_sf_syst_err size "<<m_sf_syst_err.size()<<endl;
  cout<<"m_sf size "<<m_sf.size()<<endl;

  cout<<endl<<"printing values"<<endl;

  for(unsigned int i=0; i<m_region.size(); i++) cout<<"region "<<i<<" "<<m_region[i]<<endl;
  for(unsigned int i=0; i<m_eta_min.size(); i++) cout<<"eta_min "<<i<<" "<<m_eta_min[i]<<endl;
  for(unsigned int i=0; i<m_eta_max.size(); i++) cout<<"eta_max "<<i<<" "<<m_eta_max[i]<<endl;
  for(unsigned int i=0; i<m_phi_min.size(); i++) cout<<"phi_min "<<i<<" "<<m_phi_min[i]<<endl;
  for(unsigned int i=0; i<m_phi_max.size(); i++) cout<<"phi_max "<<i<<" "<<m_phi_max[i]<<endl;
  for(unsigned int i=0; i<m_pt_min.size(); i++) cout<<"pt_min "<<i<<" "<<m_pt_min[i]<<endl;
  for(unsigned int i=0; i<m_pt_max.size(); i++) cout<<"pt_max "<<i<<" "<<m_pt_max[i]<<endl;
  for(unsigned int i=0; i<m_sf_syst_err.size(); i++) cout<<"sf_syst "<<i<<" "<<m_sf_syst_err[i]<<endl;
  for(unsigned int i=0; i<m_sf_stat_err.size(); i++) cout<<"sf stat "<<i<<" "<<m_sf_stat_err[i]<<endl;
  for(unsigned int i=0; i<m_sf.size(); i++) cout<<"sf "<<i<<" "<<m_sf[i]<<endl; 

  return;
}
//*****************************************************************************

////////////////////////
// Pseudo CONSTRUCTOR //
////////////////////////

void AnalysisMuonEfficiencyScaleFactors::init(
                                    const std::string & alg,
                                    const std::vector<double> int_lum,
                                    const std::string & unit,
                                    std::string dir) {

/////////////////
// CHECK INPUT //
/////////////////

    if (alg!="STACO_CB" && alg!="STACO_CB_plus_ST" && alg!="Muid_CB" &&
        alg!="Muid_CB_plus_ST" && alg != "CaloTag" ) {
        cerr << "Class AnalysisMuonEfficiencyScaleFactors, constructor: ERROR!"
             << endl
             << "Unknown algorithm type " << alg <<"!\n"
             << "Allowed names are STACO_CB, STACO_CB_plus_ST, "
             << "Muid_CB, Muid_CB_plus_ST, CaloTag \n";
        exit(1);
    }

    if (unit=="MeV") {
        m_unit = 0.001;
    }
    if (unit=="GeV") {
        m_unit = 1;
    }
    if (unit!="MeV" && unit!="GeV") {
        cerr << "Class AnalysisMuonEfficiencyScaleFactors, constructor: ERROR!"
             << endl
             << "Unsupported unit " << unit <<"!\n"
             << "Allowed units are MeV and GeV.\n";
        exit(1);
    }

///////////////
// VARIABLES //
///////////////

    // default to InstallArea/share for the files if running in Athena
    if (dir == "") {
      char *m_tmparea=getenv("TestArea");
      if (m_tmparea != NULL) {
        dir =  string(m_tmparea)+"/InstallArea/share/";
        std::cout<< "AnalysisMuonEfficiencyScaleFactors: Using default dir: "<<dir <<endl;
      }
    }
    else {
      std::cout<<"AnalysisMuonEfficiencyScaleFactors:: Using user defined path!"<<std::endl;
      std::cout<<"                     "<<dir<<std::endl;
    }
      

    std::vector<std::string> ID_files(11);
    if (dir!="") {
      ID_files[0] = dir+"/"+"ID_periodB.txt";
      ID_files[1] = dir+"/"+"ID_periodD.txt";
      ID_files[2] = dir+"/"+"ID_periodE.txt";
      ID_files[3] = dir+"/"+"ID_periodF.txt";
      ID_files[4] = dir+"/"+"ID_periodG.txt";
      ID_files[5] = dir+"/"+"ID_periodH.txt";
      ID_files[6] = dir+"/"+"ID_periodI.txt";
      ID_files[7] = dir+"/"+"ID_periodJ.txt";
      ID_files[8] = dir+"/"+"ID_periodK.txt";
      ID_files[9] = dir+"/"+"ID_periodL.txt";
      ID_files[10] = dir+"/"+"ID_periodM.txt";
    }
    else {
        ID_files[0] = "ID_periodB.txt";
        ID_files[1] = "ID_periodD.txt";
        ID_files[2] = "ID_periodE.txt";
        ID_files[3] = "ID_periodF.txt";
        ID_files[4] = "ID_periodG.txt";
        ID_files[5] = "ID_periodH.txt";
        ID_files[6] = "ID_periodI.txt";
        ID_files[7] = "ID_periodJ.txt";
        ID_files[8] = "ID_periodK.txt";
	ID_files[9] = "ID_periodL.txt";
        ID_files[10] = "ID_periodM.txt";
    }

    std::vector<std::string> muon_files(11);
    if (dir!="") {
      muon_files[0] = dir+"/"+alg+"_periodB.txt";
      muon_files[1] = dir+"/"+alg+"_periodD.txt";
      muon_files[2] = dir+"/"+alg+"_periodE.txt";
      muon_files[3] = dir+"/"+alg+"_periodF.txt";
      muon_files[4] = dir+"/"+alg+"_periodG.txt";
      muon_files[5] = dir+"/"+alg+"_periodH.txt";
      muon_files[6] = dir+"/"+alg+"_periodI.txt";
      muon_files[7] = dir+"/"+alg+"_periodJ.txt";
      muon_files[8] = dir+"/"+alg+"_periodK.txt";
      muon_files[9] = dir+"/"+alg+"_periodL.txt";
      muon_files[10] = dir+"/"+alg+"_periodM.txt";
    }
    else {
      muon_files[0] = alg+"_periodB.txt";
      muon_files[1] = alg+"_periodD.txt";
      muon_files[2] = alg+"_periodE.txt";
      muon_files[3] = alg+"_periodF.txt";
      muon_files[4] = alg+"_periodG.txt";
      muon_files[5] = alg+"_periodH.txt";
      muon_files[6] = alg+"_periodI.txt";
      muon_files[7] = alg+"_periodJ.txt";
      muon_files[8] = alg+"_periodK.txt";
      muon_files[9] = alg+"_periodL.txt";
      muon_files[10] = alg+"_periodM.txt";
    }


///////////////////////////////////
// READ SCALE FACTORS FROM FILES //
///////////////////////////////////

    vector< vector<int> >    region_ID(ID_files.size());
    vector< vector<double> > eta_min_ID(ID_files.size());
    vector< vector<double> > eta_max_ID(ID_files.size());
    vector< vector<double> > phi_min_ID(ID_files.size());
    vector< vector<double> > phi_max_ID(ID_files.size());
    vector< vector<double> > pt_min_ID(ID_files.size());
    vector< vector<double> > pt_max_ID(ID_files.size());
    vector< vector<double> > sf_ID(ID_files.size());
    vector< vector<double> > sf_stat_err_ID(ID_files.size());
    vector< vector<double> > sf_syst_err_ID(ID_files.size());

    vector< vector<int> >    region_muon(muon_files.size());
    vector< vector<double> > eta_min_muon(muon_files.size());
    vector< vector<double> > eta_max_muon(muon_files.size());
    vector< vector<double> > phi_min_muon(muon_files.size());
    vector< vector<double> > phi_max_muon(muon_files.size());
    vector< vector<double> > pt_min_muon(muon_files.size());
    vector< vector<double> > pt_max_muon(muon_files.size());
    vector< vector<double> > sf_muon(muon_files.size());
    vector< vector<double> > sf_stat_err_muon(muon_files.size());
    vector< vector<double> > sf_syst_err_muon(muon_files.size());

    for (unsigned int k=0; k<ID_files.size(); k++) {
        read_file(ID_files[k], region_ID[k],
                  eta_min_ID[k], eta_max_ID[k],
                  phi_min_ID[k], phi_max_ID[k],
                  pt_min_ID[k], pt_max_ID[k],
                  sf_ID[k], sf_stat_err_ID[k], sf_syst_err_ID[k]);
        read_file(muon_files[k], region_muon[k],
                  eta_min_muon[k], eta_max_muon[k],
                  phi_min_muon[k], phi_max_muon[k],
                  pt_min_muon[k], pt_max_muon[k],
                  sf_muon[k], sf_stat_err_muon[k], sf_syst_err_muon[k]);
    }

// prepare map //
    m_region= region_muon[0];
    m_eta_min = eta_min_muon[0];
    m_eta_max = eta_max_muon[0];
    m_phi_min = phi_min_muon[0];
    m_phi_max = phi_max_muon[0];
    m_pt_min = pt_min_muon[0];
    m_pt_max = pt_max_muon[0];
    m_sf = vector<double>(eta_min_muon[0].size());
    m_sf_stat_err = vector<double>(eta_min_muon[0].size());
    m_sf_syst_err = vector<double>(eta_min_muon[0].size());

//////////////////////////////////////////////////////
// COMPUTE THE SCALE FACTORS FOR THE GIVEN ANALYSIS //
//////////////////////////////////////////////////////

// number of periods to be used //
    unsigned int nb_periods(muon_files.size());

// checks //
    if (int_lum.size()==0) {
        std::cerr << "Class AnalysisMuonEfficiencyScaleFactors, constructor!\n"
            << "ERROR! No luminosities provided!\n";
        exit(1);
    }
    if (int_lum.size()<muon_files.size()) {
        std::cerr << "Class AnalysisMuonEfficiencyScaleFactors, constructor!\n"
            << "WARNING! Luminosities only provided for part of the periods.\n"
            << "Integrated luminosities of later periods will be set to 0.\n";
        nb_periods = int_lum.size();
    }
    if (int_lum.size()>muon_files.size()) {
        std::cerr << "Class AnalysisMuonEfficiencyScaleFactors, constructor!\n"
            << "WARNING! Luminosities only provided for part of the periods.\n"
            << "Integrated luminosities of additional periods will be ignored."
            << std::endl;
    }

// calculate normalization //
    double lumi_tot(0.0);
    for (unsigned int k=0; k<nb_periods; k++) {
        lumi_tot = lumi_tot+int_lum[k];
    }

// calculate average scale factors //
    for (unsigned int k=0; k<nb_periods; k++) {
        double weight = int_lum[k]/lumi_tot;
        for (unsigned int l=0; l<sf_muon[k].size(); l++) {
            if (alg == "CaloTag") m_sf[l] = m_sf[l]+weight*sf_muon[k][l];
            else m_sf[l] = m_sf[l]+weight*sf_ID[k][l]*sf_muon[k][l];

            double variance(sf_muon[k][l]*sf_muon[k][l]*
                            sf_stat_err_ID[k][l]*sf_stat_err_ID[k][l]+
                            sf_ID[k][l]*sf_ID[k][l]*
                            sf_stat_err_muon[k][l]*sf_stat_err_muon[k][l]);
            if (alg == "CaloTag") variance = sf_stat_err_muon[k][l]*sf_stat_err_muon[k][l];

            m_sf_stat_err[l] = m_sf_stat_err[l]+weight*weight*variance;

            if (alg == "CaloTag") m_sf_syst_err[l] = m_sf_syst_err[l]+weight*
                                                        sf_syst_err_muon[k][l];
            else m_sf_syst_err[l] = m_sf_syst_err[l]+
                                            weight*(sf_syst_err_ID[k][l]+
                                                    sf_syst_err_muon[k][l]);
        }
    }

    for (unsigned int l=0; l<m_sf_stat_err.size(); l++) {
        m_sf_stat_err[l] = sqrt(m_sf_stat_err[l]);
    }

}

//*****************************************************************************

////////////////////////
// METHOD scaleFactor //
////////////////////////

double AnalysisMuonEfficiencyScaleFactors::scaleFactor(
                                            const TLorentzVector & tlv) const {

    int bin=get_pt_eta_phi_bin_index(tlv);
    if(bin<0) return 1.;
    return m_sf[bin];

}

//*****************************************************************************

///////////////////////////////////
// METHOD scaleFactorUncertainty //
///////////////////////////////////

double AnalysisMuonEfficiencyScaleFactors::scaleFactorUncertainty(
                                            const TLorentzVector & tlv) const {
    int bin=get_pt_eta_phi_bin_index(tlv);
    if(bin<0) return 0.;
    return m_sf_stat_err[bin];

}

//*****************************************************************************

/////////////////////////////////////////////
// METHOD scaleFactorSystematicUncertainty //
/////////////////////////////////////////////

double AnalysisMuonEfficiencyScaleFactors::scaleFactorSystematicUncertainty(
                                            const TLorentzVector & tlv) const {

    int bin=get_pt_eta_phi_bin_index(tlv);
    if(bin<0) return 0.;
    return m_sf_syst_err[bin]+4.2e-6*tlv.E()*m_unit;

}

//*****************************************************************************

//////////////////////
// METHOD read_file //
//////////////////////

void AnalysisMuonEfficiencyScaleFactors::read_file(
                   const std::string & name,
		   std::vector<int> & region,
                   std::vector<double> & eta_min,
                   std::vector<double> & eta_max,
                   std::vector<double> & phi_min,
                   std::vector<double> & phi_max,
                   std::vector<double> & pt_min,
                   std::vector<double> & pt_max,
                   std::vector<double> & sf,
                   std::vector<double> & sf_stat_err,
                   std::vector<double> & sf_syst_err) const {

///////////
// CLEAR //
///////////

    region.clear();
    eta_min.clear();
    eta_max.clear();
    phi_min.clear();
    phi_max.clear();
    pt_min.clear();
    pt_max.clear();
    sf.clear();
    sf_stat_err.clear();
    sf_syst_err.clear();

////////////////////
// OPEN DATA FILE //
////////////////////

    ifstream infile(name.c_str());
    if (infile.fail()) {
        std::cerr << "Class AnalysisMuonEfficiencyScaleFactors, "
             << "method read_file.\n"
             << "ERROR! Could not open file " << name
             << std::endl;
        exit(1);
    }

// read data file //
    double r_eta_min, r_eta_max, r_phi_min, r_phi_max, r_pt_min, r_pt_max, r_sf;
    double r_sf_err, r_sf_syst_err;
    int r_region;
    while (!infile.eof()) {
        infile >> r_pt_min >> r_pt_max >> r_region
               >> r_eta_min >> r_eta_max
               >> r_phi_min >> r_phi_max 
               >> r_sf >> r_sf_err >> r_sf_syst_err;
        if (infile.eof()) {
            break;
        }
	region.push_back(r_region);
        pt_min.push_back(r_pt_min);
        pt_max.push_back(r_pt_max);
        eta_min.push_back(r_eta_min);
        eta_max.push_back(r_eta_max);
        phi_min.push_back(r_phi_min);
        phi_max.push_back(r_phi_max);
        sf.push_back(r_sf);
        sf_stat_err.push_back(r_sf_err);
        sf_syst_err.push_back(r_sf_syst_err);
    }

    return;

}

//*****************************************************************************

/////////////////////////////////////
// METHOD get_pt_eta_phi_bin_index //
/////////////////////////////////////

int AnalysisMuonEfficiencyScaleFactors::get_pt_eta_phi_bin_index(
                                            const TLorentzVector & tlv) const {
					    
    int region_bin=m_EPbin.symmetricBin(&tlv);
    if(tlv.Eta()<0) region_bin=-region_bin;
    
    for (unsigned int k=0; k<m_pt_min.size(); k++) {
        if (m_unit*tlv.Pt()<m_pt_min[k] || m_unit*tlv.Pt()>m_pt_max[k]) {
            continue;
        }
	if(m_region[k]<100)
	{
	  if(region_bin==0) continue;
	  if(region_bin==m_region[k]) return k;
	  else continue;
	}
        if (tlv.Eta()<m_eta_min[k] || tlv.Eta()>m_eta_max[k]) {
            continue;
        }
        if (tlv.Phi()<m_phi_min[k] || tlv.Phi()>m_phi_max[k]) {
            continue;
        }
        return k;
    }

    std::cerr << "Class AnalysisMuonEfficiencyScaleFactors, "
             << "get_pt_eta_phi_bin_index!\n"
             << "ERROR! No bin with scale factor information could be found."
             << "\npt, eta, phi: "
             << m_unit*tlv.Pt() << "\t"
             << tlv.Eta() << "\t"
             << tlv.Phi()
             << std::endl;
             
     return -1;
    //exit(1);

}

}
}
