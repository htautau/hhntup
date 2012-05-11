#include "EnergyRescaler.h"

#include <vector>
#include <iostream>
#include <fstream>

#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <cctype>


using namespace std;

namespace AnalysisFramework
{
namespace External
{

///default constructor
EnergyRescaler::EnergyRescaler()
{

   //default seed
   SetRandomSeed();
}


//destructor
EnergyRescaler::~EnergyRescaler()
{


}

void EnergyRescaler::SetRandomSeed( unsigned seed){

   m_random3.SetSeed(seed);

}

  
bool EnergyRescaler::readCalibConstants(std::string fname)
{ 

   
   if( m_corrVec.size()) {
      std::cout<<" WARNING having already  "<<m_corrVec.size()<<"  corrections "<<std::endl;
      m_corrVec.clear();

   }

   
   std::ifstream infile; 
  
   infile.open(fname.c_str()); 


   if(!infile.good()) 
   { 
      cout<<"CANT OPEN CUT FILE " << fname <<" GOING TO EXIT"<<endl; 
      exit(1); 
   }     
   
   std::cout<<"READING FILE  "<< fname<<std::endl;

   calibMap myMap;


   while( !infile.eof() )
   {

      double eta=-10.,etaErr=-10.,  phi=-10., phiErr=-10., alpha=-100., err=-100.;

     

      infile >>eta>>etaErr>>phi>> phiErr>>alpha>>err;
      if( !infile.good()  ) break; 

      myMap.eta= eta;
      myMap.etaBinSize =etaErr;
      myMap.phi = phi;
      myMap.phiBinSize = phiErr;
      myMap.alpha = alpha;
      myMap.alphaErr = err;
     
      m_corrVec.push_back(myMap);


   }
  


   return true; 
} 
 

double EnergyRescaler::applyEnergyCorrectionGeV(double eta, double phi, double energy, double et,  int value, std::string ptype) const
{ 

   double corrEnergy=-999.0;
   
   if(m_corrVec.size()==0)
   {
      std::cout<<"NO CORRECTIONS EXISTS, PLEASE EITHER SUPPLY A CORRECTION FILE OR USE THE DEFAULT CORRECTIONS"<<std::endl;
   }

   for (unsigned int i=0; i< m_corrVec.size(); i++)
   {

      if( 
         eta>=( m_corrVec.at(i).eta - m_corrVec.at(i).etaBinSize/2.) && eta< ( m_corrVec.at(i).eta+m_corrVec.at(i).etaBinSize/2.)  &&
         phi>=( m_corrVec.at(i).phi - m_corrVec.at(i).phiBinSize/2.) && phi< ( m_corrVec.at(i).phi+m_corrVec.at(i).phiBinSize/2.) 
         ) 
      { 

         
        
         for(std::string::iterator p = ptype.begin(); ptype.end() != p; ++p)
         *p = toupper(*p);


         double er_up=-99,er_do=0; 
         double scale=0.;


         switch (value)
         {
            default:
            {
               scale=m_corrVec.at(i).alpha;
               break;
            }
            case NOMINAL:
            {
               scale=m_corrVec.at(i).alpha;
               break;
            }
            case ERR_UP:
            {
               scale=m_corrVec.at(i).alpha;
               getErrorGeV(eta,et, er_up, er_do, ptype);
               scale+=er_do;
               break;
            }
            case ERR_DOWN:
            {
               scale=m_corrVec.at(i).alpha;
               getErrorGeV(eta,et, er_up, er_do, ptype);
               scale+=er_up;
               break;
            }
         }

        

         corrEnergy =  energy/(1.+ scale);

         // std::cout<<" eta : "<<eta <<" uncorrected energy : "<<    energy <<" corr energy : "<< corrEnergy<<" scale : "<< m_corrVec.at(i).alpha <<endl;
         break;
      }
   }

   if( corrEnergy==-999.)return energy;
   else return corrEnergy;

} 



void EnergyRescaler::getErrorGeV(double cl_eta,double cl_et, double &er_up, double &er_do, std::string ptype,bool withXMAT,bool withPS) const
{
  // Quick and dirty
  // Need to optimized

  er_up=-1;
  er_do=-1;

  static const int nbins=8;
  static double boundaries[nbins+1]={0,0.6,1.00,1.37,1.52,1.8,2.47,3.2,4.9};

  //systematics 
  static double stat[nbins]               ={0.0010, 0.0020, 0.0020, 777, 0.0020, 0.0020, 0.005, 0.01};
  static double sys_mcclosure[nbins]      ={0.0010, 0.0010, 0.0010, 777, 0.0010, 0.0010, 0.002, 0.002};
  static double sys_comparison[nbins]     ={0.0010, 0.0010, 0.0010, 777, 0.0010, 0.0010, 0.01, 0.008}; 
  static double sys_pileup[nbins]         ={0.0010, 0.0010, 0.0010, 777, 0.0010, 0.0010, 0.001, 0.001}; //check forward? 
  static double sys_medium2tight_up[nbins]={0.0010, 0.0010, 0.0010, 777, 0.0020, 0.0020, 0,0};
  static double sys_loose2tight_forward[nbins]    ={0.0000, 0.0000, 0.0000, 777, 0.0000, 0.0000, 0.012, 0.01};

  static double sys_masscut[nbins]        ={0.0010, 0.0010, 0.0010, 777, 0.0020, 0.0020, 0.002, 0.006};
  static double sys_elecLin[nbins]        ={0.0010, 0.0010, 0.0010, 777, 0.0010, 0.0010, 0.001, 0.001};//forward?
  static double sys_xtalkE1[nbins]        ={0.0010, 0.0010, 0.0010, 777, 0.0010, 0.0010, 0.001, 0.001};//forward?
  static double sys_lowpt  [nbins]        ={0.0100, 0.0100, 0.0100, 777, 0.0100, 0.0100, 0.01,  0.01};
  static double sys_HV     [nbins]        ={0.0000, 0.0000, 0.0000, 777, 0.0000, 0.0000, 0.006, 0.008};

  static double PS_B[nbins]={-0.00025312, -0.0011569, -0.00211677, 777, -0.00175762*2, 000, 000, 000};
//  static double PS_A[nbins]={ 0.00187341, 0.00840421,  0.016034  , 777,   0.0127718*2, 000, 000, 000};

  //Material electron
 
  static double elec_XMAT_a[nbins] = {-0.0083,  -0.013,	-0.025,	777, -0.023, -0.019, -0.042, -0.014  };
  static double elec_XMAT_b[nbins] = {-0.055, -0.019, -0.014, 777, -0.026, -0.019, -0.041, -0.016};
  static double elec_XMAT_MAX[nbins]={  0.003,  0.008,  0.015,  777,   0.010,   0.01  ,0.01, 0.01};  

 //Material photon
  static double pho_XMAT_MAX[nbins]={  0.003,  0.005,  0.010,  777,   0.010,   0.01  ,0, 0.} ;
 static double pho_PS_shift[nbins] ={  0.001,  0.002,  0.003,  777,   0.002*2,   0.000  ,0, 0.} ;

  int bin =-1;
  for(int i=0;i<nbins;i++)
    {
      if(fabs(cl_eta)>=boundaries[i] && fabs(cl_eta)<boundaries[i+1])
	{
	  bin=i;
	  break;
	}
    }
  if(bin==-1) return;

  //==================================
  //crack region
  //==================================
  if(bin==3)
    { 
      er_up=0.05;
      er_do=-0.05;
      return;
    }

  //==================================
  //PS
  //==================================
  double PS_up=0,PS_do=0;
  
  if(withPS==true)
    {          
      if(abs(cl_eta)<1.8)
	{
	  double shift=0;
	  if(ptype=="UNCONVERTED_PHOTON" || ptype =="CONVERTED_PHOTON")
	    {
	  shift=pho_PS_shift[bin];
	    }
	  double PS=PS_B[bin]*(log(cl_et)-log(40))-shift;
	  if(PS>=0)
	    {
	      PS_up=PS;
	      PS_do=-PS;
	    }
	  else
	    {
	      PS_up=-PS;
	      PS_do=PS;
	    }
	}
      
    }


  //==================================
  //material
  //==================================
  double XMat_up = 0;
  double XMat_do = 0;
  
  if(withXMAT==true)
    {
      
      if(ptype=="ELECTRON")
	{
	  double xmat= 0;

	  //	  cout<<elec_XMAT_a[bin]<<" "<<elec_XMAT_b[bin] <<endl;
	  xmat= elec_XMAT_a[bin]*exp(cl_et* elec_XMAT_b[bin])+elec_XMAT_MAX[bin];
	  //alpha = a * exp(Pt * b) + c
	  
	  if(xmat>0)
	    {
	      XMat_up=xmat;
	    }
	  else  XMat_do=xmat;

// 	  if(XMat_up>elec_XMAT_MAX[bin])
// 	    XMat_up=elec_XMAT_MAX[bin];	
	  
	}



      
      else if(ptype=="UNCONVERTED_PHOTON" || ptype =="CONVERTED_PHOTON")
	//else if(ptype=="PHOTON")
	{      
	  XMat_up=pho_XMAT_MAX[bin];
	  XMat_do=0;
	  //      cout<<XMat_up<<endl;
	}                 
    }
  //==================================
  //lowpt
  //==================================
  double lowpt=0;
  if(cl_et<20)
    {
      lowpt =sys_lowpt[bin]/(10-20)*(cl_et-20);
    }

  er_up= sqrt(stat[bin]*stat[bin]+
 	      sys_mcclosure[bin]*sys_mcclosure[bin]+
 	      sys_comparison[bin]*sys_comparison[bin]+
 	      sys_pileup[bin]*sys_pileup[bin]+
  	      sys_medium2tight_up[bin]*sys_medium2tight_up[bin]+
	      sys_loose2tight_forward[bin]*sys_loose2tight_forward[bin]+
  	      sys_masscut[bin]*sys_masscut[bin]+
  	      sys_HV[bin]*sys_HV[bin]+
  	      sys_elecLin[bin]*sys_elecLin[bin]+
	      sys_xtalkE1[bin]*sys_xtalkE1[bin]+
 	      PS_up*PS_up+
 	      XMat_up*XMat_up+
	      lowpt*lowpt);
  
  er_do= -sqrt(stat[bin]*stat[bin]+
	      sys_mcclosure[bin]*sys_mcclosure[bin]+
	      sys_comparison[bin]*sys_comparison[bin]+
	      sys_pileup[bin]*sys_pileup[bin]+
	      sys_loose2tight_forward[bin]*sys_loose2tight_forward[bin]+
	      sys_masscut[bin]*sys_masscut[bin]+
  	      sys_HV[bin]*sys_HV[bin]+
	      sys_elecLin[bin]*sys_elecLin[bin]+
	      sys_xtalkE1[bin]*sys_xtalkE1[bin]+
	      PS_do*PS_do+
	      XMat_do*XMat_do+
	      lowpt*lowpt);
//   er_up=XMat_up;
//   er_do=XMat_do;

  return;
}





/////default constants, for 60 eta bins
bool  EnergyRescaler::useDefaultCalibConstants( std::string corr_version)
{
   
   //if(corr_version!="2010")return true;

   const int netaBins=58;
   const int nphiBins=58;
   
   ///for 2011
   const int netaBinsFor2011=32;
   const int nphiBinsFor2011=32;
   
  
   double m_2pi = 2.*acos(-1.); 
   //cout<<" pi : "<<acos(-1.)<<"  M_PI : "<<M_PI<<" 2pi : "<<m_2pi<<endl;
   

  

   if( m_corrVec.size()) {
      std::cout<<" WARNING having already  "<<m_corrVec.size()<<"  corrections "<<std::endl;
      m_corrVec.clear();

   }
   

   

  
   double eta_tmp[netaBins]={   
      -4.45, -3.60, -3.00, -2.65, -2.435, -2.35, -2.25, -2.15, -2.05, -1.95, -1.85,	-1.75,	-1.65,
      -1.56, -1.445,	-1.335,	-1.25, -1.15,	-1.05,	-0.95,	-0.85,	-0.75,	-0.65,	-0.55,	-0.45,	-0.35,
      -0.25, -0.15,	-0.05,	0.05,  0.15,	0.25,	0.35,	0.45,	0.55,	0.65,	0.75,	0.85,	0.95,
       1.05, 1.15,	1.25,	1.335, 1.445,	1.56,	1.65,	1.75,	1.85,	1.95,	2.05,	2.15,	2.25,
       2.35,	2.435,   2.65,	3.00,	3.60,	4.45
         };
   

   


    double etaBinWidth_tmp[netaBins]= {
       0.9, 0.8, 0.4, 0.3, 0.07,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.08,	0.15,	0.07,	0.1,	0.1,	0.1,	        0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	        0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.1,	0.07,	0.15,	0.08,	0.1,	0.1,	0.1,	0.1,	        0.1,	0.1,	0.1,	0.1,	0.07,    0.3,	0.4,	0.8,	0.9 };

   
    

   
    double phi_tmp[nphiBins]={
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 	0, 0, 0, 0, 0, 0, 0, 0,	0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0,	0, 0 };

    double phiBinWidth_tmp[nphiBins]={
       m_2pi, m_2pi, m_2pi, m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,
       m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi, m_2pi,m_2pi,m_2pi,m_2pi,
       m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi, m_2pi,
       m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,
       m_2pi,m_2pi, m_2pi,m_2pi,m_2pi,m_2pi,m_2pi, m_2pi, m_2pi, m_2pi, m_2pi, m_2pi };


    
    

   double alpha_tmp[netaBins]= {
    0.046962, 0.044332,	-0.007777,-0.045341, -0.053683,-0.005040,0.003357,0.003830,0.016857,0.012810,
    0.008882,-0.016864,-0.021280,-0.011215,-0.008958,-0.010248, 0.005199,0.007815,-0.000701,0.003837,0.003222,
    -0.004017,-0.005908,-0.009560,-0.005382,-0.012158,-0.005060,-0.008274, -0.006701,-0.002734,-0.012920, 
    -0.010972,-0.006823,-0.007234,-0.002612,-0.004301,0.001580,-0.001986,-0.001306,0.005748, 0.002906, 
    0.001381,-0.001584,0.000799,-0.002511,-0.030902,-0.016416,-0.004976,0.002408,0.018706,-0.004309,-0.002673,
    -0.001735,-0.050173, -0.034703,	-0.003296,	0.047351,	0.028374
   };

   

   

   double alphaError_tmp[netaBins]={ 
      0.019646, 0.011500, 0.008761, 0.009514, 0.008341,0.006345,0.005522,0.005089,0.005387,0.005861,
      0.005821,0.005119,0.006227,0.009148,0.003920,0.004309,0.002929,  0.003882,0.004054,0.003716,0.003673,
      0.003832,0.003275,0.002075,0.004004,0.002497,0.003182,0.002512,0.003974,0.002302, 0.003670,0.003322,
      0.003978,0.002164,0.001984,0.002093,0.002372,0.003843,0.004138,0.004277,0.003003,0.004690,0.005480,
      0.006306,	0.007338,0.005939,0.004472,0.004535,0.005453,0.008538,0.004554,0.003382,0.005504,0.007577,
      0.010095,	0.009122,	0.013400,	0.02588
   };


   //////NUMBERS for 2011 data
   double etaFor2011_tmp[netaBinsFor2011]={ 
      -4.05,	-3.00,	-2.65,	-2.385,	-2.2,	-2.00,	-1.8,	-1.61,	-1.445,	-1.285,	-1.1,	-0.9,	-0.7,	-0.5,	-0.3,	-0.1,	0.1,	0.3,	0.5,	0.7,	0.9,	1.1,	1.285,	1.445,	1.61,	1.8,	2.0,	2.2,	2.385,	2.65,	3.0,	4.05
   };

    double etaBinWidthFor2011_tmp[netaBinsFor2011]= {
       1.7,	0.4,	0.3,  0.17,	0.2,	0.2,	0.2,	0.18,	0.15,	0.17,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.2,	0.17,	0.15,	0.18,	0.2,	0.2,	0.2,	0.17,   0.3,	0.4,	1.7
    };

    double phiFor2011_tmp[nphiBinsFor2011]={
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

    double phiBinWidthFor2011_tmp[nphiBinsFor2011]={
       m_2pi, m_2pi, m_2pi, m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,
       m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi, m_2pi, m_2pi, m_2pi, m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,
       m_2pi,m_2pi,m_2pi,m_2pi,m_2pi,m_2pi
       };


    
    double alphaFor2011_tmp[netaBinsFor2011]= {
       0.01767670,
       0.01193680,
       0.02863680,
       -0.00182836,
       0.00379822,
       0.00949851,
       0.00388223,
       -0.00630638,
       -0.01459190,
       0.00191748,
       0.00120518,
       -0.00116847,
       0.00316491,
       0.00534383,
       -0.00259050,
       -0.00303995,
       -0.00114055,
       -0.00152093,
       0.00346042,
       0.00316020,
       -6.51986e-05,
       0.00237845,
       0.00391017,
       -0.01056940,
       -0.01190480,
       -0.00177316,
       -0.00101658,
       -0.00430559,
       -0.01101630,
       0.01681920,
       0.01459770,
       0.00328200
       
    };
   

    double alphaErrorFor2011_tmp[netaBinsFor2011]={ 
      0.00365611,
      0.00269144,
      0.00253709,
      0.00124745,
      0.00094802,
      0.00094137,
      0.00097271,
      0.00130736,
      0.00170269,
      0.00091738,
      0.00076709,
      0.00072777,
      0.00067059,
      0.00063910,
      0.00064532,
      0.00063576,
      0.00066304,
      0.00064027,
      0.00066183,
      0.00068352,
      0.00075179,
      0.00080414,
      0.00093506,
      0.00167934,
      0.00131510,
      0.00097932,
      0.00095425,
      0.00097240,
      0.00127444,
      0.00258217,
      0.00265290,
      0.00360700

   };




    calibMap myMap;
    
    if(corr_version!="2010"){
       
       for(int i=0; i<netaBinsFor2011; i++){

          //for(int j=0; j<nphiBins; j++){///this should be a 2d array finally
          //}

          myMap.eta= etaFor2011_tmp[i];
          myMap.etaBinSize =etaBinWidthFor2011_tmp[i];
          myMap.alpha = alphaFor2011_tmp[i];
          myMap.phi = phiFor2011_tmp[i];
          myMap.phiBinSize = phiBinWidthFor2011_tmp[i];
          myMap.alphaErr = alphaErrorFor2011_tmp[i];

          m_corrVec.push_back(myMap);

       }
        

    }else{

       
       for(int i=0; i<netaBins; i++){

          //for(int j=0; j<nphiBins; j++){///this should be a 2d array finally
          //}

          myMap.eta= eta_tmp[i];
          myMap.etaBinSize =etaBinWidth_tmp[i];
          myMap.alpha = alpha_tmp[i];
          myMap.phi = phi_tmp[i];
          myMap.phiBinSize = phiBinWidth_tmp[i];
          myMap.alphaErr = alphaError_tmp[i];

          m_corrVec.push_back(myMap);

       }


    }//if-else


   return true;
}



bool EnergyRescaler::printMap() const
{

   for (unsigned int i=0; i< m_corrVec.size(); i++)
   {
      std::cout<<"eta :  "<< m_corrVec.at(i).eta <<
         " etaErr : " <<m_corrVec.at(i).etaBinSize<<
         " phi    :  "<<m_corrVec.at(i).phi<<
         " phiErr :  "<<m_corrVec.at(i).phiBinSize<<
         " alpha  :  "<<m_corrVec.at(i).alpha<<
         " alphaErr : "<<m_corrVec.at(i).alphaErr<<endl;
   }


   return true;

}


// NEW FUNCTIONS START HERE (MB)


// sampling term inMC, parametrization from Iro

double EnergyRescaler::mcSamplingTerm(double cl_eta) {

  double aeta = fabs( cl_eta );
  double sampling = 0.;

  if ( aeta<0.8 )
    sampling = 0.091;

  else if ( aeta<1.37 )
    sampling = 0.036 + 0.130 * aeta;

  else if ( aeta<1.52 )
    sampling = 0.27;

  else if ( aeta<2.0 )
    sampling = 0.85 - 0.36 * aeta;

  else if ( aeta<2.3 )
    sampling = 0.16;
  
  else if ( aeta<2.5 )
    sampling = -1.05 + 0.52 * aeta;

  return sampling;

}


// sampling term uncertainty

double EnergyRescaler::mcSamplingTermRelError( double cl_eta ) {

  cl_eta = cl_eta*1.;
  double relerr = 0.1;

  return relerr;

}


// noise term in MC (from Iro)

double EnergyRescaler::mcNoiseTerm( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double noise = 0.;

  double noise37[25]={ 0.27, 0.27, 0.27, 0.27, 0.27, 
		       0.26, 0.25, 0.23, 0.21, 0.19, 
		       0.17, 0.16, 0.15, 0.14, 0.27, 
		       0.23, 0.17, 0.15, 0.13, 0.10, 
		       0.07, 0.06, 0.05, 0.04, 0.03 };  

  int ieta = (int) (aeta/0.1);

  if ( ieta<25 )
    noise =  noise37[ieta];

  return noise;

}


// constant term in MC (local)

double EnergyRescaler::mcConstantTerm( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double cst = 0.;

  if ( aeta<0.6 )
    cst = 0.005;

  else if ( aeta<1.75 )
    cst = 0.003;

  else if ( aeta<2.5 )
    cst = 0.0055 * (2.69 - aeta);

  //cst = 0.005;

  return cst;

}


// constant term fitted in data (long range)

double EnergyRescaler::dataConstantTerm( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double cst = 0.;

  if ( aeta<1.37 )
    cst = 0.011;

  else if ( aeta<2.47 )
    cst = 0.018;

  else if ( aeta<3.2 )
    cst = 0.06;

  else if ( aeta<4.9 )
    cst = 0.02;

  return cst;

}


double EnergyRescaler::dataConstantTermError( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double err = 0.;

  if ( aeta<1.37 )
    err = 0.006;

  else if ( aeta<2.5 )
    err = 0.006;
  
  return err;
}


double EnergyRescaler::dataConstantTermUpError( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double err = 0.;

  if ( aeta<1.37 )
    err = 0.005;

  else if ( aeta<2.5 )
    err = 0.006;
  
  return err;
}


double EnergyRescaler::dataConstantTermDownError( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double err = 0.;

  if ( aeta<1.37 )
    err = -0.006;

  else if ( aeta<2.5 )
    err = -0.006;
  
  return err;
}


// fitted Z peak resolution, data, in GeV

double EnergyRescaler::dataZPeakResolution( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double res = 0.;

  if ( aeta<1.37 )
    res = 1.62;

  else if ( aeta<2.5 )
    res = 1.99;

  return res;

}


// fitted Z peak resolution, MC, in GeV

double EnergyRescaler::mcZPeakResolution( double cl_eta ) {

  double aeta = fabs( cl_eta );
  double res = 0.;

  if ( aeta<1.37 )
    res = 1.45;

  else if ( aeta<2.5 )
    res = 1.63;

  return res;

}


// correlated part of constant term uncertainty, in data (approx.)

double EnergyRescaler::dataConstantTermCorError( double cl_eta ) {

  double mz = 91.2;
  
  double resData = dataZPeakResolution( cl_eta );
  double resMC   = mcZPeakResolution( cl_eta );
  double cmc     = mcConstantTerm( cl_eta );

  double smpup = 1. + mcSamplingTermRelError( cl_eta );
  double smpdo = 1. - mcSamplingTermRelError( cl_eta );

  double central = sqrt( 2*(resData*resData - resMC*resMC)/mz/mz + cmc*cmc );
  double vardown = sqrt( 2*(resData*resData - resMC*resMC*smpup*smpup)/mz/mz + cmc*cmc );
  double varup   = sqrt( 2*(resData*resData - resMC*resMC*smpdo*smpdo)/mz/mz + cmc*cmc );

  double errdown = fabs( central - vardown );
  double errup   = fabs( central - varup );

  return max( errup, errdown );

}


// total resolution uncertainty (fractional)

void EnergyRescaler::resolutionError( double energy, double cl_eta, double& errUp, double& errDown ) {

  double Cdata     = dataConstantTerm( cl_eta );
  double Cdata_cor = dataConstantTermCorError( cl_eta );
  double Cdata_err = dataConstantTermError( cl_eta );
  
  double Cdata_unc = 0.;
  if( Cdata_err > Cdata_cor )
    Cdata_unc = sqrt( Cdata_err*Cdata_err - Cdata_cor*Cdata_cor );
  if( Cdata_unc<0.001 )
    Cdata_unc = 0.001; // preserve at least the stat error

  double Smc       = mcSamplingTerm( cl_eta );
  double Smc_err   = mcSamplingTermRelError( cl_eta );

  double central = fcn_sigma( energy,  Cdata,  0.,  Smc,  0.);

  double err1 = fcn_sigma( energy, Cdata,  Cdata_unc, Smc,  0.    ) - central; 
  double err2 = fcn_sigma( energy, Cdata, -Cdata_unc, Smc,  0.    ) - central; 
  double err3 = fcn_sigma( energy, Cdata, -Cdata_cor, Smc, Smc_err) - central;
  double err4 = -err3;

  errUp = 0;
  if( err1>0 ) errUp = sqrt( errUp*errUp + err1*err1);
  if( err2>0 ) errUp = sqrt( errUp*errUp + err2*err2);
  if( err3>0 ) errUp = sqrt( errUp*errUp + err3*err3);
  if( err4>0 ) errUp = sqrt( errUp*errUp + err4*err4);

  errDown   = -errUp;
}


// total resolution (fractional)

double EnergyRescaler::resolution( double energy, double cl_eta, bool withCT ) {

  double a = mcSamplingTerm( cl_eta );
  double b = mcNoiseTerm( cl_eta );
  double c = mcConstantTerm( cl_eta );
  double d = dataConstantTerm( cl_eta );

  double sig2 = a*a/energy + b*b/energy/energy + c*c;
  if( withCT )
    sig2 += d*d;

  return sqrt(sig2);

}


// internal use only

double EnergyRescaler::fcn_sigma(double energy, double Cdata, double Cdata_er, double S, double S_er) {

  double sigma2 = std::pow((Cdata+Cdata_er)*energy,2) + std::pow(S*(1+S_er)*std::sqrt(energy),2);
  
  double sigma=0; 
  if (sigma2>0) 
    sigma=sqrt(sigma2);
  
  return sigma/energy;

}


// derive smearing correction

double EnergyRescaler::getSmearingCorrectionGeV(double eta, double energy, int value, bool mc_withCT, std::string /*corr_version*/) 
{
  double resMC, resData, resVar, errUp, errDown;
  resMC   = resolution( energy, eta, false );
  resData = resolution( energy, eta, true );
  resolutionError( energy, eta, errUp, errDown );

  double Cmc = 0.007;

  //=====================================
  //choose which value to use
  //=====================================

  if( value == 1 )
    resData += errDown;
  else if( value == 2 )
    resData += errUp;
  else if( value != 0 ) {
    std::cout << "getSmearingCorrection : wrong value, return 1" << endl;
    return 1;
  }
  
  //=====================================
  //Smearing procedure
  //=====================================

  double sigma2 = std::pow( resData*energy, 2 ) - std::pow( resMC*energy, 2 );
  if (mc_withCT==true) 
    sigma2 = sigma2 - std::pow( Cmc*energy, 2 );
  
  if (sigma2<=0) return 1;
  
  double sigma = sqrt(sigma2);
  double DeltaE0 = m_random3.Gaus(0,sigma);

  double cor0=(energy+DeltaE0)/energy;
  
  return cor0;
  
}


// a calibration correction for crack electrons, to be applied to both data and MC

double EnergyRescaler::applyMCCalibrationGeV(double eta, double ET, std::string ptype) {

  if( ptype != "ELECTRON" )
    return 1.;

  double aeta = fabs(eta);

  if( aeta<1.42 || aeta>1.55 )
    return 1.;

  const int nBoundaries = 18;
  double ETBoundaries[nBoundaries]   = { 0., 5.4, 8.5, 12.9, 16., 20., 
					 25., 30., 35., 40., 45., 50., 
					 55., 60., 65., 70., 75., 99999. };
  
  double CalibFactors[nBoundaries-1] = { 0.884845, 0.898526, 0.902439, 0.91899, 0.925868, 0.929440, 
					 0.948080, 0.943788, 0.96026, 0.955709, 0.964285, 0.95762, 
					 0.970385, 0.963489, 0.968149, 0.970799, 0.961617 };

  int i0 = -1;

  for ( int i=0; i<nBoundaries; i++)
    if ( ET>ETBoundaries[i] && ET<=ETBoundaries[i+1] )
      i0 = i;

  if( i0>=0 && i0<nBoundaries )
    return 1./CalibFactors[i0];

  return 1.;

}


// AF -> G4 correction

double EnergyRescaler::applyAFtoG4(double eta) {

  double aeta = fabs(eta);
  if( aeta>2.47 )
    return 1.;
  
  const int nBoundaries = 27;
  double EtaBoundaries[nBoundaries] = { -2.47, -2.3, -2.1, -1.9, -1.7, -1.52, -1.37, -1.2, 
					-1.0, -0.8, -0.6, -0.4, -0.2,  0.0,   0.2,   0.4, 
					0.6,  0.8,  1.0,  1.2,  1.37, 1.52, 1.7, 1.9, 2.1, 2.3, 2.47 };
  
  double CalibFactors[nBoundaries-1] = { 1.01228, 1.00454, 1.00427, 1.00744, 1.00244, 1.00424, 
					 0.993664, 1.00191, 1.00366, 0.999175, 0.997696, 0.997001, 
					 0.999441, 0.999655, 0.996957, 0.997899, 0.999005, 1.00341, 
					 1.00156, 0.993236, 0.996897, 1.00218, 1.00679, 1.00396,
					 1.00353, 1.01224 };
  
  int i0 = -1;
  
  for ( int i=0; i<nBoundaries; i++)
    if ( eta>EtaBoundaries[i] && eta<=EtaBoundaries[i+1] )
      i0 = i;
  
  if( i0>=0 && i0<nBoundaries )
    return CalibFactors[i0];
  
  return 1.;
  
}  

}
}
