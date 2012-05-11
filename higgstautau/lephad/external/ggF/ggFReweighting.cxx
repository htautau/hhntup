#include "ggFReweighting.h"

namespace AnalysisFramework
{
namespace External
{

//using namespace std;
ggFReweighting::ggFReweighting(std::string generator_, int higgsMass_, std::string option_, std::string inputFilePath_, std::string mc_)
{
  _inputFilePath = inputFilePath_;
  initialize(generator_, higgsMass_, option_, mc_);
}

ggFReweighting::~ggFReweighting()
{
  ;
}

void ggFReweighting::initialize(std::string generator_, int higgsMass_, std::string option_, std::string mc_)
{
  //// generator and higgsMass
  _generator = generator_;
  _higgsMass = higgsMass_;
  _option = option_;
  _mc = mc_;
  if (mc_=="mc11")  _input = "ggFHiggsPtWeight_mc11.root";
  else _input = "ggFHiggsPtWeight_mc10.root";

  std::vector<int> _errMass;
  _errMass.push_back( 120 );
  _errMass.push_back( 130 );
  _errMass.push_back( 160 );
  _errMass.push_back( 200 );
  _errMass.push_back( 300 );
  _errMass.push_back( 400 );
  _errMass.push_back( 500 );
  _errMass.push_back( 600 );

  initializeGeneratorLib(_generatorLib);
  if(_generatorLib.find((_generator)) == _generatorLib.end())
    {
      std::cout << _generator  << " process is not available, please check !!!" << std::endl;
      exit(1);
    }
  initializeHiggsMassLib(_higgsMassLib, _generator, _input);
  if(_higgsMassLib.find(_higgsMass) == _higgsMassLib.end())
    {
      std::cout << "Higgs Mass = "
		<< _higgsMass << " is not available for "  << (_generator)
		<< " please check !!!" << std::endl;
      exit(1);
    }

  ////// check the option /////////////
  if(_option != "Mean" )
    {
      bool findPdfSet = false;
	  bool findScale = false;
      char PdfSetName[200];
	  char ScaleName[200];
      for(int i = 0; i < 40; i++)
	{
	  sprintf(PdfSetName, "Pdfset%d", i+1);
	  if(PdfSetName == _option)
	    {
	      findPdfSet = true;
	      break;
	    }
	}
	  for(int i=0; i<14; i++){
		sprintf(ScaleName, "Scale%d", i);
		if (ScaleName == _option){
			findScale = true;
			break;
		}
	  }

      if(!findPdfSet && !findScale)
	{
	  std::cout << "Option = " << _option << " is not available !" << std::endl;
	  exit(1);
	}
      if( std::find( _errMass.begin(), _errMass.end(), _higgsMass ) == _errMass.end() ){
	std::cout<< "Warning, error sets not available for this mass!" << std::endl;
	exit(1);
      }

    }

  //// initialize availabe generator and mass points
  TFile* file = new TFile((_inputFilePath + _input).c_str(), "read");
  char histName[200];

  if(_option == "Mean") sprintf(histName, "%s_H%d", _generator.c_str(), _higgsMass);
  else if(_option =="Scale0") sprintf(histName, "%s_H%d", _generator.c_str(), _higgsMass);
  else sprintf(histName, "%s_H%d_%s", _generator.c_str(), _higgsMass, _option.c_str());
  _weightHist = (TH1D*)file->Get(histName);
  _weightHist->SetDirectory(0);
  file->Close();
}

int ggFReweighting::higgsMass()
{
  return _higgsMass;
}

std::string ggFReweighting::option()
{
  return _option;
}

std::string ggFReweighting::generator()
{
  return _generator;
}

std::string ggFReweighting::mc(){
  return _mc;
}

double ggFReweighting::getWeight(double higgsPt)
{
  if ( _weightHist->GetBinContent(_weightHist->FindBin(higgsPt)) > 0 ){
    return _weightHist->GetBinContent(_weightHist->FindBin(higgsPt));
  }
  else return 0;
}

double ggFReweighting::getStatError(double higgsPt)
{
  return _weightHist->GetBinError(_weightHist->FindBin(higgsPt));
}

std::pair<double, double> ggFReweighting::getWeightAndStatError(double higgsPt)
{
  std::pair<double, double> result(getWeight(higgsPt), getStatError(higgsPt));
  return result;
}

void ggFReweighting::printInfo()
{
  std::set<int>::iterator itMass;
  std::set<std::string>::iterator itGenerator;
  std::cout << " ggFReweiting : reweight normalized Higgs Pt to NNLL+NNLL by HqT " << std::endl;
  std::cout << " Statistical error: generator statistical error " << std::endl;
  std::cout << " Systematic error : "<< std::endl;
  std::cout << "  (1) Pdf Uncertianty  MSTW2008 NNLO 90%C.L eigenvectors: PdfsetN, N={1,40} " << std::endl;
  std::cout << "  (2) Scale uncertainty:  ScaleN, N={1,13}" << std::endl;
  std::cout << " Available generators  : " << std::endl;
  std::set<std::string> generatorLib;
  std::set<int> higgsMassLib;
  initializeGeneratorLib(generatorLib);
  for(itGenerator = generatorLib.begin(); itGenerator != generatorLib.end(); itGenerator++)
    {
      std::cout << std::setw(8) << *itGenerator << std::endl;
      initializeHiggsMassLib(higgsMassLib, *itGenerator, _input);
      std::cout << " Available mass points : [";
      for(itMass = higgsMassLib.begin(); itMass != higgsMassLib.end(); itMass++)
	{
	  std::cout << std::setw(5) << *itMass;
	}
      std::cout << " ] GeV" << std::endl;

    }
  std::cout<<"Available error points (scale and pdf): [ 120, 130, 160, 200, 300, 400, 500, 600 ] GeV"<<std::endl;
  std::cout << std::endl;
}


void initializeGeneratorLib(std::set<std::string> &generatorLib)
{
  generatorLib.clear();
  generatorLib.insert("McAtNlo");
  generatorLib.insert("PowHeg");
  generatorLib.insert("PowHegMSSM");
}

void initializeHiggsMassLib(std::set<int> &higgsMassLib, std::string generator, std::string input)
{
  int i;
  int mcatnlo_mass [] = {110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600};
  int powheg_mass[] =  {100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155,
		   160, 165, 170, 175, 180, 185, 190, 195, 200, 210,
		   220, 240, 260, 280, 300, 320, 340, 360, 380, 400,
		   420, 440, 460, 480, 500, 520, 540, 560, 580, 600};
  int powhegMSSM_mass[] = { 90, 100, 110, 120, 130, 140, 150, 170,
			    200, 250, 300, 350, 400, 450, 500, 600};
  higgsMassLib.clear();
  if(generator == "McAtNlo" && input.find("mc10") != std::string::npos)
    {
      for(i = 0; i < 24; i++) higgsMassLib.insert(mcatnlo_mass[i]);
    }
  if(generator == "McAtNlo" && input.find("mc11") != std::string::npos)
    {
      for(i = 0; i < 39; i++) higgsMassLib.insert(mcatnlo_mass[i]);
    }

  if(generator == "PowHeg")
    {
      for(i = 0; i < 40; i++) higgsMassLib.insert(powheg_mass[i]);
    }

  if(generator == "PowHegMSSM")
    {
      for(i = 0; i < 16; i++) higgsMassLib.insert(powhegMSSM_mass[i]);
    }
}

//ClassImp(ggFReweighting)

}
}

