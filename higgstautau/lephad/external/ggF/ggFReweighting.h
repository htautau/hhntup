#ifndef __MYGGFREWEIGHTING_H
#define __MYGGFREWEIGHTING_H

#include <vector>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unistd.h>
#include <stdlib.h>
#include <iomanip>

#include "TMath.h"
#include "TFile.h"
#include "TH1.h"

namespace AnalysisFramework
{
namespace External
{

class ggFReweighting {
 public :  
  ggFReweighting(std::string generator, int higgsMass, std::string option = "Mean", std::string inputFilePath = "./", std::string mc="mc11");

  ~ggFReweighting();
  
  void initialize(std::string generator, int higgsMass, std::string option = "Mean", std::string mc="mc11");
  
  int higgsMass();
  
  std::string option();

  std::string generator();

  std::string mc();

  double getWeight(double higgsPt);
  
  double getStatError(double higgsPt);

  std::pair<double, double> getWeightAndStatError(double higgsPt);

  void printInfo();
  ClassDef(ggFReweighting, 0);

 private :
  std::string _generator;
  int _higgsMass;
  TH1D* _weightHist;
  std::string _inputFilePath;
  std::set<std::string> _generatorLib;
  std::set<int> _higgsMassLib;
  std::string _option;
  std::string _input;
  std::string _mc;
};

void initializeGeneratorLib(std::set<std::string> &generatorLib);

void initializeHiggsMassLib(std::set<int> &higgsMassLib, std::string generator, std::string input);

}
}

#endif
 
