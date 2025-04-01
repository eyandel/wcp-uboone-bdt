// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <set>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/tagger.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"


#include "WCPLEEANA/eval.h"

using namespace std;
using namespace LEEana;

#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"

int main( int argc, char** argv )
{
  if (argc < 2) {
    std::cout << "bdt_convert #input_file" << std::endl;
    return -1;
  }
  TString input_file = argv[1];

  TFile *file1 = new TFile(input_file);
  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_weight = (TTree*)file1->Get("wcpselection/T_weight");

  T_BDTvars->AddFriend(T_weight);
  TH1F *h1 = new TH1F("h1","h1",100,0,100);

  int ncount = T_BDTvars->GetEntries("numu_cc_flag==-1");
  T_BDTvars->Project("h1","All_UBGenie","numu_cc_flag==-1");

  std::cout << input_file << " " << ncount << " " << h1->GetMean() << std::endl;
    
    

  return 0;
}
