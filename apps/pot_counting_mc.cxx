#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "WCPLEEANA/pot.h"

using namespace std;
using namespace LEEana;



int main(int argc, char** argv){

  if (argc < 2){
    std::cout << "pot_counting_mc #mc_files " << std::endl;
    return -1;
  }
  TString mc_file = argv[1];
  TFile *file = new TFile(mc_file);
  TTree *T_eval;
  TTree *T_pot;
  if (file->Get("wcpselection/T_eval")){
    T_eval = (TTree*)file->Get("wcpselection/T_eval");
    T_pot = (TTree*)file->Get("wcpselection/T_pot");
  }else{
    T_eval = (TTree*)file->Get("wcpselection/T_eval_cv");
    T_pot = (TTree*)file->Get("wcpselection/T_pot_cv");
  }
  std::cout << mc_file << std::endl;
  std::cout << "Trigger: " << T_eval->GetEntries() << std::endl;
  double total_pot = 0;
  double pot_tor875;
  //Erin
  int r = 0, s = 0;
  T_pot->SetBranchAddress("pot_tor875",&pot_tor875);
  T_pot->SetBranchAddress("runNo",&r);
  T_pot->SetBranchAddress("subRunNo",&s);

  std::map<std::pair<int, int>, bool >  already_seen;

  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);
    auto it = already_seen.find(std::make_pair(r,s));
    //if (mc_file.Contains("single_photon")){
      if (it != already_seen.end()) continue;
    //}
    already_seen[std::make_pair(r,s)] = true;
    total_pot += pot_tor875;
  }
  std::cout << "POT: " << total_pot << std::endl;

  return 0;
}
