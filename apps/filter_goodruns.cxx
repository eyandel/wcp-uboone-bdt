// cz: code modified from tutorials/tmva/TMVAClassification.C
// nn : replicated code from bdt_convert but standalone good runs filter app


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

#include "WCPLEEANA/tree_wrangler.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  if (argc < 3) {
    std::cout << "filter_goodruns #input_file #output_file -n[flag_numi]" << std::endl;
    return -1;
  }

  TString input_file = argv[1];
  TString out_file = argv[2];

  bool flag_data = true;
  int flag_numi = 0;

  bool flag_config = false;
  std::string config_file_name="config.txt";
  char delimiter = ',';

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
      case 'n':
        flag_numi = atoi(&argv[i][2]);
        break;
      case 't':
        config_file_name = &argv[i][2];
        flag_config = true;
        break;
      case 'd':
        delimiter = argv[i][2];//In case you want to change what character you use to sperate your trees in the config
        break;
    }
  }

  tree_wrangler wrangler(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_pot(flag_config, config_file_name, delimiter,true);

  TFile *file1 = new TFile(input_file);

  TTree *T_BDTvars = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_eval = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_spacepoints = (TTree*)file1->Get("wcpselection/T_spacepoints");

  //Load other trees from directories as specified by the config file
  std::vector<TTree*>* old_trees = new std::vector<TTree*>;
  old_trees = wrangler.get_old_trees(file1);
  std::vector<TTree*>* old_trees_pot = new std::vector<TTree*>;
  old_trees_pot = wrangler_pot.get_old_trees(file1);

  if (T_eval->GetBranch("weight_cv")) flag_data =false;

  if(flag_data)
    std::cout << "Recognized input file as data! Filtering.." << std::endl;
  else{
    std::cout << "Not a data file! Exiting..." << std::endl;
    return 0;
  }


  int run = 0;
  int runNo = 0;
  T_eval->SetBranchAddress("run", &run);
  T_pot->SetBranchAddress("runNo", &runNo);

  std::vector<int>good_run_list_vec = wrangler.get_good_run_list();
  std::set<int> good_runlist_set(good_run_list_vec.begin(), good_run_list_vec.end());
  
  std::vector<int> low_lifetime_runs = wrangler.get_low_lifetime_runs();
  std::set<int> low_lifetime_set(low_lifetime_runs.begin(), low_lifetime_runs.end());
  
  std::vector<int> low_neutrino_count_numi_run2RHC = wrangler.get_low_neutrino_count_numi_run2RHC();
  std::set<int> low_neutrino_count_numi_run2RHC_set(low_neutrino_count_numi_run2RHC.begin(), low_neutrino_count_numi_run2RHC.end());

  TFile *file2 = new TFile(out_file,"RECREATE");

  //Setup the directories specified in the config file
  std::vector<TTree*>* new_trees = new std::vector<TTree*>;
  new_trees = wrangler.set_new_trees(file2);
  std::vector<TTree*>* new_trees_pot = new std::vector<TTree*>;
  new_trees_pot = wrangler_pot.set_new_trees(file2);

  file2->mkdir("wcpselection");
  file2->cd("wcpselection");
  TTree *t1 = (TTree*) T_eval->CloneTree(0);
  TTree *t2 = (TTree*) T_pot->CloneTree(0);
  TTree *t3 = (TTree*) T_PFeval->CloneTree(0);
  TTree *t4 = (TTree*) T_BDTvars->CloneTree(0);
  TTree *t5 = (TTree*) T_KINEvars->CloneTree(0);
  TTree *new_T_spacepoints = T_spacepoints->CloneTree(0);

  T_eval->SetBranchStatus("*",1);
  T_pot->SetBranchStatus("*",1);

  std::cout << "Filtering from saved goodruns list for Runs 1-5" << std::endl;

  for (int i=0;i!=T_BDTvars->GetEntries();i++){
    T_BDTvars->GetEntry(i);
    T_eval->GetEntry(i);
    T_KINEvars->GetEntry(i);
    T_PFeval->GetEntry(i);

    if (flag_data){
      if (good_runlist_set.find(run) == good_runlist_set.end()) continue;
      if (low_lifetime_set.find(run) != low_lifetime_set.end()) continue;
      if (flag_numi && low_neutrino_count_numi_run2RHC_set.find(run) != low_neutrino_count_numi_run2RHC_set.end()) continue;
      // bad run in run 1 due to beam filter bnb
      // if (run <= 5367 && run >= 5320) continue;
      // ext bnb in run 1, high rate
      if (run>=7004 && run <=7070) continue;
      // ext bnb in run 2, high rate not in good list anyway
      //      if ((run>=10287 && run <= 10304) || (run>=12277 && run <=12350)) continue;
      // ext bnb in run 2, low rate not in good list anyway
      // if ((run>=9768 && run <= 10070) || (run>=10102 && run <=10246)) continue;
      // bnb run 2 high rate
      if (run >= 8321 && run <=8404) continue;
      // bnb run 3 high rate
      if (run >=15369 && run <= 15402) continue;
    }

    t1->Fill();
    t3->Fill();
    t4->Fill();
    t5->Fill();

    T_spacepoints->GetEntry(i);
    new_T_spacepoints->Fill();

    for(auto tree_it=old_trees->begin(); tree_it!=old_trees->end(); tree_it++){
        (*tree_it)->GetEntry(i);
    }

    for(auto tree_it=new_trees->begin(); tree_it!=new_trees->end(); tree_it++){
        (*tree_it)->Fill();
    }

  }

  for (Int_t i=0;i!=T_pot->GetEntries();i++){
    T_pot->GetEntry(i);

    if (flag_data){
      if (good_runlist_set.find(runNo) == good_runlist_set.end()) continue;
      if (low_lifetime_set.find(runNo) != low_lifetime_set.end()) continue;
      if (flag_numi && low_neutrino_count_numi_run2RHC_set.find(runNo) != low_neutrino_count_numi_run2RHC_set.end()) continue;
      //bad run in run 1 due to beam filter bnb
      //if (runNo <= 5367 && runNo >= 5320) continue;
      // ext bnb in run 1, high rate
      if (runNo >=7004 && runNo <=7070) continue;
      // ext bnb in run 2, high rate not in good list anyway
      //if ((runNo>=10287 && runNo <= 10304) || (runNo>=12277 && runNo <=12350)) continue;
      // ext bnb in run 2, low rate not in good list anyway
      //if ((runNo>=9768 && runNo <= 10070) || (runNo>=10102 && runNo <=10246)) continue;
      // bnb run 2 high rate
      if (runNo >= 8321 && runNo <=8404) continue;
      // bnb run 3 high rate
      if (runNo >=15369 && runNo <= 15402) continue;
    }
    t2->Fill();
    for(auto tree_it=old_trees_pot->begin(); tree_it!=old_trees_pot->end(); tree_it++){
      (*tree_it)->GetEntry(i);
    }

    for(auto tree_it=new_trees_pot->begin(); tree_it!=new_trees_pot->end(); tree_it++){
      (*tree_it)->Fill();
    }
  }

  file2->Write("",TFile::kOverwrite);
  file2->Close();

  return 0;
}
