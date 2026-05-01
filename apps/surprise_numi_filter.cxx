// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
#include <iomanip>
#include <iostream>
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
#include "WCPLEEANA/cuts.h"

#include "WCPLEEANA/tree_wrangler.h"


int main( int argc, char** argv )
{
  if(argc==2 && argv[1][1]=='h'){
    std::cout<<"TODO"<<std::endl;
    return 0;
  }
  else if(argc==2 && argv[1][1]=='H'){
    print_help_wrangler_config(true);
    return 0;
  }
  else if (argc < 5) {
    std::cout << "numi_filter #input_file #prefix_outfile #config.txt #filter_level" << std::endl; 
    return -1;
  }

  TString input_file = argv[1];
  TString prefix_out = argv[2];

  std::string config_file_name=argv[3];
  bool flag_config = true;

  Int_t filter_level = atoi(argv[4]);

  int flag_overwrite = 0;

  char delimiter = ',';

  for (Int_t i = 3; i < argc; ++i) {

    // Skip non-flags
    if (argv[i][0] != '-') continue;

    char flag = argv[i][1];
    char* value_ptr = nullptr;

    // Case 1: attached value (-xVALUE)
    if (argv[i][2] != '\0') {
      value_ptr = &argv[i][2];
    }
    // Case 2: separate value (-x VALUE)
    else if (i + 1 < argc && argv[i+1][0] != '-') {
      value_ptr = argv[i + 1];
      ++i; // consume next argument
    }

    // Guard against missing values
    if (!value_ptr) {
      std::cerr << "Missing value for -" << flag << std::endl;
      continue;
    }

    switch(flag){

    case 'o':
      flag_overwrite = atoi(value_ptr);
      break;

    case 'd':
      delimiter = value_ptr[0];
      break;

    }

  }

  std::string  outfile_name;

  if (filter_level == 1){
    std::cout << "FHC mode" << std::endl;
    outfile_name = prefix_out + "_FHC.root";
  }else{
    std::cout << "RHC mode" << std::endl;
    outfile_name = prefix_out + "_RHC.root";
  }


  // Check if the output file exists if overwrite is not set.
  if(flag_overwrite!=1){
    TFile *f = nullptr;
    if (!gSystem->AccessPathName(outfile_name.c_str())) f = TFile::Open(outfile_name.c_str(), "READ");
    if (f && !f->IsZombie()) {
        std::cout<<'\n'<< "File exists. Exiting." << std::endl;
        f->Close();
        std::cout<<"Outputfile file "<<outfile_name<<" already exists and overwrite not set."<<std::endl;
        std::cout<<"Pick a new file name or use -o to force and overwrite of existing file."<<'\n'<<std::endl;
        return 1;
    }
  }


  TFile *file1 = new TFile(input_file);


  // Initiate the tree wranglers
  tree_wrangler wrangler(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_ex(flag_config, config_file_name, delimiter,2);
  tree_wrangler wrangler_pot(flag_config, config_file_name, delimiter,1);
  

  //Load other trees from directories as specified by the config file
  wrangler.get_old_trees(file1);
  wrangler_ex.get_old_trees(file1);
  wrangler_pot.get_old_trees(file1);


  // Figure out which tree you can load for the rse map
  int run;
  int subrun;
  int event;
  TTree *T_rse = nullptr;
  int found_rse_tree = get_T_rse(file1, T_rse, run, subrun, event);
  if(!found_rse_tree){
    std::cout<<'\n'<<"Could not find RSE tree. Exiting."<<std::endl;
    return 1;
  }


  TFile *file2 = new TFile(outfile_name.c_str(),"RECREATE");


  //Setup the directories specified in the config file
  wrangler.set_new_trees(file2);
  wrangler_ex.set_new_trees(file2);
  wrangler_pot.set_new_trees(file2);


  // Build the pairs of pot trees
  wrangler_pot.grow_pot_arboretum();


  // Begin filling the event-level trees

  int nentries = T_rse->GetEntries();
  std::cout<<"Begin looping over "<<nentries<<" events"<<std::endl; 
  for (int i=0;i!=nentries;i++){

    if (i%10000 == 0) std::cout << i/1000 << " k " << std::setprecision(3) << double(i)/nentries*100. << " %"<< std::endl;

    T_rse->GetEntry(i);

    if (filter_level==1 && (run > 6748 && run <=7001 || // RHC
			    run >=10140 && run <= 11949 ||
			    run >= 13697 && run <=17566 ||
			    run >= 19668 && run <=21410 
			    ) ) continue;
    if (filter_level!=1 && (run <=6748 || // FHC
			    run >=8784 && run <=10139 ||
			    run >= 21411 && run <= 23259 ||
			    run >= 24256 && run <= 25763
			    )) continue;
    
    for(auto tree_it=wrangler.old_trees->begin(); tree_it!=wrangler.old_trees->end(); tree_it++){
        (*tree_it)->GetEntry(i);
    }
    for(auto tree_it=wrangler.new_trees->begin(); tree_it!=wrangler.new_trees->end(); tree_it++){
        (*tree_it)->Fill();
    }
    for(auto tree_it=wrangler_ex.old_trees->begin(); tree_it!=wrangler_ex.old_trees->end(); tree_it++){
        (*tree_it)->GetEntry(i);
    }
    for(auto tree_it=wrangler_ex.new_trees->begin(); tree_it!=wrangler_ex.new_trees->end(); tree_it++){
        (*tree_it)->Fill();
    }

  }


  // Loop over each POT tree seperatly
  for(auto pot_tree_it=wrangler_pot.pot_arboretum->begin(); pot_tree_it!=wrangler_pot.pot_arboretum->end(); pot_tree_it++){

    std::cout<<"Begin looping over "<<(*pot_tree_it)->old_pot_tree->GetName()<<" tree with "<<nentries<<" entries"<<std::endl;
    nentries = (*pot_tree_it)->old_pot_tree->GetEntries();
    for (Int_t i=0;i!=nentries;i++){

      if (i%10000 == 0) std::cout << i/1000 << " k " << std::setprecision(3) << double(i)/nentries*100. << " %"<< std::endl;

      (*pot_tree_it)->old_pot_tree->GetEntry(i);

    if (filter_level==1 && ((*pot_tree_it)->runNo > 6748 && (*pot_tree_it)->runNo <=7001 || // RHC
                            (*pot_tree_it)->runNo >=10140 && (*pot_tree_it)->runNo <= 11949 ||
                            (*pot_tree_it)->runNo >= 13697 && (*pot_tree_it)->runNo <=17566 ||
                            (*pot_tree_it)->runNo >= 19668 && (*pot_tree_it)->runNo <=21410
                            ) ) continue;
    if (filter_level!=1 && ((*pot_tree_it)->runNo <=6748 || // FHC
                            (*pot_tree_it)->runNo >=8784 && (*pot_tree_it)->runNo <=10139 ||
                            (*pot_tree_it)->runNo >= 21411 && (*pot_tree_it)->runNo <= 23259 ||
                            (*pot_tree_it)->runNo >= 24256 && (*pot_tree_it)->runNo <= 25763
                            )) continue;

      (*pot_tree_it)->new_pot_tree->Fill();

    }//i, loop over events in a given pot tree set

  }//pot_trees_it, loop over sets of pot trees

  
  file2->Write("",TFile::kOverwrite);
  file2->Close();


  return 0;

  
  
}
