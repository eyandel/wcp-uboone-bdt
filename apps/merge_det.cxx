// cz: code modified from tutorials/tmva/TMVAClassification.C

#include <cstdlib>
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

#include "WCPLEEANA/eval.h"

#include "WCPLEEANA/tree_wrangler.h"

using namespace std;
using namespace LEEana;

//#include "WCPLEEANA/bdt.h"
#include "WCPLEEANA/pot.h"
#include "WCPLEEANA/pfeval.h"
#include "WCPLEEANA/kine.h"

int main( int argc, char** argv )
{
  if (argc < 4) {
    std::cout << "merge_det #input_file_cv #input_file_det #output_file " << std::endl;
    return -1;
  }

  TString input_file_cv = argv[1];
  TString input_file_det = argv[2];
  TString out_file = argv[3];
  bool flag_config = false;
  std::string config_file_name="config.txt";
  char delimiter = ',';
  bool override_pot = false;
  for (Int_t i=4;i!=argc;i++){
    switch(argv[i][1]){
    case 't':
       config_file_name = &argv[i][2];
       flag_config = true;
      break;
    case 'd':
        delimiter = argv[i][2];//In case you want to change what character you use to sperate your trees in the config
      break;
    case 'p':
        override_pot = &argv[i][2];
        if(override_pot) std::cout<<"Overriding the pot of each subrun"<<std::endl;
      break;
    }
  }

  tree_wrangler wrangler_cv(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_det(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_pot_cv(flag_config, config_file_name, delimiter,true);
  tree_wrangler wrangler_pot_det(flag_config, config_file_name, delimiter,true);

  TFile *file1 = new TFile(input_file_cv);
  TTree *T_BDTvars_cv = (TTree*)file1->Get("wcpselection/T_BDTvars");
  TTree *T_eval_cv = (TTree*)file1->Get("wcpselection/T_eval");
  TTree *T_pot_cv = (TTree*)file1->Get("wcpselection/T_pot");
  TTree *T_PFeval_cv = (TTree*)file1->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars_cv = (TTree*)file1->Get("wcpselection/T_KINEvars");
  TTree *T_spacepoints_cv = (TTree*)file1->Get("wcpselection/T_spacepoints");

  TFile *file2 = new TFile(input_file_det);
  TTree *T_BDTvars_det = (TTree*)file2->Get("wcpselection/T_BDTvars");
  TTree *T_eval_det = (TTree*)file2->Get("wcpselection/T_eval");
  TTree *T_pot_det = (TTree*)file2->Get("wcpselection/T_pot");
  TTree *T_PFeval_det = (TTree*)file2->Get("wcpselection/T_PFeval");
  TTree *T_KINEvars_det = (TTree*)file2->Get("wcpselection/T_KINEvars");
  TTree *T_spacepoints_det = (TTree*)file2->Get("wcpselection/T_spacepoints");

  //Load other trees from directories as specified by the config file
  std::vector<TTree*>* old_trees_cv = new std::vector<TTree*>;
  old_trees_cv = wrangler_cv.get_old_trees(file1);
  std::vector<TTree*>* old_trees_det = new std::vector<TTree*>;
  old_trees_det = wrangler_det.get_old_trees(file2);
  std::vector<TTree*>* old_trees_pot_cv = new std::vector<TTree*>;
  old_trees_pot_cv = wrangler_pot_cv.get_old_trees(file1);
  std::vector<TTree*>* old_trees_pot_det = new std::vector<TTree*>;
  old_trees_pot_det = wrangler_pot_det.get_old_trees(file2);

  TFile *file3 = new TFile(out_file,"RECREATE");

  //Setup the directories specified in the config file
  std::vector<TTree*>* new_trees_cv = new std::vector<TTree*>;
  new_trees_cv = wrangler_cv.set_new_trees(file3,true,"_cv");
  std::vector<TTree*>* new_trees_det = new std::vector<TTree*>;
  new_trees_det = wrangler_det.set_new_trees(file3,true,"_det");
  std::vector<TTree*>* new_trees_pot_cv = new std::vector<TTree*>;
  new_trees_pot_cv = wrangler_pot_cv.set_new_trees(file3,true,"_cv");
  std::vector<TTree*>* new_trees_pot_det = new std::vector<TTree*>;
  new_trees_pot_det = wrangler_pot_det.set_new_trees(file3,true,"_det");

  //Always do WC
  file3->mkdir("wcpselection");
  file3->cd("wcpselection");
  TTree *t1_cv = new TTree("T_eval_cv","T_eval_cv");
  TTree *t2_cv = new TTree("T_pot_cv","T_pot_cv");
  TTree *t3_cv = new TTree("T_PFeval_cv", "T_PFeval_cv");
  TTree *t5_cv = new TTree("T_KINEvars_cv", "T_KINEvars_cv");
  TTree *t4_cv = new TTree("T_BDTvars_cv","T_BDTvars_cv");
  TTree *T_spacepoints_new_cv = T_spacepoints_cv->CloneTree(0);
  T_spacepoints_new_cv->SetObject("T_spacepoints_cv","T_spacepoints_cv");

  TTree *t1_det = new TTree("T_eval_det","T_eval_det");
  TTree *t2_det = new TTree("T_pot_det","T_pot_det");
  TTree *t3_det = new TTree("T_PFeval_det", "T_PFeval_det");
  TTree *t5_det = new TTree("T_KINEvars_det", "T_KINEvars_det");
  TTree *t4_det = new TTree("T_BDTvars_det","T_BDTvars_det");
  TTree *T_spacepoints_new_det = T_spacepoints_det->CloneTree(0);
  T_spacepoints_new_det->SetObject("T_spacepoints_det","T_spacepoints_det");

  EvalInfo eval_cv;
  eval_cv.file_type = new std::string();
  POTInfo pot_cv;
  TaggerInfo tagger_cv;
  PFevalInfo pfeval_cv;
  KineInfo kine_cv;
  kine_cv.kine_energy_particle = new std::vector<float>;
  kine_cv.kine_energy_info = new std::vector<int>;
  kine_cv.kine_particle_type = new std::vector<int>;
  kine_cv.kine_energy_included = new std::vector<int>;


  tagger_cv.pio_2_v_dis2 = new std::vector<float>;
  tagger_cv.pio_2_v_angle2 = new std::vector<float>;
  tagger_cv.pio_2_v_acc_length = new std::vector<float>;
  tagger_cv.pio_2_v_flag = new std::vector<float>;
  tagger_cv.sig_1_v_angle = new std::vector<float>;
  tagger_cv.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_1_v_energy = new std::vector<float>;
  tagger_cv.sig_1_v_energy_1 = new std::vector<float>;
  tagger_cv.sig_1_v_flag = new std::vector<float>;
  tagger_cv.sig_2_v_energy = new std::vector<float>;
  tagger_cv.sig_2_v_shower_angle = new std::vector<float>;
  tagger_cv.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_flag = new std::vector<float>;
  tagger_cv.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_energy = new std::vector<float>;
  tagger_cv.stw_2_v_angle = new std::vector<float>;
  tagger_cv.stw_2_v_dir_length = new std::vector<float>;
  tagger_cv.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_flag = new std::vector<float>;
  tagger_cv.stw_3_v_angle = new std::vector<float>;
  tagger_cv.stw_3_v_dir_length = new std::vector<float>;
  tagger_cv.stw_3_v_energy = new std::vector<float>;
  tagger_cv.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_3_v_flag = new std::vector<float>;
  tagger_cv.stw_4_v_angle = new std::vector<float>;
  tagger_cv.stw_4_v_dis = new std::vector<float>;
  tagger_cv.stw_4_v_energy = new std::vector<float>;
  tagger_cv.stw_4_v_flag = new std::vector<float>;
  tagger_cv.br3_3_v_energy = new std::vector<float>;
  tagger_cv.br3_3_v_angle = new std::vector<float>;
  tagger_cv.br3_3_v_dir_length = new std::vector<float>;
  tagger_cv.br3_3_v_length = new std::vector<float>;
  tagger_cv.br3_3_v_flag = new std::vector<float>;
  tagger_cv.br3_5_v_dir_length = new std::vector<float>;
  tagger_cv.br3_5_v_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_cv.br3_5_v_n_seg = new std::vector<float>;
  tagger_cv.br3_5_v_angle = new std::vector<float>;
  tagger_cv.br3_5_v_sg_length = new std::vector<float>;
  tagger_cv.br3_5_v_energy = new std::vector<float>;
  tagger_cv.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_cv.br3_5_v_n_segs = new std::vector<float>;
  tagger_cv.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_cv.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag = new std::vector<float>;
  tagger_cv.br3_6_v_angle = new std::vector<float>;
  tagger_cv.br3_6_v_angle1 = new std::vector<float>;
  tagger_cv.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.br3_6_v_direct_length = new std::vector<float>;
  tagger_cv.br3_6_v_length = new std::vector<float>;
  tagger_cv.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_cv.br3_6_v_energy = new std::vector<float>;
  tagger_cv.br3_6_v_flag = new std::vector<float>;
  tagger_cv.tro_1_v_particle_type = new std::vector<float>;
  tagger_cv.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.tro_1_v_min_dis = new std::vector<float>;
  tagger_cv.tro_1_v_sg1_length = new std::vector<float>;
  tagger_cv.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_1_v_tmp_length = new std::vector<float>;
  tagger_cv.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_cv.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_cv.tro_1_v_flag = new std::vector<float>;
  tagger_cv.tro_2_v_energy = new std::vector<float>;
  tagger_cv.tro_2_v_stem_length = new std::vector<float>;
  tagger_cv.tro_2_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_2_v_max_length = new std::vector<float>;
  tagger_cv.tro_2_v_angle = new std::vector<float>;
  tagger_cv.tro_2_v_flag = new std::vector<float>;
  tagger_cv.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_cv.tro_4_v_angle = new std::vector<float>;
  tagger_cv.tro_4_v_angle1 = new std::vector<float>;
  tagger_cv.tro_4_v_angle2 = new std::vector<float>;
  tagger_cv.tro_4_v_length = new std::vector<float>;
  tagger_cv.tro_4_v_length1 = new std::vector<float>;
  tagger_cv.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_energy = new std::vector<float>;
  tagger_cv.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.tro_4_v_flag = new std::vector<float>;
  tagger_cv.tro_5_v_max_angle = new std::vector<float>;
  tagger_cv.tro_5_v_min_angle = new std::vector<float>;
  tagger_cv.tro_5_v_max_length = new std::vector<float>;
  tagger_cv.tro_5_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_5_v_min_count = new std::vector<float>;
  tagger_cv.tro_5_v_max_count = new std::vector<float>;
  tagger_cv.tro_5_v_energy = new std::vector<float>;
  tagger_cv.tro_5_v_flag = new std::vector<float>;
  tagger_cv.lol_1_v_energy = new std::vector<float>;
  tagger_cv.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_1_v_nseg = new std::vector<float>;
  tagger_cv.lol_1_v_angle = new std::vector<float>;
  tagger_cv.lol_1_v_flag = new std::vector<float>;
  tagger_cv.lol_2_v_length = new std::vector<float>;
  tagger_cv.lol_2_v_angle = new std::vector<float>;
  tagger_cv.lol_2_v_type = new std::vector<float>;
  tagger_cv.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_2_v_energy = new std::vector<float>;
  tagger_cv.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_cv.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.lol_2_v_flag = new std::vector<float>;
  tagger_cv.cosmict_flag_10 = new std::vector<float>;
  tagger_cv.cosmict_10_flag_inside = new std::vector<float>;
  tagger_cv.cosmict_10_vtx_z = new std::vector<float>;
  tagger_cv.cosmict_10_flag_shower = new std::vector<float>;
  tagger_cv.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_cv.cosmict_10_angle_beam = new std::vector<float>;
  tagger_cv.cosmict_10_length = new std::vector<float>;
  tagger_cv.numu_cc_flag_1 = new std::vector<float>;
  tagger_cv.numu_cc_1_particle_type = new std::vector<float>;
  tagger_cv.numu_cc_1_length = new std::vector<float>;
  tagger_cv.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_cv.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_cv.numu_cc_1_direct_length = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_cv.numu_cc_flag_2 = new std::vector<float>;
  tagger_cv.numu_cc_2_length = new std::vector<float>;
  tagger_cv.numu_cc_2_total_length = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_cv.pio_2_v_dis2 = new std::vector<float>;
  tagger_cv.pio_2_v_angle2 = new std::vector<float>;
  tagger_cv.pio_2_v_acc_length = new std::vector<float>;
  tagger_cv.pio_2_v_flag = new std::vector<float>;
  tagger_cv.sig_1_v_angle = new std::vector<float>;
  tagger_cv.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_1_v_energy = new std::vector<float>;
  tagger_cv.sig_1_v_energy_1 = new std::vector<float>;
  tagger_cv.sig_1_v_flag = new std::vector<float>;
  tagger_cv.sig_2_v_energy = new std::vector<float>;
  tagger_cv.sig_2_v_shower_angle = new std::vector<float>;
  tagger_cv.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_cv.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_cv.sig_2_v_flag = new std::vector<float>;
  tagger_cv.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_energy = new std::vector<float>;
  tagger_cv.stw_2_v_angle = new std::vector<float>;
  tagger_cv.stw_2_v_dir_length = new std::vector<float>;
  tagger_cv.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_cv.stw_2_v_flag = new std::vector<float>;
  tagger_cv.stw_3_v_angle = new std::vector<float>;
  tagger_cv.stw_3_v_dir_length = new std::vector<float>;
  tagger_cv.stw_3_v_energy = new std::vector<float>;
  tagger_cv.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.stw_3_v_flag = new std::vector<float>;
  tagger_cv.stw_4_v_angle = new std::vector<float>;
  tagger_cv.stw_4_v_dis = new std::vector<float>;
  tagger_cv.stw_4_v_energy = new std::vector<float>;
  tagger_cv.stw_4_v_flag = new std::vector<float>;
  tagger_cv.br3_3_v_energy = new std::vector<float>;
  tagger_cv.br3_3_v_angle = new std::vector<float>;
  tagger_cv.br3_3_v_dir_length = new std::vector<float>;
  tagger_cv.br3_3_v_length = new std::vector<float>;
  tagger_cv.br3_3_v_flag = new std::vector<float>;
  tagger_cv.br3_5_v_dir_length = new std::vector<float>;
  tagger_cv.br3_5_v_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_cv.br3_5_v_n_seg = new std::vector<float>;
  tagger_cv.br3_5_v_angle = new std::vector<float>;
  tagger_cv.br3_5_v_sg_length = new std::vector<float>;
  tagger_cv.br3_5_v_energy = new std::vector<float>;
  tagger_cv.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_cv.br3_5_v_n_segs = new std::vector<float>;
  tagger_cv.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_cv.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_cv.br3_5_v_flag = new std::vector<float>;
  tagger_cv.br3_6_v_angle = new std::vector<float>;
  tagger_cv.br3_6_v_angle1 = new std::vector<float>;
  tagger_cv.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.br3_6_v_direct_length = new std::vector<float>;
  tagger_cv.br3_6_v_length = new std::vector<float>;
  tagger_cv.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_cv.br3_6_v_energy = new std::vector<float>;
  tagger_cv.br3_6_v_flag = new std::vector<float>;
  tagger_cv.tro_1_v_particle_type = new std::vector<float>;
  tagger_cv.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.tro_1_v_min_dis = new std::vector<float>;
  tagger_cv.tro_1_v_sg1_length = new std::vector<float>;
  tagger_cv.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_1_v_tmp_length = new std::vector<float>;
  tagger_cv.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_cv.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_cv.tro_1_v_flag = new std::vector<float>;
  tagger_cv.tro_2_v_energy = new std::vector<float>;
  tagger_cv.tro_2_v_stem_length = new std::vector<float>;
  tagger_cv.tro_2_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_2_v_max_length = new std::vector<float>;
  tagger_cv.tro_2_v_angle = new std::vector<float>;
  tagger_cv.tro_2_v_flag = new std::vector<float>;
  tagger_cv.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_cv.tro_4_v_angle = new std::vector<float>;
  tagger_cv.tro_4_v_angle1 = new std::vector<float>;
  tagger_cv.tro_4_v_angle2 = new std::vector<float>;
  tagger_cv.tro_4_v_length = new std::vector<float>;
  tagger_cv.tro_4_v_length1 = new std::vector<float>;
  tagger_cv.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_cv.tro_4_v_energy = new std::vector<float>;
  tagger_cv.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_cv.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_cv.tro_4_v_flag = new std::vector<float>;
  tagger_cv.tro_5_v_max_angle = new std::vector<float>;
  tagger_cv.tro_5_v_min_angle = new std::vector<float>;
  tagger_cv.tro_5_v_max_length = new std::vector<float>;
  tagger_cv.tro_5_v_iso_angle = new std::vector<float>;
  tagger_cv.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_cv.tro_5_v_min_count = new std::vector<float>;
  tagger_cv.tro_5_v_max_count = new std::vector<float>;
  tagger_cv.tro_5_v_energy = new std::vector<float>;
  tagger_cv.tro_5_v_flag = new std::vector<float>;
  tagger_cv.lol_1_v_energy = new std::vector<float>;
  tagger_cv.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_1_v_nseg = new std::vector<float>;
  tagger_cv.lol_1_v_angle = new std::vector<float>;
  tagger_cv.lol_1_v_flag = new std::vector<float>;
  tagger_cv.lol_2_v_length = new std::vector<float>;
  tagger_cv.lol_2_v_angle = new std::vector<float>;
  tagger_cv.lol_2_v_type = new std::vector<float>;
  tagger_cv.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_cv.lol_2_v_energy = new std::vector<float>;
  tagger_cv.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_cv.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_cv.lol_2_v_flag = new std::vector<float>;
  tagger_cv.cosmict_flag_10 = new std::vector<float>;
  tagger_cv.cosmict_10_flag_inside = new std::vector<float>;
  tagger_cv.cosmict_10_vtx_z = new std::vector<float>;
  tagger_cv.cosmict_10_flag_shower = new std::vector<float>;
  tagger_cv.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_cv.cosmict_10_angle_beam = new std::vector<float>;
  tagger_cv.cosmict_10_length = new std::vector<float>;
  tagger_cv.numu_cc_flag_1 = new std::vector<float>;
  tagger_cv.numu_cc_1_particle_type = new std::vector<float>;
  tagger_cv.numu_cc_1_length = new std::vector<float>;
  tagger_cv.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_cv.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_cv.numu_cc_1_direct_length = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_cv.numu_cc_flag_2 = new std::vector<float>;
  tagger_cv.numu_cc_2_length = new std::vector<float>;
  tagger_cv.numu_cc_2_total_length = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_cv.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_cv.ssm_kine_energy_particle = new std::vector<float>;
  tagger_cv.ssm_kine_energy_info = new std::vector<int>;
  tagger_cv.ssm_kine_particle_type = new std::vector<int>;
  tagger_cv.ssm_kine_energy_included = new std::vector<int>;
  tagger_cv.ssm_cosmict_flag_10 = new std::vector<float>;
  tagger_cv.WCPMTInfoPePred = new std::vector<double>;
  tagger_cv.WCPMTInfoPeMeas = new std::vector<double>;
  tagger_cv.WCPMTInfoPeMeasErr = new std::vector<double>;

  EvalInfo eval_det;
  eval_det.file_type = new std::string();
  POTInfo pot_det;
  TaggerInfo tagger_det;
  PFevalInfo pfeval_det;
  KineInfo kine_det;
  kine_det.kine_energy_particle = new std::vector<float>;
  kine_det.kine_energy_info = new std::vector<int>;
  kine_det.kine_particle_type = new std::vector<int>;
  kine_det.kine_energy_included = new std::vector<int>;


  tagger_det.pio_2_v_dis2 = new std::vector<float>;
  tagger_det.pio_2_v_angle2 = new std::vector<float>;
  tagger_det.pio_2_v_acc_length = new std::vector<float>;
  tagger_det.pio_2_v_flag = new std::vector<float>;
  tagger_det.sig_1_v_angle = new std::vector<float>;
  tagger_det.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_1_v_energy = new std::vector<float>;
  tagger_det.sig_1_v_energy_1 = new std::vector<float>;
  tagger_det.sig_1_v_flag = new std::vector<float>;
  tagger_det.sig_2_v_energy = new std::vector<float>;
  tagger_det.sig_2_v_shower_angle = new std::vector<float>;
  tagger_det.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_flag = new std::vector<float>;
  tagger_det.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_energy = new std::vector<float>;
  tagger_det.stw_2_v_angle = new std::vector<float>;
  tagger_det.stw_2_v_dir_length = new std::vector<float>;
  tagger_det.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_flag = new std::vector<float>;
  tagger_det.stw_3_v_angle = new std::vector<float>;
  tagger_det.stw_3_v_dir_length = new std::vector<float>;
  tagger_det.stw_3_v_energy = new std::vector<float>;
  tagger_det.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_3_v_flag = new std::vector<float>;
  tagger_det.stw_4_v_angle = new std::vector<float>;
  tagger_det.stw_4_v_dis = new std::vector<float>;
  tagger_det.stw_4_v_energy = new std::vector<float>;
  tagger_det.stw_4_v_flag = new std::vector<float>;
  tagger_det.br3_3_v_energy = new std::vector<float>;
  tagger_det.br3_3_v_angle = new std::vector<float>;
  tagger_det.br3_3_v_dir_length = new std::vector<float>;
  tagger_det.br3_3_v_length = new std::vector<float>;
  tagger_det.br3_3_v_flag = new std::vector<float>;
  tagger_det.br3_5_v_dir_length = new std::vector<float>;
  tagger_det.br3_5_v_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_det.br3_5_v_n_seg = new std::vector<float>;
  tagger_det.br3_5_v_angle = new std::vector<float>;
  tagger_det.br3_5_v_sg_length = new std::vector<float>;
  tagger_det.br3_5_v_energy = new std::vector<float>;
  tagger_det.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_det.br3_5_v_n_segs = new std::vector<float>;
  tagger_det.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_det.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag = new std::vector<float>;
  tagger_det.br3_6_v_angle = new std::vector<float>;
  tagger_det.br3_6_v_angle1 = new std::vector<float>;
  tagger_det.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.br3_6_v_direct_length = new std::vector<float>;
  tagger_det.br3_6_v_length = new std::vector<float>;
  tagger_det.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_det.br3_6_v_energy = new std::vector<float>;
  tagger_det.br3_6_v_flag = new std::vector<float>;
  tagger_det.tro_1_v_particle_type = new std::vector<float>;
  tagger_det.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_det.tro_1_v_min_dis = new std::vector<float>;
  tagger_det.tro_1_v_sg1_length = new std::vector<float>;
  tagger_det.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_1_v_tmp_length = new std::vector<float>;
  tagger_det.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_det.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_det.tro_1_v_flag = new std::vector<float>;
  tagger_det.tro_2_v_energy = new std::vector<float>;
  tagger_det.tro_2_v_stem_length = new std::vector<float>;
  tagger_det.tro_2_v_iso_angle = new std::vector<float>;
  tagger_det.tro_2_v_max_length = new std::vector<float>;
  tagger_det.tro_2_v_angle = new std::vector<float>;
  tagger_det.tro_2_v_flag = new std::vector<float>;
  tagger_det.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_det.tro_4_v_angle = new std::vector<float>;
  tagger_det.tro_4_v_angle1 = new std::vector<float>;
  tagger_det.tro_4_v_angle2 = new std::vector<float>;
  tagger_det.tro_4_v_length = new std::vector<float>;
  tagger_det.tro_4_v_length1 = new std::vector<float>;
  tagger_det.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_energy = new std::vector<float>;
  tagger_det.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.tro_4_v_flag = new std::vector<float>;
  tagger_det.tro_5_v_max_angle = new std::vector<float>;
  tagger_det.tro_5_v_min_angle = new std::vector<float>;
  tagger_det.tro_5_v_max_length = new std::vector<float>;
  tagger_det.tro_5_v_iso_angle = new std::vector<float>;
  tagger_det.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_5_v_min_count = new std::vector<float>;
  tagger_det.tro_5_v_max_count = new std::vector<float>;
  tagger_det.tro_5_v_energy = new std::vector<float>;
  tagger_det.tro_5_v_flag = new std::vector<float>;
  tagger_det.lol_1_v_energy = new std::vector<float>;
  tagger_det.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_1_v_nseg = new std::vector<float>;
  tagger_det.lol_1_v_angle = new std::vector<float>;
  tagger_det.lol_1_v_flag = new std::vector<float>;
  tagger_det.lol_2_v_length = new std::vector<float>;
  tagger_det.lol_2_v_angle = new std::vector<float>;
  tagger_det.lol_2_v_type = new std::vector<float>;
  tagger_det.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_2_v_energy = new std::vector<float>;
  tagger_det.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_det.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_det.lol_2_v_flag = new std::vector<float>;
  tagger_det.cosmict_flag_10 = new std::vector<float>;
  tagger_det.cosmict_10_flag_inside = new std::vector<float>;
  tagger_det.cosmict_10_vtx_z = new std::vector<float>;
  tagger_det.cosmict_10_flag_shower = new std::vector<float>;
  tagger_det.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_det.cosmict_10_angle_beam = new std::vector<float>;
  tagger_det.cosmict_10_length = new std::vector<float>;
  tagger_det.numu_cc_flag_1 = new std::vector<float>;
  tagger_det.numu_cc_1_particle_type = new std::vector<float>;
  tagger_det.numu_cc_1_length = new std::vector<float>;
  tagger_det.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_det.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_det.numu_cc_1_direct_length = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_det.numu_cc_flag_2 = new std::vector<float>;
  tagger_det.numu_cc_2_length = new std::vector<float>;
  tagger_det.numu_cc_2_total_length = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_det.pio_2_v_dis2 = new std::vector<float>;
  tagger_det.pio_2_v_angle2 = new std::vector<float>;
  tagger_det.pio_2_v_acc_length = new std::vector<float>;
  tagger_det.pio_2_v_flag = new std::vector<float>;
  tagger_det.sig_1_v_angle = new std::vector<float>;
  tagger_det.sig_1_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_1_v_energy = new std::vector<float>;
  tagger_det.sig_1_v_energy_1 = new std::vector<float>;
  tagger_det.sig_1_v_flag = new std::vector<float>;
  tagger_det.sig_2_v_energy = new std::vector<float>;
  tagger_det.sig_2_v_shower_angle = new std::vector<float>;
  tagger_det.sig_2_v_flag_single_shower = new std::vector<float>;
  tagger_det.sig_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_start_dQ_dx = new std::vector<float>;
  tagger_det.sig_2_v_flag = new std::vector<float>;
  tagger_det.stw_2_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_energy = new std::vector<float>;
  tagger_det.stw_2_v_angle = new std::vector<float>;
  tagger_det.stw_2_v_dir_length = new std::vector<float>;
  tagger_det.stw_2_v_max_dQ_dx = new std::vector<float>;
  tagger_det.stw_2_v_flag = new std::vector<float>;
  tagger_det.stw_3_v_angle = new std::vector<float>;
  tagger_det.stw_3_v_dir_length = new std::vector<float>;
  tagger_det.stw_3_v_energy = new std::vector<float>;
  tagger_det.stw_3_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.stw_3_v_flag = new std::vector<float>;
  tagger_det.stw_4_v_angle = new std::vector<float>;
  tagger_det.stw_4_v_dis = new std::vector<float>;
  tagger_det.stw_4_v_energy = new std::vector<float>;
  tagger_det.stw_4_v_flag = new std::vector<float>;
  tagger_det.br3_3_v_energy = new std::vector<float>;
  tagger_det.br3_3_v_angle = new std::vector<float>;
  tagger_det.br3_3_v_dir_length = new std::vector<float>;
  tagger_det.br3_3_v_length = new std::vector<float>;
  tagger_det.br3_3_v_flag = new std::vector<float>;
  tagger_det.br3_5_v_dir_length = new std::vector<float>;
  tagger_det.br3_5_v_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag_avoid_muon_check = new std::vector<float>;
  tagger_det.br3_5_v_n_seg = new std::vector<float>;
  tagger_det.br3_5_v_angle = new std::vector<float>;
  tagger_det.br3_5_v_sg_length = new std::vector<float>;
  tagger_det.br3_5_v_energy = new std::vector<float>;
  tagger_det.br3_5_v_n_main_segs = new std::vector<float>;
  tagger_det.br3_5_v_n_segs = new std::vector<float>;
  tagger_det.br3_5_v_shower_main_length = new std::vector<float>;
  tagger_det.br3_5_v_shower_total_length = new std::vector<float>;
  tagger_det.br3_5_v_flag = new std::vector<float>;
  tagger_det.br3_6_v_angle = new std::vector<float>;
  tagger_det.br3_6_v_angle1 = new std::vector<float>;
  tagger_det.br3_6_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.br3_6_v_direct_length = new std::vector<float>;
  tagger_det.br3_6_v_length = new std::vector<float>;
  tagger_det.br3_6_v_n_other_vtx_segs = new std::vector<float>;
  tagger_det.br3_6_v_energy = new std::vector<float>;
  tagger_det.br3_6_v_flag = new std::vector<float>;
  tagger_det.tro_1_v_particle_type = new std::vector<float>;
  tagger_det.tro_1_v_flag_dir_weak = new std::vector<float>;
  tagger_det.tro_1_v_min_dis = new std::vector<float>;
  tagger_det.tro_1_v_sg1_length = new std::vector<float>;
  tagger_det.tro_1_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_1_v_max_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_1_v_tmp_length = new std::vector<float>;
  tagger_det.tro_1_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_1_v_dQ_dx_cut = new std::vector<float>;
  tagger_det.tro_1_v_flag_shower_topology = new std::vector<float>;
  tagger_det.tro_1_v_flag = new std::vector<float>;
  tagger_det.tro_2_v_energy = new std::vector<float>;
  tagger_det.tro_2_v_stem_length = new std::vector<float>;
  tagger_det.tro_2_v_iso_angle = new std::vector<float>;
  tagger_det.tro_2_v_max_length = new std::vector<float>;
  tagger_det.tro_2_v_angle = new std::vector<float>;
  tagger_det.tro_2_v_flag = new std::vector<float>;
  tagger_det.tro_4_v_dir2_mag = new std::vector<float>;
  tagger_det.tro_4_v_angle = new std::vector<float>;
  tagger_det.tro_4_v_angle1 = new std::vector<float>;
  tagger_det.tro_4_v_angle2 = new std::vector<float>;
  tagger_det.tro_4_v_length = new std::vector<float>;
  tagger_det.tro_4_v_length1 = new std::vector<float>;
  tagger_det.tro_4_v_medium_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_end_dQ_dx = new std::vector<float>;
  tagger_det.tro_4_v_energy = new std::vector<float>;
  tagger_det.tro_4_v_shower_main_length = new std::vector<float>;
  tagger_det.tro_4_v_flag_shower_trajectory = new std::vector<float>;
  tagger_det.tro_4_v_flag = new std::vector<float>;
  tagger_det.tro_5_v_max_angle = new std::vector<float>;
  tagger_det.tro_5_v_min_angle = new std::vector<float>;
  tagger_det.tro_5_v_max_length = new std::vector<float>;
  tagger_det.tro_5_v_iso_angle = new std::vector<float>;
  tagger_det.tro_5_v_n_vtx_segs = new std::vector<float>;
  tagger_det.tro_5_v_min_count = new std::vector<float>;
  tagger_det.tro_5_v_max_count = new std::vector<float>;
  tagger_det.tro_5_v_energy = new std::vector<float>;
  tagger_det.tro_5_v_flag = new std::vector<float>;
  tagger_det.lol_1_v_energy = new std::vector<float>;
  tagger_det.lol_1_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_1_v_nseg = new std::vector<float>;
  tagger_det.lol_1_v_angle = new std::vector<float>;
  tagger_det.lol_1_v_flag = new std::vector<float>;
  tagger_det.lol_2_v_length = new std::vector<float>;
  tagger_det.lol_2_v_angle = new std::vector<float>;
  tagger_det.lol_2_v_type = new std::vector<float>;
  tagger_det.lol_2_v_vtx_n_segs = new std::vector<float>;
  tagger_det.lol_2_v_energy = new std::vector<float>;
  tagger_det.lol_2_v_shower_main_length = new std::vector<float>;
  tagger_det.lol_2_v_flag_dir_weak = new std::vector<float>;
  tagger_det.lol_2_v_flag = new std::vector<float>;
  tagger_det.cosmict_flag_10 = new std::vector<float>;
  tagger_det.cosmict_10_flag_inside = new std::vector<float>;
  tagger_det.cosmict_10_vtx_z = new std::vector<float>;
  tagger_det.cosmict_10_flag_shower = new std::vector<float>;
  tagger_det.cosmict_10_flag_dir_weak = new std::vector<float>;
  tagger_det.cosmict_10_angle_beam = new std::vector<float>;
  tagger_det.cosmict_10_length = new std::vector<float>;
  tagger_det.numu_cc_flag_1 = new std::vector<float>;
  tagger_det.numu_cc_1_particle_type = new std::vector<float>;
  tagger_det.numu_cc_1_length = new std::vector<float>;
  tagger_det.numu_cc_1_medium_dQ_dx = new std::vector<float>;
  tagger_det.numu_cc_1_dQ_dx_cut = new std::vector<float>;
  tagger_det.numu_cc_1_direct_length = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_1_n_daughter_all = new std::vector<float>;
  tagger_det.numu_cc_flag_2 = new std::vector<float>;
  tagger_det.numu_cc_2_length = new std::vector<float>;
  tagger_det.numu_cc_2_total_length = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_tracks = new std::vector<float>;
  tagger_det.numu_cc_2_n_daughter_all = new std::vector<float>;
  tagger_det.ssm_kine_energy_particle = new std::vector<float>;
  tagger_det.ssm_kine_energy_info = new std::vector<int>;
  tagger_det.ssm_kine_particle_type = new std::vector<int>;
  tagger_det.ssm_kine_energy_included = new std::vector<int>;
  tagger_det.ssm_cosmict_flag_10 = new std::vector<float>;
  tagger_det.WCPMTInfoPePred = new std::vector<double>;
  tagger_det.WCPMTInfoPeMeas = new std::vector<double>;
  tagger_det.WCPMTInfoPeMeasErr = new std::vector<double>;


  set_tree_address(T_BDTvars_cv, tagger_cv,2 );
  set_tree_address(T_eval_cv, eval_cv);
  set_tree_address(T_PFeval_cv, pfeval_cv);
  set_tree_address(T_pot_cv, pot_cv);
  set_tree_address(T_KINEvars_cv, kine_cv);


  put_tree_address(t4_cv, tagger_cv,2);
  put_tree_address(t1_cv, eval_cv);
  put_tree_address(t3_cv, pfeval_cv);
  put_tree_address(t2_cv, pot_cv);
  put_tree_address(t5_cv, kine_cv);

  set_tree_address(T_BDTvars_det, tagger_det,2 );
  set_tree_address(T_eval_det, eval_det);
  set_tree_address(T_PFeval_det, pfeval_det);
  set_tree_address(T_pot_det, pot_det);
  set_tree_address(T_KINEvars_det, kine_det);

  put_tree_address(t4_det, tagger_det,2);
  put_tree_address(t1_det, eval_det);
  put_tree_address(t3_det, pfeval_det);
  put_tree_address(t2_det, pot_det);
  put_tree_address(t5_det, kine_det);



  std::map<std::pair<int, int>, int> map_re_entry_cv;
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_cv;

  bool flag_presel = false;

  for (int i=0;i!=T_eval_cv->GetEntries();i++){
    T_eval_cv->GetEntry(i);
    T_BDTvars_cv->GetEntry(i);

    map_rs_re_cv[std::make_pair(eval_cv.run, eval_cv.subrun)].insert(std::make_pair(eval_cv.run, eval_cv.event));

    int tmp_match_found = eval_cv.match_found;
    if (eval_cv.is_match_found_int) tmp_match_found = eval_cv.match_found_asInt;

    flag_presel = false;
    if (tmp_match_found != 0 && eval_cv.stm_eventtype != 0 && eval_cv.stm_lowenergy ==0 && eval_cv.stm_LM ==0 && eval_cv.stm_TGM ==0 && eval_cv.stm_STM==0 && eval_cv.stm_FullDead == 0 && eval_cv.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }

    if (tmp_match_found == -1  || (tmp_match_found == 1 && eval_cv.stm_lowenergy == -1) || (flag_presel && tagger_cv.numu_cc_flag == -1)) continue;

    map_re_entry_cv[std::make_pair(eval_cv.run, eval_cv.event)] = i;

  }

  std::map<std::pair<int, int>, int> map_re_entry_det;
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_det;
  for (int i=0;i!=T_eval_det->GetEntries();i++){
    T_eval_det->GetEntry(i);
    T_BDTvars_det->GetEntry(i);

    map_rs_re_det[std::make_pair(eval_det.run, eval_det.subrun)].insert(std::make_pair(eval_det.run, eval_det.event));

    int tmp_match_found = eval_det.match_found;
    if (eval_det.is_match_found_int) tmp_match_found = eval_det.match_found_asInt;

    flag_presel = false;
    if (tmp_match_found != 0 && eval_det.stm_eventtype != 0 && eval_det.stm_lowenergy ==0 && eval_det.stm_LM ==0 && eval_det.stm_TGM ==0 && eval_det.stm_STM==0 && eval_det.stm_FullDead == 0 && eval_det.stm_clusterlength >0) {
      flag_presel = true; // preselection ...
    }

    if (tmp_match_found == -1  || (tmp_match_found == 1&& eval_det.stm_lowenergy == -1)  || (flag_presel && tagger_det.numu_cc_flag == -1)) continue;

    map_re_entry_det[std::make_pair(eval_det.run, eval_det.event)] = i;

  }

  std::map<std::tuple<int, int>, std::pair<int, double> > map_rs_entry_pot_cv;
  for (Int_t i=0;i!=T_pot_cv->GetEntries();i++){
    T_pot_cv->GetEntry(i);
    map_rs_entry_pot_cv[std::make_pair(pot_cv.runNo,pot_cv.subRunNo)] = std::make_pair(i, pot_cv.pot_tor875);
  }

  std::map<std::tuple<int, int>, std::pair<int, double> > map_rs_entry_pot_det;
  for (Int_t i=0;i!=T_pot_det->GetEntries();i++){
    T_pot_det->GetEntry(i);
    map_rs_entry_pot_det[std::make_pair(pot_det.runNo,pot_det.subRunNo)] = std::make_pair(i, pot_det.pot_tor875);
  }



   T_eval_cv->SetBranchStatus("*",1);
   T_PFeval_cv->SetBranchStatus("*",1);
   T_BDTvars_cv->SetBranchStatus("*",1);
   T_KINEvars_cv->SetBranchStatus("*",1);

   T_eval_det->SetBranchStatus("*",1);
   T_PFeval_det->SetBranchStatus("*",1);
   T_BDTvars_det->SetBranchStatus("*",1);
   T_KINEvars_det->SetBranchStatus("*",1);



  // fill the trees ...
  std::map<std::pair<int, int>, std::set<std::pair<int, int> > > map_rs_re_common;
  std::map<int, int> map_cv_det_index;

  for (auto it = map_re_entry_cv.begin(); it != map_re_entry_cv.end(); it++){
    auto it1 = map_re_entry_det.find(it->first);
    if (it1 != map_re_entry_det.end()){ // common ...
      map_cv_det_index[it->second] = it1->second;
    }
  }

  std::cout<<"Start Filling the trees"<<std::endl;
  for (auto it = map_cv_det_index.begin(); it != map_cv_det_index.end(); it++){
      T_BDTvars_cv->GetEntry(it->first);
      T_eval_cv->GetEntry(it->first);
      T_KINEvars_cv->GetEntry(it->first);
      T_PFeval_cv->GetEntry(it->first);

      T_BDTvars_det->GetEntry(it->second);
      T_eval_det->GetEntry(it->second);
      T_KINEvars_det->GetEntry(it->second);
      T_PFeval_det->GetEntry(it->second);

      T_spacepoints_cv->GetEntry(it->first);

      T_spacepoints_det->GetEntry(it->second);

      for(auto tree_it=old_trees_cv->begin(); tree_it!=old_trees_cv->end(); tree_it++){
          (*tree_it)->GetEntry(it->first);
      }
      for(auto tree_it=old_trees_det->begin(); tree_it!=old_trees_det->end(); tree_it++){
          (*tree_it)->GetEntry(it->second);
      }


      map_rs_re_common[std::make_pair(eval_cv.run, eval_cv.subrun)].insert(std::make_pair(eval_cv.run, eval_cv.event));

      t1_cv->Fill();
      t3_cv->Fill();
      t4_cv->Fill();
      t5_cv->Fill();

      t1_det->Fill();
      t3_det->Fill();
      t4_det->Fill();
      t5_det->Fill();

      T_spacepoints_new_cv->Fill();

      T_spacepoints_new_det->Fill();

      for(auto tree_it=new_trees_cv->begin(); tree_it!=new_trees_cv->end(); tree_it++){
          (*tree_it)->Fill();
      }
      for(auto tree_it=new_trees_det->begin(); tree_it!=new_trees_det->end(); tree_it++){
          (*tree_it)->Fill();
      }

  }
  double cv_pot=0;
  double det_pot=0;
  for (auto it = map_rs_entry_pot_cv.begin(); it != map_rs_entry_pot_cv.end(); it++){
    cv_pot += it->second.second;
  }
  for (auto it = map_rs_entry_pot_det.begin(); it != map_rs_entry_pot_det.end(); it++){
    det_pot += it->second.second;
  }

  double common_pot = 0;
  for (auto it = map_rs_re_common.begin(); it != map_rs_re_common.end(); it++){
    auto it1 = map_rs_entry_pot_cv.find(it->first);
    auto it2 = map_rs_entry_pot_det.find(it->first);

    if (it1 != map_rs_entry_pot_cv.end() && it2 != map_rs_entry_pot_det.end()){
      T_pot_cv->GetEntry(it1->second.first);
      T_pot_det->GetEntry(it2->second.first);

      double ratio = it->second.size() * 1.0/map_rs_re_cv[it->first].size();
      pot_cv.pot_tor875 *= ratio;
      pot_cv.pot_tor875good *= ratio;

      ratio = it->second.size() * 1.0 /map_rs_re_det[it->first].size();
      pot_det.pot_tor875 *= ratio;
      pot_det.pot_tor875good *= ratio;

      common_pot += pot_cv.pot_tor875;

      if(override_pot){
        pot_cv.pot_tor875=5e18;
        pot_cv.pot_tor875good=5e18;
        pot_det.pot_tor875=5e18;
        pot_det.pot_tor875good=5e18;
      }

      t2_cv->Fill();
      t2_det->Fill();
      for(auto tree_it=old_trees_pot_cv->begin(); tree_it!=old_trees_pot_cv->end(); tree_it++){
        (*tree_it)->GetEntry(it1->second.first);
      }
      for(auto tree_it=new_trees_pot_cv->begin(); tree_it!=new_trees_pot_cv->end(); tree_it++){
        (*tree_it)->Fill();
      }
      for(auto tree_it=old_trees_pot_det->begin(); tree_it!=old_trees_pot_det->end(); tree_it++){
        (*tree_it)->GetEntry(it2->second.first);
      }
      for(auto tree_it=new_trees_pot_det->begin(); tree_it!=new_trees_pot_det->end(); tree_it++){
        (*tree_it)->Fill();
      }    
    }
  }

  int nevents_cv = 0;
  int nevents_det = 0;
  for (auto it = map_rs_re_cv.begin(); it != map_rs_re_cv.end(); it++){
    nevents_cv += it->second.size();
  }
  for (auto it = map_rs_re_det.begin(); it != map_rs_re_det.end(); it++){
    nevents_det += it->second.size();
  }
  std::cout << out_file << std::endl;
  std::cout << "Events: " << t1_cv->GetEntries() << " " << nevents_cv<<"/"<<T_eval_cv->GetEntries() << " " << nevents_det << "/" << T_eval_det->GetEntries() << std::endl;
  std::cout << "POT:    " << common_pot << " " << cv_pot << " " << det_pot << std::endl;



  file3->Write("",TFile::kOverwrite);
  file3->Close();


  return 0;
}
