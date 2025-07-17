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

#include "WCPLEEANA/Util.h"

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

  TTree *NeutrinoSelectionFilter_cv;
  TTree *SubRun_cv;
  TDirectory *shrreco3d_cv;
  TDirectory *proximity_cv;
  TDirectory *nuselection_cv;
  bool has_pelee_cv = false;
  if (file1->GetDirectory("nuselection")){
    has_pelee_cv = true;
    TDirectory *topdir_cv = gDirectory;
    std::cout<<"nuselection_cv"<<std::endl;
    file1->cd("nuselection");
    nuselection_cv = gDirectory;
    std::cout<<"shrreco3d_cv"<<std::endl;
    file1->cd("shrreco3d");
    shrreco3d_cv = gDirectory;
    std::cout<<"proximity_cv"<<std::endl;
    file1->cd("proximity");
    proximity_cv = gDirectory;
    topdir_cv->cd();
    NeutrinoSelectionFilter_cv = (TTree*)file1->Get("nuselection/NeutrinoSelectionFilter");
    SubRun_cv = (TTree*)file1->Get("nuselection/SubRun");
    //TTree *_energy_tree_cv = (TTree*)file1->Get("shrreco3d/_energy_tree");
    //TTree *_dedx_tree_cv = (TTree*)file1->Get("shrreco3d/_dedx_tree");
    //TTree *_rcshr_tree_cv = (TTree*)file1->Get("shrreco3d/_rcshr_tree");
    //TTree *_clus_tree_cv = (TTree*)file1->Get("proximity/_clus_tree");
  }

  TTree *vertex_tree_cv;
  //TTree *pot_tree_cv;
  TTree *eventweight_tree_cv;
  //TTree *ncdelta_slice_tree_cv;
  //TTree *run_subrun_tree_cv;
  TTree *geant4_tree_cv;
  //TTree *true_eventweight_tree_cv;
  TDirectory *singlephotonana_cv;
  bool has_glee_cv = false;
  if (file1->GetDirectory("singlephotonana")){
    has_glee_cv = true;
    TDirectory *topdir_cv = gDirectory;
    std::cout<<"singlephotonana_cv"<<std::endl;
    file1->cd("singlephotonana");
    singlephotonana_cv = gDirectory;
    topdir_cv->cd();
    vertex_tree_cv = (TTree*)file1->Get("singlephotonana/vertex_tree");
    //pot_tree_cv = (TTree*)file1->Get("singlephotonana/pot_tree");
    eventweight_tree_cv = (TTree*)file1->Get("singlephotonana/eventweight_tree");
    //ncdelta_slice_tree_cv = (TTree*)file1->Get("singlephotonana/ncdelta_slice_tree");
    //run_subrun_tree_cv = (TTree*)file1->Get("singlephotonana/run_subrun_tree");
    geant4_tree_cv = (TTree*)file1->Get("singlephotonana/geant4_tree");
  }

  TTree *EventTree_cv;
  bool has_lantern_cv = false;
  if (file1->GetDirectory("lantern")){
    has_lantern_cv = true;
    EventTree_cv = (TTree*)file1->Get("lantern/EventTree");
  }




  TTree *NeutrinoSelectionFilter_det;
  TTree *SubRun_det;
  TDirectory *shrreco3d_det;
  TDirectory *proximity_det;
  TDirectory *nuselection_det;
  bool has_pelee_det = false;
  if (file2->GetDirectory("nuselection")){
    has_pelee_det = true;
    TDirectory *topdir_det = gDirectory;
    file2->cd("nuselection");
    nuselection_det = gDirectory;
    file2->cd("shrreco3d");
    shrreco3d_det = gDirectory;
    file2->cd("proximity");
    proximity_det = gDirectory;
    topdir_det->cd();
    NeutrinoSelectionFilter_det = (TTree*)file2->Get("nuselection/NeutrinoSelectionFilter");
    SubRun_det = (TTree*)file2->Get("nuselection/SubRun");
    //TTree *_energy_tree_det = (TTree*)file2->Get("shrreco3d/_energy_tree");
    //TTree *_dedx_tree_det = (TTree*)file2->Get("shrreco3d/_dedx_tree");
    //TTree *_rcshr_tree_det = (TTree*)file2->Get("shrreco3d/_rcshr_tree");
    //TTree *_clus_tree_det = (TTree*)file2->Get("proximity/_clus_tree");
  }

  TTree *vertex_tree_det;
  //TTree *pot_tree_det;
  TTree *eventweight_tree_det;
  //TTree *ncdelta_slice_tree_det;
  //TTree *run_subrun_tree_det;
  TTree *geant4_tree_det;
  //TTree *true_eventweight_tree_det;
  TDirectory *singlephotonana_det;
  bool has_glee_det = false;
  if (file2->GetDirectory("singlephotonana")){
    has_glee_det = true;
    TDirectory *topdir_det = gDirectory;
    file2->cd("singlephotonana");
    singlephotonana_det = gDirectory;
    topdir_det->cd();
    vertex_tree_det = (TTree*)file2->Get("singlephotonana/vertex_tree");
    //pot_tree_det = (TTree*)file2->Get("singlephotonana/pot_tree");
    eventweight_tree_det = (TTree*)file2->Get("singlephotonana/eventweight_tree");
    //ncdelta_slice_tree_det = (TTree*)file2->Get("singlephotonana/ncdelta_slice_tree");
    //run_subrun_tree_det = (TTree*)file2->Get("singlephotonana/run_subrun_tree");
    geant4_tree_det = (TTree*)file2->Get("singlephotonana/geant4_tree");
  }

  TTree *EventTree_det;
  bool has_lantern_det = false;
  if (file2->GetDirectory("lantern")){
    has_lantern_det = true;
    EventTree_det = (TTree*)file2->Get("lantern/EventTree");
  }




  TFile *file3 = new TFile(out_file,"RECREATE");

  TTree *new_NeutrinoSelectionFilter_cv;// = new TTree("NeutrinoSelectionFilter","NeutrinoSelectionFilter");
  TTree *new_SubRun_cv;// = new TTree("SubRun","SubRun");
  TTree *new_NeutrinoSelectionFilter_det;// = new TTree("NeutrinoSelectionFilter","NeutrinoSelectionFilter");
  TTree *new_SubRun_det;// = new TTree("SubRun","SubRun");
  if (has_pelee_cv && has_pelee_det){
    TDirectory *topdirout = gDirectory;
    file3->mkdir("nuselection");
    file3->cd("nuselection");
    //new_NeutrinoSelectionFilter_cv = NeutrinoSelectionFilter_cv->CloneTree(0);
    //new_SubRun_cv = SubRun_cv->CloneTree(0);
    //new_NeutrinoSelectionFilter_det = NeutrinoSelectionFilter_det->CloneTree(0);
    //new_SubRun_det = SubRun_det->CloneTree(0);
    CopyDir(nuselection_cv,"_cv",true);
    CopyDir(nuselection_det,"_det",true);
    topdirout->cd();
    file3->mkdir("shrreco3d");
    file3->cd("shrreco3d");
    CopyDir(shrreco3d_cv,"_cv",true);
    CopyDir(shrreco3d_det,"_det",true);
    topdirout->cd();
    file3->mkdir("proximity");
    file3->cd("proximity");
    CopyDir(proximity_cv,"_cv",true);
    CopyDir(proximity_det,"_det",true);
    topdirout->cd();
    new_NeutrinoSelectionFilter_cv = (TTree*)file3->Get("nuselection/NeutrinoSelectionFilter_cv");
    //new_NeutrinoSelectionFilter_cv->SetBranchStatus("*",1);
    new_NeutrinoSelectionFilter_det = (TTree*)file3->Get("nuselection/NeutrinoSelectionFilter_det");
    //new_NeutrinoSelectionFilter_det->SetBranchStatus("*",1);
  }

  TTree *new_vertex_tree_cv;
  //TTree *pot_tree_cv;
  TTree *new_eventweight_tree_cv;
  //TTree *ncdelta_slice_tree_cv;
  //TTree *new_run_subrun_tree_cv;
  TTree *new_geant4_tree_cv;
  //TTree *true_eventweight_tree_cv;
  TTree *new_vertex_tree_det;
  //TTree *pot_tree_det;
  TTree *new_eventweight_tree_det;
  //TTree *ncdelta_slice_tree_det;
  //TTree *new_run_subrun_tree_det;
  TTree *new_geant4_tree_det;
  //TTree *true_eventweight_tree_det;
  std::cout<<"Checking glee"<<std::endl;
  if (has_glee_cv && has_glee_det){
  std::cout<<"Starting glee"<<std::endl;
    TDirectory *topdirout = gDirectory;
    file3->mkdir("singlephotonana");
    file3->cd("singlephotonana");
    CopyDir(singlephotonana_cv,"_cv",true);
    CopyDir(singlephotonana_det,"_det",true);
    topdirout->cd();
    new_vertex_tree_cv = (TTree*)file3->Get("singlephotonana/vertex_tree_cv");
    //new_vertex_tree_cv->SetBranchStatus("*",1);
    new_vertex_tree_det = (TTree*)file3->Get("singlephotonana/vertex_tree_det");
    //new_vertex_tree_det->SetBranchStatus("*",1);
    
  }else{std::cout<<"No glee"<<std::endl;}

  TTree *new_EventTree_cv;
  TTree *new_EventTree_det;
  std::cout<<"Checking lantern"<<std::endl;

  if (has_lantern_cv && has_lantern_det){
  std::cout<<"Starting lantern"<<std::endl;
    TDirectory *topdirout = gDirectory;
    file3->mkdir("lantern");
    file3->cd("lantern");
    new_EventTree_cv = EventTree_cv->CloneTree(0);
    new_EventTree_det = EventTree_det->CloneTree(0);
    topdirout->cd();
  }else{std::cout<<"No lantern"<<std::endl;}


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

    if (has_pelee_cv && has_pelee_det){
      NeutrinoSelectionFilter_cv->GetEntry(it->first);
      NeutrinoSelectionFilter_det->GetEntry(it->second);
    }

    if (has_glee_cv && has_glee_det){

      vertex_tree_cv->GetEntry(it->first);
      //eventweight_tree_cv->GetEntry(it->first);
      //geant4_tree_cv->GetEntry(it->first);

      vertex_tree_det->GetEntry(it->second);
      //eventweight_tree_det->GetEntry(it->second);
      //geant4_tree_det->GetEntry(it->second);

    }

    if (has_lantern_cv && has_lantern_det){
      EventTree_cv->GetEntry(it->first);
      EventTree_det->GetEntry(it->second);
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

      if (has_pelee_cv && has_pelee_det){
        new_NeutrinoSelectionFilter_cv->Fill();
        new_NeutrinoSelectionFilter_det->Fill();
      }

      if (has_glee_cv && has_glee_det){
        new_vertex_tree_cv->Fill();
        //new_eventweight_tree_cv->Fill();
        //new_geant4_tree_cv->Fill();
        new_vertex_tree_det->Fill();
        //new_eventweight_tree_det->Fill();
        //new_geant4_tree_det->Fill();
      }

      if (has_lantern_cv && has_lantern_det){
        new_EventTree_cv->Fill();
        new_EventTree_det->Fill();
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

      t2_cv->Fill();
      t2_det->Fill();
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



  file3->Write();
  file3->Close();


  return 0;
}
