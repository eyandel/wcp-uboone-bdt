#ifndef UBOONE_LEE_PANDORA
#define UBOONE_LEE_PANDORA

namespace LEEana{
struct PandoraInfo{
  Int_t run;
  Int_t subrun;
  Int_t event;

  Int_t nslice;
  Int_t slice_orig_pass_id;
  Int_t n_pfps;

  std::vector<float> *trk_llr_pid_score_v;
  std::vector<unsigned int>* pfp_generation_v;
  std::vector<float>* trk_score_v;
  std::vector<float>* trk_energy_proton_v;
  std::vector<int>* pfpdg;  
};

 void set_tree_address(TTree *tree0, PandoraInfo& pandora_info);
 void put_tree_address(TTree *tree0, PandoraInfo& pandora_info);
 void clear_pandora_info(PandoraInfo& pandora_info);
 void init_pointers(PandoraInfo& pandora_info);
 void del_pointers(PandoraInfo& pandora_info);

}

void LEEana::init_pointers(PandoraInfo& pandora_info) {
  pandora_info.trk_llr_pid_score_v = new std::vector<float>;
  pandora_info.pfp_generation_v = new std::vector<unsigned int>;
  pandora_info.trk_score_v = new std::vector<float>;
  pandora_info.trk_energy_proton_v = new std::vector<float>;
  pandora_info.pfpdg = new std::vector<int>;
}

void LEEana::del_pointers(PandoraInfo& pandora_info) {
  delete pandora_info.trk_llr_pid_score_v;
  delete pandora_info.pfp_generation_v;
  delete pandora_info.trk_score_v;
  delete pandora_info.trk_energy_proton_v;
  delete pandora_info.pfpdg;
}

void LEEana::clear_pandora_info(PandoraInfo& pandora_info) {
  pandora_info.trk_llr_pid_score_v->clear();
  pandora_info.pfp_generation_v->clear();
  pandora_info.trk_score_v->clear();
  pandora_info.trk_energy_proton_v->clear();
  pandora_info.pfpdg->clear();
}

void LEEana::set_tree_address(TTree *tree0, PandoraInfo& pandora_info) {

  init_pointers(pandora_info);

  tree0->SetBranchAddress("run", &pandora_info.run);
  tree0->SetBranchAddress("sub", &pandora_info.subrun);
  tree0->SetBranchAddress("evt", &pandora_info.event);

  tree0->SetBranchAddress("nslice", &pandora_info.nslice);
  tree0->SetBranchAddress("slice_orig_pass_id", &pandora_info.slice_orig_pass_id);
  tree0->SetBranchAddress("n_pfps", &pandora_info.n_pfps);

  tree0->SetBranchAddress("trk_llr_pid_score_v", &pandora_info.trk_llr_pid_score_v);
  tree0->SetBranchAddress("pfp_generation_v", &pandora_info.pfp_generation_v);
  tree0->SetBranchAddress("trk_score_v", &pandora_info.trk_score_v);
  tree0->SetBranchAddress("trk_energy_proton_v", &pandora_info.trk_energy_proton_v);
  tree0->SetBranchAddress("pfpdg", &pandora_info.pfpdg);
}

void LEEana::put_tree_address(TTree *tree0, PandoraInfo& pandora_info) {
  tree0->Branch("run", &pandora_info.run, "run/I");
  tree0->Branch("subrun", &pandora_info.subrun, "subrun/I");
  tree0->Branch("event", &pandora_info.event, "event/I");

  tree0->Branch("nslice", &pandora_info.nslice);
  tree0->Branch("slice_orig_pass_id", &pandora_info.slice_orig_pass_id);
  tree0->Branch("n_pfps", &pandora_info.n_pfps);

  tree0->Branch("trk_llr_pid_score_v", &pandora_info.trk_llr_pid_score_v);
  tree0->Branch("pfp_generation_v", &pandora_info.pfp_generation_v);
  tree0->Branch("trk_score_v", &pandora_info.trk_score_v);
  tree0->Branch("trk_energy_proton_v", &pandora_info.trk_energy_proton_v);
  tree0->Branch("pfpdg", &pandora_info.pfpdg);
}

#endif
