#ifndef UBOONE_LEE_PFEVAL
#define UBOONE_LEE_PFEVAL

namespace LEEana{
struct PFevalInfo{
  bool flag_NCDelta;
  bool flag_showerMomentum;
  bool flag_recoprotonMomentum;
  bool flag_pf_truth;
  bool flag_pf_reco;
  bool flag_init_pointers;
  //Erin
  bool flag_single_photon;
  bool flag_nsbeam;
  //
  bool flag_mcs;
  bool flag_larpid;
  bool flag_backtracking;
  bool flag_ns_time_cor;
  bool flag_Phtot;
  bool flag_PMT;

  Int_t run;
  Int_t subrun;
  Int_t event;
  Int_t neutrino_type;
  Float_t reco_nuvtxX;
  Float_t reco_nuvtxY;
  Float_t reco_nuvtxZ;
  Float_t reco_showervtxX;
  Float_t reco_showervtxY;
  Float_t reco_showervtxZ;
  Float_t reco_showerKE;
  Float_t reco_muonvtxX;
  Float_t reco_muonvtxY;
  Float_t reco_muonvtxZ;
  Float_t reco_muonMomentum[4];

  Float_t reco_showerMomentum[4];
  // new variables added, the will go with reco_showerMomentum ...
  Int_t reco_Nproton;
  Int_t mcflux_run;
  Int_t mcflux_evtno;
  Int_t mcflux_ndecay;
  Int_t mcflux_ntype;
  Float_t mcflux_nuEnergy;
  Float_t mcflux_vx;
  Float_t mcflux_vy;
  Float_t mcflux_vz;
  Float_t mcflux_genx;
  Float_t mcflux_geny;
  Float_t mcflux_genz;
  Float_t mcflux_dk2gen;
  Float_t mcflux_gen2vtx;
  Int_t truth_nuScatType;
  Float_t truth_nu_pos[4];
  Float_t truth_showerMomentum[4];
  Float_t truth_nu_momentum[4];

  Float_t nuvtx_diff;
  Float_t showervtx_diff;
  Float_t muonvtx_diff;

  Float_t truth_corr_nuvtxX;
  Float_t truth_corr_nuvtxY;
  Float_t truth_corr_nuvtxZ;
  Float_t truth_corr_showervtxX;
  Float_t truth_corr_showervtxY;
  Float_t truth_corr_showervtxZ;
  Float_t truth_showerKE;
  //Erin
  Int_t   truth_showerPdg;
  Int_t   truth_showerMother;
  //
  Float_t truth_corr_muonvtxX;
  Float_t truth_corr_muonvtxY;
  Float_t truth_corr_muonvtxZ;
  Float_t truth_muonvtxX;
  Float_t truth_muonvtxY;
  Float_t truth_muonvtxZ;
  Float_t truth_muonendX;
  Float_t truth_muonendY;
  Float_t truth_muonendZ;
  Float_t truth_muonMomentum[4];
  Float_t truth_nuEnergy;
  Float_t truth_energyInside;
  Float_t truth_electronInside;
  Int_t truth_nuPdg;
  Bool_t truth_isCC;
  Float_t truth_vtxX;
  Float_t truth_vtxY;
  Float_t truth_vtxZ;
  Float_t truth_nuTime;
  Int_t truth_nuIntType;

  //
  Int_t truth_NprimPio;
  Float_t truth_pio_energy_1;
  Float_t truth_pio_energy_2;
  Float_t truth_pio_angle;
  Int_t truth_NCDelta;
  //Erin
  Int_t truth_Npi0;
  Int_t   truth_single_photon;
  Float_t truth_photon_angle;
  Float_t truth_photon_dis;
  Float_t evtDeltaTimeNS;
  Float_t evtTimeNS;
  //
  Float_t reco_protonMomentum[4];


  // PF particle
  Int_t truth_Ntrack;  // number of tracks in MC
  Int_t truth_id[10000];  // track id; size == truth_Ntrack
  Int_t truth_pdg[10000];  // track particle pdg; size == truth_Ntrack
  std::vector<std::string > *truth_process;
  Int_t truth_mother[10000];  // mother id of this track; size == truth_Ntrack
  Float_t truth_startXYZT[10000][4];  // start position of this track; size == truth_Ntrack
  Float_t truth_endXYZT[10000][4];  // end position of this track; size == truth_Ntrack
  Float_t truth_startMomentum[10000][4];  // start momentum of this track; size == truth_Ntrack
  Float_t truth_endMomentum[10000][4];  // end momentum of this track; size == truth_Ntrack
  std::vector<std::vector<Int_t> > *truth_daughters;  // daughters id of this track; vector

  TObjArray *fMC_trackPosition;

  Int_t reco_Ntrack;  // number of tracks in MC
  Int_t reco_id[10000];  // track id; size == reco_Ntrack
  Int_t reco_pdg[10000];  // track particle pdg; size == reco_Ntrack
  std::vector<std::string > *reco_process;
  Int_t reco_mother[10000];  // mother id of this track; size == reco_Ntrack
  Float_t reco_startXYZT[10000][4];  // start position of this track; size == reco_Ntrack
  Float_t reco_endXYZT[10000][4];  // end position of this track; size == reco_Ntrack
  Float_t reco_startMomentum[10000][4];  // start momentum of this track; size == reco_Ntrack
  Float_t reco_endMomentum[10000][4];  // end momentum of this track; size == reco_Ntrack

  std::vector<std::vector<Int_t> > *reco_daughters;  // daughters id of this track; vector

  Int_t mc_isnu; // is neutrino Int_teraction
  Int_t mc_nGeniePrimaries; // number of Genie primaries
  Int_t mc_nu_pdg; // pdg code of neutrino
  Int_t mc_nu_ccnc; // cc or nc
  Int_t mc_nu_mode; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
  Int_t mc_nu_intType; // interaction type
  Int_t mc_nu_target; // target interaction
  Int_t mc_hitnuc; // hit nucleon
  Int_t mc_hitquark; // hit quark

  Double_t mc_nu_Q2; // Q^2
  Double_t mc_nu_W; // W
  Double_t mc_nu_X; // X
  Double_t mc_nu_Y; // Y
  Double_t mc_nu_Pt; // Pt
  Double_t mc_nu_Theta; // angle relative to lepton
  Float_t mc_nu_pos[4];  // interaction position of nu
  Float_t mc_nu_mom[4];  // interaction momentum of nu

  //new mcs variables
  Double_t mcs_mu_tracklen; // Muon length in cm
  Double_t mcs_emu_tracklen; // Pure range based, in GeV
  Double_t mcs_emu_MCS; // MCS based, in GeV
  Double_t mcs_ambiguity_MCS; // "Goodness" of MCS estimation

  //new parpid and backtrackng variables
  Int_t reco_larpid_classified[10000]; //LArPID pdg label
  Int_t reco_larpid_pdg[10000]; //LArPID pdg label
  Int_t reco_larpid_proccess[10000]; //LArPID primary vs. secondary production process class
  Float_t reco_larpid_completeness[10000]; //LArPID reconstruction completeness estimate
  Float_t reco_larpid_purity[10000]; //LArPID reconstruction purity estimate
  Float_t reco_larpid_pidScore_el[10000]; //LArPID electron score from particle classifier
  Float_t reco_larpid_pidScore_ph[10000]; //LArPID photon score from particle classifier
  Float_t reco_larpid_pidScore_mu[10000]; //LArPID muon score from particle classifier
  Float_t reco_larpid_pidScore_pi[10000]; //LArPID pion score from particle classifier
  Float_t reco_larpid_pidScore_pr[10000]; //LArPID proton score from particle classifier
  Float_t reco_larpid_procScore_prim[10000]; //LArPID primary score from production process classifier
  Float_t reco_larpid_procScore_ntrl[10000]; //LArPID secondary with neutral parent score from process classifier
  Float_t reco_larpid_procScore_chgd[10000]; //LArPID secondary with charged parent score from process classifier
  Int_t reco_truthMatch_pdg[10000]; //MC back tracking: pdg of best match simulated particle
  Int_t reco_truthMatch_id[10000]; //MC back tracking: track id (truth_id variable) of best match sim particle
  Float_t reco_truthMatch_purity[10000]; //MC back tracking: fraction of reco particle from best match sim particle
  Float_t reco_truthMatch_completeness[10000]; //MC backtracking: fraction of best match sim particle reconstructed
  Int_t reco_truthMatch_nSimParts[10000]; //MC backtracking: number of different sim particle types matched
  std::vector<std::vector<int>> *reco_truthMatch_simPart_pdg; //MC backtracking: pdg of contributing particle type
  std::vector<std::vector<float>> *reco_truthMatch_simPart_purity;//MC backtracking: purity contributing particle type

  //new ns timing vars
  Float_t evtTimeNS_cor;
  Float_t Ph_Tot;
  std::vector<int> *PMT_ID;
  std::vector<float> *PMT_Time;
  std::vector<float> *PMT_Amp;
  std::vector<float> *PMT_TimeProp;
  std::vector<float> *PMT_TimeDP;
  std::vector<float> *PMT_TimeDL;
  std::vector<bool> *PMT_Sat;
  Float_t RWM_Time;
  Float_t cor_nu_time;
  Float_t cor_nu_time_nospill;
  Float_t cor_nu_time_spill;
  Float_t cor_nu_deltatime;
  

};


 void set_tree_address(TTree *tree0, PFevalInfo& tagger_info, int flag = 1);
 void put_tree_address(TTree *tree0, PFevalInfo& tagger_info, int flag = 1);
 void clear_pfeval_info(PFevalInfo& tagger_info);
 void init_pointers(PFevalInfo& tagger_info);
 void del_pointers(PFevalInfo& tagger_info);
}

void LEEana::init_pointers(PFevalInfo& tagger_info){

  tagger_info.truth_process = new std::vector<std::string >;
  tagger_info.truth_daughters = new   std::vector<std::vector<Int_t> >;
  tagger_info.fMC_trackPosition = new TObjArray();
  tagger_info.fMC_trackPosition->SetOwner(kTRUE);
  tagger_info.reco_process = new   std::vector<std::string >;
  tagger_info.reco_daughters = new   std::vector<std::vector<Int_t> > ;

  tagger_info.reco_truthMatch_simPart_pdg = new std::vector<std::vector<int>>; 
  tagger_info.reco_truthMatch_simPart_purity = new std::vector<std::vector<float>>;

  tagger_info.PMT_ID = new std::vector<int>;
  tagger_info.PMT_Time = new std::vector<float>;
  tagger_info.PMT_Amp = new std::vector<float>;
  tagger_info.PMT_TimeProp = new std::vector<float>;
  tagger_info.PMT_TimeDP = new std::vector<float>;
  tagger_info.PMT_TimeDL = new std::vector<float>;
  tagger_info.PMT_Sat = new std::vector<bool>;
}

void LEEana::del_pointers(PFevalInfo& tagger_info){
  delete tagger_info.truth_process;
  delete tagger_info.truth_daughters;
  delete tagger_info.fMC_trackPosition;
  delete tagger_info.reco_process;
  delete tagger_info.reco_daughters;

  delete tagger_info.reco_truthMatch_simPart_pdg; 
  delete tagger_info.reco_truthMatch_simPart_purity;

  delete tagger_info.PMT_ID;
  delete tagger_info.PMT_Time;
  delete tagger_info.PMT_Amp;
  delete tagger_info.PMT_TimeProp;
  delete tagger_info.PMT_TimeDP;
  delete tagger_info.PMT_TimeDL;
  delete tagger_info.PMT_Sat;
}


void LEEana::clear_pfeval_info(PFevalInfo& tagger_info){
  tagger_info.flag_NCDelta = false;
  tagger_info.flag_showerMomentum = false;
  tagger_info.flag_recoprotonMomentum = false;
  tagger_info.flag_pf_truth = false;
  tagger_info.flag_pf_reco = false;

  if (!tagger_info.flag_init_pointers){
    init_pointers(tagger_info);
    tagger_info.flag_init_pointers = true;
  }

  tagger_info.run=0;
  tagger_info.subrun=0;
  tagger_info.event=0;
  tagger_info.neutrino_type=0;
  tagger_info.reco_nuvtxX=0;
  tagger_info.reco_nuvtxY=0;
  tagger_info.reco_nuvtxZ=0;
  tagger_info.reco_showervtxX=0;
  tagger_info.reco_showervtxY=0;
  tagger_info.reco_showervtxZ=0;
  tagger_info.reco_showerKE=0;
  tagger_info.reco_muonvtxX=0;
  tagger_info.reco_muonvtxY=0;
  tagger_info.reco_muonvtxZ=0;

  tagger_info.nuvtx_diff=0;
  tagger_info.showervtx_diff=0;
  tagger_info.muonvtx_diff=0;

  tagger_info.truth_corr_nuvtxX=0;
  tagger_info.truth_corr_nuvtxY=0;
  tagger_info.truth_corr_nuvtxZ=0;
  tagger_info.truth_corr_showervtxX=0;
  tagger_info.truth_corr_showervtxY=0;
  tagger_info.truth_corr_showervtxZ=0;
  tagger_info.truth_showerKE=0;
  //Erin
  tagger_info.truth_showerPdg=0;
  tagger_info.truth_showerMother=0;
  //
  tagger_info.truth_corr_muonvtxX=0;
  tagger_info.truth_corr_muonvtxY=0;
  tagger_info.truth_corr_muonvtxZ=0;
  tagger_info.truth_muonvtxX=0;
  tagger_info.truth_muonvtxY=0;
  tagger_info.truth_muonvtxZ=0;
  tagger_info.truth_muonendX=0;
  tagger_info.truth_muonendY=0;
  tagger_info.truth_muonendZ=0;

  tagger_info.truth_nuEnergy=0;
  tagger_info.truth_energyInside=0;
  tagger_info.truth_electronInside=0;
  tagger_info.truth_nuPdg=0;
  tagger_info.truth_isCC=0;
  tagger_info.truth_vtxX=0;
  tagger_info.truth_vtxY=0;
  tagger_info.truth_vtxZ=0;
  tagger_info.truth_nuTime=0;
  tagger_info.truth_nuIntType=0;

  //
  tagger_info.truth_NprimPio=0;
  tagger_info.truth_pio_energy_1=0;
  tagger_info.truth_pio_energy_2=0;
  tagger_info.truth_pio_angle=0;
  tagger_info.truth_NCDelta=0;
  //Erin
  tagger_info.truth_Npi0=0;
  tagger_info.truth_single_photon=0;
  tagger_info.truth_photon_angle=0;
  tagger_info.truth_photon_dis=0;
  tagger_info.evtDeltaTimeNS=-9999.;
  tagger_info.evtTimeNS=-9999.;
  //


  tagger_info.reco_Nproton=0;
  tagger_info.mcflux_run=0;
  tagger_info.mcflux_evtno=0;
  tagger_info.mcflux_ndecay=0;
  tagger_info.mcflux_ntype=0;
  tagger_info.mcflux_nuEnergy=0;
  tagger_info.mcflux_vx=0;
  tagger_info.mcflux_vy=0;
  tagger_info.mcflux_vz=0;
  tagger_info.mcflux_genx=0;
  tagger_info.mcflux_geny=0;
  tagger_info.mcflux_genz=0;
  tagger_info.mcflux_dk2gen=0;
  tagger_info.mcflux_gen2vtx=0;
  tagger_info.truth_nuScatType=0;


  for (Int_t i=0;i!=4;i++){
    tagger_info.reco_muonMomentum[i]=0;
    tagger_info.reco_showerMomentum[i]=0;
    tagger_info.truth_muonMomentum[i]=0;
    tagger_info.reco_protonMomentum[i]=0;

    tagger_info.truth_showerMomentum[i] = 0;
    tagger_info.truth_nu_pos[i] =0;
    tagger_info.truth_nu_momentum[i]=0;
  }

  //PF ...
  tagger_info.truth_Ntrack = 0;  // number of tracks in MC
  tagger_info.truth_process->clear();
  tagger_info.truth_daughters->clear();  // daughters id of this track; vector
  tagger_info.fMC_trackPosition->Clear("");

  tagger_info.reco_Ntrack = 0;  // number of tracks in MC
  tagger_info.reco_process->clear();

  tagger_info.reco_daughters->clear();  // daughters id of this track; vector

  tagger_info.mc_isnu=0; // is neutrino eraction
  tagger_info.mc_nGeniePrimaries=0; // number of Genie primaries
  tagger_info.mc_nu_pdg=0; // pdg code of neutrino
  tagger_info.mc_nu_ccnc=0; // cc or nc
  tagger_info.mc_nu_mode=0; // mode: http://nusoft.fnal.gov/larsoft/doxsvn/html/MCNeutrino_8h_source.html
  tagger_info.mc_nu_intType=0; // interaction type
  tagger_info.mc_nu_target=0; // target interaction
  tagger_info.mc_hitnuc=0; // hit nucleon
  tagger_info.mc_hitquark=0; // hit quark

  tagger_info.mc_nu_Q2=0; // Q^2
  tagger_info.mc_nu_W=0; // W
  tagger_info.mc_nu_X=0; // X
  tagger_info.mc_nu_Y=0; // Y
  tagger_info.mc_nu_Pt=0; // Pt
  tagger_info.mc_nu_Theta=0; // angle relative to lepton
  for (Int_t i=0;i!=4;i++){
    tagger_info.mc_nu_pos[i] = 0;  // interaction position of nu
    tagger_info.mc_nu_mom[i] = 0;  // interaction momentum of nu
  }

  //new mcs variables
  tagger_info.mcs_mu_tracklen=0; // Muon length in cm
  tagger_info.mcs_emu_tracklen=0; // Pure range based, in GeV
  tagger_info.mcs_emu_MCS=0; // MCS based, in GeV
  tagger_info.mcs_ambiguity_MCS=0; // "Goodness" of MCS estimation

  //new parpid and backtrackng variables
  tagger_info.reco_truthMatch_simPart_pdg->clear();
  tagger_info.reco_truthMatch_simPart_purity->clear();

  //new ns timing vars
  tagger_info.evtTimeNS_cor = 0;
  tagger_info.Ph_Tot = 0;
  tagger_info.PMT_ID->clear();
  tagger_info.PMT_Time->clear();
  tagger_info.PMT_Amp->clear();
  tagger_info.PMT_TimeProp->clear();
  tagger_info.PMT_TimeDP->clear();
  tagger_info.PMT_TimeDL->clear();
  tagger_info.PMT_Sat->clear();
  tagger_info.RWM_Time = 0;
  tagger_info.cor_nu_time = -9999.;
  tagger_info.cor_nu_time_nospill = -9999.;
  tagger_info.cor_nu_time_spill = -9999.;
  tagger_info.cor_nu_deltatime = -9999.;
  

}

void LEEana::set_tree_address(TTree *tree0, PFevalInfo& tagger_info, int flag){
  //  std::cout << "test" << std::endl;

  tagger_info.flag_NCDelta = false;
  //Erin
  tagger_info.flag_single_photon = false;
  tagger_info.flag_nsbeam = false;
  //
  tagger_info.flag_showerMomentum = false;
  tagger_info.flag_recoprotonMomentum = false;
  tagger_info.flag_pf_truth = false;
  tagger_info.flag_pf_reco = false;
  tagger_info.flag_init_pointers = false;

  tagger_info.flag_mcs = false;
  tagger_info.flag_larpid = false;
  tagger_info.flag_backtracking = false;
  tagger_info.flag_ns_time_cor = false;
  tagger_info.flag_Phtot = false;
  tagger_info.flag_PMT = false;



  tree0->SetBranchAddress("run", &tagger_info.run);
  tree0->SetBranchAddress("subrun", &tagger_info.subrun);
  tree0->SetBranchAddress("event", &tagger_info.event);
  tree0->SetBranchAddress("neutrino_type", &tagger_info.neutrino_type);
  tree0->SetBranchAddress("reco_nuvtxX", &tagger_info.reco_nuvtxX);
  tree0->SetBranchAddress("reco_nuvtxY", &tagger_info.reco_nuvtxY);
  tree0->SetBranchAddress("reco_nuvtxZ", &tagger_info.reco_nuvtxZ);
  tree0->SetBranchAddress("reco_showervtxX", &tagger_info.reco_showervtxX);
  tree0->SetBranchAddress("reco_showervtxY", &tagger_info.reco_showervtxY);
  tree0->SetBranchAddress("reco_showervtxZ", &tagger_info.reco_showervtxZ);
  tree0->SetBranchAddress("reco_showerKE", &tagger_info.reco_showerKE);
  tree0->SetBranchAddress("reco_muonvtxX", &tagger_info.reco_muonvtxX);
  tree0->SetBranchAddress("reco_muonvtxY", &tagger_info.reco_muonvtxY);
  tree0->SetBranchAddress("reco_muonvtxZ", &tagger_info.reco_muonvtxZ);
  tree0->SetBranchAddress("reco_muonMomentum", &tagger_info.reco_muonMomentum[0]);

  if (tree0->GetBranch("reco_showerMomentum")){
    tagger_info.flag_showerMomentum = true;
    tree0->SetBranchAddress("reco_showerMomentum",&tagger_info.reco_showerMomentum[0]);
    tree0->SetBranchAddress("reco_Nproton",&tagger_info.reco_Nproton);

    if (flag==1){
      tree0->SetBranchAddress("truth_showerMomentum",&tagger_info.truth_showerMomentum[0]);

      tree0->SetBranchAddress("mcflux_run",&tagger_info.mcflux_run);
      tree0->SetBranchAddress("mcflux_evtno",&tagger_info.mcflux_evtno);
      tree0->SetBranchAddress("mcflux_ndecay",&tagger_info.mcflux_ndecay);
      tree0->SetBranchAddress("mcflux_ntype",&tagger_info.mcflux_ntype);
      tree0->SetBranchAddress("mcflux_nuEnergy",&tagger_info.mcflux_nuEnergy);
      tree0->SetBranchAddress("mcflux_vx",&tagger_info.mcflux_vx);
      tree0->SetBranchAddress("mcflux_vy",&tagger_info.mcflux_vy);
      tree0->SetBranchAddress("mcflux_vz",&tagger_info.mcflux_vz);
      tree0->SetBranchAddress("mcflux_genx",&tagger_info.mcflux_genx);
      tree0->SetBranchAddress("mcflux_geny",&tagger_info.mcflux_geny);
      tree0->SetBranchAddress("mcflux_genz",&tagger_info.mcflux_genz);
      tree0->SetBranchAddress("mcflux_dk2gen",&tagger_info.mcflux_dk2gen);
      tree0->SetBranchAddress("mcflux_gen2vtx",&tagger_info.mcflux_gen2vtx);

      tree0->SetBranchAddress("truth_nuScatType",&tagger_info.truth_nuScatType);

      tree0->SetBranchAddress("truth_nu_pos", &tagger_info.truth_nu_pos[0]);
      tree0->SetBranchAddress("truth_nu_momentum", &tagger_info.truth_nu_momentum[0]);
    }
  }

  if (flag==1){
    tree0->SetBranchAddress("nuvtx_diff", &tagger_info.nuvtx_diff);
    tree0->SetBranchAddress("showervtx_diff", &tagger_info.showervtx_diff);
    tree0->SetBranchAddress("muonvtx_diff", &tagger_info.muonvtx_diff);
    tree0->SetBranchAddress("truth_corr_nuvtxX", &tagger_info.truth_corr_nuvtxX);
    tree0->SetBranchAddress("truth_corr_nuvtxY", &tagger_info.truth_corr_nuvtxY);
    tree0->SetBranchAddress("truth_corr_nuvtxZ", &tagger_info.truth_corr_nuvtxZ);
    tree0->SetBranchAddress("truth_corr_showervtxX", &tagger_info.truth_corr_showervtxX);
    tree0->SetBranchAddress("truth_corr_showervtxY", &tagger_info.truth_corr_showervtxY);
    tree0->SetBranchAddress("truth_corr_showervtxZ", &tagger_info.truth_corr_showervtxZ);
    tree0->SetBranchAddress("truth_showerKE", &tagger_info.truth_showerKE);
    //Erin
    if (tree0->GetBranch("truth_single_photon")){
      tree0->SetBranchAddress("truth_showerPdg", &tagger_info.truth_showerPdg);
      tree0->SetBranchAddress("truth_showerMother", &tagger_info.truth_showerMother);
    }
    //
    tree0->SetBranchAddress("truth_corr_muonvtxX", &tagger_info.truth_corr_muonvtxX);
    tree0->SetBranchAddress("truth_corr_muonvtxY", &tagger_info.truth_corr_muonvtxY);
    tree0->SetBranchAddress("truth_corr_muonvtxZ", &tagger_info.truth_corr_muonvtxZ);
    tree0->SetBranchAddress("truth_muonvtxX", &tagger_info.truth_muonvtxX);
    tree0->SetBranchAddress("truth_muonvtxY", &tagger_info.truth_muonvtxY);
    tree0->SetBranchAddress("truth_muonvtxZ", &tagger_info.truth_muonvtxZ);
    tree0->SetBranchAddress("truth_muonendX", &tagger_info.truth_muonendX);
    tree0->SetBranchAddress("truth_muonendY", &tagger_info.truth_muonendY);
    tree0->SetBranchAddress("truth_muonendZ", &tagger_info.truth_muonendZ);
    tree0->SetBranchAddress("truth_muonMomentum", &tagger_info.truth_muonMomentum[0]);
    tree0->SetBranchAddress("truth_nuEnergy", &tagger_info.truth_nuEnergy);
    tree0->SetBranchAddress("truth_energyInside", &tagger_info.truth_energyInside);
    tree0->SetBranchAddress("truth_electronInside", &tagger_info.truth_electronInside);
    tree0->SetBranchAddress("truth_nuPdg", &tagger_info.truth_nuPdg);
    tree0->SetBranchAddress("truth_isCC", &tagger_info.truth_isCC);
    tree0->SetBranchAddress("truth_vtxX", &tagger_info.truth_vtxX);
    tree0->SetBranchAddress("truth_vtxY", &tagger_info.truth_vtxY);
    tree0->SetBranchAddress("truth_vtxZ", &tagger_info.truth_vtxZ);
    tree0->SetBranchAddress("truth_nuTime", &tagger_info.truth_nuTime);
    tree0->SetBranchAddress("truth_nuIntType", &tagger_info.truth_nuIntType);

    if (tree0->GetBranch("truth_NCDelta")){
      tagger_info.flag_NCDelta = true;
      tree0->SetBranchAddress("truth_NCDelta",&tagger_info.truth_NCDelta);
      tree0->SetBranchAddress("truth_NprimPio",&tagger_info.truth_NprimPio);
      tree0->SetBranchAddress("truth_pio_energy_1",&tagger_info.truth_pio_energy_1);
      tree0->SetBranchAddress("truth_pio_energy_2",&tagger_info.truth_pio_energy_2);
      tree0->SetBranchAddress("truth_pio_angle",&tagger_info.truth_pio_angle);
      //tree0->SetBranchAddress("reco_protonMomentum",&tagger_info.reco_protonMomentum[0]);
    }
    //Erin
    if (tree0->GetBranch("truth_single_photon")){
      tagger_info.flag_single_photon = true;
      tree0->SetBranchAddress("truth_Npi0",&tagger_info.truth_Npi0);
      tree0->SetBranchAddress("truth_single_photon",&tagger_info.truth_single_photon);
      tree0->SetBranchAddress("truth_photon_angle",&tagger_info.truth_photon_angle);
      tree0->SetBranchAddress("truth_photon_dis",&tagger_info.truth_photon_dis);
    }
  }

  if (tree0->GetBranch("evtTimeNS")){
    tagger_info.flag_nsbeam = true;
    tree0->SetBranchAddress("evtDeltaTimeNS",&tagger_info.evtDeltaTimeNS);
    tree0->SetBranchAddress("evtTimeNS",&tagger_info.evtTimeNS);
  }
  //

  if (tree0->GetBranch("reco_protonMomentum")){
    tagger_info.flag_recoprotonMomentum = true;
    tree0->SetBranchAddress("reco_protonMomentum",&tagger_info.reco_protonMomentum[0]);
  }


  if (tree0->GetBranch("truth_Ntrack")){
    tagger_info.flag_pf_truth = true;
    if (!tagger_info.flag_init_pointers){
      init_pointers(tagger_info);
      tagger_info.flag_init_pointers = true;
    }

    tree0->SetBranchAddress("truth_Ntrack", &tagger_info.truth_Ntrack);
    tree0->SetBranchAddress("truth_id", &tagger_info.truth_id);
    tree0->SetBranchAddress("truth_pdg", &tagger_info.truth_pdg);
    tree0->SetBranchAddress("truth_process", &tagger_info.truth_process);
    tree0->SetBranchAddress("truth_mother", &tagger_info.truth_mother);
    tree0->SetBranchAddress("truth_startXYZT", &tagger_info.truth_startXYZT);
    tree0->SetBranchAddress("truth_endXYZT", &tagger_info.truth_endXYZT);
    tree0->SetBranchAddress("truth_startMomentum", &tagger_info.truth_startMomentum);
    tree0->SetBranchAddress("truth_endMomentum", &tagger_info.truth_endMomentum);
    tree0->SetBranchAddress("truth_daughters", &tagger_info.truth_daughters);
    tree0->SetBranchAddress("fMC_trackPosition", &tagger_info.fMC_trackPosition);

    tree0->SetBranchAddress("mc_isnu", &tagger_info.mc_isnu);
    tree0->SetBranchAddress("mc_nGeniePrimaries", &tagger_info.mc_nGeniePrimaries);
    tree0->SetBranchAddress("mc_nu_pdg", &tagger_info.mc_nu_pdg);
    tree0->SetBranchAddress("mc_nu_ccnc", &tagger_info.mc_nu_ccnc);
    tree0->SetBranchAddress("mc_nu_mode", &tagger_info.mc_nu_mode);
    tree0->SetBranchAddress("mc_nu_intType", &tagger_info.mc_nu_intType);
    tree0->SetBranchAddress("mc_nu_target", &tagger_info.mc_nu_target);
    tree0->SetBranchAddress("mc_hitnuc", &tagger_info.mc_hitnuc);
    tree0->SetBranchAddress("mc_hitquark", &tagger_info.mc_hitquark);
    tree0->SetBranchAddress("mc_nu_Q2", &tagger_info.mc_nu_Q2);
    tree0->SetBranchAddress("mc_nu_W", &tagger_info.mc_nu_W);
    tree0->SetBranchAddress("mc_nu_X", &tagger_info.mc_nu_X);
    tree0->SetBranchAddress("mc_nu_Y", &tagger_info.mc_nu_Y);
    tree0->SetBranchAddress("mc_nu_Pt", &tagger_info.mc_nu_Pt);
    tree0->SetBranchAddress("mc_nu_Theta", &tagger_info.mc_nu_Theta);
    tree0->SetBranchAddress("mc_nu_pos", &tagger_info.mc_nu_pos);
    tree0->SetBranchAddress("mc_nu_mom", &tagger_info.mc_nu_mom);

  }

  if (tree0->GetBranch("reco_Ntrack")){
    tagger_info.flag_pf_reco = true;
    if (!tagger_info.flag_init_pointers){
      init_pointers(tagger_info);
      tagger_info.flag_init_pointers = true;
    }

    tree0->SetBranchAddress("reco_Ntrack", &tagger_info.reco_Ntrack);
    tree0->SetBranchAddress("reco_id", &tagger_info.reco_id);
    tree0->SetBranchAddress("reco_pdg", &tagger_info.reco_pdg);
    tree0->SetBranchAddress("reco_process", &tagger_info.reco_process);
    tree0->SetBranchAddress("reco_mother", &tagger_info.reco_mother);
    tree0->SetBranchAddress("reco_startXYZT", &tagger_info.reco_startXYZT);
    tree0->SetBranchAddress("reco_endXYZT", &tagger_info.reco_endXYZT);
    tree0->SetBranchAddress("reco_startMomentum", &tagger_info.reco_startMomentum);
    tree0->SetBranchAddress("reco_endMomentum", &tagger_info.reco_endMomentum);
    tree0->SetBranchAddress("reco_daughters", &tagger_info.reco_daughters);
  }



  //new mcs vars
  if (tree0->GetBranch("mcs_mu_tracklen")){
    tagger_info.flag_mcs = true;
    tree0->SetBranchAddress("mcs_mu_tracklen",&tagger_info.mcs_mu_tracklen);
    tree0->SetBranchAddress("mcs_emu_tracklen",&tagger_info.mcs_emu_tracklen);
    tree0->SetBranchAddress("mcs_emu_MCS",&tagger_info.mcs_emu_MCS);
    tree0->SetBranchAddress("mcs_ambiguity_MCS",&tagger_info.mcs_ambiguity_MCS);
  }

  //new larpid and backtracking variables
  if (tree0->GetBranch("reco_larpid_classified")){
    tagger_info.flag_larpid = true;
    tree0->SetBranchAddress("reco_larpid_classified",&tagger_info.reco_larpid_classified);
    tree0->SetBranchAddress("reco_larpid_pdg",&tagger_info.reco_larpid_pdg);
    tree0->SetBranchAddress("reco_larpid_proccess",&tagger_info.reco_larpid_proccess);
    tree0->SetBranchAddress("reco_larpid_completeness",&tagger_info.reco_larpid_completeness);
    tree0->SetBranchAddress("reco_larpid_purity",&tagger_info.reco_larpid_purity);
    tree0->SetBranchAddress("reco_larpid_pidScore_el",&tagger_info.reco_larpid_pidScore_el);
    tree0->SetBranchAddress("reco_larpid_pidScore_ph",&tagger_info.reco_larpid_pidScore_ph);
    tree0->SetBranchAddress("reco_larpid_pidScore_mu",&tagger_info.reco_larpid_pidScore_mu);
    tree0->SetBranchAddress("reco_larpid_pidScore_pi",&tagger_info.reco_larpid_pidScore_pi);
    tree0->SetBranchAddress("reco_larpid_pidScore_pr",&tagger_info.reco_larpid_pidScore_pr);
    tree0->SetBranchAddress("reco_larpid_procScore_prim",&tagger_info.reco_larpid_procScore_prim);
    tree0->SetBranchAddress("reco_larpid_procScore_ntrl",&tagger_info.reco_larpid_procScore_ntrl);
    tree0->SetBranchAddress("reco_larpid_procScore_chgd",&tagger_info.reco_larpid_procScore_chgd);
  }

  if (tree0->GetBranch("reco_truthMatch_pdg")){
    tagger_info.flag_backtracking = true;
    tree0->SetBranchAddress("reco_truthMatch_pdg",&tagger_info.reco_truthMatch_pdg);
    tree0->SetBranchAddress("reco_truthMatch_id",&tagger_info.reco_truthMatch_id);
    tree0->SetBranchAddress("reco_truthMatch_purity",&tagger_info.reco_truthMatch_purity);
    tree0->SetBranchAddress("reco_truthMatch_completeness",&tagger_info.reco_truthMatch_completeness);
    tree0->SetBranchAddress("reco_truthMatch_nSimParts",&tagger_info.reco_truthMatch_nSimParts);
    tree0->SetBranchAddress("reco_truthMatch_simPart_pdg",&tagger_info.reco_truthMatch_simPart_pdg);
    tree0->SetBranchAddress("reco_truthMatch_simPart_purity",&tagger_info.reco_truthMatch_simPart_purity);
  }

  //new ns timing vars
  if (tree0->GetBranch("evtTimeNS_cor")){
    tagger_info.flag_ns_time_cor = true;
    tree0->SetBranchAddress("evtTimeNS_cor",&tagger_info.evtTimeNS_cor);
    tree0->SetBranchAddress("cor_nu_time",&tagger_info.cor_nu_time);
    tree0->SetBranchAddress("cor_nu_time_nospill",&tagger_info.cor_nu_time_nospill);
    tree0->SetBranchAddress("cor_nu_time_spill",&tagger_info.cor_nu_time_spill);
    tree0->SetBranchAddress("cor_nu_deltatime",&tagger_info.cor_nu_deltatime);
  }
  if (tree0->GetBranch("Ph_Tot")){
    tagger_info.flag_Phtot = true;
    tree0->SetBranchAddress("Ph_Tot",&tagger_info.Ph_Tot);
  }
  if (tree0->GetBranch("PMT_ID")){
    tagger_info.flag_PMT = true;
    tree0->SetBranchAddress("PMT_ID",&tagger_info.PMT_ID);
    tree0->SetBranchAddress("PMT_Time",&tagger_info.PMT_Time);
    tree0->SetBranchAddress("PMT_Amp",&tagger_info.PMT_Amp);
    tree0->SetBranchAddress("PMT_TimeProp",&tagger_info.PMT_TimeProp);
    tree0->SetBranchAddress("PMT_TimeDP",&tagger_info.PMT_TimeDP);
    tree0->SetBranchAddress("PMT_TimeDL",&tagger_info.PMT_TimeDL);
    tree0->SetBranchAddress("PMT_Sat",&tagger_info.PMT_Sat);
    tree0->SetBranchAddress("RWM_Time",&tagger_info.RWM_Time);
  }


}

void LEEana::put_tree_address(TTree *tree0, PFevalInfo& tagger_info, int flag){
  tree0->Branch("run", &tagger_info.run, "run/I");
  tree0->Branch("subrun", &tagger_info.subrun,"subrun/I");
  tree0->Branch("event", &tagger_info.event,"event/I");
  tree0->Branch("neutrino_type", &tagger_info.neutrino_type,"neutrino_type/I");
  tree0->Branch("reco_nuvtxX", &tagger_info.reco_nuvtxX,"reco_nuvtxX/F");
  tree0->Branch("reco_nuvtxY", &tagger_info.reco_nuvtxY,"reco_nuvtxY/F");
  tree0->Branch("reco_nuvtxZ", &tagger_info.reco_nuvtxZ,"reco_nuvtxZ/F");
  tree0->Branch("reco_showervtxX", &tagger_info.reco_showervtxX,"reco_showervtxX/F");
  tree0->Branch("reco_showervtxY", &tagger_info.reco_showervtxY,"reco_showervtxY/F");
  tree0->Branch("reco_showervtxZ", &tagger_info.reco_showervtxZ,"reco_showervtxZ/F");
  tree0->Branch("reco_showerKE", &tagger_info.reco_showerKE,"reco_showerKE/F");
  tree0->Branch("reco_muonvtxX", &tagger_info.reco_muonvtxX,"reco_muonvtxX/F");
  tree0->Branch("reco_muonvtxY", &tagger_info.reco_muonvtxY,"reco_muonvtxY/F");
  tree0->Branch("reco_muonvtxZ", &tagger_info.reco_muonvtxZ,"reco_muonvtxZ/F");
  tree0->Branch("reco_muonMomentum", &tagger_info.reco_muonMomentum[0],"reco_muonMomentum[4]/F");

  if (flag==1){
    tree0->Branch("nuvtx_diff", &tagger_info.nuvtx_diff,"nuvtx_diff/F");
    tree0->Branch("showervtx_diff", &tagger_info.showervtx_diff,"showervtx_diff/F");
    tree0->Branch("muonvtx_diff", &tagger_info.muonvtx_diff,"data/F");
    tree0->Branch("truth_corr_nuvtxX", &tagger_info.truth_corr_nuvtxX,"data/F");
    tree0->Branch("truth_corr_nuvtxY", &tagger_info.truth_corr_nuvtxY,"data/F");
    tree0->Branch("truth_corr_nuvtxZ", &tagger_info.truth_corr_nuvtxZ,"data/F");
    tree0->Branch("truth_corr_showervtxX", &tagger_info.truth_corr_showervtxX,"data/F");
    tree0->Branch("truth_corr_showervtxY", &tagger_info.truth_corr_showervtxY,"data/F");
    tree0->Branch("truth_corr_showervtxZ", &tagger_info.truth_corr_showervtxZ,"data/F");
    tree0->Branch("truth_showerKE", &tagger_info.truth_showerKE,"data/F");
    //Erin
    if (tagger_info.flag_single_photon){
      tree0->Branch("truth_showerPdg", &tagger_info.truth_showerPdg,"data/I");
      tree0->Branch("truth_showerMother", &tagger_info.truth_showerMother,"data/I");
    }
    //
    tree0->Branch("truth_corr_muonvtxX", &tagger_info.truth_corr_muonvtxX,"data/F");
    tree0->Branch("truth_corr_muonvtxY", &tagger_info.truth_corr_muonvtxY,"data/F");
    tree0->Branch("truth_corr_muonvtxZ", &tagger_info.truth_corr_muonvtxZ,"data/F");
    tree0->Branch("truth_muonvtxX", &tagger_info.truth_muonvtxX,"data/F");
    tree0->Branch("truth_muonvtxY", &tagger_info.truth_muonvtxY,"data/F");
    tree0->Branch("truth_muonvtxZ", &tagger_info.truth_muonvtxZ,"data/F");
    tree0->Branch("truth_muonendX", &tagger_info.truth_muonendX,"data/F");
    tree0->Branch("truth_muonendY", &tagger_info.truth_muonendY,"data/F");
    tree0->Branch("truth_muonendZ", &tagger_info.truth_muonendZ,"data/F");
    tree0->Branch("truth_muonMomentum", &tagger_info.truth_muonMomentum[0],"truth_muonMomentum[4]/F");
    tree0->Branch("truth_nuEnergy", &tagger_info.truth_nuEnergy,"data/F");
    tree0->Branch("truth_energyInside", &tagger_info.truth_energyInside,"data/F");
    tree0->Branch("truth_electronInside", &tagger_info.truth_electronInside,"data/F");
    tree0->Branch("truth_nuPdg", &tagger_info.truth_nuPdg,"data/I");
    tree0->Branch("truth_isCC", &tagger_info.truth_isCC,"data/O");
    tree0->Branch("truth_vtxX", &tagger_info.truth_vtxX,"data/F");
    tree0->Branch("truth_vtxY", &tagger_info.truth_vtxY,"data/F");
    tree0->Branch("truth_vtxZ", &tagger_info.truth_vtxZ,"data/F");
    tree0->Branch("truth_nuTime", &tagger_info.truth_nuTime,"data/F");
    tree0->Branch("truth_nuIntType", &tagger_info.truth_nuIntType,"data/I");

    if (tagger_info.flag_NCDelta){
      tree0->Branch("truth_NCDelta",&tagger_info.truth_NCDelta,"truth_NCDelta/I");
      tree0->Branch("truth_NprimPio",&tagger_info.truth_NprimPio,"truth_NprimPio/I");
      tree0->Branch("truth_pio_energy_1",&tagger_info.truth_pio_energy_1,"truth_pio_energy_1/F");
      tree0->Branch("truth_pio_energy_2",&tagger_info.truth_pio_energy_2,"truth_pio_energy_2/F");
      tree0->Branch("truth_pio_angle",&tagger_info.truth_pio_angle,"truth_pio_angle/F");
      //tree0->Branch("reco_protonMomentum",&tagger_info.reco_protonMomentum[0],"reco_protonMomentum[4]/F");
    }

    //Erin
    if (tagger_info.flag_single_photon){
      tree0->Branch("truth_Npi0",&tagger_info.truth_Npi0,"truth_Npi0/I");
      tree0->Branch("truth_single_photon",&tagger_info.truth_single_photon, "truth_single_photon/I");
      tree0->Branch("truth_photon_angle",&tagger_info.truth_photon_angle, "truth_photon_angle/F");
      tree0->Branch("truth_photon_dis",&tagger_info.truth_photon_dis, "truth_photon_dis/F");
    }

  }

  if (tagger_info.flag_nsbeam){
    tree0->Branch("evtDeltaTimeNS",&tagger_info.evtDeltaTimeNS,"evtDeltaTimeNS/F");
    tree0->Branch("evtTimeNS",&tagger_info.evtTimeNS,"evtTimeNS/F");
  }
  //

  if (tagger_info.flag_recoprotonMomentum){
    tree0->Branch("reco_protonMomentum",&tagger_info.reco_protonMomentum[0],"reco_protonMomentum[4]/F");
  }

  if (tagger_info.flag_showerMomentum){

    tree0->Branch("reco_showerMomentum",&tagger_info.reco_showerMomentum[0],"reco_showerMomentum[4]/F");
    tree0->Branch("reco_Nproton",&tagger_info.reco_Nproton,"reco_Nproton/I");


    if (flag==1){
      tree0->Branch("truth_showerMomentum",&tagger_info.truth_showerMomentum[0],"truth_showerMomentum[4]/F");
      tree0->Branch("mcflux_run",&tagger_info.mcflux_run,"mcflux_run/I");
      tree0->Branch("mcflux_evtno",&tagger_info.mcflux_evtno,"mcflux_evtno/I");
      tree0->Branch("mcflux_ndecay",&tagger_info.mcflux_ndecay,"mcflux_ndecay/I");
      tree0->Branch("mcflux_ntype",&tagger_info.mcflux_ntype,"mcflux_ntype/I");
      tree0->Branch("truth_nuScatType",&tagger_info.truth_nuScatType,"truth_nuScatType/I");

      tree0->Branch("mcflux_nuEnergy",&tagger_info.mcflux_nuEnergy,"mcflux_nuEnergy/F");
      tree0->Branch("mcflux_vx", &tagger_info.mcflux_vx,"mcflux_vx/F");
      tree0->Branch("mcflux_vy", &tagger_info.mcflux_vy,"mcflux_vy/F");
      tree0->Branch("mcflux_vz", &tagger_info.mcflux_vz,"mcflux_vz/F");
      tree0->Branch("mcflux_genx", &tagger_info.mcflux_genx, "mcflux_genx/F");
      tree0->Branch("mcflux_geny", &tagger_info.mcflux_geny, "mcflux_geny/F");
      tree0->Branch("mcflux_genz", &tagger_info.mcflux_genz, "mcflux_genz/F");
      tree0->Branch("mcflux_dk2gen", &tagger_info.mcflux_dk2gen, "mcflux_dk2gen/F");
      tree0->Branch("mcflux_gen2vtx", &tagger_info.mcflux_gen2vtx, "mcflux_gen2vtx/F");
      tree0->Branch("truth_nu_pos", &tagger_info.truth_nu_pos[0],"truth_nu_pos[4]/F");
      tree0->Branch("truth_nu_momentum", &tagger_info.truth_nu_momentum[0], "truth_nu_momentum[4]/F");
    }
  }


  if (tagger_info.flag_pf_truth){

    tree0->Branch("truth_Ntrack", &tagger_info.truth_Ntrack);
    tree0->Branch("truth_id", &tagger_info.truth_id, "truth_id[truth_Ntrack]/I");
    tree0->Branch("truth_pdg", &tagger_info.truth_pdg, "truth_pdg[truth_Ntrack]/I");
    tree0->Branch("truth_process", &tagger_info.truth_process);
    tree0->Branch("truth_mother", &tagger_info.truth_mother, "truth_mother[truth_Ntrack]/I");
    tree0->Branch("truth_startXYZT", &tagger_info.truth_startXYZT, "truth_startXYZT[truth_Ntrack][4]/F");
    tree0->Branch("truth_endXYZT", &tagger_info.truth_endXYZT, "truth_endXYZT[truth_Ntrack][4]/F");
    tree0->Branch("truth_startMomentum", &tagger_info.truth_startMomentum, "truth_startMomentum[truth_Ntrack][4]/F");
    tree0->Branch("truth_endMomentum", &tagger_info.truth_endMomentum, "truth_endMomentum[truth_Ntrack][4]/F");
    tree0->Branch("truth_daughters", tagger_info.truth_daughters);
    tree0->Branch("fMC_trackPosition", &tagger_info.fMC_trackPosition);



    tree0->Branch("mc_isnu", &tagger_info.mc_isnu);
    tree0->Branch("mc_nGeniePrimaries", &tagger_info.mc_nGeniePrimaries);
    tree0->Branch("mc_nu_pdg", &tagger_info.mc_nu_pdg);
    tree0->Branch("mc_nu_ccnc", &tagger_info.mc_nu_ccnc);
    tree0->Branch("mc_nu_mode", &tagger_info.mc_nu_mode);
    tree0->Branch("mc_nu_intType", &tagger_info.mc_nu_intType);
    tree0->Branch("mc_nu_target", &tagger_info.mc_nu_target);
    tree0->Branch("mc_hitnuc", &tagger_info.mc_hitnuc);
    tree0->Branch("mc_hitquark", &tagger_info.mc_hitquark);
    tree0->Branch("mc_nu_Q2", &tagger_info.mc_nu_Q2);
    tree0->Branch("mc_nu_W", &tagger_info.mc_nu_W);
    tree0->Branch("mc_nu_X", &tagger_info.mc_nu_X);
    tree0->Branch("mc_nu_Y", &tagger_info.mc_nu_Y);
    tree0->Branch("mc_nu_Pt", &tagger_info.mc_nu_Pt);
    tree0->Branch("mc_nu_Theta", &tagger_info.mc_nu_Theta);
    tree0->Branch("mc_nu_pos", &tagger_info.mc_nu_pos, "mc_nu_pos[4]/F");
    tree0->Branch("mc_nu_mom", &tagger_info.mc_nu_mom, "mc_nu_mom[4]/F");
  }

  if (tagger_info.flag_pf_reco){
    tree0->Branch("reco_Ntrack", &tagger_info.reco_Ntrack);
    tree0->Branch("reco_id", &tagger_info.reco_id, "reco_id[reco_Ntrack]/I");
    tree0->Branch("reco_pdg", &tagger_info.reco_pdg, "reco_pdg[reco_Ntrack]/I");
    tree0->Branch("reco_process", &tagger_info.reco_process);
    tree0->Branch("reco_mother", &tagger_info.reco_mother, "reco_mother[reco_Ntrack]/I");
    tree0->Branch("reco_startXYZT", &tagger_info.reco_startXYZT, "reco_startXYZT[reco_Ntrack][4]/F");
    tree0->Branch("reco_endXYZT", &tagger_info.reco_endXYZT, "reco_endXYZT[reco_Ntrack][4]/F");
    tree0->Branch("reco_startMomentum", &tagger_info.reco_startMomentum, "reco_startMomentum[reco_Ntrack][4]/F");
    tree0->Branch("reco_endMomentum", &tagger_info.reco_endMomentum, "reco_endMomentum[reco_Ntrack][4]/F");
    tree0->Branch("reco_daughters", tagger_info.reco_daughters);
  }

  
  //new mcs vars
  
  if (tagger_info.flag_mcs){
    tree0->Branch("mcs_mu_tracklen",&tagger_info.mcs_mu_tracklen,"mcs_mu_tracklen/D");
    tree0->Branch("mcs_emu_tracklen",&tagger_info.mcs_emu_tracklen,"mcs_emu_tracklen/D");
    tree0->Branch("mcs_emu_MCS",&tagger_info.mcs_emu_MCS,"mcs_emu_MCS/D");
    tree0->Branch("mcs_ambiguity_MCS",&tagger_info.mcs_ambiguity_MCS,"mcs_ambiguity_MCS/D");
  }

  //new larpid and backtracking variables
  if (tagger_info.flag_larpid){
    tree0->Branch("reco_larpid_classified", &tagger_info.reco_larpid_classified, "reco_larpid_classified[reco_Ntrack]/I");
    tree0->Branch("reco_larpid_pdg", &tagger_info.reco_larpid_pdg, "reco_larpid_pdg[reco_Ntrack]/I");
    tree0->Branch("reco_larpid_proccess", &tagger_info.reco_larpid_proccess, "reco_larpid_proccess[reco_Ntrack]/I");
    tree0->Branch("reco_larpid_completeness", &tagger_info.reco_larpid_completeness, "reco_larpid_completeness[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_purity", &tagger_info.reco_larpid_purity, "reco_larpid_purity[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_pidScore_el", &tagger_info.reco_larpid_pidScore_el, "reco_larpid_pidScore_el[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_pidScore_ph", &tagger_info.reco_larpid_pidScore_ph, "reco_larpid_pidScore_ph[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_pidScore_mu", &tagger_info.reco_larpid_pidScore_mu, "reco_larpid_pidScore_mu[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_pidScore_pi", &tagger_info.reco_larpid_pidScore_pi, "reco_larpid_pidScore_pi[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_pidScore_pr", &tagger_info.reco_larpid_pidScore_pr, "reco_larpid_pidScore_pr[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_procScore_prim", &tagger_info.reco_larpid_procScore_prim, "reco_larpid_procScore_prim[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_procScore_ntrl", &tagger_info.reco_larpid_procScore_ntrl, "reco_larpid_procScore_ntrl[reco_Ntrack]/F");
    tree0->Branch("reco_larpid_procScore_chgd", &tagger_info.reco_larpid_procScore_chgd, "reco_larpid_procScore_chgd[reco_Ntrack]/F");
  }

  if (tagger_info.flag_backtracking){
    tree0->Branch("reco_truthMatch_pdg", &tagger_info.reco_truthMatch_pdg, "reco_truthMatch_pdg[reco_Ntrack]/I");
    tree0->Branch("reco_truthMatch_id", &tagger_info.reco_truthMatch_id, "reco_truthMatch_id[reco_Ntrack]/I");
    tree0->Branch("reco_truthMatch_purity", &tagger_info.reco_truthMatch_purity, "reco_truthMatch_purity[reco_Ntrack]/F");
    tree0->Branch("reco_truthMatch_completeness", &tagger_info.reco_truthMatch_completeness, "reco_truthMatch_completeness[reco_Ntrack]/F");
    tree0->Branch("reco_truthMatch_nSimParts", &tagger_info.reco_truthMatch_nSimParts, "reco_truthMatch_nSimParts[reco_Ntrack]/I");
    tree0->Branch("reco_truthMatch_simPart_pdg", &tagger_info.reco_truthMatch_simPart_pdg);
    tree0->Branch("reco_truthMatch_simPart_purity", &tagger_info.reco_truthMatch_simPart_purity);
  }

  //new ns timing vars
  if (tagger_info.flag_ns_time_cor){
    tree0->Branch("evtTimeNS_cor",&tagger_info.evtTimeNS_cor,"evtTimeNS_cor/F");
    tree0->Branch("cor_nu_time",&tagger_info.cor_nu_time,"cor_nu_time/F");
    tree0->Branch("cor_nu_time_nospill",&tagger_info.cor_nu_time_nospill,"cor_nu_time_nospill/F");
    tree0->Branch("cor_nu_time_spill",&tagger_info.cor_nu_time_spill,"cor_nu_time_spill/F");
    tree0->Branch("cor_nu_deltatime",&tagger_info.cor_nu_deltatime,"cor_nu_deltatime/F");
  }
  if (tagger_info.flag_Phtot){
    tree0->Branch("Ph_Tot",&tagger_info.Ph_Tot,"Ph_Tot/F");
  }
  if (tagger_info.flag_PMT){
    tree0->Branch("PMT_ID",&tagger_info.PMT_ID);
    tree0->Branch("PMT_Time",&tagger_info.PMT_Time);
    tree0->Branch("PMT_Amp",&tagger_info.PMT_Amp);
    tree0->Branch("PMT_TimeProp",&tagger_info.PMT_TimeProp);
    tree0->Branch("PMT_TimeDP",&tagger_info.PMT_TimeDP);
    tree0->Branch("PMT_TimeDL",&tagger_info.PMT_TimeDL);
    tree0->Branch("PMT_Sat",&tagger_info.PMT_Sat);
    tree0->Branch("RWM_Time",&tagger_info.RWM_Time,"RWM_Time/F");
  }
  


  

}



#endif
