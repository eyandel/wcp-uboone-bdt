#ifndef UBOONE_LEE_LANTERN
#define UBOONE_LEE_LANTERN

namespace LEEana{
struct LanternInfo{
  Int_t run;
  Int_t subrun;
  Int_t event;
 
  Int_t nTracks;

  Int_t trackIsSecondary[10000];
  Int_t trackPID[10000];
  Float_t trackMuScore[10000];
  Float_t trackPrScore[10000];
  Float_t trackPiScore[10000];
  Float_t trackElScore[10000];
  Float_t trackPhScore[10000];
  Float_t trackRecoE[10000];
  Float_t trackDistToVtx[10000];

};

 void set_tree_address(TTree *tree0, LanternInfo& lantern_info);
 void put_tree_address(TTree *tree0, LanternInfo& lantern_info);

}


void LEEana::set_tree_address(TTree *tree0, LanternInfo& lantern_info) {
  tree0->SetBranchAddress("run", &lantern_info.run);
  tree0->SetBranchAddress("subrun", &lantern_info.subrun);
  tree0->SetBranchAddress("event", &lantern_info.event);

  tree0->SetBranchAddress("nTracks", &lantern_info.nTracks);

  tree0->SetBranchAddress("trackIsSecondary", &lantern_info.trackIsSecondary);
  tree0->SetBranchAddress("trackPID", &lantern_info.trackPID);
  tree0->SetBranchAddress("trackMuScore", &lantern_info.trackMuScore);
  tree0->SetBranchAddress("trackPrScore", &lantern_info.trackPrScore);
  tree0->SetBranchAddress("trackPiScore", &lantern_info.trackPiScore);
  tree0->SetBranchAddress("trackElScore", &lantern_info.trackElScore);
  tree0->SetBranchAddress("trackPhScore", &lantern_info.trackPhScore);
  tree0->SetBranchAddress("trackRecoE", &lantern_info.trackRecoE);
  tree0->SetBranchAddress("trackDistToVtx", &lantern_info.trackDistToVtx);
}

void LEEana::put_tree_address(TTree *tree0, LanternInfo& lantern_info) {
  tree0->Branch("run", &lantern_info.run, "run/I");
  tree0->Branch("subrun", &lantern_info.subrun, "subrun/I");
  tree0->Branch("event", &lantern_info.event, "event/I");

  tree0->Branch("nTracks", &lantern_info.nTracks);

  tree0->Branch("trackIsSecondary", &lantern_info.trackIsSecondary);
  tree0->Branch("trackPID", &lantern_info.trackPID);
  tree0->Branch("trackMuScore", &lantern_info.trackMuScore);
  tree0->Branch("trackPrScore", &lantern_info.trackPrScore);
  tree0->Branch("trackPiScore", &lantern_info.trackPiScore);
  tree0->Branch("trackElScore", &lantern_info.trackElScore);
  tree0->Branch("trackPhScore", &lantern_info.trackPhScore);
  tree0->Branch("trackRecoE", &lantern_info.trackRecoE);
  tree0->Branch("trackDistToVtx", &lantern_info.trackDistToVtx);
}

#endif
