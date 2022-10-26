void print_knob(){
  auto file = TFile::Open("/home/erin/Documents/MicroBoone/spbdt/data/run3_checkout/checkout_bnboverlay_run3_sp_1.root");
  auto t = (TTree*)file->Get("wcpweights/T_wgt");
  // t->Show(0);
  // t->Print();
  gInterpreter->GenerateDictionary("map<string,vector<float> >", "map");



  auto mcweight = new map<string,vector<float> >;
  t->SetBranchAddress("mcweight", &mcweight);

  int n=0;
  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if (mcweight->size()>0 && n==0) {
      cout << mcweight->size() << endl;
      n ++;

      for(auto const& e: *mcweight){
        cout << e.first << " " << e.second.size() << endl;
      }

    }
  }




}
