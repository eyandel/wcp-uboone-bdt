#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "WCPLEEANA/TLee.h"

#include "WCPLEEANA/Configure_Lee.h"

#include "TApplication.h"

#include <chrono>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
  usage:

  make clean
  make
  ./read_TLee_v20 -f 1 -p 1

  ---> README:
  ---> Makefile: comment the line "ROOTSYS=/home/xji/data0/software/root_build", if you have your own "ROOTSYS"
  ---> minuit2 is in the ROOT
*/

int main(int argc, char** argv)
{
  TString roostr = "";

  double scaleF_POT = 1;
  int ifile = 1;

  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-p")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
  }

  cout<<endl<<" ---> check, scaleF_POT "<<scaleF_POT<<", ifile "<<ifile<<endl<<endl;

  //////////////////////////////////////////////////////////////////////////////////////// Draw style

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  if( !config_Lee::flag_display_graphics ) {
    gROOT->SetBatch( 1 );
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  TApplication theApp("theApp",&argc,argv);

  TLee *Lee_test = new TLee();

  ////////// just do it one time in the whole procedure

  Lee_test->channels_observation   = config_Lee::channels_observation;
  Lee_test->syst_cov_flux_Xs_begin = config_Lee::syst_cov_flux_Xs_begin;
  Lee_test->syst_cov_flux_Xs_end   = config_Lee::syst_cov_flux_Xs_end;
  Lee_test->syst_cov_mc_stat_begin = config_Lee::syst_cov_mc_stat_begin;
  Lee_test->syst_cov_mc_stat_end   = config_Lee::syst_cov_mc_stat_end;

  Lee_test->flag_lookelsewhere     = config_Lee::flag_lookelsewhere;

  Lee_test->Set_config_file_directory(config_Lee::spectra_file, config_Lee::flux_Xs_directory,
                                      config_Lee::detector_directory, config_Lee::mc_directory);

  int size_array_LEE_ch = sizeof(config_Lee::array_LEE_ch)/sizeof(config_Lee::array_LEE_ch[0]);
  for(int idx=0; idx<size_array_LEE_ch; idx++) {
    if( config_Lee::array_LEE_ch[idx]!=0 ) Lee_test->map_Lee_ch[config_Lee::array_LEE_ch[idx]] = 1;
  }

  int size_array_LEE_Np_ch = sizeof(config_Lee::array_LEE_Np_ch)/sizeof(config_Lee::array_LEE_Np_ch[0]);
  for(int idx=0; idx<size_array_LEE_Np_ch; idx++) {
    if( config_Lee::array_LEE_Np_ch[idx]!=0 ) Lee_test->map_Lee_Np_ch[config_Lee::array_LEE_Np_ch[idx]] = 1;
  }

  int size_array_LEE_0p_ch = sizeof(config_Lee::array_LEE_0p_ch)/sizeof(config_Lee::array_LEE_0p_ch[0]);
  for(int idx=0; idx<size_array_LEE_0p_ch; idx++) {
    if( config_Lee::array_LEE_0p_ch[idx]!=0 ) Lee_test->map_Lee_0p_ch[config_Lee::array_LEE_0p_ch[idx]] = 1;
  }

  Lee_test->scaleF_POT = scaleF_POT;

  Lee_test->flag_syst_flux_Xs    = config_Lee::flag_syst_flux_Xs;
  Lee_test->flag_syst_detector   = config_Lee::flag_syst_detector;
  Lee_test->flag_syst_additional = config_Lee::flag_syst_additional;
  Lee_test->flag_syst_mc_stat    = config_Lee::flag_syst_mc_stat;

  Lee_test->Set_Spectra_MatrixCov();
  Lee_test->Set_POT_implement();
  Lee_test->Set_TransformMatrix();

  ////////// can do any times

  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
  Lee_test->scaleF_Lee_Np = config_Lee::Lee_Np_strength_for_outputfile_covariance_matrix;
  Lee_test->scaleF_Lee_0p = config_Lee::Lee_0p_strength_for_outputfile_covariance_matrix;

  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  Lee_test->scaleF_Lee_Np = config_Lee::Lee_Np_strength_for_GoF;
  Lee_test->scaleF_Lee_0p = config_Lee::Lee_0p_strength_for_GoF;

  Lee_test->Set_Collapse();

  //////////

  if( 0 ) {// shape only covariance

    int nbins = Lee_test->bins_newworld;
    TMatrixD matrix_pred = Lee_test->matrix_pred_newworld;
    TMatrixD matrix_syst = Lee_test->matrix_absolute_flux_cov_newworld;
    TMatrixD matrix_shape(nbins, nbins);
    TMatrixD matrix_mixed(nbins, nbins);
    TMatrixD matrix_norm(nbins, nbins);

    ///
    double N_T = 0;
    for(int idx=0; idx<nbins; idx++) N_T += matrix_pred(0, idx);

    ///
    double M_kl = 0;
    for(int i=0; i<nbins; i++) {
      for(int j=0; j<nbins; j++) {
        M_kl += matrix_syst(i,j);
      }
    }

    ///
    for(int i=0; i<nbins; i++) {
      for(int j=0; j<nbins; j++) {
        double N_i = matrix_pred(0, i);
        double N_j = matrix_pred(0, j);

        double M_ij = matrix_syst(i,j);
        double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
        double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);

        matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
        matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;
        matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;
      }
    }

    Lee_test->matrix_absolute_flux_cov_newworld = matrix_shape;

  }// shape only covariance

  /////////////////////////////
  /////////////////////////////

  TFile *file_collapsed_covariance_matrix = new TFile("file_collapsed_covariance_matrix.root", "recreate");

  TTree *tree_config = new TTree("tree", "configure information");
  int flag_syst_flux_Xs = config_Lee::flag_syst_flux_Xs;
  int flag_syst_detector = config_Lee::flag_syst_detector;
  int flag_syst_additional = config_Lee::flag_syst_additional;
  int flag_syst_mc_stat = config_Lee::flag_syst_mc_stat;
  double user_Lee_strength_for_output_covariance_matrix = config_Lee::Lee_strength_for_outputfile_covariance_matrix;
  double user_scaleF_POT = scaleF_POT;
  vector<double>vc_val_GOF;
  vector<int>vc_val_GOF_NDF;
  tree_config->Branch("flag_syst_flux_Xs", &flag_syst_flux_Xs, "flag_syst_flux_Xs/I" );
  tree_config->Branch("flag_syst_detector", &flag_syst_detector, "flag_syst_detector/I" );
  tree_config->Branch("flag_syst_additional", &flag_syst_additional, "flag_syst_additional/I" );
  tree_config->Branch("flag_syst_mc_stat", &flag_syst_mc_stat, "flag_syst_mc_stat/I" );
  tree_config->Branch("user_Lee_strength_for_output_covariance_matrix", &user_Lee_strength_for_output_covariance_matrix,
                      "user_Lee_strength_for_output_covariance_matrix/D" );
  tree_config->Branch("user_scaleF_POT", &user_scaleF_POT, "user_scaleF_POT/D" );
  tree_config->Branch("vc_val_GOF", &vc_val_GOF);
  tree_config->Branch("vc_val_GOF_NDF", &vc_val_GOF_NDF);
  file_collapsed_covariance_matrix->cd();

  Lee_test->matrix_absolute_cov_newworld.Write("matrix_absolute_cov_newworld");// (bins, bins)
  Lee_test->matrix_absolute_flux_cov_newworld.Write("matrix_absolute_flux_cov_newworld");
  Lee_test->matrix_absolute_Xs_cov_newworld.Write("matrix_absolute_Xs_cov_newworld");
  Lee_test->matrix_absolute_detector_cov_newworld.Write("matrix_absolute_detector_cov_newworld");
  Lee_test->matrix_absolute_mc_stat_cov_newworld.Write("matrix_absolute_mc_stat_cov_newworld");
  Lee_test->matrix_absolute_additional_cov_newworld.Write("matrix_absolute_additional_cov_newworld");

  for(auto it=Lee_test->matrix_input_cov_detector_sub.begin(); it!=Lee_test->matrix_input_cov_detector_sub.end(); it++) {
    int idx = it->first;
    roostr = TString::Format("matrix_absolute_detector_sub_cov_newworld_%02d", idx);
    Lee_test->matrix_absolute_detector_sub_cov_newworld[idx].Write(roostr);
  }

  Lee_test->matrix_pred_newworld.Write("matrix_pred_newworld");// (1, bins)
  Lee_test->matrix_data_newworld.Write("matrix_data_newworld");// (1, bins)

  for(auto it_sub=Lee_test->matrix_sub_flux_geant4_Xs_newworld.begin(); it_sub!=Lee_test->matrix_sub_flux_geant4_Xs_newworld.end(); it_sub++) {
    int index = it_sub->first;
    roostr = TString::Format("matrix_sub_flux_geant4_Xs_newworld_%d", index);
    Lee_test->matrix_sub_flux_geant4_Xs_newworld[index].Write(roostr);
  }

  //file_collapsed_covariance_matrix->Close();

  //////////

  if( config_Lee::flag_plotting_systematics ) Lee_test->Plotting_systematics();

  //////////////////////////////////////////////////////////////////////////////////////// Goodness of fit

  //Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_GoF;
  //Lee_test->Set_Collapse();

  if( config_Lee::flag_GoF_output2file_default_0 ) {
    file_collapsed_covariance_matrix->cd();

    for(auto it=Lee_test->map_data_spectrum_ch_bin.begin(); it!=Lee_test->map_data_spectrum_ch_bin.end(); it++) {
      int val_ch = it->first;
      int size_map = it->second.size();
      int size_before = 0;
      for(int idx=1; idx<val_ch; idx++) {
        int size_current = Lee_test->map_data_spectrum_ch_bin[idx].size();
        size_before += size_current;
      }

      vector<int>vc_target_chs;
      for(int ibin=1; ibin<size_map; ibin++) {
        vc_target_chs.push_back( size_before + ibin -1 );
      }

      vector<int>vc_support_chs;

      Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 100 + val_ch );

      vc_val_GOF.push_back( Lee_test->val_GOF_noConstrain );
      vc_val_GOF_NDF.push_back( Lee_test->val_GOF_NDF );
      //cout<<" ---> "<<Lee_test->val_GOF_NDF<<"\t"<<Lee_test->val_GOF_noConstrain<<endl;
    }

    tree_config->Fill();
    tree_config->Write();
    file_collapsed_covariance_matrix->Close();
  }

  // bool flag_both_numuCC            = config_Lee::flag_both_numuCC;// 1
  // bool flag_CCpi0_FC_by_numuCC     = config_Lee::flag_CCpi0_FC_by_numuCC;// 2
  // bool flag_CCpi0_PC_by_numuCC     = config_Lee::flag_CCpi0_PC_by_numuCC;// 3
  // bool flag_NCpi0_by_numuCC        = config_Lee::flag_NCpi0_by_numuCC;// 4
  // bool flag_nueCC_PC_by_numuCC_pi0 = config_Lee::flag_nueCC_PC_by_numuCC_pi0;// 5
  // bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = config_Lee::flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC;// 6, HghE>800 MeV
  // bool flag_nueCC_LowE_FC_by_all   = config_Lee::flag_nueCC_LowE_FC_by_all;// 7
  // bool flag_nueCC_FC_by_all        = config_Lee::flag_nueCC_FC_by_all;// 8

  //////////////////////////////////////////////////////////////////////////////////////// user's goodness-of-fit test

  // Lee_test->scaleF_Lee = ?;
  // Lee_test->scaleF_Lee_Np = ?;
  // Lee_test->scaleF_Lee_0p = ?;

  // Lee_test->Set_Collapse();

  ///////// validated by 1d = 5, and 2d = [5,5], [1,5], [5,1]

  if( 0 ) {
    Lee_test->scaleF_Lee = 5;
    Lee_test->Set_Collapse();

    vector<int>vc_target_chs;
    vc_target_chs.push_back( 0 );
    vc_target_chs.push_back( 2 );

    vector<int>vc_support_chs;
    for(int ibin=1; ibin<=16*4; ibin++) vc_support_chs.push_back( 4+ibin -1 );

    Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 12301 );
  }

  if( 0 ) {
    Lee_test->scaleF_Lee_Np = 5;
    Lee_test->scaleF_Lee_0p = 5;
    Lee_test->Set_Collapse();

    vector<int>vc_target_chs;
    vc_target_chs.push_back( 0 );
    vc_target_chs.push_back( 2 );

    vector<int>vc_support_chs;
    for(int ibin=1; ibin<=16*4; ibin++) vc_support_chs.push_back( 4+ibin -1 );

    Lee_test->Exe_Goodness_of_fit_detailed( vc_target_chs, vc_support_chs, 12302 );
  }

  //////////////////////////////////////////////////////////////////////////////////////// Asimov/Data fitting

  if( 0 ) {
    Lee_test->scaleF_Lee_Np = 5;
    Lee_test->scaleF_Lee_0p = 5;
    Lee_test->Set_Collapse();

    ///////// Four options:
    //Lee_test->Set_measured_data();// (1), really data
    // Lee_test->Set_toy_Asimov();// (2), Asiomv sample
    // Lee_test->Set_Variations(int num_toy);// (3), generate many toy-MC
    // Lee_test->Set_toy_Variation(int itoy);// which toy-MC to be used
    // Lee_test->Set_fakedata(TMatrixD matrix_fakedata);// (4), user's defined fakedata

    Lee_test->Set_toy_Asimov();

    Lee_test->Minimization_Lee_Np_0p_strength_FullCov(4, 4, "");
    // Lee_test->Minimization_Lee_Np_0p_strength_FullCov(5, 4, "Np");// fix Np if the string containts "Np", e.g "NNp"
    // Lee_test->Minimization_Lee_Np_0p_strength_FullCov(4, 5, "0p");// fix 0p if the string containts "0p", e.g. "A_0p"
    // Lee_test->Minimization_Lee_Np_0p_strength_FullCov(5, 5, "Np_0p");// fit both if the string containts "Np" and "0p", e.g. "Np_abc_0p"

    ///////// the fitting results/info are saved in
    // Lee_test->minimization_status;// (int type), should be 0 for a successful fitting
    // Lee_test->minimization_chi2;
    // Lee_test->minimization_Lee_Np_strength_val;
    // Lee_test->minimization_Lee_Np_strength_err;// error from the default Minuit2 setting, which is at the points delta_chi2 = chi2_var - chi2_min = 1
    // Lee_test->minimization_Lee_0p_strength_val;
    // Lee_test->minimization_Lee_0p_strength_err;// error from the default Minuit2 setting, which is at the points delta_chi2 = chi2_var - chi2_min = 1

    cout<<endl;
    cout<<" ---> Fitting results"<<endl;
    cout<<" Lee_test->minimization_status: "<<Lee_test->minimization_status<<endl;
    cout<<" Lee_test->minimization_chi2:   "<<Lee_test->minimization_chi2<<endl;
    cout<<" Lee_test->minimization_Lee_Np_strength_val: "<<Lee_test->minimization_Lee_Np_strength_val
        <<" +/- "<<Lee_test->minimization_Lee_Np_strength_err<<endl;
    cout<<" Lee_test->minimization_Lee_0p_strength_val: "<<Lee_test->minimization_Lee_0p_strength_val
        <<" +/- "<<Lee_test->minimization_Lee_0p_strength_err<<endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////// Feldman-Cousins approach

  /// The (Np, 0p) space are divided into many grid. e.g. 10x10: Np(0, 10), 0p(0, 10) ---> TH2D *h2d_space = new TH2D("h2d_space", "", 10, 0, 10, 10, 0, 10);
  ///
  /// At each space point (f_Np, f_0p), we will cacluate the confidence level following Feldman-Cousins approach as following:
  ///
  /// (1) Generate the delta_chi2 distribution "distribution_dchi2" by many pseudo experiments, which is corresponding the the pdf of dchi2 at one space point.
  ///     * dchi2 = chi2_(f_Np, f_0p) - chi2_min
  ///     * pseduo experiments are generated considering both statistical and systematic uncertainties:
  ///
  ///     For each pseudo experiment, we will caculate the chi2_(f_Np, f_0p) and chi2_min
  ///     ...
  ///     By many pseudo experiments, we will have the distribution of dchi2
  ///
  /// (2) Calculate the dchi2 of "data", "data" can be real data, or Asimov sample, or others wanted to study
  ///     dchi2_data = chi2_(f_Np, f_0p)_data - chi2_min_data
  ///
  /// (3) Calculte the pvalue or the confidence value at one space point (f_Np, f_0p)
  ///     pvalue: N( distribution_dchi2 > dchi2_data ) / N( distribution_dchi2 )
  ///     confidence level = 1 - pvalue
  ///

  //////////////////////////////// Scripts:

  if( 1 ) {

    ///////

    int grid_Np = 0;
    int grid_0p = 0;
    double true_Np = 0;
    double true_0p = 0;
    vector<int>vec_min_status;
    vector<double>vec_chi2_var;
    vector<double>vec_min_chi2;
    vector<double>vec_dchi2;
    vector<double>vec_min_fNp_val;
    vector<double>vec_min_fNp_err;
    vector<double>vec_min_f0p_val;
    vector<double>vec_min_f0p_err;

    std::cout << "i file: " << ifile << "\n";

    if (ifile==-1) {
      std::cout << "creating asimov root file\n";
      roostr = TString("sub_fit_Asimov.root");
    } else if (ifile==0) {
      roostr = TString("sub_fit_data.root");
    } else {
      roostr = TString::Format("sub_fit_%06d.root", ifile);
    }

    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");

    tree->Branch( "grid_Np",          &grid_Np,          "grid_Np/I" );
    tree->Branch( "grid_0p",          &grid_0p,          "grid_0p/I" );
    tree->Branch( "true_Np",          &true_Np,          "true_Np/D" );
    tree->Branch( "true_0p",          &true_0p,          "true_0p/D" );
    tree->Branch( "vec_min_status",   &vec_min_status );
    tree->Branch( "vec_chi2_var",     &vec_chi2_var );
    tree->Branch( "vec_min_chi2",     &vec_min_chi2 );
    tree->Branch( "vec_dchi2",        &vec_dchi2 );
    tree->Branch( "vec_min_fNp_val",  &vec_min_fNp_val );
    tree->Branch( "vec_min_fNp_err",  &vec_min_fNp_err );
    tree->Branch( "vec_min_f0p_val",  &vec_min_f0p_val );
    tree->Branch( "vec_min_f0p_err",  &vec_min_f0p_err );


    ///////
    int Ntoys = 1; // number of toy-MC used to generate the distribution_dchi2
    // changed this from 100 to 1 since it shouldn't matter for Asimov

    if (ifile==0 or ifile==-1) { // run it with data or Asimov
      Ntoys = 1;
    }

    /////// 2d space of (Np, 0p)
    int bins_Np = 16;
    int bins_0p = 16;

    TH2D *h2d_space = new TH2D("h2d_space", "", bins_Np, -0.5, 15.5, bins_0p, -0.5, 15.5);

    double pars_2d[2] = {0};


    /*
    if (ifile==0) {
          Lee_test->Set_measured_data();
    } else if (ifile==-1) {
          Lee_test->scaleF_Lee_Np = 1.;
          Lee_test->scaleF_Lee_0p = 1.;
          Lee_test->Set_Collapse();// apply the values
          Lee_test->Set_toy_Asimov();
    } else {
          Lee_test->Set_Variations( Ntoys );
    }
    */

    /////// scan the entrie space defined
    for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++) {
      for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++) {

	cout<<TString::Format(" processing Np/0p: %3d %3d", bin_Np, bin_0p)<<endl;

	/// give values of one space point
	double grid_Np_val = h2d_space->GetXaxis()->GetBinCenter( bin_Np );
	double grid_0p_val = h2d_space->GetYaxis()->GetBinCenter( bin_0p );

	///
	grid_Np = bin_Np;
	grid_0p = bin_0p;
	true_Np = grid_Np_val;
	true_0p = grid_0p_val;
	vec_min_status.clear();
	vec_chi2_var.clear();
	vec_min_chi2.clear();
	vec_dchi2.clear();
	vec_min_fNp_val.clear();
	vec_min_fNp_err.clear();
	vec_min_f0p_val.clear();
	vec_min_f0p_err.clear();

	pars_2d[0] = grid_Np_val;
	pars_2d[1] = grid_0p_val;

	/// generate pseudo experiments

	if (ifile==0) { // make the data file
          Lee_test->Set_measured_data();
        } else if (ifile==-1) { // make the Asimov file at the Standard Model
	  Lee_test->scaleF_Lee_Np = 1.;
	  Lee_test->scaleF_Lee_0p = 1.;
	  Lee_test->Set_Collapse();// apply the values
          Lee_test->Set_toy_Asimov();
        } else { // make the variations file
	  Lee_test->scaleF_Lee_Np = grid_Np_val;
	  Lee_test->scaleF_Lee_0p = grid_0p_val;
	  Lee_test->Set_Collapse();// apply the values
          Lee_test->Set_Variations( Ntoys );
        }

	for(int itoy=1; itoy<=Ntoys; itoy++) {// scan each pseudo experiment

	  if (ifile > 0) { // only for variation files
	    Lee_test->Set_toy_Variation(itoy);
          }

	  double chi2_var = Lee_test->FCN_Np_0p( pars_2d );// calcualte the chi2 value at the point

          //cout << "lhagaman debug, pars2d, chi2: " << pars_2d[0] << ", " << pars_2d[1] << ", " << chi2_var << "\n";

	  /////// do minimization: initial value is important. May find local minimum if the values are not suitable

	  double initial_Np = grid_Np_val;
	  double initial_0p = grid_0p_val;

	  /// a simple way: the true values as the initial value
	  //if( 1 ) {
	  //  initial_Np = grid_Np_val;
	  //  initial_0p = grid_0p_val;
	  //}

	  /// a more exact way: scan the space to find suitable initial values

	  ///////

	  Lee_test->Minimization_Lee_Np_0p_strength_FullCov(initial_Np, initial_0p, "");

	  double chi2_min = Lee_test->minimization_chi2;

	  double dchi2 = chi2_var - chi2_min;

	  if (dchi2 < 0.) dchi2 = 0.; // protects against small negative values, mostly relevant for Asimov at the exact grid point

	  ///////
	  vec_min_status.push_back( Lee_test->minimization_status );
	  vec_chi2_var.push_back( chi2_var );
	  vec_min_chi2.push_back( chi2_min );
	  vec_dchi2.push_back( dchi2 );
	  vec_min_fNp_val.push_back( Lee_test->minimization_Lee_Np_strength_val );
	  vec_min_fNp_err.push_back( Lee_test->minimization_Lee_Np_strength_err);
	  vec_min_f0p_val.push_back( Lee_test->minimization_Lee_0p_strength_val );
	  vec_min_f0p_err.push_back( Lee_test->minimization_Lee_0p_strength_err);

	}// for(int itoy=1; itoy<=Ntoys; itoy++)

	////// save the information into root-file for each space point

	tree->Fill();

      }// for(int bin_0p=1; bin_0p<=bins_0p; bin_0p++)
    }// for(int bin_Np=1; bin_Np<=bins_Np; bin_Np++)

    tree->Write();
    subroofile->Close();

  }

  ////////////////////////////////////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" -----------------------------------"<<endl;
  cout<<" Check the initialization at the end"<<endl;
  cout<<" -----------------------------------"<<endl;

  cout<<endl;
  cout<<" ---> flag_syst_flux_Xs    "<<Lee_test->flag_syst_flux_Xs<<endl;
  cout<<" ---> flag_syst_detector   "<<Lee_test->flag_syst_detector<<endl;
  cout<<" ---> flag_syst_additional "<<Lee_test->flag_syst_additional<<endl;
  cout<<" ---> flag_syst_mc_stat    "<<Lee_test->flag_syst_mc_stat<<endl;
  cout<<endl;

  cout<<" ---> LEE channel size (set by array_LEE_ch in config): "<<Lee_test->map_Lee_ch.size()<<endl;
  if( (int)(Lee_test->map_Lee_ch.size()) ) {
    for(auto it_map_Lee=Lee_test->map_Lee_ch.begin(); it_map_Lee!=Lee_test->map_Lee_ch.end(); it_map_Lee++) {
      cout<<" ---> LEE channel: "<< it_map_Lee->first<<endl;
    }
  }
  cout<<endl;

  cout<<" ---> LEE_Np channel size (set by array_LEE_Np_ch in config): "<<Lee_test->map_Lee_Np_ch.size()<<endl;
  if( (int)(Lee_test->map_Lee_Np_ch.size()) ) {
    for(auto it_map_Lee_Np=Lee_test->map_Lee_Np_ch.begin(); it_map_Lee_Np!=Lee_test->map_Lee_Np_ch.end(); it_map_Lee_Np++) {
      cout<<" ---> LEE_Np channel: "<< it_map_Lee_Np->first<<endl;
    }
  }
  cout<<endl;

  cout<<" ---> LEE_0p channel size (set by array_LEE_0p_ch in config): "<<Lee_test->map_Lee_0p_ch.size()<<endl;
  if( (int)(Lee_test->map_Lee_0p_ch.size()) ) {
    for(auto it_map_Lee_0p=Lee_test->map_Lee_0p_ch.begin(); it_map_Lee_0p!=Lee_test->map_Lee_0p_ch.end(); it_map_Lee_0p++) {
      cout<<" ---> LEE_0p channel: "<< it_map_Lee_0p->first<<endl;
    }
  }
  cout<<endl;

  cout<<" ---> MC stat file: "<<Lee_test->syst_cov_mc_stat_begin<<".log - "<<Lee_test->syst_cov_mc_stat_end<<".log"<<endl;
  cout<<endl;

  cout<<" ---> check: scaleF_POT "<<scaleF_POT<<", file_flag "<<ifile<<endl<<endl;

  cout<<endl<<endl;
  cout<<" ---> Complete all the program"<<endl;
  cout<<endl<<endl;

  if( config_Lee::flag_display_graphics ) {
    cout<<endl<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl<<endl;

    theApp.Run();
  }

  return 0;
}
