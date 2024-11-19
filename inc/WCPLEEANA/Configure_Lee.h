namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes

  TString spectra_file = "./merge.root";
  //Erin
  TString time_file = "./merge_time.root";
  //
  TString flux_Xs_directory = "/exp/uboone/data/users/eyandel/hist_rootfiles/XsFlux/";
  TString detector_directory = "/exp/uboone/data/users/eyandel/hist_rootfiles/det/";
  TString mc_directory = "/exp/uboone/data/users/eyandel/hist_rootfiles/mc_stat/";



  int channels_observation = 0;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse
                               // NOTE: This value is not used in the lastest version

  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 17;//cov_18.root is uncorrelated reweighting and cov_19.root is correlated

  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  int syst_cov_mc_stat_end   = 1;


  ///////////////////////////////

  //int array_LEE_ch[1] = {0};
  int array_LEE_ch[3] = {4,5,6}; // for 1d strength. Pay attention to the setting.
                             // Confict will happen if both the 1d and 2d(defined below) work.
                             // Three options: 1d works, 2d works, or neither works

  //int array_LEE_ch[12] = {2, 5, 8, 11, 14, 17, 3, 6, 9, 12, 15, 18};

  /////////////////////////////// separated Np and 0p strengthes: 2d strengthes

  int array_LEE_Np_ch[1] = {0};// element value "0" will not be set to the LEE_ch
  int array_LEE_0p_ch[1] = {0};// element value "0" will not be set to the LEE_ch

  //int array_LEE_Np_ch[1] = {0};
  //int array_LEE_0p_ch[1] = {0};

  ///////////////////////////////

  /// some places may need to be changed when use different file-formats
  /// void TLee::Set_Spectra_MatrixCov()
  /// (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  /// (*) map_Lee_ch                     -----> tag LEE channels
  /// (*) map_detectorfile_str           -----> detector files

  ////////// display graphics flag

  bool flag_display_graphics = false;

  ////////// systematics flag

  bool flag_syst_flux_Xs    = true;
  bool flag_syst_detector   = true;
  bool flag_syst_additional = true;
  bool flag_syst_mc_stat    = true;
  bool flag_syst_reweight        = false;
  bool flag_syst_reweight_cor    = false;
  //ns timing error
  bool flag_syst_time = false;

  double Lee_strength_for_outputfile_covariance_matrix = 0;
  //double Lee_strength_for_outputfile_covariance_matrix = 1.0;
  //double Lee_strength_for_outputfile_covariance_matrix = 0.2123;
  //double Lee_strength_for_outputfile_covariance_matrix = 0.1328;

  double Lee_Np_strength_for_outputfile_covariance_matrix = 0;
  double Lee_0p_strength_for_outputfile_covariance_matrix = 0;

  bool flag_plotting_systematics   = true;

  ////////// goodness of fit

  double Lee_strength_for_GoF         = 0;
  //double Lee_strength_for_GoF         = 1.0;
  //double Lee_strength_for_GoF         = 0.2123;
  //double Lee_strength_for_GoF         = 0.1328;

  double Lee_Np_strength_for_GoF      = 0;
  double Lee_0p_strength_for_GoF      = 0;

  bool flag_GoF_output2file_default_0 = true;

  bool flag_lookelsewhere             = false;

  bool flag_both_numuCC            = false;// 1
  bool flag_CCpi0_FC_by_numuCC     = false;// 2
  bool flag_CCpi0_PC_by_numuCC     = false;// 3
  bool flag_NCpi0_by_numuCC        = false;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = false;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = false;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = false;// 7
  bool flag_nueCC_FC_by_all        = false;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = true;
  bool flag_Lee_scan_data     = true;

  bool flag_GOF = true;

  ////////// MicroBooNE suggested

  bool flag_chi2_data_H0 = 1;
  bool flag_dchi2_H0toH1 = 1;

  ////////// Advanced tools

  ///// void TLee::Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed)
  ///// do the fitting on the spectra and cov_total after constraint ?
  bool flag_Lee_minimization_after_constraint = 0;// hardcoded, only for the standard 7-ch fitting

}
