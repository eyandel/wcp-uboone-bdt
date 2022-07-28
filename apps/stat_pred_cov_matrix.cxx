#include <iostream>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h" 
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TVectorD.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{
  
  if (argc < 2){
    std::cout << "./stat_cov_matrix -r[run period]" << std::endl;
  }
  int run = 1; // run 1 ...
  bool flag_osc = false;
  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':
      run = atoi(&argv[i][2]); // which run period
      break;
    case 'o':
      flag_osc = atoi(&argv[i][2]); // which run period
      break;
    }
  }
  
  CovMatrix cov;
  if (flag_osc) cov.add_osc_config();

  cov.print_rw(cov.get_rw_info()); 
 
  // Get the file based on runno ...
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();
  // Construct the histogram ...

  // outfilename ...
  TString outfile_name = "./hist_rootfiles/run_pred_stat.root";

  TH1F *htemp;
  std::map<TString, TH1F*> map_histoname_hist;
  std::map<int, TH1F*> map_covch_hist;
  std::map<int, TH1F*> map_obsch_hist;
  
  for (auto it = map_inputfile_info.begin(); it!=map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);

    if (filetype == 5 || filetype == 15) continue;// skip data ...
    //    std::cout << filetype << " " << period << " " << out_filename << " " << file_no << std::endl;    
    if (period == run || run == 0){ // for a period ...
      //outfile_name = out_filename;
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename, 0);
      
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	TString histoname = std::get<0>(*it1);
	Int_t nbin = std::get<1>(*it1);
	float llimit = std::get<2>(*it1);
	float hlimit = std::get<3>(*it1);
	TString var_name = std::get<4>(*it1);
	TString ch_name = std::get<5>(*it1);
	TString add_cut = std::get<6>(*it1);
	TString weight = std::get<7>(*it1);
	
	//std::cout << std::get<0>(*it1)  << " " << std::get<1>(*it1) << " " << std::get<4>(*it1) << " " << std::get<5>(*it1) << " " << std::get<6>(*it1) << " " << std::get<7>(*it1) << std::endl;
	htemp = new TH1F(histoname, histoname, nbin, llimit, hlimit);
	map_histoname_hist[histoname] = htemp;

	
	
	int covch = cov.get_covch_name(ch_name);

	if (map_covch_hist.find(covch) == map_covch_hist.end()){
	  TH1F *htemp1 = (TH1F*)htemp->Clone(Form("pred_covch_%d",covch));
	  map_covch_hist[covch] = htemp1;
	}

	int obsch = cov.get_obsch_name(ch_name);
	if (map_obsch_hist.find(obsch) == map_obsch_hist.end()){
	  TH1F *htemp1 = (TH1F*)htemp->Clone(Form("pred_obsch_%d",obsch));
	  map_obsch_hist[obsch] = htemp1;
	}
      }
      //  std::cout << input_filename << " " << filetype << " " << out_filename << std::endl; 
    }
  }
  std::cout << outfile_name << std::endl;

  TMatrixD* cov_add_mat = cov.get_add_cov_matrix();
  //  TMatrixD* mat_collapse = cov.get_mat_collapse();
  //  std::cout << mat_collapse->GetNrows() << " " << mat_collapse->GetNcols() << std::endl;
  
  // create a covariance matrix for bootstrapping ...
  // TMatrixD* cov_mat = new TMatrixD(mat_collapse->GetNcols(), mat_collapse->GetNcols());
  TMatrixD* cov_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  
  // create a covariance matrix for det systematics ...
  TVectorD* vec_mean = new TVectorD(cov_add_mat->GetNrows());
  
  cov.gen_pred_stat_cov_matrix(run, map_covch_hist, map_histoname_hist, vec_mean, cov_mat);

  
  TMatrixD* frac_cov_mat = new TMatrixD(cov_add_mat->GetNrows(), cov_add_mat->GetNcols());
  for (size_t i=0; i!= frac_cov_mat->GetNrows(); i++){
    double val_1 = (*vec_mean)(i);
    for (size_t j=0; j!=frac_cov_mat->GetNrows();j++){
      double val_2 = (*vec_mean)(j);
      double val = (*cov_mat)(i,j);
      if (val_1 ==0 && val_2 == 0){
	(*frac_cov_mat)(i,j) = 0;
      }else if (val_1 ==0 || val_2 ==0){
	if (val !=0){
	  if (i==j){
	    (*frac_cov_mat)(i,j) = 1.; //
	  }else{
	    (*frac_cov_mat)(i,j) = 0;
	  }
	}else{
	  (*frac_cov_mat)(i,j) = 0;
	}
      }else{
	(*frac_cov_mat)(i,j) = val/val_1/val_2;
      }
    }
  }
  TMatrixD* mat_collapse = cov.get_mat_collapse();
  TMatrixD mat_collapse_T(mat_collapse->GetNcols(), mat_collapse->GetNrows()); 
  mat_collapse_T.Transpose(*mat_collapse);

  TMatrixD cov_mat_collapse = mat_collapse_T * (*cov_mat) * (*mat_collapse);
  TVectorD vec_mean_collapse = (mat_collapse_T)* (*vec_mean);
  
  TFile *file = new TFile(outfile_name,"RECREATE");
  vec_mean->Write(Form("full_vec_mean_%d",run));  
  cov_mat->Write(Form("full_cov_mat_%d",run));
  frac_cov_mat->Write(Form("frac_full_cov_mat_%d",run));
  cov_mat_collapse.Write(Form("cov_mat_%d",run));
  vec_mean_collapse.Write(Form("vec_mean_%d",run));  
  

  for (auto it = map_covch_hist.begin(); it != map_covch_hist.end(); it++){
    int covch = it->first;
    TH1F *htemp = it->second;
    int obsch = cov.get_obsch_fcov(covch);
    TH1F *htemp1 = map_obsch_hist[obsch];
    htemp1->Add(htemp);
    
    
    ((TH1F*)it->second)->SetDirectory(file);
  }
  for (auto it = map_obsch_hist.begin(); it != map_obsch_hist.end(); it++){
     ((TH1F*)it->second)->SetDirectory(file);
  }
  
  file->Write();
  file->Close();

  
}
