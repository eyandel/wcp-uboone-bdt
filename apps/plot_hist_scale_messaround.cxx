#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <utility>

#include "WCPLEEANA/master_cov_matrix.h"
#include "WCPLEEANA/bayes.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TString.h"
#include "TLine.h"
#include "TPDF.h"
#include "TF1.h"
#include "TMatrixD.h"

using namespace std;
using namespace LEEana;

int main( int argc, char** argv )
{

  if (argc < 2){
    std::cout << "./plot_hist -r[#run, 0 for all] -e[1 for standard, 2 for Bayesian, 3 for total uncertainty] -l[LEE strength, check against total uncertainty config] -s[path-to-external covmatrix file] -c[check MC and DATA spectra against the ones from external file] -x[lee strength to scale]" << std::endl;
  }

  int run = 1; // run 1 ...
  int flag_err = 1;// 1 for standard, 2 for Bayesian, 3 for breakdown
  std::map<int, double> sumtotalcov;
  std::map<int, std::pair<double, int>> GOF;
  float lee_strength = 0; // no LEE strength ...
  int flag_display = 1;
  int flag_breakdown = 1;
  TString cov_inputfile = "";
  int flag_check = 0;
  int flag_truthlabel = 0;
  float fit_scale = 0;

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'r':{
      run = atoi(&argv[i][2]); // which run period
    }break;
    case 'e':{
      flag_err = atoi(&argv[i][2]); // error for plotting
    }break;
    case 'l':{
      lee_strength = atof(&argv[i][2]);
    }break;
    case 'd':{
      flag_display = atoi(&argv[i][2]);
    }break;
    case 'b':{
      flag_breakdown = atoi(&argv[i][2]);
    }break;
    case 's':{
      TString sss = argv[i];
      cov_inputfile = sss(2, sss.Length()-2);
    }break;
    case 'c':{
      flag_check = atoi(&argv[i][2]);
    }break;
    case 't':{
      flag_truthlabel = atoi(&argv[i][2]);
    }break;
    case 'x':{
      fit_scale = atof(&argv[i][2]);
    }break;
  }
  }

  CovMatrix cov;

  // get data histograms ...
  // filetype, period, outfilename, external pot, fileno
  std::map<TString, std::tuple<int, int, TString, float, int, double, int> > map_inputfile_info = cov.get_map_inputfile_info();

  TFile *temp_file;
  TH1F *htemp;
  TTree *T;
  Double_t pot;
  std::map<TString, std::pair<TH1F*, double> > map_name_histogram;

  // data POT ...
  std::map<int, double> map_data_period_pot;
  //  std::vector<TH1F*> temp_histograms;

  // open all the histograms ...
  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    temp_file = new TFile(out_filename);
    T = (TTree*)temp_file->Get("T");
    T->SetBranchAddress("pot",&pot);
    T->GetEntry(0);


    if (filetype==5 || filetype == 15){
      map_data_period_pot[period] = pot;
    }

    // std::cout << pot << std::endl;

    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > all_histo_infos;
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
    std::copy(histo_infos.begin(), histo_infos.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_err2 = cov.get_histograms(input_filename,1);
    std::copy(histo_infos_err2.begin(), histo_infos_err2.end(), std::back_inserter(all_histo_infos));
    std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos_cros = cov.get_histograms(input_filename,2);
    std::copy(histo_infos_cros.begin(), histo_infos_cros.end(), std::back_inserter(all_histo_infos));

    for (size_t i=0;i!=all_histo_infos.size();i++){
      htemp = (TH1F*)temp_file->Get(std::get<0>(all_histo_infos.at(i)));
      //      std::cout << std::get<0>(all_histo_infos.at(i)) << " " << htemp->GetSum() << std::endl;
      //      temp_histograms.push_back(htemp);
      map_name_histogram[std::get<0>(all_histo_infos.at(i))] = std::make_pair(htemp, pot);
    }
  }




  // create histograms for data, create histograms for predictions
  // obsch --> histograms (data, prediction, prediction_err2
  std::map<int, std::vector<TH1F*> > map_obsch_histos;
  // obsch --> break down histograms (truth label: add_cut) prediction
  std::map<int, std::vector<TH1F*> > map_obsch_subhistos;
  // Bayesian error needed ...
  // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_bayes;
  std::map<int, std::vector< std::vector< std::tuple<double, double, double, int, double> > > > map_obsch_infos;

  for (auto it = map_inputfile_info.begin(); it != map_inputfile_info.end(); it++){
    TString input_filename = it->first;
    int filetype = std::get<0>(it->second);
    int period = std::get<1>(it->second);
    TString out_filename = std::get<2>(it->second);
    int file_no = std::get<4>(it->second);
    if (filetype == 5 || filetype == 15){
      // name, nbin, lowlimit, highlimit, variable, channel cut, additional cut, weight
      std::vector< std::tuple<TString,  int, float, float, TString, TString, TString, TString > > histo_infos = cov.get_histograms(input_filename,0);
      for (auto it1 = histo_infos.begin(); it1 != histo_infos.end(); it1++){
	int obsch = cov.get_obsch_name(std::get<5>(*it1));
	htemp = map_name_histogram[std::get<0>(*it1)].first;

	std::vector<TH1F*> vec_histos;

	TH1F *hdata = (TH1F*)htemp->Clone(Form("data_%d",obsch));
	hdata->Reset();
	TH1F *hpred = (TH1F*)htemp->Clone(Form("pred_%d",obsch));
	hpred->Reset();
	TH1F *hpred_err2 = (TH1F*)htemp->Clone(Form("pred_err2_%d",obsch));
	hpred_err2->Reset();

	vec_histos.push_back(hdata);
	vec_histos.push_back(hpred);
	vec_histos.push_back(hpred_err2);

	// get histograms ...
	map_obsch_histos[obsch] = vec_histos;
	//map_obsch_bayes[obsch].resize(htemp->GetNbinsX()+1);
	// for (Int_t i=0;i!=htemp->GetNbinsX()+1;i++){
	//   std::vector< std::tuple<double, double, double> > temp;
	//   map_obsch_bayes[obsch].push_back(temp);
	// }


	//std::cout << std::get<5>(*it1) << " " << obsch << " " << htemp->GetSum() << std::endl;
      }

      // break;
    }
  }

  // get data histograms ...
  cov.fill_data_histograms(run, map_obsch_histos, map_name_histogram);

  // get predictions and its uncertainties ...,
  cov.fill_pred_histograms(run, map_obsch_histos, map_obsch_bayes, map_obsch_infos, map_name_histogram, lee_strength, map_data_period_pot, flag_breakdown, map_obsch_subhistos);

    /* for (auto it = map_obsch_subhistos.begin(); it!= map_obsch_subhistos.end(); it++){ */
    /*         for(size_t i=0; i<it->second.size(); i++){ */
    /*             std::cout<<"DEBUG2: "<<it->first<<": "<<map_obsch_subhistos[it->first].at(i)->GetName()<<std::endl; */
    /*         } */
    /* } */

  // get Bayesian errrors ...

  if (flag_err==2){
    std::cout << lee_strength << " " << run << std::endl;

    for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
      TH1F *h1 = it->second.at(1);
      TH1F *h2 = it->second.at(2);
      int obsch = it->first;

      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = map_obsch_bayes[obsch];

      //std::cout << obsch << " " << bayes_inputs.size() << " " << bayes_inputs.at(0).size() << " " << h1->GetNbinsX() << std::endl;

      //if (obsch !=1) continue;
     double temp_sumtotalcov=0;
      for (int i=0;i!=h1->GetNbinsX()+1;i++){
	Bayes bayes;
	//	if (i!=0) continue;
	//double temp = 0, temp1=0;
	double zero_weight = 0, zero_count = 0;
	double nonzero_meas = 0, nonzero_sigma2 = 0, nonzero_weight = 0;
	for (auto it1 = bayes_inputs.begin(); it1!=bayes_inputs.end(); it1++){
	  // bayes.add_meas_component(std::get<0>((*it1).at(i)), std::get<1>((*it1).at(i)), std::get<2>((*it1).at(i)));
	  // temp += std::get<0>((*it1).at(i));
	  // temp1 += std::get<1>((*it1).at(i));
	  //std::cout << i << " " << std::get<0>((*it1).at(i)) << " " << std::get<1>((*it1).at(i)) << " " << std::get<2>((*it1).at(i)) << " " << std::endl;
	  if (std::get<0>((*it1).at(i)) == 0){
	    zero_weight += pow(std::get<2>((*it1).at(i)),2);
	    zero_count += std::get<2>((*it1).at(i));
	  }else{
	    nonzero_meas += std::get<0>((*it1).at(i));
	    nonzero_sigma2 += std::get<1>((*it1).at(i));
	    nonzero_weight += std::get<2>((*it1).at(i));
	  }
	}
	if (zero_count != 0)	bayes.add_meas_component(0,0,zero_weight/zero_count);
	if (nonzero_meas != 0)  bayes.add_meas_component(nonzero_meas, nonzero_sigma2, nonzero_weight);

	bayes.do_convolution();

	double cov = bayes.get_covariance();
	//double cov1 = bayes.get_covariance_mc();
	std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << cov  << " "  << h2->GetBinContent(i+1) << std::endl;
	// std::cout << temp << " " << temp1 << " " << h1->GetBinContent(i+1) << " " << h2->GetBinContent(i+1) << std::endl;
        if(isnan(cov) || isinf(cov)) {
            std::cout<<"Abnormal Bayesian variance.\n";
            cov = bayes.get_covariance_mc();
            //cov = h1->SetBinError(i+1, h1->GetBinError(i));
        }
	h1->SetBinError(i+1,sqrt(cov));
    if(i!=h1->GetNbinsX()) temp_sumtotalcov += cov;
      }
      sumtotalcov[obsch] = temp_sumtotalcov;
      // obsch --> bin with overflow bin --> vector of all channels (merge certain channels) --> mean and err2
      //std::map<int, std::vector< std::vector< std::tuple<double, double, double> > > > map_obsch_bayes;
    }
  }

  if (flag_err==3){
    std::cout <<"Total uncertainty from external covariance matrix: "<< cov_inputfile << std::endl;
    TFile* f_cov = new TFile(cov_inputfile, "READ");
    int flag_syst_flux_Xs = 0;
    int flag_syst_detector = 0;
    int flag_syst_additional = 0;
    int flag_syst_mc_stat = 0;
    double cov_lee_strength = 0;
    std::vector<double> *vc_val_GOF = new std::vector<double>;
    std::vector<int> *vc_val_GOF_NDF = new std::vector<int>;
    TTree* t_covconfig = (TTree*)f_cov->Get("tree");
    t_covconfig->SetBranchAddress("flag_syst_flux_Xs", &flag_syst_flux_Xs);
    t_covconfig->SetBranchAddress("flag_syst_detector", &flag_syst_detector);
    t_covconfig->SetBranchAddress("flag_syst_additional", &flag_syst_additional);
    t_covconfig->SetBranchAddress("flag_syst_mc_stat", &flag_syst_mc_stat);
    t_covconfig->SetBranchAddress("user_Lee_strength_for_output_covariance_matrix", &cov_lee_strength);
    t_covconfig->SetBranchAddress("vc_val_GOF", &vc_val_GOF);
    t_covconfig->SetBranchAddress("vc_val_GOF_NDF", &vc_val_GOF_NDF);
    t_covconfig->GetEntry(0);

    // fill GOF map
    for(size_t vv=0; vv<vc_val_GOF->size(); vv++)
    {
        GOF[vv] = std::make_pair(vc_val_GOF->at(vv), vc_val_GOF_NDF->at(vv));
    }

    // absolute cov matrix
    TMatrixD* matrix_absolute_cov = (TMatrixD*)f_cov->Get("matrix_absolute_cov_newworld");
    TMatrixD* matrix_absolute_detector_cov = (TMatrixD*)f_cov->Get("matrix_absolute_detector_cov_newworld");

    std::cout <<"User: "<< "chosen LEE strength: "<< lee_strength << " run option: " << run << std::endl;
    std::cout <<"Cov matrix config: \n"
                <<"\t syst_flux_Xs: "<< flag_syst_flux_Xs << std::endl
                <<"\t syst_detector: "<< flag_syst_detector << std::endl
                <<"\t syst_additional: "<< flag_syst_additional << std::endl
                <<"\t syst_mc_stat: "<< flag_syst_mc_stat << std::endl
                <<"\t LEE_strength: "<< cov_lee_strength << std::endl;
    if( abs(lee_strength-cov_lee_strength)>1e-6 ) {
        std::cout<<"ERROR: plot_hist -l option is inconsistent with external cov matrix LEE strength!\n";
        return 1;
    }

    // construct a map from (obsch, bin) to cov index
    std::map< std::pair<int, int>, int> obsch_bin_index;
    int index = 0;
    for(size_t i=0; i<map_obsch_histos.size(); i++)
    {
        // + overflow bin
        for(int j=1; j<=map_obsch_histos[i+1].at(1)->GetNbinsX()+1; j++)
        {
            index += 1;
            // cov index starts from 0
            obsch_bin_index.insert({{i+1, j}, index-1});
        }
    }

    TMatrixD* matrix_pred = (TMatrixD*)f_cov->Get("matrix_pred_newworld");
    TMatrixD* matrix_data = (TMatrixD*)f_cov->Get("matrix_data_newworld");
    for (auto it = map_obsch_histos.begin(); it!= map_obsch_histos.end(); it++){
      TH1F *h1 = it->second.at(1); // error --> total uncertainty
      TH1F *h2 = it->second.at(2); // bonus: error --> detector systematic uncertainty
      int obsch = it->first;

      TH1F* htemp_data = (TH1F*)h1->Clone("htemp_data");
      TH1F* htemp_pred = (TH1F*)h1->Clone("htemp_pred");
      htemp_data->Reset();
      htemp_pred->Reset();
      double temp_sumtotalcov=0;
      for (int i=0;i!=h1->GetNbinsX()+1;i++){
        int index = obsch_bin_index.find(std::make_pair(obsch, i+1))->second;
        double total_uncertainty = (*matrix_absolute_cov)(index, index); // only diagonal term
        //summation of total cov
        if(i!=h1->GetNbinsX()){ //no overflow bin in this calculation
        for(int j=0; j!=h1->GetNbinsX();j++){
            int jndex = obsch_bin_index.find(std::make_pair(obsch, j+1))->second;
            temp_sumtotalcov += (*matrix_absolute_cov)(index, jndex);
        }
        }
        //
        double detector_uncertainty = (*matrix_absolute_detector_cov)(index, index); // only diagonal term
        std::cout << obsch << " " << i << " "	  << h1->GetBinContent(i+1) << " " << total_uncertainty << " " <<sqrt(total_uncertainty)/h1->GetBinContent(i+1) << " "  << h2->GetBinContent(i+1) << " " << detector_uncertainty << std::endl;
	    h1->SetBinError(i+1,sqrt(total_uncertainty));
        h2->SetBinError(i+1,sqrt(detector_uncertainty));

        if(flag_check == 1){
            htemp_data->SetBinContent(i+1, (*matrix_data)(0, index));
            htemp_pred->SetBinContent(i+1, (*matrix_pred)(0, index));
        }
      }

      sumtotalcov[obsch] = temp_sumtotalcov;
      if(flag_check == 1){
        it->second.push_back(htemp_data);
        it->second.push_back(htemp_pred);
      }
    }
  }



  if (flag_display == 1){
    // plotting ...
    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    if(flag_breakdown == 0){
      cout<<"flag_breakdown = 0!"<<endl;
  }
  else // flag_breakdown == true
  {
    gROOT->ProcessLine(".x DrawOption.cc");
    /*gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(132);
    gStyle->SetLegendTextSize(0.04);*/

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(22); //132);
    gStyle->SetLegendTextSize(0.035);
    gStyle->SetLabelFont(22); //132);
    gStyle->SetLabelSize(0.06);
    gStyle->SetPalette(kCubehelix);
    

    auto colors = TColor::GetPalette();

    int nchannels = map_obsch_subhistos.size();
    TCanvas *canvas[nchannels];
    TCanvas *canvas3[nchannels];
    TGraphAsymmErrors *gr[nchannels];
    TGraphAsymmErrors *gratio_mc[nchannels];
    TGraphAsymmErrors *gratio_mc2[nchannels]; // to plot uncertainty
    TGraphAsymmErrors *gratio_data[nchannels];
    TGraphAsymmErrors *gratio_data2[nchannels]; // to normalize data points and compare to pred
    THStack *hstack[nchannels];
    TLegend *legend[nchannels];
    TLegend *legend2[nchannels];
    TLegend *legend3[nchannels];

    float scale_amount = 1.0;
    float signal_amount = 0.0;

    float scale_amount_NCpi1g   = 1.0;
    float scale_amount_NCdel    = 1.0;
    float scale_amount_NCother  = 1.0;
    float scale_amount_numuCC1g = 1.0;
    float scale_amount_out1g    = 1.0;
    float scale_amount_ext      = 1.0;
    float scale_amount_LEE      = 1.0;

    float signal_amount_NCpi1g   = 0.0;
    float signal_amount_NCdel    = 0.0;
    float signal_amount_NCother  = 0.0;
    float signal_amount_numuCC1g = 0.0;
    float signal_amount_out1g    = 0.0;
    float signal_amount_ext      = 0.0;
    float signal_amount_LEE      = 0.0;

    for(auto it = map_obsch_subhistos.begin(); it!= map_obsch_subhistos.end(); it++){
        int obschannel = it->first;
        std::cout<<"Channel: "<<obschannel<<std::endl;
        if (obschannel != 1) continue;
        canvas[obschannel-1] = new TCanvas(Form("canvas%d", obschannel), Form("channel%d", obschannel), 1200, 900);
        canvas3[obschannel-1] = new TCanvas(Form("canvas3%d", obschannel), Form("channel3%d", obschannel), 1200, 900);
        canvas[obschannel-1]->SetLeftMargin(0.15);
        canvas[obschannel-1]->SetRightMargin(0.1);
        canvas3[obschannel-1]->SetLeftMargin(0.15);
        canvas3[obschannel-1]->SetRightMargin(0.1);
        canvas[obschannel-1]->cd(); 
        
        TPad *pad1 = new TPad("pad1", "", 0.01,0.01,0.99,0.99,0,0,0);
        TPad *pad2 = new TPad("pad2", "", 0.01,0.01,0.99,0.3,0,0,0);
        pad1->SetBottomMargin(0.20);
        pad1->SetLeftMargin(0.12);
        pad1->SetRightMargin(0.1);
        pad2->SetTopMargin(0.05);
        pad2->SetLeftMargin(0.12);
        pad2->SetRightMargin(0.1);
        pad2->SetBottomMargin(0.20);
        //pad1->Draw();
        //pad2->Draw();
        hstack[obschannel-1] = new THStack(Form("hs%d", obschannel),"");
        legend[obschannel-1] = new TLegend(0.33, 0.59, 0.87, 0.94);
        legend3[obschannel-1] = new TLegend(0.5, 0.59, 0.85, 0.94);
        TH1F* hdata = (TH1F*)map_obsch_histos[obschannel].at(0)->Clone("hdata");
        TH1F* hdata_orignal = (TH1F*)map_obsch_histos[obschannel].at(0)->Clone("hdata_original");
        TH1F* hbadmatch = (TH1F*)hdata->Clone("hbadmatch");
        TH1F* hnumuCCinFV = (TH1F*)hdata->Clone("hnumuCCinFV");
        TH1F* hRnumuCCinFV = (TH1F*)hdata->Clone("hRnumuCCinFV");
        TH1F* hAnumuCCinFV = (TH1F*)hdata->Clone("hAnumuCCinFV");
        TH1F* hnueCCinFV = (TH1F*)hdata->Clone("hnueCCinFV");
        TH1F* hNCinFV = (TH1F*)hdata->Clone("hNCinFV");
        TH1F* hCCpi0inFV = (TH1F*)hdata->Clone("hCCpi0inFV");
        TH1F* hNCpi0inFV = (TH1F*)hdata->Clone("hNCpi0inFV");
        TH1F* houtFV = (TH1F*)hdata->Clone("houtFV");
        TH1F* hext = (TH1F*)hdata->Clone("hext");
        TH1F* hdirt = (TH1F*)hdata->Clone("hdirt");
        TH1F* hNCpi1g   = (TH1F*)hdata->Clone("hNCpi1g");
        TH1F* hNCdel    = (TH1F*)hdata->Clone("hNCdel");
        TH1F* hNCother  = (TH1F*)hdata->Clone("hNCother");
        TH1F* hnumuCC1g = (TH1F*)hdata->Clone("hnumuCC1g");
        TH1F* hout1g    = (TH1F*)hdata->Clone("hout1g");
        //scale
        TH1F* h1gscale    = (TH1F*)hdata->Clone("h1gscale");
        TH1F* hNCpi1gscale    = (TH1F*)hdata->Clone("hNCpi1gscale");
        TH1F* hNCdelscale    = (TH1F*)hdata->Clone("hNCdelscale");
        TH1F* hNCotherscale    = (TH1F*)hdata->Clone("hNCotherscale");
        TH1F* hnumuCC1gscale    = (TH1F*)hdata->Clone("hnumuCC1gscale");
        TH1F* hout1gscale    = (TH1F*)hdata->Clone("hout1gscale");
        TH1F* hextscale    = (TH1F*)hdata->Clone("hextscale");
        TH1F* hLEEscale = (TH1F*)hdata->Clone("hLEEscale");
        //
        TH1F* hLEE = (TH1F*)hdata->Clone("hLEE");
        TH1F* hCCQE = (TH1F*)hdata->Clone("hCCQE");
        TH1F* hNCQE = (TH1F*)hdata->Clone("hNCQE");
        TH1F* hCCRES = (TH1F*)hdata->Clone("hCCRES");
        TH1F* hNCRES = (TH1F*)hdata->Clone("hNCRES");
        TH1F* hCCDIS = (TH1F*)hdata->Clone("hCCDIS");
        TH1F* hNCDIS = (TH1F*)hdata->Clone("hNCDIS");
        TH1F* hCCMEC = (TH1F*)hdata->Clone("hCCMEC");
        TH1F* hNCMEC = (TH1F*)hdata->Clone("hNCMEC");
        TH1F* hOTHER = (TH1F*)hdata->Clone("hOTHER");
        TH1F* hXsecCosmic = (TH1F*)hdata->Clone("hXsecCosmic");
        TH1F* hXsecNumuCCinFV = (TH1F*)hdata->Clone("hXsecNumuCCinFV");
        TH1F* hXsecNC = (TH1F*)hdata->Clone("hXsecNC");
        TH1F* hXsecBkgCC = (TH1F*)hdata->Clone("hXsecBkgCC");
        hbadmatch->Reset();
        hnumuCCinFV->Reset();
        hRnumuCCinFV->Reset();
        hAnumuCCinFV->Reset();
        hnueCCinFV->Reset();
        hNCinFV->Reset();
        hCCpi0inFV->Reset();
        hNCpi0inFV->Reset();
        houtFV->Reset();
        hext->Reset();
        hdirt->Reset();
        hNCpi1g->Reset();
        hNCdel->Reset();
        hNCother->Reset();
        hnumuCC1g->Reset();
        hout1g->Reset();
        h1gscale->Reset();
        hNCpi1gscale->Reset();
        hNCdelscale->Reset();
        hNCotherscale->Reset();
        hnumuCC1gscale->Reset();
        hout1gscale->Reset();
        hextscale->Reset();
        hLEEscale->Reset();
        hLEE->Reset();
        hCCQE->Reset();
        hNCQE->Reset();
        hCCRES->Reset();
        hNCRES->Reset();
        hCCDIS->Reset();
        hNCDIS->Reset();
        hCCMEC->Reset();
        hNCMEC->Reset();
        hOTHER->Reset();
        hXsecCosmic->Reset();
        hXsecNumuCCinFV->Reset();
        hXsecNC->Reset();
        hXsecBkgCC->Reset();
        bool flag_leeexist = false;
        //hack
        double scalePOT = 1.0; // overall POT scaling
        double normalization = 1.0; // just for gratio_data2 normalization
        /* if(obschannel == 1 || obschannel == 3 || obschannel == 5) normalization = 1.35; */
        /* if(obschannel == 2 || obschannel == 4 || obschannel == 6) normalization = 1.28; */
        /* if(obschannel == 7 || obschannel == 9 || obschannel == 11) normalization = 0.95; */
        /* if(obschannel == 8 || obschannel == 10 || obschannel == 12) normalization = 0.96; */
        //scalePOT = 5.0/5.327;
        //end
        for(size_t i=0; i<it->second.size(); i++){
            TH1F* htemp = map_obsch_subhistos[obschannel].at(i);
            std::string histname = htemp->GetName();
            std::istringstream sss(histname);
            htemp->Scale(scalePOT);
            for(std::string line; std::getline(sss, line, '_');){
                if(line == "badmatch") {
                    std::cout<<"badmatch"<<" "<<histname<<std::endl;
                    hbadmatch->Add(htemp);
                    break;
                }
                if(line == "SPnumuCCBkg") {
                    std::cout<<"SPnumuCCBkg"<<" "<<histname<<std::endl;
                    hnumuCCinFV->Add(htemp);
                    break;
                }
                if(line == "RnumuCCinFV") {
                    std::cout<<"RnumuCCinFV"<<" "<<histname<<std::endl;
                    hRnumuCCinFV->Add(htemp);
                    break;
                }
                if(line == "AnumuCCinFV") {
                    std::cout<<"AnumuCCinFV"<<" "<<histname<<std::endl;
                    hAnumuCCinFV->Add(htemp);
                    break;
                }
                if(line == "nueCCinFV") {
                    std::cout<<"nueCCinFV"<<" "<<histname<<std::endl;
                    hnueCCinFV->Add(htemp);
                    break;
                }
                if(line == "SPNCBkg") {
                    std::cout<<"SPNCBkg"<<" "<<histname<<std::endl;
                    hNCinFV->Add(htemp);
                    break;
                }
                if(line == "SPnumuCCpi0Bkg") {
                    std::cout<<"SPnumuCCpi0Bkg"<<" "<<histname<<std::endl;
                    hCCpi0inFV->Add(htemp);
                    break;
                }
                if(line == "SPNCpi0Bkg") {
                    std::cout<<"SPNCpi0Bkg"<<" "<<histname<<std::endl;
                    hNCpi0inFV->Add(htemp);
                    break;
                }
                if(line == "SPoutFVBkg") {
                    std::cout<<"SPoutFVBkg"<<" "<<histname<<std::endl;
                    houtFV->Add(htemp);
                    break;
                }
                if(line == "ext") {
                    std::cout<<"ext"<<" "<<histname<<std::endl;
                    hext->Add(htemp);
                    hextscale->Add(htemp);
                    break;
                }
                if(line == "SPdirtBkg") {
                    std::cout<<"SPdirtBkg"<<" "<<histname<<std::endl;
                    hdirt->Add(htemp);
                    break;
                }
                if(line == "SPOutFVSig" || line == "SPNumuCCSig" || line == "SPNCOtherSig" || 
                    line == "SPNCDeltaSig" || line == "SPNCPi0Sig") {
                    std::cout<<"1gSigScale"<<" "<<histname<<std::endl;
                    h1gscale->Add(htemp);
                    //break;
                }
                if(line == "SPNCPi0Sig") {
                    std::cout<<"SPNCPi0Sig"<<" "<<histname<<std::endl;
                    hNCpi1g->Add(htemp);
                    hNCpi1gscale->Add(htemp);
                    break;
                }
                if(line == "SPNCDeltaSig") {
                    std::cout<<"SPNCDeltaSig"<<" "<<histname<<std::endl;
                    hNCdel->Add(htemp);
                    hNCdelscale->Add(htemp);
                    break;
                }
                if(line == "SPNCOtherSig") {
                    std::cout<<"SPNCOtherSig"<<" "<<histname<<std::endl;
                    hNCother->Add(htemp);
                    hNCotherscale->Add(htemp);
                    break;
                }
                if(line == "SPNumuCCSig") {
                    std::cout<<"SPNumuCCSig"<<" "<<histname<<std::endl;
                    hnumuCC1g->Add(htemp);
                    hnumuCC1gscale->Add(htemp);
                    break;
                }
                if(line == "SPOutFVSig") {
                    std::cout<<"SPOutFVSig"<<" "<<histname<<std::endl;
                    hout1g->Add(htemp);
                    hout1gscale->Add(htemp);
                    break;
                }
                if(line == "LEE") {
                    std::cout<<"LEE"<<" "<<histname<<std::endl;
                    flag_leeexist = true;
                    hLEE->Add(htemp);
                    hLEEscale->Add(htemp);
                    break;
                }
                if(line == "CCQE") {
                    std::cout<<"CCQE"<<" "<<histname<<std::endl;
                    hCCQE->Add(htemp);
                    break;
                }
                if(line == "NCQE") {
                    std::cout<<"NCQE"<<" "<<histname<<std::endl;
                    hNCQE->Add(htemp);
                    break;
                }
                if(line == "CCRES") {
                    std::cout<<"CCRES"<<" "<<histname<<std::endl;
                    hCCRES->Add(htemp);
                    break;
                }
                if(line == "NCRES") {
                    std::cout<<"NCRES"<<" "<<histname<<std::endl;
                    hNCRES->Add(htemp);
                    break;
                }
                if(line == "CCDIS") {
                    std::cout<<"CCDIS"<<" "<<histname<<std::endl;
                    hCCDIS->Add(htemp);
                    break;
                }
                if(line == "NCDIS") {
                    std::cout<<"NCDIS"<<" "<<histname<<std::endl;
                    hNCDIS->Add(htemp);
                    break;
                }
                if(line == "CCMEC") {
                    std::cout<<"CCMEC"<<" "<<histname<<std::endl;
                    hCCMEC->Add(htemp);
                    break;
                }
                if(line == "NCMEC") {
                    std::cout<<"NCMEC"<<" "<<histname<<std::endl;
                    hNCMEC->Add(htemp);
                    break;
                }
                if(line == "OTHER") {
                    std::cout<<"OTHER"<<" "<<histname<<std::endl;
                    hOTHER->Add(htemp);
                    break;
                }
                if(line == "XsecCosmic") {
                    std::cout<<"XsecCosmic"<<" "<<histname<<std::endl;
                    hXsecCosmic->Add(htemp);
                    break;
                }
                if(line == "XsecNumuCCinFV") {
                    std::cout<<"XsecNumuCCinFV"<<" "<<histname<<std::endl;
                    hXsecNumuCCinFV->Add(htemp);
                    break;
                }
                if(line == "XsecNC") {
                    std::cout<<"XsecNC"<<" "<<histname<<std::endl;
                    hXsecNC->Add(htemp);
                    break;
                }
                if(line == "XsecBkgCC") {
                    std::cout<<"XsecBkgCC"<<" "<<histname<<std::endl;
                    hXsecBkgCC->Add(htemp);
                    break;
                }
            }
        }
        //pad1->cd();
        canvas[obschannel-1]->cd();
        gStyle->SetOptTitle(kFALSE);
        float datapot = 0;
        float datapot_numi = 0;
        if(run == 0){
            for(auto it=map_data_period_pot.begin(); it!=map_data_period_pot.end(); it++)
            {
                if(it->first<10) datapot += it->second;
                if(it->first>10) datapot_numi += it->second;
            }
        }else datapot = map_data_period_pot[run];

        gr[obschannel-1] = new TGraphAsymmErrors();
        legend[obschannel-1]->SetNColumns(2);
        //legend3[obschannel-1]->SetNColumns(2);
        // numi channels
        if(obschannel==999) legend[obschannel-1]->AddEntry((TObject*)0, Form("NuMI POT: %.3e", datapot_numi*scalePOT), "");
        else {
          legend[obschannel-1]->AddEntry((TObject*)0, Form("Data POT: %.3e", datapot*scalePOT), "");
          legend3[obschannel-1]->AddEntry((TObject*)0, Form("Data POT: %.3e", datapot*scalePOT), "");
        }
        //else legend[obschannel-1]->AddEntry((TObject*)0, Form("Scaled to POT: %.3e", datapot*scalePOT), "");
        //remove for style
        //legend[obschannel-1]->AddEntry((TObject*)0, Form("#chi^{2}/ndf=%.2f/%d", GOF[obschannel-1].first, GOF[obschannel-1].second), "");
      //legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("BNB data, %.1f", hdata->Integral()*scalePOT), "lp");
        //legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("Fake data, %.1f", hdata->Integral()*scalePOT), "lp");
        //legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("Scaled BNB data, %.1f", hdata->Integral()*scalePOT), "lp");

        TH1F* hmc = (TH1F*)map_obsch_histos[obschannel].at(1)->Clone("hmc");
        TH1F* hmc2 = (TH1F*)map_obsch_histos[obschannel].at(2)->Clone("hmc2");
        TH1F* hmcerror = (TH1F*)hmc->Clone("hmcerror");
        //remove for style
        //legend[obschannel-1]->AddEntry(hmcerror, "Pred. uncertainty", "lf");

        if(flag_truthlabel==0){
        // truth labels start
        hstack[obschannel-1]->Add(hbadmatch);
      //  legend[obschannel-1]->AddEntry(hbadmatch, Form("Cosmic, %.1f", hbadmatch->Integral()), "F");
        hbadmatch->SetFillStyle(3004);
        hbadmatch->SetFillColorAlpha(kRed+2, 0.5);
        hbadmatch->SetLineColor(kRed+2);
        hbadmatch->SetLineWidth(4);

        hstack[obschannel-1]->Add(hext);
      //  legend[obschannel-1]->AddEntry(hext, Form("EXT, %.1f", hext->Integral()), "F");
        hext->SetFillStyle(3004);
        hext->SetFillColorAlpha(kOrange+3, 0.5);
        hext->SetLineColor(kOrange+3);
        hext->SetLineWidth(4);

        hstack[obschannel-1]->Add(hdirt);
      //  legend[obschannel-1]->AddEntry(hdirt, Form("Dirt, %.1f", hdirt->Integral()), "F");
        hdirt->SetFillStyle(3224);
        hdirt->SetFillColorAlpha(kGray, 0.5);
        hdirt->SetLineColor(kGray+2);
        hdirt->SetLineWidth(4);

        hstack[obschannel-1]->Add(houtFV);
      //  legend[obschannel-1]->AddEntry(houtFV, Form("out FV, %.1f", houtFV->Integral()), "F");
        houtFV->SetFillStyle(3224);
        houtFV->SetFillColorAlpha(kOrange+1, 0.5);
        houtFV->SetLineColor(kOrange+1);
        houtFV->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCpi0inFV);
      //  legend[obschannel-1]->AddEntry(hNCpi0inFV, Form("NC #pi^{0} in FV,  %.1f", hNCpi0inFV->Integral()), "F");
        hNCpi0inFV->SetFillStyle(1001);
        hNCpi0inFV->SetFillColorAlpha(38, 0.5);
        hNCpi0inFV->SetLineColor(38);
        hNCpi0inFV->SetLineWidth(4);

        hstack[obschannel-1]->Add(hCCpi0inFV);
      //  legend[obschannel-1]->AddEntry(hCCpi0inFV, Form("CC #pi^{0} in FV, %.1f", hCCpi0inFV->Integral()), "F");
        hCCpi0inFV->SetFillStyle(1001);
        hCCpi0inFV->SetFillColorAlpha(30, 0.5);
        hCCpi0inFV->SetLineColor(30);
        hCCpi0inFV->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCinFV);
       // legend[obschannel-1]->AddEntry(hNCinFV, Form("NC in FV, %.1f", hNCinFV->Integral()), "F");
        hNCinFV->SetFillStyle(1001);
        hNCinFV->SetFillColorAlpha(kOrange+1, 0.5);
        hNCinFV->SetLineColor(kOrange+1);
        hNCinFV->SetLineWidth(4);

        hstack[obschannel-1]->Add(hnumuCCinFV);
      //  legend[obschannel-1]->AddEntry(hnumuCCinFV, Form("#nu_{#mu} CC in FV, %.1f", hnumuCCinFV->Integral()), "F");
        hnumuCCinFV->SetFillStyle(1001);
        hnumuCCinFV->SetFillColorAlpha(kAzure+6, 0.5);
        hnumuCCinFV->SetLineColor(kAzure+6);
        hnumuCCinFV->SetLineWidth(4);

      /*  hstack[obschannel-1]->Add(hRnumuCCinFV);
        legend[obschannel-1]->AddEntry(hRnumuCCinFV, Form("#nu_{#mu} CC in FV, %.1f", hRnumuCCinFV->Integral()), "F");
        hRnumuCCinFV->SetFillStyle(1001);
        hRnumuCCinFV->SetFillColorAlpha(kAzure+6, 0.5);
        hRnumuCCinFV->SetLineColor(kAzure+6);
        hRnumuCCinFV->SetLineWidth(4);

        hstack[obschannel-1]->Add(hAnumuCCinFV);
        legend[obschannel-1]->AddEntry(hAnumuCCinFV, Form("Anti #nu_{#mu} CC in FV, %.1f", hAnumuCCinFV->Integral()), "F");
        hAnumuCCinFV->SetFillStyle(1001);
        hAnumuCCinFV->SetFillColorAlpha(kAzure-6, 0.5);
        hAnumuCCinFV->SetLineColor(kAzure+6);
        hAnumuCCinFV->SetLineWidth(4);*/

        hstack[obschannel-1]->Add(hnueCCinFV);
      //  legend[obschannel-1]->AddEntry(hnueCCinFV, Form("#nu_{e} CC in FV, %.1f", hnueCCinFV->Integral()), "F");
        hnueCCinFV->SetFillStyle(1001);
        hnueCCinFV->SetFillColorAlpha(kGreen+1, 0.5);
        hnueCCinFV->SetLineColor(kGreen+1);
        hnueCCinFV->SetLineWidth(4);

        TList *stackHists = hstack[obschannel-1]->GetHists();
        if (obschannel == 1){
            float num_bkg = 0.0;
            for (int i=0;i<stackHists->GetSize();++i) {
                num_bkg+=((TH1*)stackHists->At(i))->Integral();
            }
            //cout<< "h1gscale = " << h1gscale->Integral() << endl;
            if (fit_scale == 0){
              scale_amount = abs(hdata->Integral() - num_bkg) / h1gscale->Integral();
            }else{
              scale_amount = fit_scale + 1.0;
            }

            float num_bkg_NCpi1g = num_bkg + hNCdel->Integral() + hNCother->Integral() + hnumuCC1g->Integral() + hout1g->Integral();
            float num_bkg_NCdel = num_bkg + hNCpi1g->Integral() + hNCother->Integral() + hnumuCC1g->Integral() + hout1g->Integral();
            float num_bkg_NCother = num_bkg + hNCpi1g->Integral() + hNCdel->Integral() + hnumuCC1g->Integral() + hout1g->Integral();
            float num_bkg_numuCC1g = num_bkg + hNCpi1g->Integral() + hNCdel->Integral() + hNCother->Integral() + hout1g->Integral();
            float num_bkg_out1g = num_bkg + hNCpi1g->Integral() + hNCdel->Integral() + hNCother->Integral() + hnumuCC1g->Integral();
            float num_bkg_ext = num_bkg - hext->Integral() + hNCpi1g->Integral() + hNCdel->Integral() + hNCother->Integral() + hnumuCC1g->Integral() + hout1g->Integral();
            float num_bkg_LEE = num_bkg + h1gscale->Integral();

            if (fit_scale == 0){
              scale_amount_NCpi1g   = abs(hdata->Integral() - num_bkg_NCpi1g) /   hNCpi1gscale->Integral();
              scale_amount_NCdel    = abs(hdata->Integral() - num_bkg_NCdel) /    hNCdelscale->Integral();
              scale_amount_NCother  = abs(hdata->Integral() - num_bkg_NCother) /  hNCotherscale->Integral();
              scale_amount_numuCC1g = abs(hdata->Integral() - num_bkg_numuCC1g) / hnumuCC1gscale->Integral();
              scale_amount_out1g    = abs(hdata->Integral() - num_bkg_out1g) /    hout1gscale->Integral();
              scale_amount_ext      = abs(hdata->Integral() - num_bkg_ext) /      hextscale->Integral();
              scale_amount_LEE      = abs(hdata->Integral() - num_bkg_LEE) /      hLEEscale->Integral();
            }else{
              scale_amount_NCpi1g   = ((h1gscale->Integral() / hNCpi1gscale->Integral()) * fit_scale) + 1.0;
              scale_amount_NCdel    = ((h1gscale->Integral() / hNCdelscale->Integral()) * fit_scale) + 1.0;
              scale_amount_NCother  = ((h1gscale->Integral() / hNCotherscale->Integral()) * fit_scale) + 1.0;
              scale_amount_numuCC1g = ((h1gscale->Integral() / hnumuCC1gscale->Integral()) * fit_scale) + 1.0;
              scale_amount_out1g    = ((h1gscale->Integral() / hout1gscale->Integral()) * fit_scale) + 1.0;
              scale_amount_ext      = ((h1gscale->Integral() / hextscale->Integral()) * fit_scale) + 1.0;
              scale_amount_LEE      = fit_scale;
            }
        }
            //cout<< "scale_amount = " << scale_amount << endl;
            h1gscale->Scale(abs(scale_amount - 1.0));
            signal_amount = h1gscale->Integral();
            //TList *stackHists = hstack[obschannel-1]->GetHists();
            /*for (int i=0;i<stackHists->GetSize();++i) {
                h1gscale->Add((TH1*)stackHists->At(i));
            }*/
            //legend[obschannel-1]->AddEntry(h1gscale, Form("Signal #times %.1f, %.1f", scale_amount, signal_amount), "F");
            h1gscale->SetFillStyle(0);
            //h1gscale->SetFillColorAlpha(kMagenta, 0.5);
            h1gscale->SetLineColor(colors.At(0));
            h1gscale->SetLineStyle(4);//2);
            h1gscale->SetLineWidth(8);

            hNCpi1gscale->Scale(abs(scale_amount_NCpi1g - 1.0));
            signal_amount_NCpi1g = hNCpi1gscale->Integral();
            /*for (int i=0;i<stackHists->GetSize();++i) {
                hNCpi1gscale->Add((TH1*)stackHists->At(i));
            }
            hNCpi1gscale->Add(hNCdel);
            hNCpi1gscale->Add(hNCother);
            hNCpi1gscale->Add(hnumuCC1g);
            hNCpi1gscale->Add(hout1g);*/
            hNCpi1gscale->SetFillStyle(0);
            //hNCpi1gscale->SetFillColorAlpha(kPink+5, 0.5);
            hNCpi1gscale->SetLineColor(colors.At(42));
            hNCpi1gscale->SetLineStyle(1);//4);
            hNCpi1gscale->SetLineWidth(8);

            hNCdelscale->Scale(abs(scale_amount_NCdel - 1.0));
            signal_amount_NCdel = hNCdelscale->Integral();
           /* for (int i=0;i<stackHists->GetSize();++i) {
                hNCdelscale->Add((TH1*)stackHists->At(i));
            }
            hNCdelscale->Add(hNCpi1g);
            hNCdelscale->Add(hNCother);
            hNCdelscale->Add(hnumuCC1g);
            hNCdelscale->Add(hout1g); */
            hNCdelscale->SetFillStyle(0);
            //hNCdelscale->SetFillColorAlpha(kPink-6, 0.5);
            hNCdelscale->SetLineColor(colors.At(84));
            hNCdelscale->SetLineStyle(6);
            hNCdelscale->SetLineWidth(8);

            hNCotherscale->Scale(abs(scale_amount_NCother - 1.0));
            signal_amount_NCother = hNCotherscale->Integral();
            /*for (int i=0;i<stackHists->GetSize();++i) {
                hNCotherscale->Add((TH1*)stackHists->At(i));
            }
            hNCotherscale->Add(hNCpi1g);
            hNCotherscale->Add(hNCdel);
            hNCotherscale->Add(hnumuCC1g);
            hNCotherscale->Add(hout1g);*/
            hNCotherscale->SetFillStyle(0);
            //hNCotherscale->SetFillColorAlpha(kPink-8, 0.5);
            hNCotherscale->SetLineColor(colors.At(126));
            hNCotherscale->SetLineStyle(8);
            hNCotherscale->SetLineWidth(8);

            hnumuCC1gscale->Scale(abs(scale_amount_numuCC1g - 1.0));
            signal_amount_numuCC1g = hnumuCC1gscale->Integral();
            /*for (int i=0;i<stackHists->GetSize();++i) {
                hnumuCC1gscale->Add((TH1*)stackHists->At(i));
            }
            hnumuCC1gscale->Add(hNCpi1g);
            hnumuCC1gscale->Add(hNCdel);
            hnumuCC1gscale->Add(hNCother);
            hnumuCC1gscale->Add(hout1g);*/
            hnumuCC1gscale->SetFillStyle(0);
            //hnumuCC1gscale->SetFillColorAlpha(kPink-7, 0.5);
            hnumuCC1gscale->SetLineColor(colors.At(168));
            hnumuCC1gscale->SetLineStyle(9);
            hnumuCC1gscale->SetLineWidth(8);

            hout1gscale->Scale(abs(scale_amount_out1g - 1.0));
            signal_amount_out1g = hout1gscale->Integral();
            /*for (int i=0;i<stackHists->GetSize();++i) {
                hout1gscale->Add((TH1*)stackHists->At(i));
            }
            hout1gscale->Add(hNCpi1g);
            hout1gscale->Add(hNCdel);
            hout1gscale->Add(hNCother);
            hout1gscale->Add(hnumuCC1g);*/
            hout1gscale->SetFillStyle(0);
            //hout1gscale->SetFillColorAlpha(kPink, 0.5);
            hout1gscale->SetLineColor(colors.At(210));
            hout1gscale->SetLineStyle(10);
            hout1gscale->SetLineWidth(8);

            //EXT scale
            /*hextscale->Scale(abs(scale_amount_ext - 1.0));
            signal_amount_ext = hextscale->Integral();
            hextscale->SetFillStyle(0);
            hextscale->SetLineColor(kGray+2);
            hextscale->SetLineStyle(1);
            hextscale->SetLineWidth(4);*/

            if (flag_leeexist){
              hLEEscale->Scale(abs(scale_amount_LEE));
              signal_amount_LEE = hLEEscale->Integral();
              hLEEscale->SetFillStyle(0);
              //hLEEscale->SetFillColorAlpha(kMagenta, 0.5);
              hLEEscale->SetLineColor(kViolet);//colors.At(252));
              hLEEscale->SetLineStyle(2);//1);
              hLEEscale->SetLineWidth(8);
            }
        //}

        if (0){

        hstack[obschannel-1]->Add(hNCpi1g);
        legend[obschannel-1]->AddEntry(hNCpi1g, Form("NC #pi^{0} 1#gamma, %.1f", hNCpi1g->Integral()), "F");
        hNCpi1g->SetFillStyle(1001);
        hNCpi1g->SetFillColorAlpha(kPink+5, 0.5);
        hNCpi1g->SetLineColor(kPink+5);
        hNCpi1g->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCdel);
        legend[obschannel-1]->AddEntry(hNCdel, Form("NC #Delta 1#gamma, %.1f", hNCdel->Integral()), "F");
        hNCdel->SetFillStyle(1001);
        hNCdel->SetFillColorAlpha(kPink-6, 0.5);
        hNCdel->SetLineColor(kPink-6);
        hNCdel->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCother);
        legend[obschannel-1]->AddEntry(hNCother, Form("NC Other 1#gamma, %.1f", hNCother->Integral()), "F");
        hNCother->SetFillStyle(1001);
        hNCother->SetFillColorAlpha(kPink-8, 0.5);
        hNCother->SetLineColor(kPink-8);
        hNCother->SetLineWidth(4);

        hstack[obschannel-1]->Add(hnumuCC1g);
        legend[obschannel-1]->AddEntry(hnumuCC1g, Form("#nu_{#mu}CC 1#gamma #mu<100MeV, %.1f", hnumuCC1g->Integral()), "F");
        hnumuCC1g->SetFillStyle(1001);
        hnumuCC1g->SetFillColorAlpha(kPink-7, 0.5);
        hnumuCC1g->SetLineColor(kPink-7);
        hnumuCC1g->SetLineWidth(4);

        hstack[obschannel-1]->Add(hout1g);
        legend[obschannel-1]->AddEntry(hout1g, Form("out of FV 1#gamma, %.1f", hout1g->Integral()), "F");
        hout1g->SetFillStyle(1001);
        hout1g->SetFillColorAlpha(kPink, 0.5);
        hout1g->SetLineColor(kPink);
        hout1g->SetLineWidth(4);

        

        if(flag_leeexist){
        hstack[obschannel-1]->Add(hLEE);
        legend[obschannel-1]->AddEntry(hLEE, Form("LEE, %.1f", hLEE->Integral()), "F");
        hLEE->SetFillStyle(1001);
        hLEE->SetFillColorAlpha(kMagenta, 0.5);
        hLEE->SetLineColor(kMagenta);
        hLEE->SetLineWidth(4);
        }
        // truth labels end
        }//turn off breakdown
        }
        if(flag_truthlabel==1){
        hstack[obschannel-1]->Add(hbadmatch);
        legend[obschannel-1]->AddEntry(hbadmatch, Form("Cosmic, %.1f", hbadmatch->Integral()), "F");
        hbadmatch->SetFillStyle(3004);
        hbadmatch->SetFillColorAlpha(kRed+2, 0.5);
        hbadmatch->SetLineColor(kRed+2);
        hbadmatch->SetLineWidth(4);

        hstack[obschannel-1]->Add(hext);
        legend[obschannel-1]->AddEntry(hext, Form("EXT, %.1f", hext->Integral()), "F");
        hext->SetFillStyle(3004);
        hext->SetFillColorAlpha(kOrange+3, 0.5);
        hext->SetLineColor(kOrange+3);
        hext->SetLineWidth(4);

        hstack[obschannel-1]->Add(hdirt);
        legend[obschannel-1]->AddEntry(hdirt, Form("Dirt, %.1f", hdirt->Integral()), "F");
        hdirt->SetFillStyle(3224);
        hdirt->SetFillColorAlpha(kGray, 0.5);
        hdirt->SetLineColor(kGray+2);
        hdirt->SetLineWidth(4);

        hstack[obschannel-1]->Add(hOTHER);
        legend[obschannel-1]->AddEntry(hOTHER, Form("OTHER, %.1f", hOTHER->Integral()), "F");
        hOTHER->SetFillStyle(3224);
        hOTHER->SetFillColorAlpha(kOrange+1, 0.5);
        hOTHER->SetLineColor(kOrange+1);
        hOTHER->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCDIS);
        legend[obschannel-1]->AddEntry(hNCDIS, Form("NCDIS,  %.1f", hNCDIS->Integral()), "F");
        hNCDIS->SetFillStyle(1001);
        hNCDIS->SetFillColorAlpha(38, 0.5);
        hNCDIS->SetLineColor(38);
        hNCDIS->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCMEC);
        legend[obschannel-1]->AddEntry(hNCMEC, Form("NCMEC,  %.1f", hNCMEC->Integral()), "F");
        hNCMEC->SetFillStyle(1001);
        hNCMEC->SetFillColorAlpha(30, 0.5);
        hNCMEC->SetLineColor(30);
        hNCMEC->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCRES);
        legend[obschannel-1]->AddEntry(hNCRES, Form("NCRES, %.1f", hNCRES->Integral()), "F");
        hNCRES->SetFillStyle(1001);
        hNCRES->SetFillColorAlpha(kOrange+1, 0.5);
        hNCRES->SetLineColor(kOrange+1);
        hNCRES->SetLineWidth(4);

        hstack[obschannel-1]->Add(hNCQE);
        legend[obschannel-1]->AddEntry(hNCQE, Form("NCQE, %.1f", hNCQE->Integral()), "F");
        hNCQE->SetFillStyle(1001);
        hNCQE->SetFillColorAlpha(kAzure+6, 0.5);
        hNCQE->SetLineColor(kAzure+6);
        hNCQE->SetLineWidth(4);

        hstack[obschannel-1]->Add(hCCDIS);
        legend[obschannel-1]->AddEntry(hCCDIS, Form("CCDIS,  %.1f", hCCDIS->Integral()), "F");
        hCCDIS->SetFillStyle(1001);
        hCCDIS->SetFillColorAlpha(kGray+5, 0.5);
        hCCDIS->SetLineColor(kGray+5);
        hCCDIS->SetLineWidth(4);

        hstack[obschannel-1]->Add(hCCMEC);
        legend[obschannel-1]->AddEntry(hCCMEC, Form("CCMEC,  %.1f", hCCMEC->Integral()), "F");
        hCCMEC->SetFillStyle(1001);
        hCCMEC->SetFillColorAlpha(kAzure, 0.5);
        hCCMEC->SetLineColor(kAzure);
        hCCMEC->SetLineWidth(4);

        hstack[obschannel-1]->Add(hCCRES);
        legend[obschannel-1]->AddEntry(hCCRES, Form("CCRES, %.1f", hCCRES->Integral()), "F");
        hCCRES->SetFillStyle(1001);
        hCCRES->SetFillColorAlpha(kMagenta-5, 0.5);
        hCCRES->SetLineColor(kMagenta-5);
        hCCRES->SetLineWidth(4);

        hstack[obschannel-1]->Add(hCCQE);
        legend[obschannel-1]->AddEntry(hCCQE, Form("CCQE, %.1f", hCCQE->Integral()), "F");
        hCCQE->SetFillStyle(1001);
        hCCQE->SetFillColorAlpha(kGreen+1, 0.5);
        hCCQE->SetLineColor(kGreen+1);
        hCCQE->SetLineWidth(4);

        if(flag_leeexist){
        hstack[obschannel-1]->Add(hLEE);
        legend[obschannel-1]->AddEntry(hLEE, Form("LEE, %.1f", hLEE->Integral()), "F");
        hLEE->SetFillStyle(1001);
        hLEE->SetFillColorAlpha(kMagenta, 0.5);
        hLEE->SetLineColor(kMagenta);
        hLEE->SetLineWidth(4);
        }
        }
        if(flag_truthlabel==2){ // xsec style
        // truth labels start
        hstack[obschannel-1]->Add(hext);
        legend[obschannel-1]->AddEntry(hext, Form("Beam-off, %.1f", hext->Integral()), "F");
        hext->SetFillStyle(3004);
        hext->SetFillColorAlpha(kOrange+3, 0.5);
        hext->SetLineColor(kOrange+3);
        hext->SetLineWidth(4);

        hstack[obschannel-1]->Add(hdirt);
        legend[obschannel-1]->AddEntry(hdirt, Form("Dirt, %.1f", hdirt->Integral()), "F");
        hdirt->SetFillStyle(3224);
        hdirt->SetFillColorAlpha(kGray, 0.5);
        hdirt->SetLineColor(kGray+2);
        hdirt->SetLineWidth(4);

        hstack[obschannel-1]->Add(hXsecCosmic);
        legend[obschannel-1]->AddEntry(hXsecCosmic, Form("Cosmic, %.1f", hXsecCosmic->Integral()), "F");
        hXsecCosmic->SetFillStyle(3224);
        hXsecCosmic->SetFillColorAlpha(kRed+2, 0.5);
        hXsecCosmic->SetLineColor(kRed+2);
        hXsecCosmic->SetLineWidth(4);

        hstack[obschannel-1]->Add(hXsecNC);
        legend[obschannel-1]->AddEntry(hXsecNC, Form("NC,  %.1f", hXsecNC->Integral()), "F");
        hXsecNC->SetFillStyle(1001);
        hXsecNC->SetFillColorAlpha(kOrange+1, 0.5);
        hXsecNC->SetLineColor(38);
        hXsecNC->SetLineWidth(4);

        hstack[obschannel-1]->Add(hXsecBkgCC);
        legend[obschannel-1]->AddEntry(hXsecBkgCC, Form("Other CC, %.1f", hXsecBkgCC->Integral()), "F");
        hXsecBkgCC->SetFillStyle(1001);
        hXsecBkgCC->SetFillColorAlpha(kGreen+1, 0.5);
        hXsecBkgCC->SetLineColor(30);
        hXsecBkgCC->SetLineWidth(4);

        hstack[obschannel-1]->Add(hXsecNumuCCinFV);
        legend[obschannel-1]->AddEntry(hXsecNumuCCinFV, Form("#nu_{#mu} CC in FV, %.1f", hXsecNumuCCinFV->Integral()), "F");
        hXsecNumuCCinFV->SetFillStyle(1001);
        hXsecNumuCCinFV->SetFillColorAlpha(kAzure+6, 0.5);
        hXsecNumuCCinFV->SetLineColor(kAzure+6);
        hXsecNumuCCinFV->SetLineWidth(4);
        // truth labels end
        }

        hmc->Sumw2();
        hmc->Scale(scalePOT);
        //hmc->Draw("hist");
        hmc->GetYaxis()->SetTitle("Event counts");
        hmc->GetYaxis()->SetTitleSize(0.05);
        hmc->GetYaxis()->SetTitleFont(22);//132);
        hmc->GetYaxis()->SetTitleOffset(0.73);
        hmc->GetYaxis()->SetLabelFont(22);//132);
        hmc->GetYaxis()->SetLabelSize(0.04);
        //if(obschannel==9) hmc->GetXaxis()->SetRangeUser(0.5,1);
        float mcymax = hmc->GetBinContent(hmc->GetMaximumBin())*scalePOT;
        float dataymax = hdata->GetBinContent(hdata->GetMaximumBin())*scalePOT/normalization;
        if(dataymax>mcymax) mcymax = dataymax;
        hmc->SetMaximum(2.0*mcymax);
        hmc->GetYaxis()->SetRangeUser(-0.02*mcymax, 1.6*mcymax);
        hmc->SetLineColor(kBlack);
        hmc->SetLineWidth(5);

        //Erin, h1gscale comment
        //hstack[obschannel-1]->Draw("hist same");
        hmcerror->Sumw2();
        hmcerror->Scale(scalePOT);
        //Erin, h1gscale comment
        //hmcerror->Draw("same E2");
        hmcerror->SetFillColor(kGray+2);
        hmcerror->SetFillStyle(3002);
        hmcerror->SetLineWidth(0);
        hmcerror->SetLineColor(12);
        hmcerror->SetMarkerColor(0);
        hmcerror->SetMarkerSize(0);

        //Erin , plot data as just excess
        if (flag_leeexist){
          //REMOVE: looking at scaling to excess - LEE
          //hmc->Add(hLEE,-1);
        }
        hdata->Add(hmc,-1);

        legend[obschannel-1]->AddEntry(gr[obschannel-1], Form("BNB data, %.1f", hdata->Integral()*scalePOT), "lp");
        legend3[obschannel-1]->AddEntry(gr[obschannel-1], Form("BNB data, %.1f", hdata->Integral()*scalePOT), "lp");

        float excessymax = hdata->GetBinContent(hdata->GetMaximumBin())*scalePOT/normalization;
        float excessymin = hdata->GetBinContent(hdata->GetMinimumBin())*scalePOT/normalization;
        hdata->SetMaximum(3.0*excessymax);
        //hdata->GetYaxis()->SetRangeUser(-0.02*excessymax, 2.7*excessymax);
        if (excessymin > 0.0) excessymin = 0.0;
        hdata->GetYaxis()->SetRangeUser(17.0*excessymin - 0.02*excessymax, 2.7*excessymax);
        //hdata->GetYaxis()->SetRangeUser(30.0*excessymin - 1.0*excessymax, 3.5*excessymax);
        //hdata->GetXaxis()->SetRange(1, hdata->GetNbinsX() + 1);

        hdata->SetLineColor(kWhite);
        hdata->Draw("hist");

        canvas3[obschannel-1]->cd();
        hdata->Draw("hist");
        canvas[obschannel-1]->cd();

        //if (obschannel==1){
            //legend[obschannel-1]->SetTextSize(0.03);

            double scalechi2 = hdata->Chi2Test(h1gscale,"CHI2");
            double scaleKS = hdata->KolmogorovTest(h1gscale);
            //legend[obschannel-1]->AddEntry(h1gscale, Form("Signal #times %.1f, %.1f, #chi^{2} = %.2f", scale_amount, signal_amount, scalechi2), "F");
            //legend[obschannel-1]->AddEntry(h1gscale, Form("Signal #times %.1f, #chi^{2} = %.2f", scale_amount - 1.0, scalechi2), "l");
            //h1gscale->Draw("hist same");

            double scalechi2_NCpi1g = hdata->Chi2Test(hNCpi1gscale,"CHI2");
            double scaleKS_NCpi1g = hdata->KolmogorovTest(hNCpi1gscale);
            //legend[obschannel-1]->AddEntry(hNCpi1gscale, Form("NC #pi^{0} 1#gamma #times %.1f, #chi^{2} = %.2f", scale_amount_NCpi1g - 1.0, scalechi2_NCpi1g), "l");
            legend[obschannel-1]->AddEntry(hNCpi1gscale, Form("#splitline{NC #pi^{0} 1#gamma #times %.1f, }{KS = %.2f}", scale_amount_NCpi1g - 1.0, scaleKS_NCpi1g), "l");
            legend3[obschannel-1]->AddEntry(hNCpi1gscale, Form("NC #pi^{0} 1#gamma #times %.1f, KS = %.2f", scale_amount_NCpi1g - 1.0, scaleKS_NCpi1g), "l");
            //legend3[obschannel-1]->AddEntry(hNCpi1gscale, Form("NC #pi^{0} 1#gamma #times %.1f", scale_amount_NCpi1g - 1.0), "l");
            //hNCpi1gscale->GetXaxis()->SetRange(1, hNCpi1gscale->GetNbinsX() + 1);
            hNCpi1gscale->Draw("hist same");
            canvas3[obschannel-1]->cd();
            hNCpi1gscale->Draw("hist same");
            canvas[obschannel-1]->cd();

            double scalechi2_NCdel = hdata->Chi2Test(hNCdelscale,"CHI2");
            double scaleKS_NCdel = hdata->KolmogorovTest(hNCdelscale);
            //legend[obschannel-1]->AddEntry(hNCdelscale, Form("NC #Delta 1#gamma #times %.1f, #chi^{2} = %.2f", scale_amount_NCdel - 1.0, scalechi2_NCdel), "l");
            legend[obschannel-1]->AddEntry(hNCdelscale, Form("#splitline{NC #Delta 1#gamma #times %.1f, }{KS = %.2f}", scale_amount_NCdel - 1.0, scaleKS_NCdel), "l");
            legend3[obschannel-1]->AddEntry(hNCdelscale, Form("NC #Delta 1#gamma #times %.1f, KS = %.2f", scale_amount_NCdel - 1.0, scaleKS_NCdel), "l");
            //legend3[obschannel-1]->AddEntry(hNCdelscale, Form("NC #Delta 1#gamma #times %.1f", scale_amount_NCdel - 1.0), "l");
            //hNCdelscale->GetXaxis()->SetRange(1, hNCdelscale->GetNbinsX() + 1);
            hNCdelscale->Draw("hist same");    
            canvas3[obschannel-1]->cd();
            hNCdelscale->Draw("hist same"); 
            canvas[obschannel-1]->cd(); 

            double scalechi2_NCother = hdata->Chi2Test(hNCotherscale,"CHI2");
            double scaleKS_NCother = hdata->KolmogorovTest(hNCotherscale);
            //legend[obschannel-1]->AddEntry(hNCotherscale, Form("NC Other 1#gamma #times %.1f, #chi^{2} = %.2f", scale_amount_NCother - 1.0, scalechi2_NCother), "l");
            legend[obschannel-1]->AddEntry(hNCotherscale, Form("#splitline{NC Other 1#gamma #times %.1f, }{KS = %.2f}", scale_amount_NCother - 1.0, scaleKS_NCother), "l");
            //hNCotherscale->GetXaxis()->SetRange(1, hNCotherscale->GetNbinsX() + 1);
            hNCotherscale->Draw("hist same"); 

            double scalechi2_numuCC1g = hdata->Chi2Test(hnumuCC1gscale,"CHI2");
            double scaleKS_numuCC1g = hdata->KolmogorovTest(hnumuCC1gscale);
            //legend[obschannel-1]->AddEntry(hnumuCC1gscale, Form("#nu_{#mu}CC 1#gamma #times %.1f, #chi^{2} = %.2f", scale_amount_numuCC1g - 1.0, scalechi2_numuCC1g), "l");
            legend[obschannel-1]->AddEntry(hnumuCC1gscale, Form("#splitline{#nu_{#mu}CC 1#gamma #times %.1f, }{KS = %.2f}", scale_amount_numuCC1g - 1.0, scaleKS_numuCC1g), "l");
            //hnumuCC1gscale->GetXaxis()->SetRange(1, hnumuCC1gscale->GetNbinsX() + 1);
            hnumuCC1gscale->Draw("hist same");

            double scalechi2_out1g = hdata->Chi2Test(hout1gscale,"CHI2");
            double scaleKS_out1g = hdata->KolmogorovTest(hout1gscale);
            //legend[obschannel-1]->AddEntry(hout1gscale, Form("out of FV 1#gamma #times %.1f, #chi^{2} = %.2f", scale_amount_out1g - 1.0, scalechi2_out1g), "l");
            legend[obschannel-1]->AddEntry(hout1gscale, Form("#splitline{out of FV 1#gamma #times %.1f, }{KS = %.2f}", scale_amount_out1g - 1.0, scaleKS_out1g), "l");
            legend3[obschannel-1]->AddEntry(hout1gscale, Form("out of FV 1#gamma #times %.1f, KS = %.2f", scale_amount_out1g - 1.0, scaleKS_out1g), "l");
            //legend3[obschannel-1]->AddEntry(hout1gscale, Form("out of FV 1#gamma #times %.1f", scale_amount_out1g - 1.0), "l");
            //hout1gscale->GetXaxis()->SetRange(1, hout1gscale->GetNbinsX() + 1);
            hout1gscale->Draw("hist same");
            canvas3[obschannel-1]->cd();
            hout1gscale->Draw("hist same");
            canvas[obschannel-1]->cd(); 

            double scalechi2_ext = hdata->Chi2Test(hextscale,"CHI2");
            //legend[obschannel-1]->AddEntry(hextscale, Form("EXT #times %.1f, #chi^{2} = %.2f", scale_amount_ext - 1.0, scalechi2_ext), "l");
            //hextscale->Draw("hist same");

            if (flag_leeexist){
              double scalechi2_LEE = hdata->Chi2Test(hLEEscale,"CHI2");
              double scaleKS_LEE = hdata->KolmogorovTest(hLEEscale);
              //legend[obschannel-1]->AddEntry(hLEEscale, Form("LEE #times %.2f, #chi^{2} = %.1f", scale_amount_LEE, scalechi2_LEE), "l");
              legend[obschannel-1]->AddEntry(hLEEscale, Form("#splitline{LEE #times %.1f, }{KS = %.2f}", scale_amount_LEE, scaleKS_LEE), "l");
              //legend[obschannel-1]->AddEntry(hLEEscale, Form("LEE #times %.1f", scale_amount_LEE), "l");
              //hLEEscale->GetXaxis()->SetRange(1, hLEEscale->GetNbinsX() + 1);
              hLEEscale->Draw("hist same");
            }
        //}

        //get the constrained error
        TFile fileconstr("file_hists_3002.root");
        
        if (!fileconstr.IsOpen()){
            std::cout << "Error opening the constained ROOT file." << std::endl;
            break;
        }

        TMatrixD* pred_ptr = (TMatrixD*)fileconstr.Get("matrix_Y_under_X");
        TMatrixD* data_ptr = (TMatrixD*)fileconstr.Get("matrix_data_Y");
        TMatrixD* cov_ptr = (TMatrixD*)fileconstr.Get("matrix_YY_under_XX");

        TMatrixD pred_con = *pred_ptr;
        TMatrixD data_con = *data_ptr;
        TMatrixD cov_con = *cov_ptr;

        pred_con.T();
        data_con.T();
        cov_con.T();

        TH1D	*h1_pred_Y_wiConstraint = (TH1D*)fileconstr.Get("h1_pred_Y_wiConstraint_3002");

        fileconstr.Close();

        //assign error to the data (excess) points
        gratio_mc[obschannel-1] = new TGraphAsymmErrors();
        gratio_mc2[obschannel-1] = new TGraphAsymmErrors();
        gratio_data[obschannel-1] = new TGraphAsymmErrors();
        gratio_data2[obschannel-1] = new TGraphAsymmErrors();
        float maxratio = 1.5;
        for(int i=0; i<hdata->GetNbinsX()+1; i++)
        {
            double x = hdata->GetBinCenter(i+1);
            double x_err = hdata->GetBinWidth(1)*0.5;
            double y = hdata->GetBinContent(i+1);
            double y_original = hdata_orignal->GetBinContent(i+1);
            double y_err = hdata->GetBinError(i+1);
            auto bayesError = cov.get_bayes_errors(y_original);//cov.get_bayes_errors(y);
            double bayesError_low = bayesError.first;
            double bayesError_up = bayesError.second;
            double ymc = hmc->GetBinContent(i+1);
            double ymc_err = hmc->GetBinError(i+1);
            double ymc2_err = hmc2->GetBinError(i+1);
            double ymc_err_const = ymc_err_const = sqrt(cov_con(i,i));//h1_pred_Y_wiConstraint->GetBinError(i+1);
            double total_error_const_low = sqrt( pow(bayesError_low*scalePOT/normalization, 2) + pow(ymc_err_const/2., 2)  );
            double total_error_const_up = sqrt( pow(bayesError_up*scalePOT/normalization, 2) + pow(ymc_err_const/2., 2)  );
            gr[obschannel-1]->SetPoint(i,x,y*scalePOT/normalization);
            gratio_mc[obschannel-1]->SetPoint(i,x,1);
            gratio_mc2[obschannel-1]->SetPoint(i,x,1);
            if(ymc!=0){
                gratio_data[obschannel-1]->SetPoint(i,x,y/ymc);
                gratio_data2[obschannel-1]->SetPoint(i,x,y/ymc/normalization);
                if(maxratio<y/ymc) maxratio = y/ymc;
                gratio_mc[obschannel-1]->SetPointError(i, x_err, x_err, ymc_err/ymc, ymc_err/ymc);
                gratio_mc2[obschannel-1]->SetPointError(i, x_err, x_err, sqrt(ymc_err*ymc_err-ymc2_err*ymc2_err)/ymc, sqrt(ymc_err*ymc_err-ymc2_err*ymc2_err)/ymc);
            }
            else {
                gratio_data[obschannel-1]->SetPoint(i, x, 10); // invalid value
                gratio_data2[obschannel-1]->SetPoint(i, x, 10); // invalid value
                gratio_mc[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
                gratio_mc2[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
            if(flag_err==2 || flag_err==3){ //update data point errors
                //gr[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low*scalePOT/normalization, bayesError_up*scalePOT/normalization);
                gr[obschannel-1]->SetPointError(i, x_err, x_err, total_error_const_low, total_error_const_up);
                if(ymc!=0) gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low/ymc, bayesError_up/ymc);
                else gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
                if(ymc!=0) gratio_data2[obschannel-1]->SetPointError(i, x_err, x_err, bayesError_low/ymc/normalization, bayesError_up/ymc/normalization);
                else gratio_data2[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
            if(flag_err==1){
                gr[obschannel-1]->SetPointError(i, x_err, x_err, y_err*scalePOT/normalization, y_err*scalePOT/normalization);
                if(ymc!=0) gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, y_err/ymc, y_err/ymc);
                else gratio_data[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
                if(ymc!=0) gratio_data2[obschannel-1]->SetPointError(i, x_err, x_err, y_err/ymc/normalization, y_err/ymc/normalization);
                else gratio_data2[obschannel-1]->SetPointError(i, x_err, x_err, 0, 0);
            }
        }
        //gr[obschannel-1]->GetXaxis()->SetRange(1, gr[obschannel-1]->GetNbinsX() + 1);
        gr[obschannel-1]->Draw("P same");
        gr[obschannel-1]->SetLineWidth(2);
        gr[obschannel-1]->SetMarkerStyle(20);
        gr[obschannel-1]->SetMarkerSize(1.5);
        gr[obschannel-1]->SetLineColor(kBlack);
        canvas3[obschannel-1]->cd();
        gr[obschannel-1]->Draw("P same");
        canvas[obschannel-1]->cd();
        if(flag_check == 1){
            TH1F* hcheck_data = (TH1F*)map_obsch_histos[obschannel].at(3);
            TH1F* hcheck_pred = (TH1F*)map_obsch_histos[obschannel].at(4);
            hcheck_data->Draw("hist same");
            hcheck_data->SetLineColor(kBlack);
            hcheck_data->SetLineWidth(2);
            hcheck_data->SetLineStyle(kDashed);
            hcheck_pred->Draw("hist same");
            hcheck_pred->SetLineColor(kRed);
            hcheck_pred->SetLineWidth(2);
        }


        //legend[obschannel-1]->SetFillStyle(0);
        double relerr_data = 1./TMath::Sqrt(hdata->Integral());
        double relerr_pred = TMath::Sqrt(sumtotalcov[obschannel])/hmc->Integral();
        double data_pred_ratio = hdata->Integral()/normalization/hmc->Integral();
        //remove for style
        //legend[obschannel-1]->SetHeader(Form("#SigmaDATA/#Sigma(MC+EXT)=%.2f#pm%.2f(data err)#pm%.2f(pred err)", data_pred_ratio, relerr_data*data_pred_ratio, relerr_pred*data_pred_ratio), "C");
        legend[obschannel-1]->Draw();
        canvas[obschannel-1]->Modified();
        canvas3[obschannel-1]->cd();
        legend3[obschannel-1]->Draw();
        canvas3[obschannel-1]->Modified();
        canvas[obschannel-1]->cd();         

        //pad1->Modified();
        /*pad2->cd();
        gratio_mc[obschannel-1]->Draw("a2");
        gratio_mc[obschannel-1]->SetFillColor(kBlue-10);
        //gratio_mc[obschannel-1]->SetFillColor(kRed-10);
        gratio_mc[obschannel-1]->GetYaxis()->SetRangeUser(0,int(1.5*maxratio)<2?int(1.5*maxratio):2);
        gratio_mc[obschannel-1]->GetXaxis()->SetRangeUser(hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
        gratio_mc[obschannel-1]->GetYaxis()->SetNdivisions(210);
        //if(obschannel==5 || obschannel==6) gratio_mc[obschannel-1]->GetXaxis()->SetRangeUser(0,1200);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitle("Data/Pred");
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleOffset(0.5);
        if(flag_err==3){
        gratio_mc2[obschannel-1]->Draw("2 same");
        gratio_mc2[obschannel-1]->SetFillColor(kRed-10);
        gratio_mc2[obschannel-1]->GetYaxis()->SetRangeUser(0,int(1.5*maxratio)<2?int(1.5*maxratio):2);
        gratio_mc2[obschannel-1]->GetXaxis()->SetRangeUser(hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
        }*/
        /*if(obschannel>=5) //hard coded at this moment
        {
            gratio_mc[obschannel-1]->GetXaxis()->SetTitle("Reco #pi^{0} energy [MeV]");
        }
        else*/ 
        //Erin 
        hdata->GetXaxis()->SetTitle("Reconstructed Shower Energy (MeV)");
        //hdata->GetXaxis()->SetTitle("Reco Shower Angle [degrees]");
        //hdata->GetXaxis()->SetTitle("Reconstructed Shower Cosine Angle");
        //hdata->GetXaxis()->SetTitle("Reconstructed Shower Backwards Projected Distance (cm)");
        //hdata->GetXaxis()->SetTitle("Reco Shower Forwards Projected Distance [cm]");
        //hdata->GetXaxis()->SetTitle("Reco Shower Minumum Distance to Wall [cm]");
        //hdata->GetXaxis()->SetTitle("Reco Shower Minimum Projected Distance [cm]");

        hdata->GetYaxis()->SetTitle("Excess Events - LEE");

        //hdata->SetTitle("");

        hdata->GetXaxis()->SetTitleSize(0.055);
        hdata->GetXaxis()->SetLabelSize(0.055);
        hdata->GetYaxis()->SetTitleSize(0.055);
        hdata->GetYaxis()->SetLabelSize(0.055);
        hdata->GetXaxis()->SetTitleFont(22);//132);
        hdata->GetXaxis()->SetLabelFont(22);//132);
        hdata->GetYaxis()->SetTitleFont(22);//132);
        hdata->GetYaxis()->SetLabelFont(22);//132);

        canvas[obschannel-1]->Modified();
        canvas3[obschannel-1]->Modified();
        canvas[obschannel-1]->cd();
        //pad1->Modified();
      
        
        gratio_mc[obschannel-1]->GetXaxis()->SetTitleSize(0.1);
        gratio_mc[obschannel-1]->GetXaxis()->SetLabelSize(0.1);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleSize(0.1);
        gratio_mc[obschannel-1]->GetYaxis()->SetTitleOffset(0.35);
        gratio_mc[obschannel-1]->GetYaxis()->SetLabelSize(0.1);
        //change for no pad 2
        //gratio_data[obschannel-1]->Draw("P same");
        gratio_data[obschannel-1]->SetLineWidth(2);
        gratio_data[obschannel-1]->SetMarkerStyle(20);
        gratio_data[obschannel-1]->SetMarkerSize(1.5);
        gratio_data[obschannel-1]->SetLineColor(kBlack);
        /* gratio_data2[obschannel-1]->Draw("P same"); */
        /* gratio_data2[obschannel-1]->SetLineWidth(2); */
        /* gratio_data2[obschannel-1]->SetMarkerStyle(24); */
        /* gratio_data2[obschannel-1]->SetMarkerColor(kBlack); */
        /* gratio_data2[obschannel-1]->SetMarkerSize(1.0); */
        /* gratio_data2[obschannel-1]->SetLineColor(kBlack); */

        TH1F* hist = (TH1F*)hdata->Clone("hist");
        hist->Reset();
        hist->Scale(scalePOT);
        hist->Draw("axis same");
        hist->GetYaxis()->SetNdivisions(405);

        canvas3[obschannel-1]->cd();
        hist->Draw("axis same");
        hist->GetYaxis()->SetNdivisions(405);
        canvas[obschannel-1]->cd();

        

        TLine* line;
        line = new TLine(hdata->GetXaxis()->GetXmin(),0,hdata->GetXaxis()->GetXmax(), 0);// + hdata->GetBinWidth(1),0);
        //if(obschannel==5 || obschannel==6) line = new TLine(0,1,1200,1);
        //change for no pad2
        //line->Draw();
        line->SetLineWidth(2);
        line->SetLineStyle(kDashed);
        line->Draw("same");
        canvas3[obschannel-1]->cd();
        line->Draw("same");
        canvas[obschannel-1]->cd();

        legend2[obschannel-1] = new TLegend(0.2, 0.7, 0.8, 0.95);
        legend2[obschannel-1]->SetNColumns(2);
        if(flag_err==1) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred stat uncertainty", "F");
        if(flag_err==2) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred stat uncertainty (Bayesian)", "F");
        if(flag_err==3) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred total uncertainty", "F");
        //if(flag_err==3) legend2[obschannel-1]->AddEntry(gratio_mc[obschannel-1],"Pred stat+xsec+flux uncertainty", "F");
        if(flag_err==3) legend2[obschannel-1]->AddEntry(gratio_mc2[obschannel-1],"Pred stat+xsec+flux uncertainty", "F");
        //legend2[obschannel-1]->AddEntry(gratio_data[obschannel-1],"Data with stat. uncertainty", "lp");
        //legend2[obschannel-1]->AddEntry(gratio_data2[obschannel-1],"Data with stat. uncertainty (normalized)", "lp");
        legend2[obschannel-1]->SetTextSize(0.08);
        legend2[obschannel-1]->SetFillStyle(0);
        //change for no pad 2
       // legend2[obschannel-1]->Draw();
       // pad2->Modified();

        canvas[obschannel-1]->Print((TString)hdata->GetTitle()+"_scale_minusLEE.png");
        canvas[obschannel-1]->Print(Form("scale_canvas%d_minusLEE.pdf", obschannel));

       /* if(obschannel==1) canvas[obschannel-1]->Print("scale_selection.pdf(");
        else if(obschannel==nchannels) canvas[obschannel-1]->Print("scale_selection.pdf)");
        else canvas[obschannel-1]->Print("scale_selection.pdf");*/

        canvas3[obschannel-1]->Print("scale_selection_3_minusLEE.pdf");

        canvas3[obschannel-1]->Print((TString)hdata->GetTitle()+"_scale_3_minusLEE.png");
        canvas3[obschannel-1]->Print(Form("scale_canvas%d_3_minusLEE.pdf", obschannel));

        canvas3[obschannel-1]->Print("scale_selection_3_minusLEE.pdf");


    }
    theApp.Run();
  } // flag_breakdown end
  }


  // legacy codes to save hist/matrix for TLee
  if (0){
    // prediction ...
    TFile *file3 = new TFile("merge.root","RECREATE");

    TMatrixD* mat_collapse = cov.get_mat_collapse();
    mat_collapse->Write("mat_collapse");

    // additional covariance matrix ...
    TMatrixD* mat_add_cov = cov.get_add_cov_matrix();

    std::map<int, TH1F*> map_covch_histo;

    for (auto it = map_obsch_infos.begin(); it != map_obsch_infos.end(); it++){
      int obsch = it->first;
      std::vector<std::vector< std::tuple<double, double, double, int, double> > >  bayes_inputs = it->second;
      TH1F *htemp = (TH1F*)map_obsch_histos[obsch].at(0); // data histogram ...

      for (size_t i=0;i!=bayes_inputs.size(); i++){
	int covch = std::get<3>(bayes_inputs.at(i).front());


	// double add_sys = std::get<4>(bayes_inputs.at(i).front());
	int start_bin = cov.get_covch_startbin(covch);

	if (map_covch_histo.find(covch) == map_covch_histo.end()){
	  TH1F *hnew = (TH1F*)htemp->Clone(Form("histo_%d",covch));
	  hnew->Reset();
	  map_covch_histo[covch] = hnew;
	}
	TH1F *htemp1 = map_covch_histo[covch];

	//std::cout << covch << " " << obsch << " " << i << " " << htemp1->GetSum() << std::endl;
	for (size_t j=0;j!=bayes_inputs.at(i).size();j++){
	  htemp1->SetBinContent(j+1, htemp1->GetBinContent(j+1) + std::get<0>(bayes_inputs.at(i).at(j)));
	  htemp1->SetBinError(j+1, sqrt(pow(htemp1->GetBinError(j+1),2) + std::get<1>(bayes_inputs.at(i).at(j))));

	  (*mat_add_cov)(start_bin + j, start_bin + j) += std::get<4>(bayes_inputs.at(i).at(j));

	  //  if (covch == 8 && obsch == 1)
	  //  std::cout << "bin: " << j << " " << std::get<0>(bayes_inputs.at(i).at(j)) << " " << std::get<1>(bayes_inputs.at(i).at(j)) << " " << std::get<3>(bayes_inputs.at(i).at(j)) << std::endl;
	  //	  if (std::get<4>(bayes_inputs.at(i).at(j)) > 0)
	  //  std::cout << obsch << " " << i << " " << start_bin + j << " " <<  std::get<4>(bayes_inputs.at(i).at(j)) << std::endl;
	}

	//	std::cout << obsch << " " << bayes_inputs.size() << " " << bayes_inputs.at(0).size() << " " << i << " " << covch << " " << start_bin << std::endl;
	//  break;
      }
    }

    //for (auto it = map_name_histogram.begin(); it != map_name_histogram.end(); it++){
    //  it->second.first->SetDirectory(file3);
    // }



    //    for (auto it = map_covch_histo.begin(); it != map_covch_histo.end(); it++){
    //  it->second->SetDirectory(file3);
    // }

    for (auto it = map_obsch_histos.begin(); it != map_obsch_histos.end(); it++){
      int obsch = it->first;
      TH1F *hdata = it->second.at(0);
      TH1F *hmc = it->second.at(1);

      hdata->SetName(Form("hdata_obsch_%d",obsch));
      hmc->SetName(Form("hmc_obsch_%d",obsch));
      hdata->SetDirectory(file3);
      hmc->SetDirectory(file3);
    }

    mat_add_cov->Write("cov_mat_add");

    file3->Write();
    file3->Close();
  }

}
