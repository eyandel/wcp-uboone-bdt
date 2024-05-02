#pragma once

#include <iostream>
#include <iomanip>
#include <string>

#include <ROOT/RDataFrame.hxx>
#include <memory>
#include "glob.h"
#include "utils.h"

using namespace ROOT;
using namespace ROOT::RDF;

ColumnNames_t get_columns(TTree* t)
{
  RDataFrame df(*t);
  return df.GetColumnNames();
}

void patch_ntuples_univs(std::string inputname, std::string outname = "test/test_univs.root", bool isrhc=true, bool isfullosc=false)
{
  TFile* fin1 = new TFile(inputname.c_str(), "read");
  TTree* t1bdt = (TTree*)fin1->Get("wcpselection/T_BDTvars");
  TTree* t1kine = (TTree*)fin1->Get("wcpselection/T_KINEvars");
  TTree* t1eval = (TTree*)fin1->Get("wcpselection/T_eval");
  TTree* t1weight = (TTree*)fin1->Get("wcpselection/T_weight");
  TTree* t1pfeval = (TTree*)fin1->Get("wcpselection/T_PFeval");
  TTree* t1pot = (TTree*)fin1->Get("wcpselection/T_pot");

  t1bdt->SetBranchStatus("*", 1);
  t1pfeval->SetBranchStatus("*", 1);
  t1kine->SetBranchStatus("*", 1);
  t1eval->SetBranchStatus("*", 1);
  t1weight->SetBranchStatus("*", 1);
  t1pot->SetBranchStatus("*", 1);

  auto evalcols = get_columns(t1eval);
  auto weightcols = get_columns(t1weight);

  // now add friend and create dataframe
  // DataFrame.Snapshot is not very well supported if the TTree is not flat (has a branch of some container type)
  // but thankfully T_eval is a flat TTree so I can do it this way, otherwise one would have to loop and then add a new branch etc
  t1eval->AddFriend(t1pfeval, "pfeval");
  RDataFrame df1eval(*t1eval);
  t1weight->AddFriend(t1pfeval, "pfeval");
  t1weight->AddFriend(t1eval, "eval");
  RDataFrame df1weight(*t1weight);

  auto filter_g3Chasevol = [](float nvx, float nvy, float nvz){
    return !((std::abs(nvx) > 58.42 || ((nvy > 35.56 && nvz > 350) || (nvy > 65 && nvz < 350))) && nvz < 4569.9);
  };
  auto baseline = [](float dk2gen, float gen2vtx){ return (float)(dk2gen+gen2vtx); };
  auto f_dfeval = df1eval.Define("truth_nuLength", baseline, {"pfeval.mcflux_dk2gen", "pfeval.mcflux_gen2vtx"});
  auto f_dfweight = df1weight.Define("truth_nuLength", baseline, {"pfeval.mcflux_dk2gen", "pfeval.mcflux_gen2vtx"});

  std::string cv_str = "ratio_cv_fhc_enu_length";
  if(isrhc)
    cv_str = "ratio_cv_rhc_enu_length";

  HistReader<TH2D> kNumuReweighter("/home/nitish/uboone/flux/ppfx_univs/ratios/prod/set1/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set1.root", cv_str+"_numu");
  HistReader<TH2D> kNueReweighter("/home/nitish/uboone/flux/ppfx_univs/ratios/prod/set1/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set1.root", cv_str+"_nue");
  HistReader<TH2D> kNumuBarReweighter("/home/nitish/uboone/flux/ppfx_univs/ratios/prod/set1/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set1.root", cv_str+"_numubar");
  HistReader<TH2D> kNueBarReweighter("/home/nitish/uboone/flux/ppfx_univs/ratios/prod/set1/ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set1.root", cv_str+"_nuebar");

  auto newWgt = [&](float wgt_spline, float nuE, float nuL, int nuPdg, float nvx, float nvy, float nvz){
    if(!filter_g3Chasevol(nvx, nvy, nvz))
      return 0.f;
    if(nuPdg == 14)
      return !isfullosc ? (float)(wgt_spline*kNumuReweighter(nuE/1000., nuL)) : (float)(wgt_spline*kNueReweighter(nuE/1000., nuL));
    if(nuPdg == 12)
      return !isfullosc ? (float)(wgt_spline*kNueReweighter(nuE/1000., nuL)) : (float)(wgt_spline*kNumuReweighter(nuE/1000., nuL));
    if(nuPdg == -14)
      return !isfullosc ? (float)(wgt_spline*kNumuBarReweighter(nuE/1000., nuL)) : (float)(wgt_spline*kNueBarReweighter(nuE/1000., nuL));
    if(nuPdg == -12)
      return !isfullosc ? (float)(wgt_spline*kNueBarReweighter(nuE/1000., nuL)) : (float)(wgt_spline*kNumuBarReweighter(nuE/1000., nuL));
    return wgt_spline;
  };

  std::vector<HistReaderFile<TH2D>> kNumuUnivReweighter;
  std::vector<HistReaderFile<TH2D>> kNueUnivReweighter;
  std::vector<HistReaderFile<TH2D>> kNumuBarUnivReweighter;
  std::vector<HistReaderFile<TH2D>> kNueBarUnivReweighter;
  std::vector<std::string> flavs = {"numu", "nue", "numubar", "nuebar"};
  std::vector<TFile*> files;
  for(int i = 1; i <= 6; i++){
    std::string finname = "/home/nitish/public_html/shared/uboone/flugg_study/comparisons/unisim/shared/prod/";
    finname += "ratios_g4_10_4_ppfx_le_g3Chase_standardbin_set"+std::to_string(i)+".root";
    files.push_back(new TFile(finname.c_str(), "read"));
  }
  // flux systematics are repeated after the first 600
  for(int i = 0; i < 1000; i++){
    for(auto &flav: flavs){
      int univ_localid = i % 100;
      int univ_globalid = i < 600 ? (i / 100) + 1 : (i / 100) - 5;

      std::string univ_str = "univs/"+flav+"/ratio_univ"+std::to_string(univ_localid)+"_fhc_enu_length_"+flav;
      if(isrhc)
        univ_str = "univs/"+flav+"/ratio_univ"+std::to_string(univ_localid)+"_rhc_enu_length_"+flav;

      // for fullosc fill the numu ratios into the nue container and so on
      if(flav == "numu"){
        if(!isfullosc)
          kNumuUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
        else
          kNueUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
      }
      if(flav == "nue"){
        if(!isfullosc)
          kNueUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
        else
          kNumuUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
      }
      if(flav == "numubar"){
        if(!isfullosc)
          kNumuBarUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
        else
          kNueBarUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
      }
      if(flav == "nuebar"){
        if(!isfullosc)
          kNueBarUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
        else
          kNumuBarUnivReweighter.push_back(HistReaderFile<TH2D>(files[univ_globalid-1], univ_str));
      }
    } // end flav
  } // end univs
  std::cout << kNumuUnivReweighter.size() << std::endl;


  // store the old weights as well
  f_dfeval = f_dfeval.Define("weight_spline_old", [](float wgt){return (float)wgt; }, {"weight_spline"});
  f_dfeval = f_dfeval.Redefine("weight_spline", newWgt, {"weight_spline", "truth_nuEnergy", "truth_nuLength", "truth_nuPdg", "pfeval.mcflux_vx", "pfeval.mcflux_vy", "pfeval.mcflux_vz"});
  f_dfweight = f_dfweight.Define("weight_spline_old", [](float wgt){return (float)wgt; }, {"weight_spline"});
  f_dfweight = f_dfweight.Redefine("weight_spline", newWgt, {"weight_spline", "eval.truth_nuEnergy", "truth_nuLength", "eval.truth_nuPdg", "pfeval.mcflux_vx", "pfeval.mcflux_vy", "pfeval.mcflux_vz"});

  auto newUnivWgt = [&](RVecF &wgt_univs, float wgt_spline, float nuE, float nuL, int nuPdg, float nvx, float nvy, float nvz){

    float new_cv_wgt = newWgt(wgt_spline, nuE, nuL, nuPdg, nvx, nvy, nvz);
    RVecF new_univ_wgts = wgt_univs*0.f;
    if(new_cv_wgt)
      new_univ_wgts = (wgt_univs*wgt_spline)/new_cv_wgt;

    RVecF new_ratio_wgts;
    if(nuPdg == 14)
      for(auto &reweighter: kNumuUnivReweighter) { new_ratio_wgts.push_back((float)reweighter(nuE/1000., nuL)); }
    if(nuPdg == 12)
      for(auto &reweighter: kNueUnivReweighter) { new_ratio_wgts.push_back((float)reweighter(nuE/1000., nuL)); }
    if(nuPdg == -14)
      for(auto &reweighter: kNumuBarUnivReweighter) { new_ratio_wgts.push_back((float)reweighter(nuE/1000., nuL)); }
    if(nuPdg == -12)
      for(auto &reweighter: kNueBarUnivReweighter) { new_ratio_wgts.push_back((float)reweighter(nuE/1000., nuL)); }
    new_univ_wgts *= new_ratio_wgts;

    std::vector<float> ret(new_univ_wgts.begin(), new_univ_wgts.end());
    return ret;
  };
  auto convertType = [](RVecF &wgt_univs){
    std::vector<float> ret(wgt_univs.begin(), wgt_univs.end());
    return ret;
  };
  f_dfweight = f_dfweight.Define("kminus_PrimaryHadronNormalization_old", convertType, {"kminus_PrimaryHadronNormalization"});
  f_dfweight = f_dfweight.Redefine("kminus_PrimaryHadronNormalization", newUnivWgt, {"kminus_PrimaryHadronNormalization", "weight_spline_old", "eval.truth_nuEnergy", "truth_nuLength", "eval.truth_nuPdg", "pfeval.mcflux_vx", "pfeval.mcflux_vy", "pfeval.mcflux_vz"});

  evalcols.emplace_back("weight_spline_old");
  weightcols.emplace_back("weight_spline_old");

  auto opt = RSnapshotOptions();
  opt.fMode = "UPDATE";

  auto c = f_dfeval.Count();
  c.OnPartialResult(200000, [](unsigned int n){ std::cout << "===== " << n << " =====" << std::endl; });

  f_dfeval.Snapshot("wcpselection/T_eval", outname, evalcols, opt);
  f_dfweight.Snapshot("wcpselection/T_weight", outname, weightcols, opt);
  TFile* fout = new TFile(outname.c_str(), "update");
  fout->cd("wcpselection");
  (t1bdt->CloneTree())->Write("T_BDTvars");
  (t1kine->CloneTree())->Write("T_KINEvars");
  (t1pfeval->CloneTree())->Write("T_PFeval");
  (t1pot->CloneTree())->Write("T_pot");
  fout->Close();

}
