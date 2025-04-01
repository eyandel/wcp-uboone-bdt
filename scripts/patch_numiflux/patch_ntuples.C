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

void patch_ntuples(std::string inputname, std::string outname = "test/test.root", bool isrhc=true, bool isfullosc=false)
{
  TFile* fin1 = new TFile(inputname.c_str(), "read");
  TTree* t1bdt = (TTree*)fin1->Get("wcpselection/T_BDTvars");
  TTree* t1kine = (TTree*)fin1->Get("wcpselection/T_KINEvars");
  TTree* t1eval = (TTree*)fin1->Get("wcpselection/T_eval");
  TTree* t1pfeval = (TTree*)fin1->Get("wcpselection/T_PFeval");
  TTree* t1pot = (TTree*)fin1->Get("wcpselection/T_pot");

  t1bdt->SetBranchStatus("*", 1);
  t1pfeval->SetBranchStatus("*", 1);
  t1kine->SetBranchStatus("*", 1);
  t1eval->SetBranchStatus("*", 1);
  t1pot->SetBranchStatus("*", 1);

  auto evalcols = get_columns(t1eval);

  // now add friend and create dataframe
  // DataFrame.Snapshot is not very well supported if the TTree is not flat (has a branch of some container type)
  // but thankfully T_eval is a flat TTree so I can do it this way, otherwise one would have to loop and then add a new branch etc
  t1eval->AddFriend(t1pfeval, "pfeval");
  RDataFrame df1eval(*t1eval);

  auto filter_g3Chasevol = [](float nvx, float nvy, float nvz){
    return !((std::abs(nvx) > 58.42 || ((nvy > 35.56 && nvz > 350) || (nvy > 65 && nvz < 350))) && nvz < 4569.9);
  };
  auto baseline = [](float dk2gen, float gen2vtx){ return (float)(dk2gen+gen2vtx); };
  auto f_dfeval = df1eval.Define("truth_nuLength", baseline, {"pfeval.mcflux_dk2gen", "pfeval.mcflux_gen2vtx"});

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
  // store the old weights as well
  f_dfeval = f_dfeval.Define("weight_spline_old", [](float wgt){return (float)wgt; }, {"weight_spline"});
  f_dfeval = f_dfeval.Redefine("weight_spline", newWgt, {"weight_spline", "truth_nuEnergy", "truth_nuLength", "truth_nuPdg", "pfeval.mcflux_vx", "pfeval.mcflux_vy", "pfeval.mcflux_vz"});
  evalcols.emplace_back("weight_spline_old");

  auto opt = RSnapshotOptions();
  opt.fMode = "RECREATE";

  auto c = f_dfeval.Count();
  c.OnPartialResult(200000, [](unsigned int n){ std::cout << "===== " << n << " =====" << std::endl; });

  f_dfeval.Snapshot("wcpselection/T_eval", outname, evalcols, opt);
  TFile* fout = new TFile(outname.c_str(), "update");
  fout->cd("wcpselection");
  (t1bdt->CloneTree())->Write("T_BDTvars");
  (t1kine->CloneTree())->Write("T_KINEvars");
  (t1pfeval->CloneTree())->Write("T_PFeval");
  (t1pot->CloneTree())->Write("T_pot");
  fout->Close();

}
