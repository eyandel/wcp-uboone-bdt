#include <cstdlib>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TH1F.h"

#include "WCPLEEANA/tree_wrangler.h"

using namespace std;
using namespace LEEana;

//One can use this program to generate the NuMI flux geometry weights root file.
//It takes three arguments.
//1. path to the processed checkout CV file. it can be either overolay or intrinsic nue
//2. path to NuMI pre-generated geometry weight histograms. The file can be found on github under the same directory as this current file. File name is NuMI_Geometry_Weights_Histograms.root.
//3. specify mode of the horn, run1 is "FHC" , run3 is "RHC"
//This program produces an output file name "nucleoninexsec_FluxUnisim.root": navigate to, for example, ../processed_checkout_rootfiles/prodgenie_numi_intrinsic_nue_overlay_run1 and replace the original nucleoninexsec_FluxUnisim.root with this new output. Done!
//Example to call this program Example root -l 'extractNuMIGeometryWeights.C("checkout_prodgenie_numi_overlay_run1.root","NuMI_Geometry_Weights_Histograms.root","FHC")'

//void extractNuMIGeometryWeights(string path_to_CV, string path_to_weightHistos,string hornMode){

int main( int argc, char** argv )
{
  string path_to_CV="/data1/xqian/MicroBooNE/processed_checkout_rootfiles/checkout_prodgenie_numi_intrinsic_nue_overlay_run1.root";
  string path_to_weightHistos = "/home/xqian/wire-cell/wcp-uboone-bdt/scripts/NuMI_Geometry_Weights_Histograms.root";
  string hornMode = "FHC";
  string outfile_name = "/data1/xqian/MicroBooNE/processed_checkout_rootfiles/prodgenie_numi_intrinsic_nue_overlay_run1/nucleoninexsec_FluxUnisim.root";
  bool flag_config = false;
  std::string config_file_name="config.txt";
  char delimiter = ',';
  if (argc < 3){
    std::cout << "applyNuMIGeometryWeights -i[#input_file] -o[#out_putfile] -w[#weight_file] -m[mode]" << std::endl;
    return -1;
  }

  for (Int_t i=1;i!=argc;i++){
    switch(argv[i][1]){
    case 'i':
      path_to_CV = &argv[i][2];
      break;
    case 'o':
      outfile_name = &argv[i][2];
      break;
    case 'w':
      path_to_weightHistos = &argv[i][2];
      break;
    case 'm':
      hornMode = &argv[i][2];
      break;
    case 't':
       config_file_name = &argv[i][2];
       flag_config = true;
      break;
    case 'd':
        delimiter = argv[i][2];//In case you want to change what character you use to sperate your trees in the config
      break;
    }
  }
  
  tree_wrangler wrangler(flag_config, config_file_name, delimiter);
  tree_wrangler wrangler_pot(flag_config, config_file_name, delimiter,true);

  auto CVfile = TFile::Open(path_to_CV.c_str());
  auto weightHistosFile = TFile::Open(path_to_weightHistos.c_str());
  auto T_KINEvars_copy = (TTree*)CVfile->Get("wcpselection/T_KINEvars");//kine_reco_enu
  auto T_eval_copy = (TTree*)CVfile->Get("wcpselection/T_eval");//truth_nuPdg, truth_nuEnergy,run,subrun,event
  auto T_pot_copy = (TTree*)CVfile->Get("wcpselection/T_pot");//truth_nuPdg, truth_nuEnergy,run,subrun,event
  auto T_PFeval_copy = (TTree*)CVfile->Get("wcpselection/T_PFeval");//truth_nuPdg, truth_nuEnergy,run,subrun,event
  auto T_BDTvars_copy = (TTree*)CVfile->Get("wcpselection/T_BDTvars");//truth_nuPdg, truth_nuEnergy,run,subrun,event
  TTree *T_spacepoints_copy = (TTree*)CVfile->Get("wcpselection/T_spacepoints");

  //Load other trees from directories as specified by the config file
  std::vector<TTree*>* old_trees = new std::vector<TTree*>;
  old_trees = wrangler.get_old_trees(CVfile);
  std::vector<TTree*>* old_trees_pot = new std::vector<TTree*>;
  old_trees_pot = wrangler_pot.get_old_trees(CVfile);

  auto h_nue_FHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_FHC_nue_CV_AV_TPC");
  auto h_nue_FHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_FHC_nue_CV_AV_TPC");

  auto h_nuebar_FHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_FHC_nuebar_CV_AV_TPC");
  auto h_nuebar_FHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_FHC_nuebar_CV_AV_TPC");

  auto h_numu_FHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_FHC_numu_CV_AV_TPC");
  auto h_numu_FHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_FHC_numu_CV_AV_TPC");

  auto h_numubar_FHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_FHC_numubar_CV_AV_TPC");
  auto h_numubar_FHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_FHC_numubar_CV_AV_TPC");

  auto h_nue_RHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_RHC_nue_CV_AV_TPC");
  auto h_nue_RHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_RHC_nue_CV_AV_TPC");

  auto h_nuebar_RHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_RHC_nuebar_CV_AV_TPC");
  auto h_nuebar_RHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_RHC_nuebar_CV_AV_TPC");

  auto h_numu_RHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_RHC_numu_CV_AV_TPC");
  auto h_numu_RHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_RHC_numu_CV_AV_TPC");

  auto h_numubar_RHC_variation1=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run1_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation2=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run2_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation3=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run3_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation4=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run4_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation5=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run5_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation6=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run6_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation7=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run7_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation8=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run8_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation9=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run9_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation10=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run10_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation11=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run11_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation12=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run12_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation13=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run13_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation14=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run14_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation15=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run15_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation16=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run16_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation17=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run17_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation18=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run18_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation19=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run19_RHC_numubar_CV_AV_TPC");
  auto h_numubar_RHC_variation20=(TH1D*)weightHistosFile->Get("EnergyVarBin/ratio_run20_RHC_numubar_CV_AV_TPC");


  Int_t run, subrun, event, truth_nuPdg;
  Float_t kine_reco_Enu, truth_nuEnergy;
  std::string* file_type = new std::string;
  Float_t weight_cv, weight_spline, weight_lee;
  weight_lee=-1;
  bool mcweight_filled = true;
  //auto kine_reco_enu = new map<string,vector<float> >;
  T_eval_copy->SetBranchAddress("run", &run);
  T_eval_copy->SetBranchAddress("subrun", &subrun);
  T_eval_copy->SetBranchAddress("event", &event);
  T_eval_copy->SetBranchAddress("truth_nuPdg", &truth_nuPdg);
  T_KINEvars_copy->SetBranchAddress("kine_reco_Enu", &kine_reco_Enu);
  T_eval_copy->SetBranchAddress("truth_nuEnergy", &truth_nuEnergy);
  T_eval_copy->SetBranchAddress("file_type", &file_type);
  T_eval_copy->SetBranchAddress("weight_cv", &weight_cv);
  T_eval_copy->SetBranchAddress("weight_spline", &weight_spline);
  T_eval_copy->SetBranchAddress("weight_lee", &weight_lee);

  //new tree for beamline variation weights
  auto mcweight = new std::vector<float>;
  //  auto ofile = new TFile("nucleoninexsec_FluxUnisim.root","RECREATE");
  auto ofile = new TFile(outfile_name.c_str(),"RECREATE");

  //Setup the directories specified in the config file
  std::vector<TTree*>* new_trees = new std::vector<TTree*>;
  new_trees = wrangler.set_new_trees(ofile);
  std::vector<TTree*>* new_trees_pot = new std::vector<TTree*>;
  new_trees_pot = wrangler_pot.set_new_trees(ofile);

  ofile->mkdir("wcpselection")->cd();
  // ofile->cd();
  TTree* UBTree = nullptr;
  UBTree = new TTree("T_weight","T_weight");
  UBTree->Branch("run", &run, "run/I");
  UBTree->Branch("subrun", &subrun, "subrun/I");
  UBTree->Branch("event", &event, "event/I");
  UBTree->Branch("file_type", &file_type);
  UBTree->Branch("weight_cv", &weight_cv, "weight_cv/F");
  UBTree->Branch("weight_spline", &weight_spline, "weight_spline/F");
  UBTree->Branch("weight_lee", &weight_lee, "weight_lee/F");
  UBTree->Branch("mcweight_filled", &mcweight_filled);
  UBTree->Branch("nucleoninexsec_FluxUnisim", &mcweight);

  auto T_eval = T_eval_copy->CloneTree(0);
  auto T_PFeval = T_PFeval_copy->CloneTree(0);
  auto T_KINEvars = T_KINEvars_copy->CloneTree(0);
  auto T_BDTvars = T_BDTvars_copy->CloneTree(0);
  auto T_pot = T_pot_copy->CloneTree(0);
  TTree *T_spacepoints = T_spacepoints_copy->CloneTree(0);

  // ofile->cd("wcpselection");



  // float h0 = 0;
  // float h1 = 0;
  // float h2 = 0;
  // float h3 = 0;
  // float h4 = 0;
  // float h5 = 0;
  // float h6 = 0;
  // float h7 = 0;
  // float h8 = 0;
  // float h9 = 0;
  // float h10 = 0;
  // float h11 = 0;
  // float h12 = 0;
  // float h13 = 0;
  // float h14 = 0;
  // float h15 = 0;
  // float h16 = 0;
  // float h17 = 0;
  // float h18 = 0;
  // float h19= 0;
  // float hcv = 0;

  for (size_t i=0;i < T_pot_copy->GetEntries(); i++){
    T_pot_copy->GetEntry(i);
    T_pot->Fill();
    for(auto tree_it=old_trees_pot->begin(); tree_it!=old_trees_pot->end(); tree_it++){
        (*tree_it)->GetEntry(i);
    }

    for(auto tree_it=new_trees_pot->begin(); tree_it!=new_trees_pot->end(); tree_it++){
        (*tree_it)->Fill();
    }
  }
  

  for (size_t i=0; i<T_eval_copy->GetEntries(); i++) {
    T_eval_copy->GetEntry(i);
    T_PFeval_copy->GetEntry(i);
    T_KINEvars_copy->GetEntry(i);
    T_BDTvars_copy->GetEntry(i);
    //

    // if (truth_nuEnergy < 1200){
    //   continue;
    // }

    // if (truth_nuEnergy > 3000){
    //   continue;
    // }


    // if (kine_reco_Enu < 3000){
    //   continue;
    // }

    // if (kine_reco_Enu > 6000){
    //   continue;
    // }

    // hcv = hcv + 1;

    // weight_spline=weight_spline*sqrt(10);

    T_eval->Fill();
    T_PFeval->Fill();
    T_KINEvars->Fill();
    T_BDTvars->Fill();

    T_spacepoints_copy->GetEntry(i);
    T_spacepoints->Fill();

    for(auto tree_it=old_trees->begin(); tree_it!=old_trees->end(); tree_it++){
      (*tree_it)->GetEntry(i);
    }

    for(auto tree_it=new_trees->begin(); tree_it!=new_trees->end(); tree_it++){
      (*tree_it)->Fill();
    }



    if(truth_nuPdg==12 && hornMode=="FHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation1->GetBinContent(h_nue_FHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation2->GetBinContent(h_nue_FHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation3->GetBinContent(h_nue_FHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation4->GetBinContent(h_nue_FHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation5->GetBinContent(h_nue_FHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation6->GetBinContent(h_nue_FHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation7->GetBinContent(h_nue_FHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation8->GetBinContent(h_nue_FHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation9->GetBinContent(h_nue_FHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation10->GetBinContent(h_nue_FHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation11->GetBinContent(h_nue_FHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation12->GetBinContent(h_nue_FHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation13->GetBinContent(h_nue_FHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation14->GetBinContent(h_nue_FHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation15->GetBinContent(h_nue_FHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation16->GetBinContent(h_nue_FHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation17->GetBinContent(h_nue_FHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation18->GetBinContent(h_nue_FHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation19->GetBinContent(h_nue_FHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_FHC_variation20->GetBinContent(h_nue_FHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==-12 && hornMode=="FHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation1->GetBinContent(h_nuebar_FHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation2->GetBinContent(h_nuebar_FHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation3->GetBinContent(h_nuebar_FHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation4->GetBinContent(h_nuebar_FHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation5->GetBinContent(h_nuebar_FHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation6->GetBinContent(h_nuebar_FHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation7->GetBinContent(h_nuebar_FHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation8->GetBinContent(h_nuebar_FHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation9->GetBinContent(h_nuebar_FHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation10->GetBinContent(h_nuebar_FHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation11->GetBinContent(h_nuebar_FHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation12->GetBinContent(h_nuebar_FHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation13->GetBinContent(h_nuebar_FHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation14->GetBinContent(h_nuebar_FHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation15->GetBinContent(h_nuebar_FHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation16->GetBinContent(h_nuebar_FHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation17->GetBinContent(h_nuebar_FHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation18->GetBinContent(h_nuebar_FHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation19->GetBinContent(h_nuebar_FHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_FHC_variation20->GetBinContent(h_nuebar_FHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==14 && hornMode=="FHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation1->GetBinContent(h_numu_FHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation2->GetBinContent(h_numu_FHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation3->GetBinContent(h_numu_FHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation4->GetBinContent(h_numu_FHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation5->GetBinContent(h_numu_FHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation6->GetBinContent(h_numu_FHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation7->GetBinContent(h_numu_FHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation8->GetBinContent(h_numu_FHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation9->GetBinContent(h_numu_FHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation10->GetBinContent(h_numu_FHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation11->GetBinContent(h_numu_FHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation12->GetBinContent(h_numu_FHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation13->GetBinContent(h_numu_FHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation14->GetBinContent(h_numu_FHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation15->GetBinContent(h_numu_FHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation16->GetBinContent(h_numu_FHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation17->GetBinContent(h_numu_FHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation18->GetBinContent(h_numu_FHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation19->GetBinContent(h_numu_FHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_FHC_variation20->GetBinContent(h_numu_FHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==-14 && hornMode=="FHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation1->GetBinContent(h_numubar_FHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation2->GetBinContent(h_numubar_FHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation3->GetBinContent(h_numubar_FHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation4->GetBinContent(h_numubar_FHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation5->GetBinContent(h_numubar_FHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation6->GetBinContent(h_numubar_FHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation7->GetBinContent(h_numubar_FHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation8->GetBinContent(h_numubar_FHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation9->GetBinContent(h_numubar_FHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation10->GetBinContent(h_numubar_FHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation11->GetBinContent(h_numubar_FHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation12->GetBinContent(h_numubar_FHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation13->GetBinContent(h_numubar_FHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation14->GetBinContent(h_numubar_FHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation15->GetBinContent(h_numubar_FHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation16->GetBinContent(h_numubar_FHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation17->GetBinContent(h_numubar_FHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation18->GetBinContent(h_numubar_FHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation19->GetBinContent(h_numubar_FHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_FHC_variation20->GetBinContent(h_numubar_FHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==12 && hornMode=="RHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation1->GetBinContent(h_nue_RHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation2->GetBinContent(h_nue_RHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation3->GetBinContent(h_nue_RHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation4->GetBinContent(h_nue_RHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation5->GetBinContent(h_nue_RHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation6->GetBinContent(h_nue_RHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation7->GetBinContent(h_nue_RHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation8->GetBinContent(h_nue_RHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation9->GetBinContent(h_nue_RHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation10->GetBinContent(h_nue_RHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation11->GetBinContent(h_nue_RHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation12->GetBinContent(h_nue_RHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation13->GetBinContent(h_nue_RHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation14->GetBinContent(h_nue_RHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation15->GetBinContent(h_nue_RHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation16->GetBinContent(h_nue_RHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation17->GetBinContent(h_nue_RHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation18->GetBinContent(h_nue_RHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation19->GetBinContent(h_nue_RHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nue_RHC_variation20->GetBinContent(h_nue_RHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==-12 && hornMode=="RHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation1->GetBinContent(h_nuebar_RHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation2->GetBinContent(h_nuebar_RHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation3->GetBinContent(h_nuebar_RHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation4->GetBinContent(h_nuebar_RHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation5->GetBinContent(h_nuebar_RHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation6->GetBinContent(h_nuebar_RHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation7->GetBinContent(h_nuebar_RHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation8->GetBinContent(h_nuebar_RHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation9->GetBinContent(h_nuebar_RHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation10->GetBinContent(h_nuebar_RHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation11->GetBinContent(h_nuebar_RHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation12->GetBinContent(h_nuebar_RHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation13->GetBinContent(h_nuebar_RHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation14->GetBinContent(h_nuebar_RHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation15->GetBinContent(h_nuebar_RHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation16->GetBinContent(h_nuebar_RHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation17->GetBinContent(h_nuebar_RHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation18->GetBinContent(h_nuebar_RHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation19->GetBinContent(h_nuebar_RHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_nuebar_RHC_variation20->GetBinContent(h_nuebar_RHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==14 && hornMode=="RHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation1->GetBinContent(h_numu_RHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation2->GetBinContent(h_numu_RHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation3->GetBinContent(h_numu_RHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation4->GetBinContent(h_numu_RHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation5->GetBinContent(h_numu_RHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation6->GetBinContent(h_numu_RHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation7->GetBinContent(h_numu_RHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation8->GetBinContent(h_numu_RHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation9->GetBinContent(h_numu_RHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation10->GetBinContent(h_numu_RHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation11->GetBinContent(h_numu_RHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation12->GetBinContent(h_numu_RHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation13->GetBinContent(h_numu_RHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation14->GetBinContent(h_numu_RHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation15->GetBinContent(h_numu_RHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation16->GetBinContent(h_numu_RHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation17->GetBinContent(h_numu_RHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation18->GetBinContent(h_numu_RHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation19->GetBinContent(h_numu_RHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numu_RHC_variation20->GetBinContent(h_numu_RHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }

    if(truth_nuPdg==-14 && hornMode=="RHC"){
      for (size_t j=0; j<50; j++) {
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation1->GetBinContent(h_numubar_RHC_variation1->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation2->GetBinContent(h_numubar_RHC_variation2->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation3->GetBinContent(h_numubar_RHC_variation3->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation4->GetBinContent(h_numubar_RHC_variation4->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation5->GetBinContent(h_numubar_RHC_variation5->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation6->GetBinContent(h_numubar_RHC_variation6->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation7->GetBinContent(h_numubar_RHC_variation7->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation8->GetBinContent(h_numubar_RHC_variation8->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation9->GetBinContent(h_numubar_RHC_variation9->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation10->GetBinContent(h_numubar_RHC_variation10->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation11->GetBinContent(h_numubar_RHC_variation11->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation12->GetBinContent(h_numubar_RHC_variation12->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation13->GetBinContent(h_numubar_RHC_variation13->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation14->GetBinContent(h_numubar_RHC_variation14->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation15->GetBinContent(h_numubar_RHC_variation15->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation16->GetBinContent(h_numubar_RHC_variation16->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation17->GetBinContent(h_numubar_RHC_variation17->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation18->GetBinContent(h_numubar_RHC_variation18->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation19->GetBinContent(h_numubar_RHC_variation19->FindBin(truth_nuEnergy/1000))));
        mcweight->push_back(1+sqrt(10)*(-1+h_numubar_RHC_variation20->GetBinContent(h_numubar_RHC_variation20->FindBin(truth_nuEnergy/1000))));
      }
    }


    // h0 = h0 + mcweight->at(0);
    // h1 = h1 + mcweight->at(1);
    // h2 = h2 + mcweight->at(2);
    // h3 = h3 + mcweight->at(3);
    // h4 = h4 + mcweight->at(4);
    // h5 = h5 + mcweight->at(5);
    // h6 = h6 + mcweight->at(6);
    // h7 = h7 + mcweight->at(7);
    // h8 = h8 + mcweight->at(8);
    // h9 = h9 + mcweight->at(9);
    // h10 = h10 + mcweight->at(10);
    // h11 = h11 + mcweight->at(11);
    // h12 = h12 + mcweight->at(12);
    // h13 = h13 + mcweight->at(13);
    // h14 = h14 + mcweight->at(14);
    // h15 = h15 + mcweight->at(15);
    // h16 = h16 + mcweight->at(16);
    // h17 = h17 + mcweight->at(17);
    // h18 = h18 + mcweight->at(18);
    // h19 = h19 + mcweight->at(19);


    UBTree->Fill();
    mcweight->clear();
   }

  // float var = ((h0 - hcv)*(h0 - hcv) + (h1 - hcv)*(h1 - hcv) + (h2 - hcv)*(h2 - hcv) + (h3 - hcv)*(h3 - hcv) + (h4 - hcv)*(h4 - hcv) + (h5 - hcv)*(h5 - hcv) + (h6 - hcv)*(h6 - hcv) + (h7 - hcv)*(h7 - hcv) + (h8 - hcv)*(h8 - hcv) + (h9 - hcv)*(h9 - hcv) + (h10 - hcv)*(h10 - hcv) + (h11 - hcv)*(h11 - hcv) + (h12 - hcv)*(h12 - hcv) + (h13 - hcv)*(h13 - hcv) + (h14 - hcv)*(h14 - hcv) + (h15 - hcv)*(h15 - hcv) + (h16 - hcv)*(h16 - hcv) +  (h17 - hcv)*(h17 - hcv) + (h18 - hcv)*(h18 - hcv) + (h19 - hcv)*(h19 - hcv))/20;
  // std::cout<<"variance is"<<var<<std::endl;
  // std::cout<<"std is"<<sqrt(var)<<std::endl;
  // std::cout<<"rel uncer is"<<sqrt(var)/hcv<<std::endl;

  // ofile->cd("wcpselection");
/*
  UBTree->Write();


  T_pot->Write();
  T_eval->Write();
  T_PFeval->Write();
  T_BDTvars->Write();
  T_KINEvars->Write();
*/
  ofile->Write("",TFile::kOverwrite);
  ofile->Close();
  CVfile->Close();
  weightHistosFile->Close();

}
