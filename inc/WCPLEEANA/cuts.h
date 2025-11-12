#ifndef UBOONE_LEE_CUTS
#define UBOONE_LEE_CUTS

// define cuts here ...
#include "TCut.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TH1F.h"

#include "tagger.h"
#include "kine.h"
#include "eval.h"
#include "pfeval.h"
#include "space.h"
#include "pandora.h"
#include "lantern.h"

#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

namespace LEEana{
  // this is for the real data, for fake data this should be 1 ...
  double em_charge_scale = 0.95;
  //double em_charge_scale = 1.0;

  // correct reco neutrino energy and reco shower energy
  double get_reco_Enu_corr(KineInfo& kine, bool flag_data);
  double get_reco_showerKE_corr(PFevalInfo& pfeval, bool flag_data);

  double get_reco_Eproton(KineInfo& kine);
  double get_reco_Epion(KineInfo& kine);

  double get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern);
  double get_truth_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, TString var_name);

  bool get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern);
  bool get_rw_cut_pass(TString cut, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);
  double get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval, KineInfo& kine, TaggerInfo& tagger, std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool,  std::vector<double>, std::vector<double>  > > > rw_info, std::map<int, std::tuple< double, double, double, double > > time_info, bool flag_data=false);
  int get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine);

  bool wayToSort(int i, int j) { return i > j; };
  // WC prim, WC all, LArPID prim, LArPID all, each is sorted in order of KE
  std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> get_range_proton_KE(PFevalInfo& pfeval, SpaceInfo& space, bool return_MeV);
  std::vector<double> get_pandora_proton_KE(PandoraInfo& pandora, double TRACK_SCORE_CUT, bool return_MeV);
  std::vector<double> get_lantern_KE(LanternInfo& lantern, int pdg, double vtx_cut, bool return_MeV);

  // generic neutrino cuts
  // TCut generic_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length >15";
  bool is_generic(EvalInfo& info);

  // preselection cuts
  // TCut preselect_cut = "match_found == 1 && stm_eventtype != 0 &&stm_lowenergy ==0 && stm_LM ==0 && stm_TGM ==0 && stm_STM==0 && stm_FullDead == 0 && stm_cluster_length > 0";
  bool is_preselection(EvalInfo& info);

  // nueCC cuts
  // TCut nueCC_cut = "numu_cc_flag >=0 && nue_score > 7.0";
  bool is_nueCC(TaggerInfo& tagger_info);
  bool is_loosenueCC(TaggerInfo& tagger_info);

  bool is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data);
  bool is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data);

  // numuCC cuts
  // TCut numuCC_cut = "numu_cc_flag >=0 && numu_score > 0.9";
  bool is_numuCC(TaggerInfo& tagger_info);
  bool is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_1p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);
  bool is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval);

  bool is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data);
  bool is_numuCC_cutbased(TaggerInfo& tagger_info);

  // pio cuts (with and without vertex)
  // TCut pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_pi0(KineInfo& kine, bool flag_data);

  // must be with vertex ...
  // TCut cc_pi0_cut = "(kine_pio_flag==1 && kine_pio_vtx_dis < 9 || kine_pio_flag ==2) && kine_pio_energy_1 > 40 && kine_pio_energy_2 > 25 && kine_pio_dis_1 < 110 && kine_pio_dis_2 < 120 && kine_pio_angle > 0  && kine_pio_angle < 174 && kine_pio_mass > 22 && kine_pio_mass < 300";
  bool is_cc_pi0(KineInfo& kine, bool flag_data);




  // NC cuts
  // TCut NC_cut = "(!cosmict_flag) && numu_score < 0.0";
  bool is_NC(TaggerInfo& tagger_info);
  bool is_NCpio_sel(TaggerInfo& tagger_info, KineInfo& kine);
  bool is_NCdelta_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);

  //Erin
  bool is_singlephoton_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singleshower_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singleshower_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_numu_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_other_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_ncpi0_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_nue_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_nue_sel_allshw(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_nsbeam(PFevalInfo& pfeval, EvalInfo& eval);
  bool is_nsbeam_photon(PFevalInfo& pfeval, EvalInfo& eval);
  bool is_singlephoton_pre(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_other(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_other(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_eff_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval);
  bool is_singlephoton_oneshw(TaggerInfo& tagger_info, PFevalInfo& pfeval);


  // TCut FC_cut = "match_isFC==1";
  // TCut PC_cut = "match_isFC==0";
  bool is_FC(EvalInfo& eval);


  // TCut truth_nueCC_inside = "abs(truth_nuPdg)==12 && truth_isCC==1 && truth_vtxInside==1";
  // TCut truth_numuCC_inside = "abs(truth_nuPdg)==14 && truth_isCC==1 && truth_vtxInside==1";
  bool is_truth_nueCC_inside(EvalInfo& eval);
  bool is_truth_numuCC_inside(EvalInfo& eval);

  bool is_true_0p(PFevalInfo& pfeval);

  int mcc8_pmuon_costheta_bin(float pmuon, float costh);
  int alt_var_index(std::string var1, float val1, std::string var2, float val2, std::string config="./configurations/alt_var_xbins.txt");
  std::map<std::string, TH1F> map_var_hist; // variable name and binning


  std::vector<double> get_proton_length_bins(){
    std::vector<double> proton_length_bins = {7.65721e-06, 0.003103, 0.0091284, 0.0176375, 0.0283912, 0.0412807, 0.0562127, 0.0731138, 0.0919263, 0.112603, 0.135101, 0.159272, 0.185117, 0.212807, 0.24211, 0.273135, 0.305831, 0.340111, 0.376083, 0.413583, 0.452689, 0.493293, 0.535259, 0.57868, 0.62366, 0.670316, 0.718547, 0.768222, 0.819402, 0.871993, 0.926046, 0.981485, 1.03822, 1.09633, 1.15586, 1.21689, 1.27932, 1.34303, 1.40806, 1.47449, 1.54236, 1.61158, 1.68204, 1.75379, 1.82687, 1.90134, 1.97711, 2.05409, 2.13233, 2.21187, 2.29275, 2.37489, 2.45822, 2.54277, 2.62858, 2.71569, 2.80402, 2.89352, 2.98421, 3.07612, 3.16929, 3.26367, 3.35917, 3.45585, 3.55372, 3.65281, 3.75307, 3.85444, 3.95696, 4.06063, 4.1655, 4.2715, 4.37859, 4.48679, 4.59613, 4.70662, 4.81823, 4.9309, 5.04466, 5.15953, 5.27553, 5.39262, 5.51076, 5.62996, 5.75025, 5.87164, 5.9941, 6.11758, 6.24211, 6.3677, 6.49437, 6.62208, 6.75079, 6.88053, 7.0113, 7.14313, 7.27598, 7.40982, 7.54466, 7.68052, 7.81742, 7.95524, 8.09388, 8.23334, 8.37365, 8.51481, 8.65683, 8.79972, 8.94348, 9.08814, 9.23371, 9.38019, 9.52759, 9.67593, 9.82522, 9.97548, 10.1267, 10.2789, 10.4321, 10.5864, 10.7416, 10.8979, 11.0553, 11.2137, 11.3732, 11.5338, 11.6953, 11.8576, 12.0207, 12.1846, 12.3493, 12.5147, 12.681, 12.8481, 13.016, 13.1847, 13.3543, 13.5248, 13.696, 13.8682, 14.0412, 14.2152, 14.39, 14.5657, 14.7423, 14.9198, 15.0983, 15.2777, 15.4581, 15.6394, 15.8217, 16.0049, 16.1887, 16.3733, 16.5586, 16.7447, 16.9315, 17.119, 17.3073, 17.4963, 17.6862, 17.8767, 18.0681, 18.2603, 18.4532, 18.6469, 18.8415, 19.0368, 19.233, 19.43, 19.6278, 19.8265, 20.026, 20.2264, 20.4276, 20.6297, 20.8325, 21.0361, 21.2403, 21.4452, 21.6507, 21.857, 22.064, 22.2716, 22.48, 22.689, 22.8988, 23.1093, 23.3205, 23.5325, 23.7451, 23.9585, 24.1727, 24.3876, 24.6033, 24.8197, 25.0368, 25.2548, 25.4735, 25.693, 25.9133, 26.1343, 26.3559, 26.5782, 26.801, 27.0245, 27.2487, 27.4735, 27.6989, 27.925, 28.1517, 28.3791, 28.6072, 28.8359, 29.0653, 29.2953, 29.5261, 29.7575, 29.9896, 30.2223, 30.4558, 30.69, 30.9249, 31.1604, 31.3967, 31.6337, 31.8714, 32.1096, 32.3484, 32.5878, 32.8278, 33.0684, 33.3095, 33.5513, 33.7937, 34.0367, 34.2803, 34.5245, 34.7693, 35.0147, 35.2608, 35.5075, 35.7548, 36.0027, 36.2513, 36.5005, 36.7503, 37.0008, 37.2519, 37.5037, 37.7562, 38.0092, 38.2628, 38.5169, 38.7716, 39.0268, 39.2825, 39.5388, 39.7957, 40.0531, 40.3111, 40.5697, 40.8288, 41.0885, 41.3487, 41.6096, 41.871, 42.1329, 42.3955, 42.6586, 42.9224, 43.1867, 43.4516, 43.7171, 43.9832, 44.2498, 44.5171, 44.7848, 45.0531, 45.3218, 45.5911, 45.8609, 46.1312, 46.402, 46.6733, 46.9451, 47.2175, 47.4904, 47.7638, 48.0377, 48.3122, 48.5872, 48.8627, 49.1388, 49.4154, 49.6925, 49.9702, 50.2485, 50.5272, 50.8065, 51.0864, 51.3668, 51.6476, 51.9288, 52.2105, 52.4926, 52.7752, 53.0582, 53.3417, 53.6256, 53.9099, 54.1947, 54.48, 54.7657, 55.0518, 55.3384, 55.6255, 55.913, 56.201, 56.4895, 56.7784, 57.0678, 57.3576, 57.6479, 57.9387, 58.23, 58.5217, 58.8139, 59.1066, 59.3997, 59.6934, 59.9875, 60.2821, 60.5771, 60.8727, 61.1687, 61.4653, 61.7623, 62.0598, 62.3578, 62.6563, 62.9553, 63.2548, 63.5548, 63.8553, 64.1563, 64.4578, 64.7598, 65.0623, 65.3653, 65.6688, 65.9728, 66.2772, 66.5819, 66.8871, 67.1926, 67.4985, 67.8048, 68.1114, 68.4185, 68.726, 69.0338, 69.342, 69.6507, 69.9597, 70.2691, 70.5789, 70.8891, 71.1998, 71.5108, 71.8222, 72.134, 72.4462, 72.7588, 73.0719, 73.3853, 73.6991, 74.0134, 74.328, 74.6431, 74.9585, 75.2744, 75.5907, 75.9074, 76.2246, 76.5421, 76.8601, 77.1784, 77.4972, 77.8164, 78.1361, 78.4561, 78.7766, 79.0975, 79.4189, 79.7406, 80.0628, 80.3854, 80.7085, 81.032, 81.3559, 81.6802, 82.0048, 82.3298, 82.6551, 82.9808, 83.3068, 83.6331, 83.9598, 84.2868, 84.6142, 84.9419, 85.27, 85.5984, 85.9271, 86.2562, 86.5857, 86.9155, 87.2456, 87.5761, 87.907, 88.2382, 88.5697, 88.9016, 89.2339, 89.5665, 89.8995, 90.2328, 90.5665, 90.9005, 91.2349, 91.5697, 91.9048, 92.2403, 92.5761, 92.9123, 93.2489, 93.5858, 93.9231, 94.2608, 94.5988, 94.9372, 95.276, 95.6151, 95.9546, 96.2945, 96.6347, 96.9753, 97.3163, 97.6577, 97.9994, 98.3415, 98.6839, 99.0265, 99.3695, 99.7128, 100.056, 100.4, 100.744, 101.089, 101.434, 101.779, 102.124, 102.47, 102.816, 103.162, 103.509, 103.856, 104.203, 104.55, 104.898, 105.246, 105.595, 105.943, 106.292, 106.642, 106.991, 107.341, 107.692, 108.042, 108.393, 108.744, 109.096, 109.448, 109.8, 110.153, 110.505, 110.858, 111.212, 111.566, 111.92, 112.274, 112.629, 112.984, 113.339, 113.695, 114.051, 114.407, 114.764, 115.121, 115.478, 115.836, 116.194, 116.552, 116.91, 117.269, 117.628, 117.987, 118.346, 118.706, 119.066, 119.426, 119.786, 120.147, 120.508, 120.869, 121.231, 121.593, 121.955, 122.317, 122.68, 123.043, 123.406, 123.769, 124.133, 124.497, 124.861, 125.225, 125.59, 125.955, 126.321, 126.686, 127.052, 127.418, 127.785, 128.151, 128.518, 128.885, 129.253, 129.621, 129.989, 130.357, 130.726, 131.094, 131.464, 131.833, 132.203, 132.573, 132.943, 133.314, 133.684, 134.055, 134.427, 134.798, 135.17, 135.542, 135.914, 136.287, 136.659, 137.032, 137.406, 137.779, 138.153, 138.526, 138.901, 139.275, 139.649, 140.024, 140.399, 140.775, 141.15, 141.526, 141.902, 142.278, 142.655, 143.031, 143.408, 143.785, 144.163, 144.54, 144.918, 145.296, 145.675, 146.053, 146.432, 146.811, 147.191, 147.57, 147.95, 148.33, 148.71, 149.091, 149.472, 149.853, 150.234, 150.616, 150.997, 151.379, 151.762, 152.144, 152.527, 152.91, 153.293, 153.676, 154.06, 154.444, 154.828, 155.212, 155.596, 155.981, 156.366, 156.751, 157.136, 157.521, 157.907, 158.293, 158.679, 159.065, 159.452, 159.838, 160.225, 160.612, 161.0, 161.387, 161.775, 162.163, 162.551, 162.939, 163.328, 163.717, 164.106, 164.495, 164.884, 165.274, 165.664, 166.054, 166.444, 166.835, 167.225, 167.616, 168.007, 168.399, 168.79, 169.182, 169.574, 169.966, 170.359, 170.751, 171.144, 171.537, 171.93, 172.324, 172.717, 173.111, 173.505, 173.899, 174.294, 174.688, 175.083, 175.478, 175.873, 176.268, 176.664, 177.06, 177.455, 177.851, 178.248, 178.644, 179.041, 179.437, 179.834, 180.231, 180.629, 181.026, 181.424, 181.822, 182.22, 182.618, 183.017, 183.415, 183.814, 184.213, 184.612, 185.011, 185.411, 185.811, 186.211, 186.611, 187.011, 187.412, 187.812, 188.213, 188.614, 189.015, 189.417, 189.818, 190.22, 190.622, 191.024, 191.427, 191.829, 192.232, 192.635, 193.038, 193.441, 193.845, 194.248, 194.652, 195.056, 195.46, 195.864, 196.268, 196.673, 197.077, 197.482, 197.887, 198.292, 198.698, 199.103, 199.509, 199.914, 200.32, 200.726, 201.133, 201.539, 201.946, 202.352, 202.759, 203.166, 203.574, 203.981, 204.389, 204.796, 205.204, 205.612, 206.02, 206.429, 206.837, 207.246, 207.655, 208.064, 208.473, 208.882, 209.292, 209.702, 210.111, 210.521, 210.932, 211.342, 211.752, 212.163, 212.574, 212.985, 213.396, 213.807, 214.219, 214.63, 215.042, 215.454, 215.865, 216.278, 216.69, 217.102, 217.515, 217.927, 218.34, 218.753, 219.166, 219.579, 219.993, 220.406, 220.82, 221.234, 221.648, 222.062, 222.476, 222.89, 223.305, 223.72, 224.134, 224.549, 224.964, 225.38, 225.795, 226.211, 226.626, 227.042, 227.458, 227.874, 228.29, 228.707, 229.123, 229.54, 229.957, 230.374, 230.791, 231.208, 231.626, 232.043, 232.461, 232.879, 233.297, 233.715, 234.133, 234.552, 234.97, 235.389, 235.807, 236.226, 236.645, 237.064, 237.484, 237.903, 238.323, 238.742, 239.162, 239.582, 240.002, 240.422, 240.842, 241.263, 241.683, 242.104, 242.525, 242.946, 243.367, 243.788, 244.209, 244.631, 245.052, 245.474, 245.896, 246.318, 246.74, 247.162, 247.585, 248.007, 248.43, 248.852, 249.275, 249.698, 250.121, 250.545, 250.968, 251.392, 251.815, 252.239, 252.663, 253.087, 253.511, 253.935, 254.36, 254.784, 255.209, 255.634, 256.059, 256.484, 256.909, 257.334, 257.759, 258.185, 258.611, 259.036, 259.462, 259.888, 260.314, 260.74, 261.166, 261.593, 262.019, 262.446, 262.873, 263.3, 263.727, 264.154, 264.581, 265.008, 265.436, 265.863, 266.291, 266.719, 267.147, 267.575, 268.003, 268.431, 268.86, 269.288, 269.717, 270.146, 270.574, 271.003, 271.432, 271.862, 272.291, 272.72, 273.15, 273.58, 274.01, 274.44, 274.87, 275.3, 275.73, 276.16, 276.591, 277.021, 277.452, 277.883, 278.314, 278.745, 279.176, 279.607, 280.038, 280.47, 280.901, 281.333, 281.764, 282.196, 282.628, 283.06, 283.492, 283.924, 284.356, 284.789, 285.221, 285.654, 286.086, 286.519, 286.952, 287.385, 287.818, 288.251, 288.684, 289.118, 289.551, 289.985, 290.418, 290.852, 291.286, 291.72, 292.154, 292.588, 293.022, 293.457, 293.891, 294.326, 294.76, 295.195, 295.63, 296.065, 296.5, 296.935, 297.37, 297.805, 298.241, 298.676, 299.112, 299.548, 299.983, 300.419, 300.855, 301.291, 301.728, 302.164, 302.6, 303.036, 303.473, 303.91, 304.346, 304.783, 305.22, 305.657, 306.094, 306.531, 306.968, 307.406, 307.843, 308.28, 308.718, 309.156, 309.594, 310.031, 310.469, 310.907, 311.346, 311.784, 312.222, 312.661, 313.099, 313.538, 313.976, 314.415, 314.854, 315.293, 315.732, 316.171, 316.61, 317.05, 317.489, 317.929, 318.368, 328.368, 338.368, 348.368, 358.368, 368.368, 378.368, 388.368, 398.368, 408.368, 418.368, 428.368, 438.368, 448.368, 458.368, 468.368, 478.368, 488.368, 498.368, 508.368, 518.3679999999999, 528.3679999999999, 538.3679999999999, 548.3679999999999, 558.3679999999999, 568.3679999999999, 578.3679999999999, 588.3679999999999, 598.3679999999999, 608.3679999999999, 618.3679999999999, 628.3679999999999, 638.3679999999999, 648.3679999999999, 658.3679999999999, 668.3679999999999, 678.3679999999999, 688.3679999999999, 698.3679999999999, 708.3679999999999, 718.3679999999999, 728.3679999999999, 738.3679999999999, 748.3679999999999, 758.3679999999999, 768.3679999999999, 778.3679999999999, 788.3679999999999, 798.3679999999999, 808.3679999999999, 818.3679999999999, 828.3679999999999, 838.3679999999999, 848.3679999999999, 858.3679999999999, 868.3679999999999, 878.3679999999999, 888.3679999999999, 898.3679999999999, 908.3679999999999, 918.3679999999999, 928.3679999999999, 938.3679999999999, 948.3679999999999, 958.3679999999999, 968.3679999999999, 978.3679999999999, 988.3679999999999, 998.3679999999999, 1008.3679999999999, 1018.3679999999999, 1028.368, 1038.368, 1048.368, 1058.368, 1068.368, 1078.368, 1088.368, 1098.368, 1108.368};
  return proton_length_bins;
}
  std::vector<double> get_proton_energy_bins(){ 
    std::vector<double> proton_energy_bins = {0.001, 1.001, 2.001, 3.001, 4.001, 5.001, 6.001, 7.001, 8.001, 9.001, 10.001, 11.001, 12.001, 13.001, 14.001, 15.001, 16.001, 17.001, 18.001, 19.001, 20.001, 21.001, 22.001, 23.001, 24.001, 25.001, 26.001, 27.001, 28.001, 29.001, 30.001, 31.001, 32.001, 33.001, 34.001, 35.001, 36.001, 37.001, 38.001, 39.001, 40.001, 41.001, 42.001, 43.001, 44.001, 45.001, 46.001, 47.001, 48.001, 49.001, 50.001, 51.001, 52.001, 53.001, 54.001, 55.001, 56.001, 57.001, 58.001, 59.001, 60.001, 61.001, 62.001, 63.001, 64.001, 65.001, 66.001, 67.001, 68.001, 69.001, 70.001, 71.001, 72.001, 73.001, 74.001, 75.001, 76.001, 77.001, 78.001, 79.001, 80.001, 81.001, 82.001, 83.001, 84.001, 85.001, 86.001, 87.001, 88.001, 89.001, 90.001, 91.001, 92.001, 93.001, 94.001, 95.001, 96.001, 97.001, 98.001, 99.001, 100.001, 101.001, 102.001, 103.001, 104.001, 105.001, 106.001, 107.001, 108.001, 109.001, 110.001, 111.001, 112.001, 113.001, 114.001, 115.001, 116.001, 117.001, 118.001, 119.001, 120.001, 121.001, 122.001, 123.001, 124.001, 125.001, 126.001, 127.001, 128.001, 129.001, 130.001, 131.001, 132.001, 133.001, 134.001, 135.001, 136.001, 137.001, 138.001, 139.001, 140.001, 141.001, 142.001, 143.001, 144.001, 145.001, 146.001, 147.001, 148.001, 149.001, 150.001, 151.001, 152.001, 153.001, 154.001, 155.001, 156.001, 157.001, 158.001, 159.001, 160.001, 161.001, 162.001, 163.001, 164.001, 165.001, 166.001, 167.001, 168.001, 169.001, 170.001, 171.001, 172.001, 173.001, 174.001, 175.001, 176.001, 177.001, 178.001, 179.001, 180.001, 181.001, 182.001, 183.001, 184.001, 185.001, 186.001, 187.001, 188.001, 189.001, 190.001, 191.001, 192.001, 193.001, 194.001, 195.001, 196.001, 197.001, 198.001, 199.001, 200.001, 201.001, 202.001, 203.001, 204.001, 205.001, 206.001, 207.001, 208.001, 209.001, 210.001, 211.001, 212.001, 213.001, 214.001, 215.001, 216.001, 217.001, 218.001, 219.001, 220.001, 221.001, 222.001, 223.001, 224.001, 225.001, 226.001, 227.001, 228.001, 229.001, 230.001, 231.001, 232.001, 233.001, 234.001, 235.001, 236.001, 237.001, 238.001, 239.001, 240.001, 241.001, 242.001, 243.001, 244.001, 245.001, 246.001, 247.001, 248.001, 249.001, 250.001, 251.001, 252.001, 253.001, 254.001, 255.001, 256.001, 257.001, 258.001, 259.001, 260.001, 261.001, 262.001, 263.001, 264.001, 265.001, 266.001, 267.001, 268.001, 269.001, 270.001, 271.001, 272.001, 273.001, 274.001, 275.001, 276.001, 277.001, 278.001, 279.001, 280.001, 281.001, 282.001, 283.001, 284.001, 285.001, 286.001, 287.001, 288.001, 289.001, 290.001, 291.001, 292.001, 293.001, 294.001, 295.001, 296.001, 297.001, 298.001, 299.001, 300.001, 301.001, 302.001, 303.001, 304.001, 305.001, 306.001, 307.001, 308.001, 309.001, 310.001, 311.001, 312.001, 313.001, 314.001, 315.001, 316.001, 317.001, 318.001, 319.001, 320.001, 321.001, 322.001, 323.001, 324.001, 325.001, 326.001, 327.001, 328.001, 329.001, 330.001, 331.001, 332.001, 333.001, 334.001, 335.001, 336.001, 337.001, 338.001, 339.001, 340.001, 341.001, 342.001, 343.001, 344.001, 345.001, 346.001, 347.001, 348.001, 349.001, 350.001, 351.001, 352.001, 353.001, 354.001, 355.001, 356.001, 357.001, 358.001, 359.001, 360.001, 361.001, 362.001, 363.001, 364.001, 365.001, 366.001, 367.001, 368.001, 369.001, 370.001, 371.001, 372.001, 373.001, 374.001, 375.001, 376.001, 377.001, 378.001, 379.001, 380.001, 381.001, 382.001, 383.001, 384.001, 385.001, 386.001, 387.001, 388.001, 389.001, 390.001, 391.001, 392.001, 393.001, 394.001, 395.001, 396.001, 397.001, 398.001, 399.001, 400.001, 401.001, 402.001, 403.001, 404.001, 405.001, 406.001, 407.001, 408.001, 409.001, 410.001, 411.001, 412.001, 413.001, 414.001, 415.001, 416.001, 417.001, 418.001, 419.001, 420.001, 421.001, 422.001, 423.001, 424.001, 425.001, 426.001, 427.001, 428.001, 429.001, 430.001, 431.001, 432.001, 433.001, 434.001, 435.001, 436.001, 437.001, 438.001, 439.001, 440.001, 441.001, 442.001, 443.001, 444.001, 445.001, 446.001, 447.001, 448.001, 449.001, 450.001, 451.001, 452.001, 453.001, 454.001, 455.001, 456.001, 457.001, 458.001, 459.001, 460.001, 461.001, 462.001, 463.001, 464.001, 465.001, 466.001, 467.001, 468.001, 469.001, 470.001, 471.001, 472.001, 473.001, 474.001, 475.001, 476.001, 477.001, 478.001, 479.001, 480.001, 481.001, 482.001, 483.001, 484.001, 485.001, 486.001, 487.001, 488.001, 489.001, 490.001, 491.001, 492.001, 493.001, 494.001, 495.001, 496.001, 497.001, 498.001, 499.001, 500.001, 501.001, 502.001, 503.001, 504.001, 505.001, 506.001, 507.001, 508.001, 509.001, 510.001, 511.001, 512.001, 513.001, 514.001, 515.001, 516.001, 517.001, 518.001, 519.001, 520.001, 521.001, 522.001, 523.001, 524.001, 525.001, 526.001, 527.001, 528.001, 529.001, 530.001, 531.001, 532.001, 533.001, 534.001, 535.001, 536.001, 537.001, 538.001, 539.001, 540.001, 541.001, 542.001, 543.001, 544.001, 545.001, 546.001, 547.001, 548.001, 549.001, 550.001, 551.001, 552.001, 553.001, 554.001, 555.001, 556.001, 557.001, 558.001, 559.001, 560.001, 561.001, 562.001, 563.001, 564.001, 565.001, 566.001, 567.001, 568.001, 569.001, 570.001, 571.001, 572.001, 573.001, 574.001, 575.001, 576.001, 577.001, 578.001, 579.001, 580.001, 581.001, 582.001, 583.001, 584.001, 585.001, 586.001, 587.001, 588.001, 589.001, 590.001, 591.001, 592.001, 593.001, 594.001, 595.001, 596.001, 597.001, 598.001, 599.001, 600.001, 601.001, 602.001, 603.001, 604.001, 605.001, 606.001, 607.001, 608.001, 609.001, 610.001, 611.001, 612.001, 613.001, 614.001, 615.001, 616.001, 617.001, 618.001, 619.001, 620.001, 621.001, 622.001, 623.001, 624.001, 625.001, 626.001, 627.001, 628.001, 629.001, 630.001, 631.001, 632.001, 633.001, 634.001, 635.001, 636.001, 637.001, 638.001, 639.001, 640.001, 641.001, 642.001, 643.001, 644.001, 645.001, 646.001, 647.001, 648.001, 649.001, 650.001, 651.001, 652.001, 653.001, 654.001, 655.001, 656.001, 657.001, 658.001, 659.001, 660.001, 661.001, 662.001, 663.001, 664.001, 665.001, 666.001, 667.001, 668.001, 669.001, 670.001, 671.001, 672.001, 673.001, 674.001, 675.001, 676.001, 677.001, 678.001, 679.001, 680.001, 681.001, 682.001, 683.001, 684.001, 685.001, 686.001, 687.001, 688.001, 689.001, 690.001, 691.001, 692.001, 693.001, 694.001, 695.001, 696.001, 697.001, 698.001, 699.001, 700.001, 701.001, 702.001, 703.001, 704.001, 705.001, 706.001, 707.001, 708.001, 709.001, 710.001, 711.001, 712.001, 713.001, 714.001, 715.001, 716.001, 717.001, 718.001, 719.001, 720.001, 721.001, 722.001, 723.001, 724.001, 725.001, 726.001, 727.001, 728.001, 729.001, 730.001, 731.001, 732.001, 733.001, 734.001, 735.001, 736.001, 737.001, 738.001, 739.001, 740.001, 741.001, 742.001, 743.001, 744.001, 745.001, 746.001, 747.001, 748.001, 749.001, 750.001, 751.001, 752.001, 753.001, 754.001, 755.001, 756.001, 757.001, 758.001, 759.001, 760.001, 761.001, 762.001, 763.001, 764.001, 765.001, 766.001, 767.001, 768.001, 769.001, 770.001, 771.001, 772.001, 773.001, 774.001, 775.001, 776.001, 777.001, 778.001, 779.001, 780.001, 781.001, 782.001, 783.001, 784.001, 785.001, 786.001, 787.001, 788.001, 789.001, 790.001, 791.001, 792.001, 793.001, 794.001, 795.001, 796.001, 797.001, 798.001, 799.001, 800.001, 801.001, 802.001, 803.001, 804.001, 805.001, 806.001, 807.001, 808.001, 809.001, 810.001, 811.001, 812.001, 813.001, 814.001, 815.001, 816.001, 817.001, 818.001, 819.001, 820.001, 821.001, 822.001, 823.001, 824.001, 825.001, 826.001, 827.001, 828.001, 829.001, 830.001, 831.001, 832.001, 833.001, 834.001, 835.001, 836.001, 837.001, 838.001, 839.001, 840.001, 841.001, 842.001, 843.001, 844.001, 845.001, 846.001, 847.001, 848.001, 849.001, 850.001, 851.001, 852.001, 853.001, 854.001, 855.001, 856.001, 857.001, 858.001, 859.001, 860.001, 861.001, 862.001, 863.001, 864.001, 865.001, 866.001, 867.001, 868.001, 869.001, 870.001, 871.001, 872.001, 873.001, 874.001, 875.001, 876.001, 877.001, 878.001, 879.001, 880.001, 881.001, 882.001, 883.001, 884.001, 885.001, 886.001, 887.001, 888.001, 889.001, 890.001, 891.001, 892.001, 893.001, 894.001, 895.001, 896.001, 897.001, 898.001, 899.001, 900.001, 901.001, 902.001, 903.001, 904.001, 905.001, 906.001, 907.001, 908.001, 909.001, 910.001, 911.001, 912.001, 913.001, 914.001, 915.001, 916.001, 917.001, 918.001, 919.001, 920.001, 921.001, 922.001, 923.001, 924.001, 925.001, 926.001, 927.001, 928.001, 929.001, 930.001, 931.001, 932.001, 933.001, 934.001, 935.001, 936.001, 937.001, 938.001, 939.001, 940.001, 941.001, 942.001, 943.001, 944.001, 945.001, 946.001, 947.001, 948.001, 949.001, 950.001, 951.001, 952.001, 953.001, 954.001, 955.001, 956.001, 957.001, 958.001, 959.001, 960.001, 961.001, 962.001, 963.001, 964.001, 965.001, 966.001, 967.001, 968.001, 969.001, 970.001, 971.001, 972.001, 973.001, 974.001, 975.001, 976.001, 977.001, 978.001, 979.001, 980.001, 981.001, 982.001, 983.001, 984.001, 985.001, 986.001, 987.001, 988.001, 989.001, 990.001, 991.001, 992.001, 993.001, 994.001, 995.001, 996.001, 997.001, 998.001, 999.001, 1020.001, 1041.001, 1062.001, 1083.001, 1104.001, 1125.001, 1146.001, 1167.001, 1188.001, 1209.001, 1230.001, 1251.001, 1272.001, 1293.001, 1314.001, 1335.001, 1356.001, 1377.001, 1398.001, 1419.001, 1440.001, 1461.001, 1482.001, 1503.001, 1524.001, 1545.001, 1566.001, 1587.001, 1608.001, 1629.001, 1650.001, 1671.001, 1692.001, 1713.001, 1734.001, 1755.001, 1776.001, 1797.001, 1818.001, 1839.001, 1860.001, 1881.001, 1902.001, 1923.001, 1944.001, 1965.001, 1986.001, 2007.001, 2028.001, 2049.001, 2070.001, 2091.001, 2112.001, 2133.001, 2154.001, 2175.001, 2196.001, 2217.001, 2238.001, 2259.001, 2280.001, 2301.001, 2322.001, 2343.001, 2364.001, 2385.001, 2406.001, 2427.001, 2448.001, 2469.001, 2490.001, 2511.001, 2532.001, 2553.001, 2574.001, 2595.001, 2616.001, 2637.001, 2658.001};
    return proton_energy_bins;
  }


}


double LEEana::get_reco_Enu_corr(KineInfo& kine, bool flag_data){
  double reco_Enu_corr = 0;
  if (kine.kine_reco_Enu > 0){
    if (flag_data){
      for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
  	if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
  	  reco_Enu_corr +=  kine.kine_energy_particle->at(j) * em_charge_scale;
  	}else{
  	  reco_Enu_corr +=  kine.kine_energy_particle->at(j);
  	}
  	//	std::cout << "p: " << kine.kine_energy_particle->at(j) << " " << kine.kine_energy_info->at(j) << " " << kine.kine_particle_type->at(j) << " " << kine.kine_energy_included->at(j) << std::endl;
      }
      reco_Enu_corr += kine.kine_reco_add_energy;
      return reco_Enu_corr;
    }
  }
  return kine.kine_reco_Enu;
}

double LEEana::get_reco_showerKE_corr(PFevalInfo& pfeval, bool flag_data){
    if (flag_data){
        return pfeval.reco_showerKE * em_charge_scale;
    } else {
        return pfeval.reco_showerKE;
    }
}


double LEEana::get_reco_Eproton(KineInfo& kine){
  double reco_Eproton=0;
  for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
    {
      int pdgcode = kine.kine_particle_type->at(i);
      if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) // proton threshold of 35 MeV
      //if(abs(pdgcode)==2212) // no proton threshold
        reco_Eproton+=kine.kine_energy_particle->at(i);

    }
  return reco_Eproton;
}

double LEEana::get_reco_Epion(KineInfo& kine){
  double reco_Epion=0;
  for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
    {
      int pdgcode = kine.kine_particle_type->at(i);
      //if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10)  // KE threshold: 10 keV
      if(abs(pdgcode)==211)  // no threshold
        reco_Epion+=kine.kine_energy_particle->at(i);
    }
  return reco_Epion;
}


bool LEEana::is_true_0p(PFevalInfo& pfeval){
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue; //Erin: CHANGE, no proton threshold
      return false;
    }
  return true;
}

std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> LEEana::get_range_proton_KE(PFevalInfo& pfeval, SpaceInfo& space, bool return_MeV){

  std::vector<double> prim_proton_KEs;
  std::vector<double> proton_KEs;
  std::vector<double> prim_larpid_proton_KEs;
  std::vector<double> larpid_proton_KEs;

  std::vector<double> proton_length_bins = get_proton_length_bins(); 
  std::vector<double> proton_energy_bins = get_proton_energy_bins();

  int n_spacepoints = space.Trecchargeblob_spacepoints_real_cluster_id->size();

  for(size_t i=0; i<pfeval.reco_Ntrack; i++){

    bool flag_wc=false;
    bool flag_larpid=false;
    if (pfeval.reco_pdg[i]!=2212) flag_wc=true;
    if (pfeval.reco_larpid_pdg[i]!=2212) flag_larpid=true;
    if(!flag_wc && !flag_larpid) continue;

    double range_proton = 0; 
    double proton_KE = 0;
    std::vector<float> *temp_spacepoints_x = new std::vector<float>;
    std::vector<float> *temp_spacepoints_y = new std::vector<float>;
    std::vector<float> *temp_spacepoints_z = new std::vector<float>;
    std::vector<float> *temp_spacepoints_q = new std::vector<float>;
    for(size_t sp=0; sp<n_spacepoints; sp++){
      if(space.Trecchargeblob_spacepoints_real_cluster_id->at(sp)==pfeval.reco_id[i]){
        temp_spacepoints_x->push_back(space.Trecchargeblob_spacepoints_x->at(sp));
        temp_spacepoints_y->push_back(space.Trecchargeblob_spacepoints_y->at(sp));
        temp_spacepoints_z->push_back(space.Trecchargeblob_spacepoints_z->at(sp));
        temp_spacepoints_q->push_back(space.Trecchargeblob_spacepoints_q->at(sp));
      }
    }
    int n_spacepoints_part = temp_spacepoints_x->size();
    if(n_spacepoints_part==0) continue;
      
    for(int sp=0; sp<n_spacepoints_part-1; sp++){
      double dx = temp_spacepoints_x->at(sp) - temp_spacepoints_x->at(sp+1);
      double dy = temp_spacepoints_y->at(sp) - temp_spacepoints_y->at(sp+1);
      double dz = temp_spacepoints_z->at(sp) - temp_spacepoints_z->at(sp+1);
      double dist = sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
      range_proton+=dist;
    }

    int nbb = proton_length_bins.size();
    for(int bb=0; bb<nbb-1; bb++){
      if(proton_length_bins.at(bb)<=range_proton && proton_length_bins.at(bb+1)>range_proton){
        double m = (proton_energy_bins.at(bb+1)-proton_energy_bins.at(bb))/(proton_length_bins.at(bb+1)-proton_length_bins.at(bb));
        proton_KE = proton_energy_bins.at(bb) + m * (range_proton-proton_length_bins.at(bb));
        break;
      }
    }
    
    if(flag_wc) proton_KEs.push_back(proton_KE);
    if(flag_larpid) larpid_proton_KEs.push_back(proton_KE);
    if(pfeval.reco_mother[i]==0){
      if(flag_wc) prim_proton_KEs.push_back(proton_KE);
      if(flag_larpid) prim_larpid_proton_KEs.push_back(proton_KE);
    }
  }

  if(proton_KEs.size()==0) proton_KEs.push_back(0);
  if(larpid_proton_KEs.size()==0) larpid_proton_KEs.push_back(0);
  if(prim_proton_KEs.size()==0) prim_proton_KEs.push_back(0);
  if(prim_larpid_proton_KEs.size()==0) prim_larpid_proton_KEs.push_back(0);
       
  std::sort(prim_proton_KEs.begin(), prim_proton_KEs.end(), wayToSort);
  std::sort(proton_KEs.begin(), proton_KEs.end(), wayToSort);
  std::sort(prim_larpid_proton_KEs.begin(), prim_larpid_proton_KEs.end(), wayToSort);
  std::sort(larpid_proton_KEs.begin(), larpid_proton_KEs.end(), wayToSort);
  std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = std::make_tuple(prim_proton_KEs,proton_KEs,prim_larpid_proton_KEs,larpid_proton_KEs);
  return result;

}


std::vector<double> LEEana::get_pandora_proton_KE(PandoraInfo& pandora, double TRACK_SCORE_CUT, bool return_MeV){

  std::vector<double> proton_KEs;
    
  if (pandora.slice_orig_pass_id != 1){
    proton_KEs.push_back(0);
    return proton_KEs;
  }

  for(size_t part=0; part<pandora.n_pfps; part++){
        
    int generation = pandora.pfp_generation_v->at(part);
    if(generation!=2) continue;

    if(pandora.pfpdg->at(part)!=13) continue;

    if(pandora.trk_llr_pid_score_v->at(part) < TRACK_SCORE_CUT){ proton_KEs.push_back(pandora.trk_energy_proton_v->at(part)); }
  }

  if(return_MeV){
    for(size_t part=0; part<proton_KEs.size(); part++){
      proton_KEs.at(part) = proton_KEs.at(part)*1000;
    }
  }
               
  if(proton_KEs.size()==0) proton_KEs.push_back(0); -  std::sort(proton_KEs.begin(), proton_KEs.end(), wayToSort);
  std::sort(proton_KEs.begin(), proton_KEs.end(), wayToSort);
  return proton_KEs;

}

std::vector<double> LEEana::get_lantern_KE(LanternInfo& lantern, int pdg, double vtx_cut, bool return_MeV){

  std::vector<double> p_KEs;

  for(size_t part=0; part<lantern.nTracks; part++){

    if(lantern.trackPID[part]!=pdg) continue;
    if(lantern.trackIsSecondary[part]!=0) continue;
    if(lantern.trackDistToVtx[part]>vtx_cut) continue;
    p_KEs.push_back(lantern.trackRecoE[part]);
  }

  if(!return_MeV){
    for(size_t part=0; part<p_KEs.size(); part++){
      p_KEs.at(part) = p_KEs.at(part)/1000;
    }
  }

  if(p_KEs.size()==0) p_KEs.push_back(0);
  std::sort(p_KEs.begin(), p_KEs.end(), wayToSort);
  return p_KEs;

}


double LEEana::get_weight(TString weight_name, EvalInfo& eval, PFevalInfo& pfeval, KineInfo& kine, TaggerInfo& tagger, std::tuple< bool, std::vector< std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double>  > > > rw_info, std::map<int, std::tuple< double, double, double, double > > time_info, bool flag_data){
  double addtl_weight = 1.0;

  // CV correction from numuCC cross section data
  // if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1){
  //   if (eval.truth_nuEnergy>200 && eval.truth_nuEnergy<=540) addtl_weight = 1.28043;
  //   else if (eval.truth_nuEnergy>540 && eval.truth_nuEnergy<=705) addtl_weight = 1.21158;
  //   else if (eval.truth_nuEnergy>705 && eval.truth_nuEnergy<=805) addtl_weight = 1.19091;
  //   else if (eval.truth_nuEnergy>805 && eval.truth_nuEnergy<=920) addtl_weight = 1.17733;
  //   else if (eval.truth_nuEnergy>920 && eval.truth_nuEnergy<=1050) addtl_weight = 1.13983;
  //   else if (eval.truth_nuEnergy>1050 && eval.truth_nuEnergy<=1200) addtl_weight = 1.07864;
  //   else if (eval.truth_nuEnergy>1200 && eval.truth_nuEnergy<=1375) addtl_weight = 1.00722;
  //   else if (eval.truth_nuEnergy>1375 && eval.truth_nuEnergy<=1570) addtl_weight = 0.93857;
  //   else if (eval.truth_nuEnergy>1570 && eval.truth_nuEnergy<=2050) addtl_weight = 0.886241;
  //   else if (eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) addtl_weight = 0.858724;
  //   else if (eval.truth_nuEnergy>4000) addtl_weight = 0.858724;
  // }
  // std::cout << "energy: " << eval.truth_nuEnergy << " addtl_weight: " << addtl_weight << std::endl;
  // end of data correction

  //Begin reweighting
  std::tuple<bool, TString, TString, double, double, bool, bool, bool, std::vector<double>, std::vector<double> > rw_info_i;

  if(std::get<0>(rw_info) && !(flag_data)){//Are you applying any reweighting?
    for(size_t rw=0; rw<std::get<1>(rw_info).size(); rw++){
      rw_info_i = std::get<1>(rw_info)[rw];
      if(std::get<0>(rw_info_i)){//Are you reweighting this channel and cut?
        TString cut_str = std::get<1>(rw_info_i);
        TString var_str = std::get<2>(rw_info_i);
        double var = get_truth_var(kine, eval, pfeval, tagger, var_str);
        double min_var = std::get<3>(rw_info_i);
        double max_var = std::get<4>(rw_info_i);
        bool underflow = std::get<5>(rw_info_i);
        bool overflow = std::get<6>(rw_info_i);
        bool equal_binning = std::get<7>(rw_info_i);
        std::vector<double> reweight = std::get<8>(rw_info_i);

        int wbin;
        bool flag_pass = get_rw_cut_pass(cut_str, eval, pfeval, tagger, kine);
        if (flag_pass){
          if (var>max_var && overflow) addtl_weight = reweight.back();
          else if(var>max_var) addtl_weight = 1;
          else if (var<min_var && underflow) addtl_weight = reweight[0];
          else if (var>min_var){
            if(equal_binning){
              double bin_len = (max_var-min_var)/reweight.size();
              if(underflow && overflow) bin_len = (max_var-min_var)/(reweight.size()-2);
              else if (underflow || overflow) bin_len = (max_var-min_var)/(reweight.size()-1);
              wbin = floor((var-min_var)/bin_len);
            }else{
              std::vector<double> bins = std::get<9>(rw_info_i);
              for(int b=0; b<bins.size()-1; b++){
                if(var<=bins[b+1] && var>bins[b]){
                  wbin = b;
                  break;
                }
              }
            }
            if(underflow) wbin++;
            addtl_weight *= reweight[wbin];
          }
        }
      }
    }
  }

  if (weight_name == "cv_spline"){
    return addtl_weight*eval.weight_cv * eval.weight_spline;
  //Erin - ns beam time scaling
  }else if (weight_name == "cv_spline_nsbeam"){
    float beam_scale = 0.86;

    bool has_muon = false; //set to true to turn off
    if (pfeval.reco_muonMomentum[3] > 0){has_muon = true;}

    if (pfeval.run >= 13697){ 
      beam_scale=std::get<0>(time_info[3]);//0.845;}//0.931503;}//0.913671 - 0.0812331; }
      if (!has_muon){ beam_scale = 0.82;}
    }
    else if (pfeval.run >= 8321){ 
      beam_scale=std::get<0>( time_info[2]);
      if (!has_muon){ beam_scale = 0.86;}
    }//0.88;}//0.919618;}//0.900644 - 0.044328;}
    else if (pfeval.run > 0 ){ 
      beam_scale=std::get<0>( time_info[1]);//0.912832;}//0.885887 - 0.0315298;}
      if (!has_muon){ beam_scale = 0.78;}
    }
    //beam_scale = beam_scale - 0.03;
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}
    //0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}
    //0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}
    //0.527815;}//0.471911 + 0.0315298;}
    if(eval.match_completeness_energy<=0.1*eval.truth_energyInside){beam_scale = 1.0-ext_rej;}
    return addtl_weight*eval.weight_cv * eval.weight_spline * beam_scale;
  }else if (weight_name == "dirt_nsbeam"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    return addtl_weight*eval.weight_cv * eval.weight_spline * (1.0-ext_rej);
  }else if (weight_name == "nsbeam_ext"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}
    //0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}
    //0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}
    //0.527815;}//0.471911 + 0.0315298;}
    //ext_rej = ext_rej + 0.03;
    float ext_scale = 1.0 - ext_rej;
    return ext_scale;
  }else if (weight_name == "cv_spline_nsbeam_cv_spline_nsbeam"){
    float beam_scale = 0.86;

    bool has_muon = false; //set to true to turn off
    if (pfeval.reco_muonMomentum[3] > 0){has_muon = true;}

    if (pfeval.run >= 13697){ 
      beam_scale=std::get<0>(time_info[3]);//0.845;}//0.931503;}//0.913671 - 0.0812331; }
      if (!has_muon){ beam_scale = 0.82;}
    }
    else if (pfeval.run >= 8321){ 
      beam_scale=std::get<0>( time_info[2]);//0.88;}//0.919618;}//0.900644 - 0.044328;}
      if (!has_muon){ beam_scale = 0.86;}
    }
    else if (pfeval.run > 0 ){ 
      beam_scale=std::get<0>( time_info[1]);//0.912832;}//0.885887 - 0.0315298;}
      if (!has_muon){ beam_scale = 0.78;}
    }
    //beam_scale = beam_scale - 0.03;
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    if(eval.match_completeness_energy<=0.1*eval.truth_energyInside){beam_scale = 1.0-ext_rej;}
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * beam_scale,2);
  }else if (weight_name == "dirt_nsbeam_dirt_nsbeam"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * (1.0-ext_rej),2);
  }else if (weight_name == "nsbeam_ext_nsbeam_ext"){
    float ext_rej = 0.47;
    if (pfeval.run >= 13697){ ext_rej = std::get<2>( time_info[3]);}//0.68;}//0.535783;}//0.471911 + 0.0812331; }
    else if (pfeval.run >= 8321){ ext_rej = std::get<2>( time_info[2]);}//0.66;}//0.532919;}//0.471911 + 0.044328;}
    else if (pfeval.run > 0 ){ ext_rej = std::get<2>( time_info[1]);}//0.527815;}//0.471911 + 0.0315298;}
    //ext_rej = ext_rej + 0.03;
    float ext_scale = 1.0 - ext_rej;
    return pow(ext_scale,2);
//cex bug fix weights
  }else if (weight_name == "cv_spline_cexbugfix"){
    double ratio_weight = 1.0;
    if (eval.truth_isCC == 0 && pfeval.truth_NprimPio==1)
    {
      //get the pi0 costheta and KE
      double truth_pi0_costheta = -1000.;
      double truth_pi0_KE = -1000.;
      int true_num_protons_35_MeV = 0;
      for(int jth=0; jth<pfeval.truth_Ntrack; jth++){
        int mother = pfeval.truth_mother[jth];
        if(mother != 0) continue;
        int pdgcode = pfeval.truth_pdg[jth];
        if(abs(pdgcode)==111){
          //N_th_pi0++;
          double px = pfeval.truth_startMomentum[jth][0]*1000.; // MeV
          double py = pfeval.truth_startMomentum[jth][1]*1000.; // MeV
          double pz = pfeval.truth_startMomentum[jth][2]*1000.; // MeV
          truth_pi0_costheta = pz / sqrt(px*px + py*py + pz*pz);
          truth_pi0_KE = pfeval.truth_startMomentum[jth][3]*1000. - 134.9768;
        }
        if (pdgcode==2212 && pfeval.truth_startMomentum[jth][3]*1000. - 938.272089 > 35.){
          true_num_protons_35_MeV++;
        }
      }
      //pick 0p or Np csv file
      std::string ratiofilename;
      if (true_num_protons_35_MeV>0){
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_Np1pi0.csv";
      }else{
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_0p1pi0.csv";
      }
      std::ifstream ratiofile(ratiofilename);
      if (!ratiofile.is_open()) {
          throw std::runtime_error("Could not open file");
      }
      std::string ratioline;
      std::getline(ratiofile, ratioline); // Skip header line
      while (std::getline(ratiofile, ratioline)) {
          std::istringstream ss(ratioline);
          char comma;
          double pi0_KE;
          double pi0_cos;
          double ratio;
          ss >> pi0_KE >> comma >> pi0_cos >> comma >> ratio;
          double pi0_cos_low = pi0_cos - 0.02;
          double pi0_cos_high = pi0_cos + 0.02;
          double pi0_KE_low = pi0_KE - 5.0;
          double pi0_KE_high = pi0_KE + 5.0;
          if (pi0_cos_low <= truth_pi0_costheta && truth_pi0_costheta <= pi0_cos_high && pi0_KE_low <= truth_pi0_KE && truth_pi0_KE <= pi0_KE_high){
            ratio_weight = ratio;
            break;
          }  
      }
    }
    if (ratio_weight > 10.0) ratio_weight = 10.0;
    if (ratio_weight < 0.0) ratio_weight = 0.0;
    return addtl_weight*eval.weight_cv * eval.weight_spline * ratio_weight;
  }else if (weight_name == "cv_spline_cexbugfix_cv_spline_cexbugfix"){
    double ratio_weight = 1.0;
    if (eval.truth_isCC == 0 && pfeval.truth_NprimPio==1)
    {
      //get the pi0 costheta and KE
      double truth_pi0_costheta = -1000.;
      double truth_pi0_KE = -1000.;
      int true_num_protons_35_MeV = 0;
      for(int jth=0; jth<pfeval.truth_Ntrack; jth++){
        int mother = pfeval.truth_mother[jth];
        if(mother != 0) continue;
        int pdgcode = pfeval.truth_pdg[jth];
        if(abs(pdgcode)==111){
          //N_th_pi0++;
          double px = pfeval.truth_startMomentum[jth][0]*1000.; // MeV
          double py = pfeval.truth_startMomentum[jth][1]*1000.; // MeV
          double pz = pfeval.truth_startMomentum[jth][2]*1000.; // MeV
          truth_pi0_costheta = pz / sqrt(px*px + py*py + pz*pz);
          truth_pi0_KE = pfeval.truth_startMomentum[jth][3]*1000. - 134.9768;
        }
        if (pdgcode==2212 && pfeval.truth_startMomentum[jth][3]*1000. - 938.272089 > 35.){
          true_num_protons_35_MeV++;
        }
      }
      //pick 0p or Np csv file
      std::string ratiofilename;
      if (true_num_protons_35_MeV>0){
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_Np1pi0.csv";
      }else{
        ratiofilename = "/exp/uboone/data/users/mismail/pi0-fsi/get_ratios/2d_ratio_0p1pi0.csv";
      }
      std::ifstream ratiofile(ratiofilename);
      if (!ratiofile.is_open()) {
          throw std::runtime_error("Could not open file");
      }
      std::string ratioline;
      std::getline(ratiofile, ratioline); // Skip header line
      while (std::getline(ratiofile, ratioline)) {
          std::istringstream ss(ratioline);
          char comma;
          double pi0_KE;
          double pi0_cos;
          double ratio;
          ss >> pi0_KE >> comma >> pi0_cos >> comma >> ratio;
          double pi0_cos_low = pi0_cos - 0.02;
          double pi0_cos_high = pi0_cos + 0.02;
          double pi0_KE_low = pi0_KE - 5.0;
          double pi0_KE_high = pi0_KE + 5.0;
          if (pi0_cos_low <= truth_pi0_costheta && truth_pi0_costheta <= pi0_cos_high && pi0_KE_low <= truth_pi0_KE && truth_pi0_KE <= pi0_KE_high){
            ratio_weight = ratio;
            break;
          }
      }
    }
    if (ratio_weight > 10.0) ratio_weight = 10.0;
    if (ratio_weight < 0.0) ratio_weight = 0.0;
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline * ratio_weight,2);
  }else if (weight_name == "cv_spline_cv_spline"){
    return pow(addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "unity" || weight_name == "unity_unity"){
    return 1;
  }else if (weight_name == "lee_cv_spline"){
    if (eval.weight_lee <= 0){
      eval.weight_lee = 1.0;
    }
    return (eval.weight_lee * addtl_weight*eval.weight_cv * eval.weight_spline);
  }else if (weight_name == "lee_cv_spline_lee_cv_spline"){
    if (eval.weight_lee <= 0){
      eval.weight_lee = 1.0;
    }
    return pow(eval.weight_lee * addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "lee_cv_spline_cv_spline" || weight_name == "cv_spline_lee_cv_spline"){
    if (eval.weight_lee <= 0){
      eval.weight_lee = 1.0;
    }
    return eval.weight_lee * pow(addtl_weight*eval.weight_cv * eval.weight_spline,2);
  }else if (weight_name == "spline"){
    return eval.weight_spline;
  }else if (weight_name == "spline_spline"){
    return pow(eval.weight_spline,2);
  }else if (weight_name == "lee_spline"){
    return (eval.weight_lee * eval.weight_spline);
  }else if (weight_name == "lee_spline_lee_spline"){
    return pow(eval.weight_lee * eval.weight_spline,2);
  }else if (weight_name == "lee_spline_spline" || weight_name == "spline_lee_spline"){
    return eval.weight_lee * pow( eval.weight_spline,2);
  }else if (weight_name == "add_weight"){//for systematics
    return addtl_weight;
  }else{
    std::cout <<"Unknown weights: " << weight_name << std::endl;
  }


  return 1;
}

double LEEana::get_truth_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger , TString var_name){
  if(var_name == "truth_energyInside"){
    return eval.truth_energyInside;
  }else {std::cout<<"Unknown truth var, check configurations"<<std::endl;}
  return 0;
}

double LEEana::get_kine_var(KineInfo& kine, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, bool flag_data, TString var_name, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern){
  //  if (var_name == "kine_reco_Enu"){
  //  return kine.kine_reco_Enu;
  //  }else
  if (var_name == "kine_reco_Enu"){
    return get_reco_Enu_corr(kine, flag_data);

  }else if (var_name == "lantern_nTracks"){
    return lantern.nTracks;

  }else if (var_name=="range_leading_prim_proton_KE"){
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = get_range_proton_KE(pfeval, space, true);
    return std::get<0>(result).at(0);
  }else if (var_name=="range_subleading_prim_proton_KE"){
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = get_range_proton_KE(pfeval, space, true);
    if(std::get<0>(result).size()<2) return 0;
    return std::get<0>(result).at(1);
  }else if (var_name=="range_leading_proton_KE"){
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = get_range_proton_KE(pfeval, space, true);
    return std::get<1>(result).at(0);
  }else if (var_name=="range_leading_prim_LArPID_proton_KE"){
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = get_range_proton_KE(pfeval, space, true);
    return std::get<2>(result).at(0);
  }else if (var_name=="range_leading_LArPID_proton_KE"){
    std::tuple<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>> result = get_range_proton_KE(pfeval, space, true);
    return std::get<3>(result).at(0);
  std::vector<double> get_lantern_KE(LanternInfo& lantern, int pdg, double vtx_cut, bool return_MeV);

  }else if (var_name == "pandora_leading_prim_proton_KE"){
    return get_pandora_proton_KE(pandora, 0.05, true).at(0);
  }else if (var_name == "panora_subleading_prim_proton_KE"){
    std::vector<double> result = get_pandora_proton_KE(pandora, 0.05, true);
    if(result.size()<2) return 0;
    return result.at(1);

  }else if (var_name == "lantern_leading_prim_proton_KE"){
    return get_lantern_KE(lantern, 2212, 10, true).at(0);
  }else if (var_name == "lantern_subleading_prim_proton_KE"){
    std::vector<double> result = get_lantern_KE(lantern, 2212, 10, true);
    if(result.size()<2) return 0;
    return result.at(1);

  }else if (var_name == "reco_showerKE"){
    return get_reco_showerKE_corr(pfeval, flag_data) * 1000.;
  }else if (var_name == "kine_reco_Eproton"){
    return get_reco_Eproton(kine);
  }else if (var_name == "kine_reco_Eproton_nothreshold"){
    double reco_Eproton=0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
        int pdgcode = kine.kine_particle_type->at(i);
        //if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) // proton threshold of 35 MeV
        if(abs(pdgcode)==2212) // no proton threshold
          reco_Eproton+=kine.kine_energy_particle->at(i);

      }
    return reco_Eproton;
  }else if (var_name == "kine_reco_Epion"){
    return get_reco_Epion(kine);
  }else if (var_name == "kine_pio_energy_1"){
    return kine.kine_pio_energy_1;
  }else if (var_name == "kine_pio_energy_2"){
    return kine.kine_pio_energy_2;
  }else if (var_name == "kine_pio_energy_max"){
    if(flag_data)
      return std::max(kine.kine_pio_energy_1*em_charge_scale, kine.kine_pio_energy_2*em_charge_scale);
    else
      return std::max(kine.kine_pio_energy_1, kine.kine_pio_energy_2);
  }else if (var_name == "kine_pio_energy_min"){
    if(flag_data)
      return std::min(kine.kine_pio_energy_1*em_charge_scale, kine.kine_pio_energy_2*em_charge_scale);
    else
      return std::min(kine.kine_pio_energy_1, kine.kine_pio_energy_2);
  }else if (var_name == "kine_pio_angle" || var_name == "kine_pio_costheta"){
    if (var_name == "kine_pio_angle")
      return kine.kine_pio_angle;
    else
      return TMath::Cos(kine.kine_pio_angle/180.*TMath::Pi());
  }else if (var_name == "match_energy"){
    return eval.match_energy;
  }else if (var_name == "pi0_energy"){
    double pi0_mass = 135;
    double alpha = fabs(kine.kine_pio_energy_1 - kine.kine_pio_energy_2)/(kine.kine_pio_energy_1 + kine.kine_pio_energy_2);
    return pi0_mass * (sqrt(2./(1-alpha*alpha)/(1-cos(kine.kine_pio_angle/180.*3.1415926)))-1);
  }else if (var_name == "pi0_mass"){
    if (kine.kine_pio_mass >0){
      //      TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      // TLorentzVector pio = p1 + p2;
      if (flag_data) {
	return kine.kine_pio_mass * em_charge_scale;
      }else{
	return kine.kine_pio_mass;
      }
    }else{
      return kine.kine_pio_mass;
    }
    //  }else if (var_name == "pi0_mass"){
    // return kine.kine_pio_mass;
  }else if (var_name == "nue_score"){
    return (tagger.nue_score<=15.99?tagger.nue_score:15.99);
  }else if (var_name == "nc_pio_score"){
    return tagger.nc_pio_score;
  }else if (var_name == "nc_delta_score"){
    return tagger.nc_delta_score;
  }else if (var_name == "numu_score"){
    return tagger.numu_score;
  }else if (var_name == "shower_energy"){
    if(flag_data)
      return tagger.mip_energy*em_charge_scale;
    else
      return tagger.mip_energy;
  }else if (var_name == "electron_energy"){
    if (flag_data){
      return pfeval.reco_showerMomentum[3] * em_charge_scale *1000;
    }else{
      return pfeval.reco_showerMomentum[3] * 1000;
    }
  }else if (var_name == "electron_polar_angle"){
    return pfeval.reco_showerMomentum[2]/pfeval.reco_showerMomentum[3];
  }else if (var_name == "shower_angle_beam"){
    return tagger.mip_angle_beam;
  }else if (var_name == "shower_angle_vertical"){
    return tagger.spt_angle_vertical;
  }else if (var_name == "shwvtx_nuvtx_dis"){
    return sqrt(pow(pfeval.reco_nuvtxX-pfeval.reco_showervtxX,2)+pow(pfeval.reco_nuvtxY-pfeval.reco_showervtxY,2)+pow(pfeval.reco_nuvtxZ-pfeval.reco_showervtxZ,2));
  }else if (var_name == "median_dQdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    return vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
  }else if (var_name == "median_dEdx"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.mip_vec_dQ_dx_2);
    dqdx.push_back(tagger.mip_vec_dQ_dx_3);
    dqdx.push_back(tagger.mip_vec_dQ_dx_4);
    dqdx.push_back(tagger.mip_vec_dQ_dx_5);
    dqdx.push_back(tagger.mip_vec_dQ_dx_6);
    dqdx.push_back(tagger.mip_vec_dQ_dx_7);
    dqdx.push_back(tagger.mip_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    return median_dedx; // MeV/cm
  }else if (var_name == "reco_showervtxX"){
    if(pfeval.reco_showerKE>0.){
      return pfeval.reco_showervtxX;
    }else{
      Float_t x = -9999.;
      Float_t max_KE = 0.02;
      for (Int_t i_p = 0; i_p < pfeval.reco_Ntrack; i_p++){
        if (pfeval.reco_pdg[i_p]==11 || pfeval.reco_pdg[i_p]==22 && pfeval.reco_startMomentum[i_p][3]>max_KE){
          x = pfeval.reco_startXYZT[i_p][0];
        }
      }
      return x;
    }
  }else if (var_name == "reco_showervtxY"){
    if(pfeval.reco_showerKE>0.){
      return pfeval.reco_showervtxY;
    }else{
      Float_t y = -9999.;
      Float_t max_KE = 0.02;
      for (Int_t i_p = 0; i_p < pfeval.reco_Ntrack; i_p++){
        if (pfeval.reco_pdg[i_p]==11 || pfeval.reco_pdg[i_p]==22 && pfeval.reco_startMomentum[i_p][3]>max_KE){
          y = pfeval.reco_startXYZT[i_p][1];
        }
      }
      return y;
    }
  }else if (var_name == "reco_showervtxZ"){
    if(pfeval.reco_showerKE>0.){
      return pfeval.reco_showervtxZ;
    }else{
      Float_t z = -9999.;
      Float_t max_KE = 0.02;
      for (Int_t i_p = 0; i_p < pfeval.reco_Ntrack; i_p++){
        if (pfeval.reco_pdg[i_p]==11 || pfeval.reco_pdg[i_p]==22 && pfeval.reco_startMomentum[i_p][3]>max_KE){
          z = pfeval.reco_startXYZT[i_p][2];
        }
      }
      return z;
    }
  }else if (var_name == "reco_nuvtxX"){
      return pfeval.reco_nuvtxX;
  }else if (var_name == "reco_nuvtxY"){
      return pfeval.reco_nuvtxY;
  }else if (var_name == "reco_nuvtxZ"){
      return pfeval.reco_nuvtxZ;
  }else if (var_name == "reco_nuvtxU"){
    return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) - pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "reco_nuvtxV"){
    return pfeval.reco_nuvtxZ * TMath::Cos(3.1415926/3.) + pfeval.reco_nuvtxY * TMath::Sin(3.1415926/3.);
  }else if (var_name == "mip_quality_n_tracks"){
      return tagger.mip_quality_n_tracks;
  }else if (var_name == "mip_quality_n_showers"){
      return tagger.mip_quality_n_showers;
  }else if (var_name == "gap_n_bad"){
      return tagger.gap_n_bad;
  }else if (var_name == "muon_KE"){
      return pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
  }else if (var_name == "reco_Emuon"){
      return pfeval.reco_muonMomentum[3]*1000; // GeV --> MeV
  /*}else if(var_name == "reco_Emuon_hybrid"){
      if (eval.match_isFC) {
        return 1000.0*kine.vlne_v4_numu_full_primaryE;
      }
      else {
        return 1000.0*pfeval.reco_muonMomentum[3];
      }
  }else if(var_name == "reco_Emuon_dlnew"){
      // std::cout << "vlne_v4_numu_full_primaryE: " << kine.vlne_v4_numu_full_primaryE << std::endl;
      return 1000.0*kine.vlne_v4_numu_full_primaryE;*/
  }else if (var_name == "muon_momentum"){
      if (pfeval.reco_muonMomentum[3] < 0) { return -1; }
      float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66; // GeV --> MeV
      return (TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66));
  }else if (var_name == "muon_costheta"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return TMath::Cos(muonMomentum.Theta());
      else
	return -2;
  }else if (var_name == "muon_theta"){

      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return muonMomentum.Theta()*180./TMath::Pi();
      else
	return -1000;
  }else if (var_name == "muon_phi"){
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
	return muonMomentum.Phi()/TMath::Pi()*180.;
      else
	return -1000;
   }else if (var_name == "reco_Eqe_muon" || var_name == "reco_Eqe_muon_Enu_diff" || var_name == "reco_Eqe_electron" || var_name == "reco_Eqe_electron_Enu_diff"){
      // everything is in MeV
      float neutron_mass = 939.57;
      float binding_energy = 30.0;
      float muon_mass = 105.66;
      float electron_mass = 0.511;
      float proton_mass = 938.27;

      float muon_costheta;
      TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
      if (pfeval.reco_muonMomentum[3]>0)
        muon_costheta = TMath::Cos(muonMomentum.Theta());

      float shower_costheta;
      TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);
      if (pfeval.reco_showerMomentum[3]>0)
        shower_costheta = TMath::Cos(showerMomentum.Theta());

      shower_costheta = TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi());

      float reco_Emuon =  pfeval.reco_muonMomentum[3]*1000.; // GeV --> MeV
      float reco_Eqe_muon = 0.5 * (2*(neutron_mass-binding_energy)*reco_Emuon - (pow(neutron_mass-binding_energy,2) + pow(muon_mass,2) - pow(proton_mass,2))) / ((neutron_mass-binding_energy) - reco_Emuon + sqrt(pow(reco_Emuon,2)-pow(muon_mass,2))*muon_costheta);

      float reco_Eelectron= pfeval.reco_showerMomentum[3]*1000.;

      reco_Eelectron = get_reco_showerKE_corr(pfeval, flag_data) * 1000.;


      float reco_Eqe_electron = 0.5 * (2*(neutron_mass-binding_energy)*reco_Eelectron - (pow(neutron_mass-binding_energy,2) + pow(electron_mass,2) - pow(proton_mass,2))) / ((neutron_mass-binding_energy) - reco_Eelectron + sqrt(pow(reco_Eelectron,2)-pow(electron_mass,2))*shower_costheta);

      if (shower_costheta > 0.999 || shower_costheta < -0.999) {
          reco_Eqe_electron = -1;
      }

      if(var_name=="reco_Eqe_muon") return reco_Eqe_muon;
      else if(var_name=="reco_Eqe_electron") return reco_Eqe_electron;
      else if(var_name=="reco_Eqe_muon_Enu_diff") return reco_Eqe_muon - get_reco_Enu_corr(kine, flag_data);
      else if(var_name=="reco_Eqe_electron_Enu_diff") return reco_Eqe_electron - get_reco_Enu_corr(kine, flag_data);
  }else if (var_name == "proton_KE"){
      return pfeval.reco_protonMomentum[3]*1000.-938.27; // GeV--> MeV
  }else if (var_name == "proton_theta"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Theta()/TMath::Pi()*180.;
  }else if (var_name == "proton_phi"){
      TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
      return protonMomentum.Phi()/TMath::Pi()*180.;
   }else if (var_name == "shower_theta" || var_name == "shower_costheta" || var_name == "shower_costheta" || var_name == "shower_phi"){
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);

    if(var_name == "shower_theta")
      return tagger.mip_angle_beam;
      //return showerMomentum.Theta()/TMath::Pi()*180.;

    if(var_name == "shower_costheta"){
      return TMath::Cos(tagger.mip_angle_beam/180.*TMath::Pi());
    }

    if(var_name == "shower_phi"){
      if (pfeval.reco_showerMomentum[3]>0)
        return showerMomentum.Phi()/TMath::Pi()*180.;
      else
        return -1000;
    }
   }else if (var_name=="shower_proton_angle_sum"){
    TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
    TLorentzVector showerMomentum(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1], pfeval.reco_showerMomentum[2], pfeval.reco_showerMomentum[3]);

    if(pfeval.reco_showerMomentum[3]>0 && pfeval.reco_protonMomentum[3]>0)
      return showerMomentum.Theta()/TMath::Pi()*180. + protonMomentum.Theta()/TMath::Pi()*180.;
    else
      return -1000;
  }
  else if (var_name=="muon_proton_angle_sum"){
    TLorentzVector protonMomentum(pfeval.reco_protonMomentum[0], pfeval.reco_protonMomentum[1], pfeval.reco_protonMomentum[2], pfeval.reco_protonMomentum[3]);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);

    if(pfeval.reco_muonMomentum[3]>0 && pfeval.reco_protonMomentum[3]>0)
      return muonMomentum.Theta()/TMath::Pi()*180. + protonMomentum.Theta()/TMath::Pi()*180.;
    else
      return -1000;
  }else if (var_name == "Ehadron"){
    if (pfeval.reco_muonMomentum[3]>0)
      return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    else
      return -1000;
    //  }else if (var_name == "Ehadron"){
      /* Float_t Ehadron = kine.kine_reco_Enu; */
      /* for(size_t i=0; i<kine.kine_energy_particle->size(); i++) */
      /* { */
      /*     int pdgcode = kine.kine_particle_type->at(i); */
      /*     if(abs(pdgcode)==13) Ehadron = Ehadron - kine.kine_energy_particle->at(i) - 105.658; */
      /*     //if(abs(pdgcode)==11) Ehadron = Ehadron - kine.kine_energy_particle->at(i); */
      /* } */
    // return kine.kine_reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
  /*}else if (var_name == "Ehadron_hybrid"){
    if (eval.match_isFC)
      return 1000.0*(kine.vlne_v4_numu_full_totalE - kine.vlne_v4_numu_full_primaryE);
    else {
      if (pfeval.reco_muonMomentum[3]>0)
        return get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
      else
        return -1000;
    }
  }else if (var_name == "Ehadron_dlnew"){
      return 1000.0*(kine.vlne_v4_numu_full_totalE - kine.vlne_v4_numu_full_primaryE);*/
  }else if (var_name == "Q2"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    // Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
    //  }else if (var_name == "Q2"){
    // Float_t Enu = kine.kine_reco_Enu;
    //Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    //Float_t Ehadron = Enu - Emu;
    //Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    //TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    //Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    //return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(1000.*1000.); // GeV^2
  }else if (var_name == "x_Bjorken"){
    Float_t Enu = get_reco_Enu_corr(kine, flag_data);
    Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    Float_t Ehadron = Enu - Emu;
    Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
    //  }else if (var_name == "x_Bjorken"){
    // Float_t Enu = kine.kine_reco_Enu;
    // Float_t Emu = pfeval.reco_muonMomentum[3]*1000.;
    // Float_t Ehadron = Enu - Emu;
    // Float_t Pmu = TMath::Sqrt(Emu*Emu - 105.658*105.658);
    // TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    // Float_t cosTheta = TMath::Cos(muonMomentum.Theta());
    // return (2*Enu*(Emu-Pmu*cosTheta)-105.658*105.658)/(2*938.272*Ehadron);
  }else if (var_name == "N_tracks"){
      int N_tracks = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==11) continue;
          if(kine.kine_energy_particle->at(i)<10) continue;
          if(abs(pdgcode)==13 || abs(pdgcode)==211){
            N_tracks += 1;
          }
          else if(kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_tracks += 1;
          }
      }
      return N_tracks;
   }else if (var_name == "N_other_tracks"){
        int Nothertracks = 0;
        for(size_t i=0; i<kine.kine_energy_particle->size(); i++){
                int pdgcode = kine.kine_particle_type->at(i);
                if((abs(pdgcode)==211 || abs(pdgcode)==13) && kine.kine_energy_particle->at(i)>10) Nothertracks++; // KE threshold: 10 MeV
        }
        return Nothertracks;
   }else if (var_name == "N_showers"){
      int N_showers = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)!=11) continue;
          if(kine.kine_energy_particle->at(i)>10) N_showers += 1;
      }
      return N_showers;
  }else if (var_name == "N_protons"){
      int N_protons = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)== 2212 && kine.kine_energy_particle->at(i)>35){ // proton KE threshold
              N_protons += 1;
          }
      }
      return N_protons;
  }else if (var_name == "N_true_protons"){
    if (flag_data) return -1;
    int np = 0;
    for(size_t i=0; i<pfeval.truth_Ntrack; i++){
      if(pfeval.truth_mother[i] != 0) continue;
      if(pfeval.truth_pdg[i] != 2212) continue;
      if(pfeval.truth_startMomentum[i][3] - 0.938272 < 0.035) continue;
      np++;
    }
    return np;
  }else if (var_name == "EhadShwrFrac"){
      double EhadShwr=0, EhadTot=0;

      if (pfeval.reco_muonMomentum[3]>0) {
        EhadTot = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
        for ( size_t j=0;j!= kine.kine_energy_particle->size();j++){
          if (kine.kine_energy_info->at(j) == 2 && kine.kine_particle_type->at(j) == 11){
            EhadShwr +=  kine.kine_energy_particle->at(j);
          }

        }
        return EhadShwr/ EhadTot;

      }
      else return -1;
  }else if (var_name == "reco_mcc8_pmuoncosth_Enu"){
    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon momentum
    float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66;
    float pmuon = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66) / 1000.0;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());
    // flattened Enu
    int indx = mcc8_pmuon_costheta_bin(pmuon, costh); // index of pmuon-costh
    double reco_Enu = get_reco_Enu_corr(kine, flag_data) / 100.0;
    if (reco_Enu<0) return -10000;
    else if (reco_Enu>25.0) return 10000; // overflow bin
    else return (indx-1)*25.0 + reco_Enu;
  }
  else if (var_name == "reco_concatenated_Pmuon"){

    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon momentum
    float KE_muon = pfeval.reco_muonMomentum[3]*1000.-105.66;
    float pmuon = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66) / 1000.0;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());

    float pmuon_MeV = pmuon * 1000.0;
    if (pmuon_MeV > 1500.0) return -10000;
    if (costh>=-1 and costh<-0.5) { return pmuon_MeV; }
    else if (costh>=-0.5 and costh<0){ return pmuon_MeV + 1500*1; }
    else if (costh>=0 and costh<0.27){ return pmuon_MeV + 1500*2; }
    else if (costh>=0.27 and costh<0.45){ return pmuon_MeV + 1500*3; }
    else if (costh>=0.45 and costh<0.62){ return pmuon_MeV + 1500*4; }
    else if (costh>=0.62 and costh<0.76){ return pmuon_MeV + 1500*5; }
    else if (costh>=0.76 and costh<0.86){ return pmuon_MeV + 1500*6; }
    else if (costh>=0.86 and costh<0.94){ return pmuon_MeV + 1500*7; }
    else if (costh>=0.94 and costh<=1.00){ return pmuon_MeV + 1500*8; }

    return -10000;
  }
  else if (var_name == "muon_momentum_costheta"){
    float muon_momentum = get_kine_var(kine, eval, pfeval, tagger, flag_data , "muon_momentum", space, pandora, lantern);
    float costheta = get_kine_var(kine, eval, pfeval, tagger, flag_data , "muon_costheta", space, pandora, lantern);
    int bin = alt_var_index("muon_momentum",muon_momentum, "costheta", costheta);
    return bin;
  }
  else if (var_name == "Ehad_muon_costheta"){
    float Ehadron = -1000;
    if (pfeval.reco_muonMomentum[3]>0)
      Ehadron = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.;
    float costheta = get_kine_var(kine, eval, pfeval, tagger, flag_data , "muon_costheta", space, pandora, lantern);
    int bin = alt_var_index("Ehadron",Ehadron, "costheta", costheta);
    return bin;
  }
  else if (var_name == "reco_concatenated_Ehad"){

    if (pfeval.reco_muonMomentum[3]<0) return -10000;
    // muon costheta
    TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());

    float Ehad = get_reco_Enu_corr(kine, flag_data) - pfeval.reco_muonMomentum[3]*1000.0;

    if (Ehad > 1500.0) return -10000;
    if (costh>=-1 and costh<-0.5) { return Ehad; }
    else if (costh>=-0.5 and costh<0){ return Ehad + 1500*1; }
    else if (costh>=0 and costh<0.27){ return Ehad + 1500*2; }
    else if (costh>=0.27 and costh<0.45){ return Ehad + 1500*3; }
    else if (costh>=0.45 and costh<0.62){ return Ehad + 1500*4; }
    else if (costh>=0.62 and costh<0.76){ return Ehad + 1500*5; }
    else if (costh>=0.76 and costh<0.86){ return Ehad + 1500*6; }
    else if (costh>=0.86 and costh<0.94){ return Ehad + 1500*7; }
    else if (costh>=0.94 and costh<=1.00){ return Ehad + 1500*8; }

    return -10000;
  }else if (var_name == "proton_pi0_total_momentum" || var_name == "proton_pi0_invariant_mass") {
        bool debug_pf_info = 0;
        TLorentzVector max_energy_proton_momentum(-1., -1., -1., -1.);
        TLorentzVector gamma_1_momentum(-1., -1., -1., -1.);
        TLorentzVector gamma_2_momentum(-1., -1., -1., -1.);
        float max_proton_energy = 0.;
        if (debug_pf_info) {
                std::cout << "********************************* starting event  ************************\n";
                std::cout << "(100 hard coded) looping over " << pfeval.reco_Ntrack << " reco particles (" << pfeval.truth_Ntrack << " true particles)\n";
        }
        for(size_t i=0; i<100; i++) // looping over all reconstructed particles
        {
                int pdgcode = pfeval.reco_pdg[i];
                if (debug_pf_info) std::cout << "investigating particle" << i << "//" << pfeval.reco_Ntrack << " : " << pdgcode << "\n";
                if(abs(pdgcode)==2212 && pfeval.reco_startMomentum[i][3] > max_proton_energy){ // new max energy proton
                        if (debug_pf_info) std::cout << "new max energy proton\n";
                        max_energy_proton_momentum = pfeval.reco_startMomentum[i];
                        max_proton_energy = pfeval.reco_startMomentum[i][3];
                }
                if(abs(pdgcode)==22 || abs(pdgcode)==11) { // reconstructed shower
                        float shower_energy = 1000. * pfeval.reco_startMomentum[i][3];
                        if (abs(shower_energy - kine.kine_pio_energy_1) / kine.kine_pio_energy_1 < 0.01) { // very close to gamma 1 energy
                                if (debug_pf_info) std::cout << "gamma 1 matched";
                                if (flag_data) {
                                        gamma_1_momentum = TLorentzVector(em_charge_scale * pfeval.reco_startMomentum[i][0],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][1],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][2],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][3]);
                                } else {
                                        gamma_1_momentum = pfeval.reco_startMomentum[i];
                                }
                        }
                        if (abs(shower_energy - kine.kine_pio_energy_2) / kine.kine_pio_energy_2 < 0.01) { // very close to gamma 2 energy
                                if (debug_pf_info) std::cout << "gamma 2 matched";
                                if (flag_data) {
                                        gamma_2_momentum = TLorentzVector(em_charge_scale * pfeval.reco_startMomentum[i][0],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][1],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][2],
                                                        em_charge_scale * pfeval.reco_startMomentum[i][3]);
                                } else {
                                        gamma_2_momentum = pfeval.reco_startMomentum[i];
                                }
                        }

                }
        }

        float proton_pi0_invariant_mass = -1.;
        float proton_pi0_total_momentum = -1.;

        float pi0_mass = 134.9768;
        float proton_mass = 938.272;

        if (max_energy_proton_momentum[3]>0 && gamma_1_momentum[3]>0 && gamma_2_momentum[3]>0) { // found proton and pi0 in PF tree
                if (debug_pf_info) std::cout << "matched proton and pi0!" << std::endl;
                TLorentzVector pi0_momentum = gamma_1_momentum + gamma_2_momentum;
                proton_pi0_invariant_mass = sqrt(pi0_mass * pi0_mass + proton_mass * proton_mass
                                + 2. * 1000. * 1000. * (max_energy_proton_momentum[3] * pi0_momentum[3]
                                                      - max_energy_proton_momentum[0] * pi0_momentum[0]
                                                      - max_energy_proton_momentum[1] * pi0_momentum[1]
                                                      - max_energy_proton_momentum[2] * pi0_momentum[2]));
                proton_pi0_total_momentum = 1000. * sqrt((max_energy_proton_momentum[0] + pi0_momentum[0]) * (max_energy_proton_momentum[0] + pi0_momentum[0])
                                                       + (max_energy_proton_momentum[1] + pi0_momentum[1]) * (max_energy_proton_momentum[1] + pi0_momentum[1])
                                                       + (max_energy_proton_momentum[2] + pi0_momentum[2]) * (max_energy_proton_momentum[2] + pi0_momentum[2]));
        }

        if (debug_pf_info) std::cout << "******************* ending event *************\n";

        if (var_name == "proton_pi0_total_momentum") {
                return proton_pi0_total_momentum;
        } else if (var_name == "proton_pi0_invariant_mass") {
                return proton_pi0_invariant_mass;
        } else {
                std::cout << "No such proton-pi0 variable: " << var_name << std::endl;
        }
  //Erin
  }else if (var_name == "single_photon_numu_score"){
    return tagger.single_photon_numu_score;
  }else if (var_name == "single_photon_other_score"){
    return tagger.single_photon_other_score;
  }else if (var_name == "single_photon_ncpi0_score"){
    return tagger.single_photon_ncpi0_score;
  }else if (var_name == "single_photon_nue_score"){
    return tagger.single_photon_nue_score;
  }else if (var_name == "shower_energy_sp"){
    if(flag_data)
      return tagger.shw_sp_energy*em_charge_scale;
    else
      return tagger.shw_sp_energy;
  }else if (var_name == "shower_angle_beam_sp"){
    return tagger.shw_sp_angle_beam;
  }else if (var_name == "cos_shower_angle_beam_sp"){
    return TMath::Cos(tagger.shw_sp_angle_beam/180.*TMath::Pi());
  }else if (var_name == "num_shower_sp"){
    return tagger.shw_sp_n_20br1_showers;
  }else if (var_name == "median_dEdx_sp"){
    if(flag_data)
      return tagger.shw_sp_vec_median_dedx;//*em_charge_scale;
    else
      return tagger.shw_sp_vec_median_dedx;
  }else if (var_name == "median_dEdx_sp_15"){
    std::vector<float> dqdx;
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_2);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_3);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_4);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_5);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_6);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_7);
    dqdx.push_back(tagger.shw_sp_vec_dQ_dx_8);
    std::sort(dqdx.begin(), dqdx.end());
    size_t vecsize = dqdx.size();
    size_t mid = vecsize/2;
    float median_dqdx = vecsize%2==0 ? (dqdx[mid]+dqdx[mid-1])/2:dqdx[mid];
    float alpha = 1.;
    float beta = 0.255;
    float median_dedx = (exp((median_dqdx*43e3) * 23.6e-6*beta/1.38/0.273) - alpha)/(beta/1.38/0.273);
    if(median_dedx<0) median_dedx = 0;
    if(median_dedx>50) median_dedx = 50;
    if(flag_data)
      return median_dedx;//*em_charge_scale;
    else
      return median_dedx; // MeV/cm
  }else if (var_name == "dQdx_0_sp"){
    if(flag_data)
      return tagger.shw_sp_vec_dQ_dx_0;//*em_charge_scale;
    else
      return tagger.shw_sp_vec_dQ_dx_0;
  }else if (var_name == "dQdx_1_sp"){
    if(flag_data)
      return tagger.shw_sp_vec_dQ_dx_1;//*em_charge_scale;
    else
      return tagger.shw_sp_vec_dQ_dx_1;
  }else if (var_name == "shw_vtx_dis_sp"){
    return tagger.shw_sp_shw_vtx_dis;
  }else if (var_name == "max_shw_dis_sp"){
    return tagger.shw_sp_max_shw_dis;
  }else if (var_name == "shw_sp_br3_1_n_shower_segments"){
    return tagger.shw_sp_br3_1_n_shower_segments;
  }else if (var_name == "shw_sp_E_indirect_max_energy"){
    return tagger.shw_sp_E_indirect_max_energy;
  }else if (var_name == "shw_sp_lem_n_3seg"){
    return tagger.shw_sp_lem_n_3seg;
  }else if (var_name == "shw_sp_lem_n_3seg"){
    return tagger.shw_sp_lem_n_3seg;
  }else if (var_name == "kine_pio_phi_2"){
    return kine.kine_pio_phi_2;
  }else if (var_name == "shw_sp_pio_flag_pio"){
    return tagger.shw_sp_pio_flag_pio;
  }else if (var_name == "shw_sp_length_total"){
    return tagger.shw_sp_length_total;
  }else if (var_name == "shw_sp_n_vertex"){
    return tagger.shw_sp_n_vertex;
  }else if (var_name == "shw_backwards_projected_dist"){
    float backwards_projected_dist = -99999.0;

    if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);


        float min_backwards_projected_dist = 1e9;

        //projecting to x walls
        if (shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
              min_backwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2];
        }
        if (isinf(min_backwards_projected_dist)) min_backwards_projected_dist = -99999.0;

        backwards_projected_dist = min_backwards_projected_dist;
      }

      return backwards_projected_dist;
  }else if (var_name == "shw_forwards_projected_dist"){
    float forwards_projected_dist = -99999.0;

    if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        std::vector<float> inv_shower_unit_vector_3d = {-(pfeval.reco_showerMomentum[0] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[1] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[2] / shower_momentum_total_3d)};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);


        float min_forwards_projected_dist = 1e9;

        //projecting to x walls
        if (inv_shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
              min_forwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (inv_shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / inv_shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / inv_shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (inv_shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2];
        }
        if (isinf(min_forwards_projected_dist)) min_forwards_projected_dist = -99999.0;

        forwards_projected_dist = min_forwards_projected_dist;
      }

      return forwards_projected_dist;
    }else if (var_name == "shw_min_dist"){
    float minimum_dist = -99999.0;

    if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        std::vector<float> inv_shower_unit_vector_3d = {-(pfeval.reco_showerMomentum[0] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[1] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[2] / shower_momentum_total_3d)};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);


        float min_dist = 1e9;

        //x walls
        if (sqrt(pow(pfeval.reco_showervtxX,2) - pow(-1.0,2)) < min_dist)
            min_dist =  sqrt(pow(pfeval.reco_showervtxX,2) - pow(-1.0,2));
        if (sqrt(pow(pfeval.reco_showervtxX,2) - pow(254.3,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxX,2) - pow(254.3,2));
      
        //y walls
        if (sqrt(pow(pfeval.reco_showervtxY,2) - pow(-115.0,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxY,2) - pow(-115.0,2));
        if (sqrt(pow(pfeval.reco_showervtxY,2) - pow(117.0,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxY,2) - pow(117.0,2));
        
        //projecting to z walls
        if (sqrt(pow(pfeval.reco_showervtxZ,2) - pow(0.6,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxZ,2) - pow(0.6,2));
        if (sqrt(pow(pfeval.reco_showervtxZ,2) - pow(1036.4,2)) < min_dist)
            min_dist = sqrt(pow(pfeval.reco_showervtxZ,2) - pow(1036.4,2));
        
        if (isinf(min_dist)) min_dist = -99999.0;

        minimum_dist = min_dist;
      }

      return minimum_dist;
    }else if (var_name == "shw_projected_dist"){
      float projected_dist = -99999.0;
      float forwards_projected_dist = -99999.0;
      float backwards_projected_dist = -99999.0;

      if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        std::vector<float> inv_shower_unit_vector_3d = {-(pfeval.reco_showerMomentum[0] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[1] / shower_momentum_total_3d),
                                 -(pfeval.reco_showerMomentum[2] / shower_momentum_total_3d)};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);

        //forwards projected dist
        float min_forwards_projected_dist = 1e9;

        //projecting to x walls
        if (inv_shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
              min_forwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / inv_shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / inv_shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (inv_shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / inv_shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / inv_shower_unit_vector_3d[1] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / inv_shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (inv_shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / inv_shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2] < min_forwards_projected_dist)
            min_forwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / inv_shower_unit_vector_3d[2];
        }
        if (isinf(min_forwards_projected_dist)) min_forwards_projected_dist = -99999.0;

        forwards_projected_dist = min_forwards_projected_dist;

        //backwards projected dist
        float min_backwards_projected_dist = 1e9;

        //projecting to x walls
        if (shower_unit_vector_3d[0] > 0){
            if ((pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
              min_backwards_projected_dist =  (pfeval.reco_showervtxX - (-1.0)) / shower_unit_vector_3d[0];
        }else{
          if ((pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxX - (254.3)) / shower_unit_vector_3d[0];
        }
        //projecting to y walls
        if (shower_unit_vector_3d[1] > 0){
          if ((pfeval.reco_showervtxY - (-115.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (-115.)) / shower_unit_vector_3d[1];
        }else{
          if ((pfeval.reco_showervtxY - (117.0)) / shower_unit_vector_3d[1] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxY - (117.)) / shower_unit_vector_3d[1];
        }
        //projecting to z walls
        if (shower_unit_vector_3d[2] > 0){
          if ((pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (0.6)) / shower_unit_vector_3d[2];
        }else{
          if ((pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2] < min_backwards_projected_dist)
            min_backwards_projected_dist = (pfeval.reco_showervtxZ - (1036.4)) / shower_unit_vector_3d[2];
        }
        if (isinf(min_backwards_projected_dist)) min_backwards_projected_dist = -99999.0;

        backwards_projected_dist = min_backwards_projected_dist;

        if (forwards_projected_dist < backwards_projected_dist) projected_dist = forwards_projected_dist;
        else projected_dist = backwards_projected_dist;
      }

      return projected_dist;
  }else if (var_name == "shw_inwardness_3d"){
    }else if (var_name == "shw_backwards_projected_dist"){
      float inwardness = -99999.0;
    
      if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);
        inwardness = inwardness_3d;
      }
      return inwardness;
  }else if (var_name == "shw_inwardness_2d"){
    }else if (var_name == "shw_backwards_projected_dist"){
      float inwardness = -99999.0;
    
      if (pfeval.reco_showerMomentum[3] > 0){
        float reco_shower_momentum_perp = sqrt(pow(pfeval.reco_showerMomentum[0],2) + pow(pfeval.reco_showerMomentum[1],2));
        float shower_theta = atan2(reco_shower_momentum_perp, pfeval.reco_showerMomentum[2]) * (180. / TMath::Pi());
        float shower_phis = atan2(pfeval.reco_showerMomentum[0], pfeval.reco_showerMomentum[1]) * (180. / TMath::Pi());

        float shower_momentum_total_3d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1] +
                                           pfeval.reco_showerMomentum[2] * pfeval.reco_showerMomentum[2]);
        std::vector<float> shower_unit_vector_3d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[2] / shower_momentum_total_3d};
        float center_x = 130.;
        float center_y = 0.;
        float center_z = 525.;
        float towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y) +
                                        (pfeval.reco_showervtxZ - center_z) * (pfeval.reco_showervtxZ - center_z));
        std::vector<float> towards_center_unit_vector_3d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length,
                                         (center_z - pfeval.reco_showervtxZ) / towards_center_length};
        float inwardness_3d = (shower_unit_vector_3d[0] * towards_center_unit_vector_3d[0]) +
                             (shower_unit_vector_3d[1] * towards_center_unit_vector_3d[1]) +
                             (shower_unit_vector_3d[2] * towards_center_unit_vector_3d[2]);

        float shower_momentum_total_2d = sqrt(pfeval.reco_showerMomentum[0] * pfeval.reco_showerMomentum[0] +
                                           pfeval.reco_showerMomentum[1] * pfeval.reco_showerMomentum[1]);
        std::vector<float> shower_unit_vector_2d = {pfeval.reco_showerMomentum[0] / shower_momentum_total_3d,
                                 pfeval.reco_showerMomentum[1] / shower_momentum_total_3d};
        towards_center_length = sqrt((pfeval.reco_showervtxX - center_x) * (pfeval.reco_showervtxX - center_x) +
                                        (pfeval.reco_showervtxY - center_y) * (pfeval.reco_showervtxY - center_y));
        std::vector<float> towards_center_unit_vector_2d = {(center_x - pfeval.reco_showervtxX) / towards_center_length,
                                         (center_y - pfeval.reco_showervtxY) / towards_center_length};
        float inwardness_2d = (shower_unit_vector_2d[0] * towards_center_unit_vector_2d[0]) +
                             (shower_unit_vector_2d[1] * towards_center_unit_vector_2d[1]);
        inwardness = inwardness_2d;
      }
      return inwardness;
  }else if (var_name == "ns_beam_time"){
    if(flag_data){
      double delta_time_calc = -9999.;
      //Merge Peaks
      double gap=18.936;
      double Shift=0;
      double TThelp=0;
      if (pfeval.run >= 17380){ Shift=2916.0; }
      else if (pfeval.run >= 13697){ Shift = 3147.3;}//3166.1;}
      else if (pfeval.run >= 10812){ Shift = 3568.5; }
      else if (pfeval.run >= 8321){ Shift = 3610.7;}
      else if (pfeval.run >= 5800){ Shift = 3164.4;}
      else if (pfeval.run >= 0){ Shift = 3168.9;}
      //else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}
      /*if (pfeval.run >= 13697){ Shift = 3166.9;}
      else if(pfeval.run>=10812){ Shift = 3568.5; }
      else if (pfeval.run >= 8321){ Shift = 3610.7;}
      else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}*/
      //9.43;}
      //if(run>8000 && run<10812){Shift=3610.7; }
      //if(run>=10812 && run <12500){Shift=3568.5; }
      TThelp=pfeval.evtTimeNS-Shift+gap*0.5;
      double TT_merged = -9999.;

      //merge peaks
      if(TThelp>=0 && TThelp<gap*81.0){
        TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
      }

      delta_time_calc = TT_merged;

      return delta_time_calc;
    }else{
      return -9999.;
    }
  }else if (var_name == "flag_0p"){
    if (is_0p(tagger, kine, pfeval)){
      return 1;
    }else{
      return 0;
    }
  }else if (var_name == "flag_FC"){
    if (is_FC(eval)){
      return 1;
    }else{
      return 0;
    }
  }else if (var_name == "run_period"){
    if (pfeval.run >= 13697){ 
      return 3;
    }
    else if (pfeval.run >= 8321){ 
      return 2;
    }
    else if (pfeval.run > 0 ){ 
      return 1;
    }
  }else if (var_name == "run_period_100_300"){
    float shwen = -1;
    if(flag_data)
      shwen = tagger.shw_sp_energy*em_charge_scale;
    else
      shwen = tagger.shw_sp_energy;
    
    if (shwen >= 100.0 && shwen < 300){
      if (pfeval.run >= 13697){ 
        return 3;
      }
      else if (pfeval.run >= 8321){ 
        return 2;
      }
      else if (pfeval.run > 0 ){ 
        return 1;
      }
    }
    return -1;
    
  //
  }else{
    std::cout << "No such variable: " << var_name << std::endl;
    exit(EXIT_FAILURE);
  }
  return -1;
}

int get_costheta_bin (float costh) {
  int nbins = 9;
  float costheta_binning[nbins+1] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};
  if (costh == costheta_binning[0]) { return 0; }
  for (int i=0;i<nbins;i++) { if (costh >  costheta_binning[i] && costh <= costheta_binning[i+1]) { return i; } }
  return -1;
}

int get_Pmuon_bin (float Pmuon) {
  int nbins = 6;
  float pmuon_binning[nbins+1] = {0, 180, 300, 450, 770, 1280, 2500};
  if (Pmuon == pmuon_binning[0]) { return 0; }
  for (int i=0;i<nbins;i++) { if (Pmuon >  pmuon_binning[i] && Pmuon <= pmuon_binning[i+1]) { return i; } }
  return -1;
}

int get_Enu_bin (float Enu) {
  int nbins = 4;
  float enu_binning[nbins+1] = {200, 705, 1050, 1570, 4000};
  if (Enu == enu_binning[0]) { return 0; }
  for (int i=0;i<nbins;i++) { if (Enu >  enu_binning[i] && Enu <= enu_binning[i+1]) { return i; } }
  return -1;
}

int LEEana::get_xs_signal_no(int cut_file, std::map<TString, int>& map_cut_xs_bin, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  for (auto it = map_cut_xs_bin.begin(); it != map_cut_xs_bin.end(); it++){
    TString cut_name = it->first;
    int number = it->second;

    double KE_muon = pfeval.truth_muonMomentum[3]*1000.-105.66; // MeV
    double pmuon   = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66); // MeV
    double Pmuon   = TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66); // MeV
    double Emuon   = pfeval.truth_muonMomentum[3]*1000; // MeV
    double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

    float pmuon_binning[7] = {0, 180, 300, 450, 770, 1280, 2500};

    float costheta_binning[10] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};			//fine binning
    //float costheta_binning[5]  = {-1,         .27,      .62,      .86,      1};		//coarse binning
    //float costheta_binning[3]    = {-1,                   .62,                1};		//very coarse binning
    TLorentzVector muonMomentum(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);
    float costh = TMath::Cos(muonMomentum.Theta());

    if (cut_file == 1){
      if (cut_name == "numuCC.inside.Enu.le.300"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 ) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.400.gt.300"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=400 && eval.truth_nuEnergy>300) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.500.gt.400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=500 && eval.truth_nuEnergy>400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.600.gt.500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=600 && eval.truth_nuEnergy>500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.700.gt.600"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=700 && eval.truth_nuEnergy>600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.800.gt.700"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=800 && eval.truth_nuEnergy>700) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.900.gt.800"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=900 && eval.truth_nuEnergy>800) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1000.gt.900"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1000 && eval.truth_nuEnergy>900) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1100.gt.1000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1100 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1100"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1100) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1000) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1500.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1500 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2100.gt.1500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2100 && eval.truth_nuEnergy>1500) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1400.gt.1200"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1400 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1600.gt.1400"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1600 && eval.truth_nuEnergy>1400) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2000.gt.1600"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2000 && eval.truth_nuEnergy>1600) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2500.gt.2000"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2500 && eval.truth_nuEnergy>2000) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2500) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.2100"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2100) return number;
      }else if (cut_name == "numuCC.inside.Enu.gt.1500"){
	if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>1500) return number;
      }else{
	         std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 2) {
      if (cut_name == "numuCC.inside.Emuon.le.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=100 && Emuon>0) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.200.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=200 && Emuon>100) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.300.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=300 && Emuon>200) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.400.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=400 && Emuon>300) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.500.gt.400"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=500 && Emuon>400) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.600.gt.500"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=600 && Emuon>500) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.700.gt.600"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=700 && Emuon>600) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.800.gt.700"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=800 && Emuon>700) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.900.gt.800"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=900 && Emuon>800) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1000.gt.900"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1000 && Emuon>900) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1200.gt.1000"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1200 && Emuon>1000) return number;
      }
      else if (cut_name == "numuCC.inside.Emuon.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1200) return number;
      }else{
          std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 3) {
      if (cut_name == "numuCC.inside.Ehadron.le.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.200.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=200 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.300.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=300 && Ehadron>200) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.400.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=400 && Ehadron>300) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.500.gt.400"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=500 && Ehadron>400) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.600.gt.500"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=600 && Ehadron>500) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.700.gt.600"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=700 && Ehadron>600) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.800.gt.700"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=800 && Ehadron>700) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.900.gt.800"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=900 && Ehadron>800) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1000.gt.900"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1000 && Ehadron>900) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.gt.1000"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1000) return number;
      }else{
          std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 4){
      bool pre_cut = eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && muonMomentum[3]>0 && Pmuon > 0 && Pmuon <= 2500;
      if (cut_name == "numuCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
	if (pre_cut && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){
	if (pre_cut && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.805.gt.705"){
	if (pre_cut && eval.truth_nuEnergy<=805 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.920.gt.805"){
	if (pre_cut && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>805) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1050.gt.920"){
	if (pre_cut && eval.truth_nuEnergy<=1050 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1050"){
	if (pre_cut && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1050) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1375.gt.1200"){
	if (pre_cut && eval.truth_nuEnergy<=1375 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1375"){
	if (pre_cut && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1375) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){
	if (pre_cut && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
	if (pre_cut && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 5) {
      if (cut_name == "numuCC.inside.Emuon.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984) return number;
      }else if (cut_name == "numuCC.inside.Emuon.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 6) {
      if (cut_name == "numuCC.inside.Ehadron.le.100.gt.30"){ // 30 MeV - 100 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=100 && Ehadron>30) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.150.gt.100"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=150 && Ehadron>100) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.225.gt.150"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=225 && Ehadron>150) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.275.gt.225"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=275 && Ehadron>225) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.336.gt.275"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=336 && Ehadron>275) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.411.gt.336"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=411 && Ehadron>336) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.502.gt.411"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=502 && Ehadron>411) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.614.gt.502"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=614 && Ehadron>502) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.750.gt.614"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=750 && Ehadron>614) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1120.gt.750"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1120 && Ehadron>750) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.2500.gt.1120"){ // 1120 - 2500 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1120 && Ehadron<= 2500) return number;
      }else{
	std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 7){
      if (cut_name == "numuCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.920.gt.705"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1200.gt.920"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 8){
      if (cut_name == "numuCC.inside.Enu.le.4000.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 9){
      if (cut_name == "numuCC.inside.Enu.le.1200.gt.200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "numuCC.inside.Enu.le.4000.gt.1200"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>1200) return number;
      }
      else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 10) {
      if (cut_name == "numuCC.inside.Ehadron.le.150.gt.30"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=150 && Ehadron>30) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.275.gt.150"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=275 && Ehadron>150) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.411.gt.275"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=411 && Ehadron>275) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.502.gt.411"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=502 && Ehadron>411) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.614.gt.502"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=614 && Ehadron>502) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.750.gt.614"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=750 && Ehadron>614) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.1120.gt.750"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron<=1120 && Ehadron>750) return number;
      }else if (cut_name == "numuCC.inside.Ehadron.le.2500.gt.1120"){ // 1120 - 2500 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron>1120 && Ehadron<= 2500) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 11){
      if (cut_name == "nueCC.inside.Enu.le.540.gt.200"){ // recommended range: 200 - 540
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=540 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.705.gt.540"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=705 && eval.truth_nuEnergy>540) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.920.gt.705"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=920 && eval.truth_nuEnergy>705) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.1200.gt.920"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>920) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.1570.gt.1200"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.2050.gt.1570"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.4000.gt.2050"){ // recommended range: 2050 - 4000
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy>2050 && eval.truth_nuEnergy<=4000) return number;
      }else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }
    else if (cut_file == 12){
      if (cut_name == "nueCC.inside.Enu.le.1200.gt.200"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>200) return number;
      }else if (cut_name == "nueCC.inside.Enu.le.4000.gt.1200"){
        if (eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>1200) return number;
      }
      else{
        std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }

    //MCC8 binning -> 40 bins, muon mometum
    else if (cut_file == 13){
      if  (cut_name == "numuCC.inside.Pmuon.theta0.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;

      }else{
	std::cout << "get_xs_signal_no: no cut found!   cut_name = " << cut_name << std::endl;
      }
    }

    //MCC8 binning -> 36 bins, muon mometum
    else if (cut_file==14) {

      if       (cut_name == "numuCC.inside.Pmuon.theta0.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta0.le.2500.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta1.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.180.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[1] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.300.gt.180"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon> pmuon_binning[1] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta2.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta3.le.2500.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta4.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta5.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta6.le.2500.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta7.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) return number;

      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.300.gt.0"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[2] && Pmuon>=pmuon_binning[0] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.450.gt.300"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[3] && Pmuon> pmuon_binning[2] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.770.gt.450"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[4] && Pmuon> pmuon_binning[3] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.1280.gt.770"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[5] && Pmuon> pmuon_binning[4] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;
      }else if (cut_name == "numuCC.inside.Pmuon.theta8.le.2500.gt.1280"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon<=pmuon_binning[6] && Pmuon> pmuon_binning[5] && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())> costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) return number;

      }else{
	std::cout << "get_xs_signal_no: no cut found!   cut_name = " << cut_name << std::endl;
      }

    }

    //very coarse angle binning
    else if (cut_file == 15){
      if (number==-1) { std::cout << "cut_name, number = " << cut_name << ", " << number << std::endl; }
      if       (cut_name == "numuCC.inside.Emuon.theta0.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984  && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta0.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<costheta_binning[1])) return number;

      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.226.gt.106"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=226 && Emuon>105.7 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.296.gt.226"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=296 && Emuon>226   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.386.gt.296"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=386 && Emuon>296   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.505.gt.386"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=505 && Emuon>386   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.577.gt.505"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=577 && Emuon>505   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.659.gt.577"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=659 && Emuon>577   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.753.gt.659"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=753 && Emuon>659   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.861.gt.753"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=861 && Emuon>753   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.984.gt.861"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=984 && Emuon>861   && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.1285.gt.984"){
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon<=1285 && Emuon>984  && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }else if (cut_name == "numuCC.inside.Emuon.theta1.le.2506.gt.1285"){ // 1285 - 2506, only 1% > 2506 MeV
        if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon>1285 && Emuon<=2506 && (muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) return number;
      }
    }

    else if (cut_file == 16){ // pmuon, costh

      if (cut_name == "numuCC.inside.Pmuon.theta0.le.180.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<180.00 && costh>=-1.00 && costh<-0.50 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta0.le.300.gt.180"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=180.00 && pmuon<300.00 && costh>=-1.00 && costh<-0.50 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta0.le.2500.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<2500.00 && costh>=-1.00 && costh<-0.50 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.180.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<180.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.300.gt.180"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=180.00 && pmuon<300.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta1.le.2500.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<2500.00 && costh>=-0.50 && costh<0.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.180.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<180.00 && costh>=0.00 && costh<0.27 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.300.gt.180"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=180.00 && pmuon<300.00 && costh>=0.00 && costh<0.27 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.00 && costh<0.27 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta2.le.2500.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<2500.00 && costh>=0.00 && costh<0.27 ) return number;
        }
        else if (cut_name == "numuCC.inside.Pmuon.theta3.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.27 && costh<0.45 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta3.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.27 && costh<0.45 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta3.le.2500.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<2500.00 && costh>=0.27 && costh<0.45 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta4.le.2500.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<2500.00 && costh>=0.45 && costh<0.62 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta5.le.2500.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<2500.00 && costh>=0.62 && costh<0.76 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta6.le.2500.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<2500.00 && costh>=0.76 && costh<0.86 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.1280.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<1280.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta7.le.2500.gt.1280"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=1280.00 && pmuon<2500.00 && costh>=0.86 && costh<0.94 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.300.gt.0"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=0.00 && pmuon<300.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.450.gt.300"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=300.00 && pmuon<450.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.770.gt.450"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=450.00 && pmuon<770.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.1280.gt.770"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=770.00 && pmuon<1280.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else if (cut_name == "numuCC.inside.Pmuon.theta8.le.2500.gt.1280"){
      if (eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pmuon>=1280.00 && pmuon<2500.00 && costh>=0.94 && costh<1.00 ) return number;
        }
      else{
      std::cout << "get_xs_signal_no: no cut found!" << std::endl;
      }
    }

    else if (cut_file == 17){ // Enu, costheta, Pmuon

      int Enu_bin      = get_Enu_bin(eval.truth_nuEnergy);
      int costheta_bin = get_costheta_bin(costh);
      int Pmuon_bin    = get_Pmuon_bin(Pmuon);
      bool pre_cut     = eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && (muonMomentum[3]>0);

      if      (cut_name == "numuCC.inside.Enu0.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==0 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta1.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==0 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta4.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta5.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta6.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta7.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta7.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta8.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==0 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==0 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu0.theta8.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==0 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }

      else if (cut_name == "numuCC.inside.Enu1.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==1 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta1.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==1 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta4.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==1 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta5.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==5 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta6.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==6 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta7.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==7 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu1.theta8.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==1 && costheta_bin==8 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }

      else if (cut_name == "numuCC.inside.Enu2.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==2 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta1.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==2 && costheta_bin==1 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==2 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==2 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta4.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==4 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta5.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==5 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta6.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==6 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta7.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==2 && costheta_bin==7 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu2.theta8.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==2 && costheta_bin==8 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }

      else if (cut_name == "numuCC.inside.Enu3.theta0.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==0 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta0.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==0 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta0.Pmuon.le.2500.gt.300" ) { if (pre_cut && Enu_bin==3 && costheta_bin==0 && Pmuon_bin>=2 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta1.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==3 && costheta_bin==1 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.180.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=0 && Pmuon_bin<=0) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.300.gt.180"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=1 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta2.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==3 && costheta_bin==2 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta3.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==3 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta3.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==3 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta3.Pmuon.le.2500.gt.450" ) { if (pre_cut && Enu_bin==3 && costheta_bin==3 && Pmuon_bin>=3 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta4.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==4 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta5.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==5 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.300.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=0 && Pmuon_bin<=1) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta6.Pmuon.le.2500.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==6 && Pmuon_bin>=4 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.450.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=0 && Pmuon_bin<=2) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta7.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==3 && costheta_bin==7 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.450.gt.0"    ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=0 && Pmuon_bin<=2) { return number; } }
      //else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.450.gt.300"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=2 && Pmuon_bin<=2) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.770.gt.450"  ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=3 && Pmuon_bin<=3) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.1280.gt.770" ) { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=4 && Pmuon_bin<=4) { return number; } }
      else if (cut_name == "numuCC.inside.Enu3.theta8.Pmuon.le.2500.gt.1280") { if (pre_cut && Enu_bin==3 && costheta_bin==8 && Pmuon_bin>=5 && Pmuon_bin<=5) { return number; } }
      else { std::cout << "get_xs_signal_no: no cut found!" << std::endl; }

    }

    //1D Enu truth binning using the same inclusive selection as the 3D binning
    else if (cut_file == 18){

      int Enu_bin      = get_Enu_bin(eval.truth_nuEnergy);
      int costheta_bin = get_costheta_bin(costh);
      int Pmuon_bin    = get_Pmuon_bin(Pmuon);
      bool pre_cut     = eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && (muonMomentum[3]>0) && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8 && Pmuon_bin>=0 && Pmuon_bin<=5;

      if      (cut_name == "numuCC.inside.Enu.le.540.gt.200"){   if (pre_cut && eval.truth_nuEnergy<=540  && eval.truth_nuEnergy>200)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.705.gt.540"){   if (pre_cut && eval.truth_nuEnergy<=705  && eval.truth_nuEnergy>540)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.805.gt.705"){   if (pre_cut && eval.truth_nuEnergy<=805  && eval.truth_nuEnergy>705)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.920.gt.805"){   if (pre_cut && eval.truth_nuEnergy<=920  && eval.truth_nuEnergy>805)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1050.gt.920"){  if (pre_cut && eval.truth_nuEnergy<=1050 && eval.truth_nuEnergy>920)   { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1200.gt.1050"){ if (pre_cut && eval.truth_nuEnergy<=1200 && eval.truth_nuEnergy>1050)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1375.gt.1200"){ if (pre_cut && eval.truth_nuEnergy<=1375 && eval.truth_nuEnergy>1200)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.1570.gt.1375"){ if (pre_cut && eval.truth_nuEnergy<=1570 && eval.truth_nuEnergy>1375)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.2050.gt.1570"){ if (pre_cut && eval.truth_nuEnergy<=2050 && eval.truth_nuEnergy>1570)  { return number; } }
      else if (cut_name == "numuCC.inside.Enu.le.4000.gt.2050"){ if (pre_cut && eval.truth_nuEnergy>2050  && eval.truth_nuEnergy<=4000) { return number; } }
      else{ std::cout << "get_xs_signal_no: no cut found!" << std::endl; }
    }
  }

  return -1;
}

bool LEEana::get_cut_pass(TString ch_name, TString add_cut, bool flag_data, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine, SpaceInfo& space, PandoraInfo& pandora, LanternInfo& lantern){

  double reco_Enu = get_reco_Enu_corr(kine, flag_data);

  double KE_muon = pfeval.truth_muonMomentum[3]*1000.-105.66; // MeV
  double Pmuon   = (TMath::Sqrt(pow(KE_muon,2) + 2*KE_muon*105.66));

  double Emuon = pfeval.truth_muonMomentum[3]*1000; // MeV
  double Ehadron = eval.truth_nuEnergy - pfeval.truth_muonMomentum[3]*1000.; // MeV

  TLorentzVector truth_muonMomentum(pfeval.truth_muonMomentum[0], pfeval.truth_muonMomentum[1], pfeval.truth_muonMomentum[2], pfeval.truth_muonMomentum[3]);

  bool flag_truth_inside = false; // in the active volume
  if (eval.truth_vtxX > -1 && eval.truth_vtxX <= 254.3 &&  eval.truth_vtxY >-115.0 && eval.truth_vtxY<=117.0 && eval.truth_vtxZ > 0.6 && eval.truth_vtxZ <=1036.4) flag_truth_inside = true;

  // definition of additional cuts
  std::map<std::string, bool> map_cuts_flag;
  if(is_far_sideband(kine, tagger, flag_data)) map_cuts_flag["farsideband"] = true;
  else map_cuts_flag["farsideband"] = false;

  if(is_near_sideband(kine, tagger, flag_data)) map_cuts_flag["nearsideband"] = true;
  else map_cuts_flag["nearsideband"] = false;

  if(is_nueCC(tagger)) map_cuts_flag["nueCC"] = true;
  else map_cuts_flag["nueCC"] = false;

  if(is_loosenueCC(tagger)) map_cuts_flag["loosenueCC"] = true;
  else map_cuts_flag["loosenueCC"] = false;

  if(is_generic(eval)) map_cuts_flag["generic"] = true;
  else map_cuts_flag["generic"] = false;

  if(eval.truth_nuEnergy <=400) map_cuts_flag["LowEintnueCC"] = true;
  else map_cuts_flag["LowEintnueCC"] = false;

  if (!(eval.truth_nuEnergy <=400)) map_cuts_flag["antiLowEintnueCC"] = true;
  else map_cuts_flag["antiLowEintnueCC"] = false;

  if(eval.truth_nuEnergy<=400) map_cuts_flag["LowEnu"] = true;
  else map_cuts_flag["LowEnu"] = false;

  if(!(eval.truth_nuEnergy<=400)) map_cuts_flag["antiLowEnu"] = true;
  else map_cuts_flag["antiLowEnu"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["badmatch"] = true;
  else map_cuts_flag["badmatch"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["numuCCinFV"] = true;
  else map_cuts_flag["numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["RnumuCCinFV"] = true;
  else map_cuts_flag["RnumuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==-14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["AnumuCCinFV"] = true;
  else map_cuts_flag["AnumuCCinFV"] = false;

  // Xs related cuts ...

  map_cuts_flag["XsnumuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1;

  map_cuts_flag["Xs_Enu_numuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && truth_muonMomentum[3]>0 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200 && Pmuon > 0 && Pmuon <= 2500;

  map_cuts_flag["Xs_Enu_mu_numuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && truth_muonMomentum[3]>0 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200 && Pmuon > 0 && Pmuon <= 2500;

  map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] = eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && truth_muonMomentum[3]>0 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200 && Pmuon > 0 && Pmuon <= 2500;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Emuon > 105.7 && Emuon<=2506) map_cuts_flag["Xs_Emu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Emu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Pmuon > 0 && Pmuon<=2500) map_cuts_flag["Xs_Pmu_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Pmu_numuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && Ehadron > 30 && Ehadron <=2500) map_cuts_flag["Xs_Ehad_numuCCinFV"] = true;
  else map_cuts_flag["Xs_Ehad_numuCCinFV"] = false;

  // xs breakdown mode
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200) map_cuts_flag["XsecNumuCCinFV"] = true;
  else map_cuts_flag["XsecNumuCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0) map_cuts_flag["XsecNC"] = true;
  else map_cuts_flag["XsecNC"] = false;

  if(eval.match_completeness_energy<=0.1*eval.truth_energyInside) map_cuts_flag["XsecCosmic"] = true;
  else map_cuts_flag["XsecCosmic"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==1 && !(eval.truth_nuPdg==14 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy>200)) map_cuts_flag["XsecBkgCC"] = true;
  else map_cuts_flag["XsecBkgCC"] = false;

  // finish Xs related cuts ...
  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && eval.truth_nuEnergy<=4000 && eval.truth_nuEnergy > 200) map_cuts_flag["Xs_Enu_nueCCinFV"] = true;
  else map_cuts_flag["Xs_Enu_nueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["nueCCinFV"] = true;
  else map_cuts_flag["nueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["RnueCCinFV"] = true;
  else map_cuts_flag["RnueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_nuPdg==-12 && eval.truth_isCC==1 && eval.truth_vtxInside==1) map_cuts_flag["AnueCCinFV"] = true;
  else map_cuts_flag["AnueCCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==0) map_cuts_flag["NCinFV"] = true;
  else map_cuts_flag["NCinFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==0) map_cuts_flag["outFV"] = true;
  else map_cuts_flag["outFV"] = false;

  if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==1) map_cuts_flag["CCpi0inFV"] = true;
  else map_cuts_flag["CCpi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_NprimPio==1) map_cuts_flag["NCpi0inFV"] = true;
  else map_cuts_flag["NCpi0inFV"] = false;

  // breakdown categories for NC Delta analysis
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1) map_cuts_flag["NCDeltainFV"] = true;
  else map_cuts_flag["NCDeltainFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && pfeval.truth_NCDelta==0) map_cuts_flag["NC1Pi0inFV"] = true;
  else map_cuts_flag["NC1Pi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_NprimPio==1 && pfeval.truth_NCDelta==0) map_cuts_flag["numuCC1Pi0inFV"] = true;
  else map_cuts_flag["numuCC1Pi0inFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_NprimPio!=1 && pfeval.truth_NCDelta==0) map_cuts_flag["numuCCotherinFV"] = true;
  else map_cuts_flag["numuCCotherinFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==1 && eval.truth_isCC==0 && pfeval.truth_NprimPio!=1 && pfeval.truth_NCDelta==0) map_cuts_flag["NCotherinFV"] = true;
  else map_cuts_flag["NCotherinFV"] = false;
  // done with NC Delta breakdown categories

  //Erin
  // breakdown categories for single photon analysis
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && eval.truth_vtxInside==1) map_cuts_flag["SPNCDeltaSig"] = true;
  else map_cuts_flag["SPNCDeltaSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && (pfeval.truth_showerMother==111) && eval.truth_vtxInside==1) map_cuts_flag["SPNCPi0Sig"] = true;
  else map_cuts_flag["SPNCPi0Sig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_showerMother!=111 && pfeval.truth_NCDelta==0 && eval.truth_vtxInside==1) map_cuts_flag["SPNCOtherSig"] = true;
  else map_cuts_flag["SPNCOtherSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_muonMomentum[3]-0.105658<0.1 && eval.truth_vtxInside==1) map_cuts_flag["SPNumuCCSig"] = true;
  else map_cuts_flag["SPNumuCCSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && (eval.truth_isCC==0 || (eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_muonMomentum[3]-0.105658<0.1)) && eval.truth_vtxInside==0) map_cuts_flag["SPOutFVSig"] = true;
  else map_cuts_flag["SPOutFVSig"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside &&  pfeval.reco_muonMomentum[3] > 0) map_cuts_flag["muon"] = true;
  else map_cuts_flag["muon"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside &&  !(pfeval.reco_muonMomentum[3] > 0)) map_cuts_flag["nomuon"] = true;
  else map_cuts_flag["nomuon"] = false;

  //for testing, no FV check
  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_NCDelta==1) map_cuts_flag["SPNCDeltaSigNoFV"] = true;
  else map_cuts_flag["SPNCDeltaSigNoFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && (pfeval.truth_showerMother==111)) map_cuts_flag["SPNCPi0SigNoFV"] = true;
  else map_cuts_flag["SPNCPi0SigNoFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==0 && pfeval.truth_showerMother!=111 && pfeval.truth_NCDelta==0) map_cuts_flag["SPNCOtherSigNoFV"] = true;
  else map_cuts_flag["SPNCOtherSigNoFV"] = false;

  if (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1 && eval.truth_isCC==1 && abs(eval.truth_nuPdg)==14 && pfeval.truth_muonMomentum[3]-0.105658<0.1) map_cuts_flag["SPNumuCCSigNoFV"] = true;
  else map_cuts_flag["SPNumuCCSigNoFV"] = false;
  //

  map_cuts_flag["SPdirtBkg"] = false;
  map_cuts_flag["SPoutFVBkg"] = false;
  map_cuts_flag["SPnumuCCBkg"] = false;
  map_cuts_flag["SPnumuCCpi0Bkg"] = false;
  map_cuts_flag["SPNCBkg"] = false;
  map_cuts_flag["SPNCpi0Bkg"] = false;

  if(!(map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] || map_cuts_flag["SPNumuCCSig"] || map_cuts_flag["SPOutFVSig"])){
      map_cuts_flag["SPdirtBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_vtxInside==0) map_cuts_flag["SPoutFVBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_Npi0==0) map_cuts_flag["SPnumuCCBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && abs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1 && pfeval.truth_Npi0>0) map_cuts_flag["SPnumuCCpi0Bkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_Npi0==0) map_cuts_flag["SPNCBkg"] = true;
      if(eval.match_completeness_energy>0.1*eval.truth_energyInside && eval.truth_isCC==0 && eval.truth_vtxInside==1 && pfeval.truth_Npi0>0) map_cuts_flag["SPNCpi0Bkg"] = true;
  }

  if((map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] || map_cuts_flag["SPNumuCCSig"] || map_cuts_flag["SPOutFVSig"])){
      if (is_true_0p(pfeval)==1) {
        map_cuts_flag["SP0p"] = true;
        map_cuts_flag["SPNp"] = false;
      }
      else {
        map_cuts_flag["SPNp"] = true;
        map_cuts_flag["SP0p"] = false;
      }
  }else{
      map_cuts_flag["SPNp"] = false;
      map_cuts_flag["SP0p"] = false;
  }
  // done with single photon breakdown categories
  //


  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCMEC"] = true;
  else map_cuts_flag["CCMEC"] = false;

  if(pfeval.truth_nuScatType==10 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCMEC"] = true;
  else map_cuts_flag["NCMEC"] = false;

  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCQE"] = true;
  else map_cuts_flag["CCQE"] = false;

  if(pfeval.truth_nuScatType==1 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCQE"] = true;
  else map_cuts_flag["NCQE"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCRES"] = true;
  else map_cuts_flag["CCRES"] = false;

  if(pfeval.truth_nuScatType==4 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCRES"] = true;
  else map_cuts_flag["NCRES"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==1 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["CCDIS"] = true;
  else map_cuts_flag["CCDIS"] = false;

  if(pfeval.truth_nuScatType==3 && eval.truth_isCC==0 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["NCDIS"] = true;
  else map_cuts_flag["NCDIS"] = false;

  if(pfeval.truth_nuScatType!=10 && pfeval.truth_nuScatType!=1 && pfeval.truth_nuScatType!=3 && pfeval.truth_nuScatType!=4 && eval.match_completeness_energy>0.1*eval.truth_energyInside) map_cuts_flag["OTHER"] = true;
  else map_cuts_flag["OTHER"] = false;

  map_cuts_flag["none"] = false;
  map_cuts_flag["LEE"] = true;

  // figure out additional cuts and flag_data ...
  bool flag_add = true;
  if(add_cut == "all") flag_add = true;
  else if( (flag_data && (add_cut=="none" || add_cut=="farsideband" || add_cut=="nearsideband" || add_cut=="nueCC" || add_cut=="generic" || add_cut=="loosenueCC")) || !flag_data ){
      std::istringstream sss(add_cut.Data());
      for(std::string line; std::getline(sss, line, '_');){
          if(map_cuts_flag.find(line)!=map_cuts_flag.end()){
              flag_add *= map_cuts_flag[line];
          }
          else{
              std::cout<<"ERROR: add_cut "<<line<<" not defined!\n";
              exit(EXIT_FAILURE);
          }
      }
  }
  else{
    std::cout<<"ERROR: add_cut "<<add_cut<<" of channel "<< ch_name <<" is not assigned to sample "<<flag_data<<" [1: data; 0: mc]\n";
    std::cout<<"Please modify inc/WCPLEEANA/cuts.h\n";
    exit(EXIT_FAILURE);
  }

  if (!flag_add) return false;

  bool flag_generic = is_generic(eval);
  bool flag_numuCC = is_numuCC(tagger);
  //bool flag_numuCC = is_numuCC(tagger) && (is_far_sideband(kine, tagger, flag_data) || is_near_sideband(kine, tagger, flag_data) );
  bool flag_numuCC_tight = is_numuCC_tight(tagger, pfeval);
  bool flag_numuCC_1mu0p = is_numuCC_1mu0p(tagger, kine, pfeval);
  bool flag_numuCC_lowEhad = is_numuCC_lowEhad(tagger, kine, pfeval, flag_data);
  bool flag_numuCC_cutbased = is_numuCC_cutbased(tagger);
  bool flag_nueCC = is_nueCC(tagger);
  bool flag_nueCC_loose = is_loosenueCC(tagger);

  bool flag_0p = is_0p(tagger, kine, pfeval);
  bool flag_1p = is_1p(tagger, kine, pfeval);
  bool flag_0pi = is_0pi(tagger, kine, pfeval);

  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);
  bool flag_FC = is_FC(eval);

  //bool flag_ncpio_sel = is_NCpio_bdt(tagger) && (!flag_0p);
  bool flag_ncpio_sel = is_NCpio_sel(tagger, kine);
  bool flag_ncdelta_sel = is_NCdelta_sel(tagger, pfeval);

  //Erin
  bool flag_singlephoton_sel = is_singlephoton_sel(tagger, pfeval);
  bool flag_singlephoton_eff_sel = is_singlephoton_eff_sel(tagger, pfeval);
  bool flag_singleshower_sel = is_singleshower_sel(tagger, pfeval);
  bool flag_singleshower_eff_sel = is_singleshower_eff_sel(tagger, pfeval);
  bool flag_singlephoton_numu_sel = is_singlephoton_numu_sel(tagger, pfeval);
  bool flag_singlephoton_other_sel = is_singlephoton_other_sel(tagger, pfeval);
  bool flag_singlephoton_ncpi0_sel = is_singlephoton_ncpi0_sel(tagger, pfeval);
  bool flag_singlephoton_nue_sel = is_singlephoton_nue_sel(tagger, pfeval);
  bool flag_singlephoton_nue_sel_allshw = is_singlephoton_nue_sel_allshw(tagger, pfeval);
  bool flag_nsbeam = is_nsbeam_photon(pfeval, eval); //set all cuts to shifted
  bool flag_nsbeam_photon = is_nsbeam_photon(pfeval, eval);
  bool flag_singlephoton_pre = is_singlephoton_pre(tagger, pfeval);
  bool flag_singlephoton_numu = is_singlephoton_numu(tagger, pfeval);
  bool flag_singlephoton_other = is_singlephoton_other(tagger, pfeval);
  bool flag_singlephoton_ncpi0 = is_singlephoton_ncpi0(tagger, pfeval);
  bool flag_singlephoton_nue = is_singlephoton_nue(tagger, pfeval);
  bool flag_singlephoton_eff_numu = is_singlephoton_eff_numu(tagger, pfeval);
  bool flag_singlephoton_eff_other = is_singlephoton_eff_other(tagger, pfeval);
  bool flag_singlephoton_eff_ncpi0 = is_singlephoton_eff_ncpi0(tagger, pfeval);
  bool flag_singlephoton_eff_nue = is_singlephoton_eff_nue(tagger, pfeval);
  bool flag_singlephoton_oneshw = is_singlephoton_oneshw(tagger, pfeval);
  //

  float costheta_binning[10] = {-1, -.5, 0, .27, .45, .62, .76, .86, .94, 1};		// PeLEE binning
  //float costheta_binning[7]  = {-1,         .27,      .62, .76, .86, .94, 1};		// coarse binning
  //float costheta_binning[3]    = {-1,                   .62,                1};	//very coarse binning
  TLorentzVector muonMomentum(pfeval.reco_muonMomentum[0], pfeval.reco_muonMomentum[1], pfeval.reco_muonMomentum[2], pfeval.reco_muonMomentum[3]);
  float reco_pmuon = TMath::Sqrt(pow(pfeval.reco_muonMomentum[0],2)+pow(pfeval.reco_muonMomentum[1],2)+pow(pfeval.reco_muonMomentum[2],2))*1000;

  int costheta_bin = get_costheta_bin(TMath::Cos(muonMomentum.Theta()));
  int Pmu_bin      = get_Pmuon_bin(reco_pmuon);
  int Enu_bin      = get_Enu_bin(reco_Enu);

  if (ch_name == "LEE_FC_nueoverlay"  || ch_name == "nueCC_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "nueCC_FC_nueoverlay_numi"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if ( ch_name == "nueCC_FC_numu2nueoverlay" ){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if ( ch_name == "nueCC_FC_numu2nueoverlay_numi" ){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext" || ch_name == "BG_nueCC_FC_dirt" || ch_name =="nueCC_FC_bnb"){
    //nueCC FC
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_ext_numi" || ch_name == "BG_nueCC_FC_dirt_numi" || ch_name =="nueCC_FC_numi"){
    //nueCC FC
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_FC_overlay_numi"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "LEE_PC_nueoverlay" || ch_name == "nueCC_PC_nueoverlay" ){
    // nueCC PC
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "nueCC_PC_nueoverlay_numi" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if ( ch_name == "nueCC_PC_numu2nueoverlay" ){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if ( ch_name == "nueCC_PC_numu2nueoverlay_numi" ){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext" || ch_name == "BG_nueCC_PC_dirt" || ch_name == "nueCC_PC_bnb"){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_ext_numi" || ch_name == "BG_nueCC_PC_dirt_numi" || ch_name == "nueCC_PC_numi"){
    // nueCC PC
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_PC_overlay_numi"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta0_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta0_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta0_FC_overlay"     || ch_name == "numuCC_signal_nu_theta0_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta0_FC_overlay" || ch_name == "numuCC_background_Emu_theta0_FC_overlay" || ch_name == "numuCC_background_Pmu_theta0_FC_overlay" || ch_name == "numuCC_background_nu_theta0_FC_overlay"
         || ch_name == "BG_numuCC_theta0_FC_ext"                 || ch_name =="BG_numuCC_theta0_FC_dirt"                 || ch_name == "numuCC_theta0_FC_bnb"
         || ch_name == "BG_numuCC_theta0_FC_ext_v2"              || ch_name =="BG_numuCC_theta0_FC_dirt_v2"              || ch_name == "numuCC_theta0_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) {
      if      (ch_name == "numuCC_signal_Enu_theta0_FC_overlay" || ch_name == "numuCC_background_Enu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta0_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta0_FC_overlay" || ch_name == "numuCC_background_Emu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta0_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta0_FC_overlay" || ch_name == "numuCC_background_Pmu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta0_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta0_FC_overlay" || ch_name ==  "numuCC_background_nu_theta0_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta0_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta1_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta1_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta1_FC_overlay"     || ch_name == "numuCC_signal_nu_theta1_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta1_FC_overlay" || ch_name == "numuCC_background_Emu_theta1_FC_overlay" || ch_name == "numuCC_background_Pmu_theta1_FC_overlay" || ch_name == "numuCC_background_nu_theta1_FC_overlay"
         || ch_name == "BG_numuCC_theta1_FC_ext"                 || ch_name =="BG_numuCC_theta1_FC_dirt"                 || ch_name == "numuCC_theta1_FC_bnb"
         || ch_name == "BG_numuCC_theta1_FC_ext_v2"              || ch_name =="BG_numuCC_theta1_FC_dirt_v2"              || ch_name == "numuCC_theta1_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) {
      if      (ch_name == "numuCC_signal_Enu_theta1_FC_overlay" || ch_name == "numuCC_background_Enu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta1_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta1_FC_overlay" || ch_name == "numuCC_background_Emu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta1_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta1_FC_overlay" || ch_name == "numuCC_background_Pmu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta1_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta1_FC_overlay" || ch_name ==  "numuCC_background_nu_theta1_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta1_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta2_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta2_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta2_FC_overlay"     || ch_name == "numuCC_signal_nu_theta2_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta2_FC_overlay" || ch_name == "numuCC_background_Emu_theta2_FC_overlay" || ch_name == "numuCC_background_Pmu_theta2_FC_overlay" || ch_name == "numuCC_background_nu_theta2_FC_overlay"
         || ch_name == "BG_numuCC_theta2_FC_ext"                 || ch_name =="BG_numuCC_theta2_FC_dirt"                 || ch_name == "numuCC_theta2_FC_bnb"
         || ch_name == "BG_numuCC_theta2_FC_ext_v2"              || ch_name =="BG_numuCC_theta2_FC_dirt_v2"              || ch_name == "numuCC_theta2_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) {
      if      (ch_name == "numuCC_signal_Enu_theta2_FC_overlay" || ch_name == "numuCC_background_Enu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta2_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta2_FC_overlay" || ch_name == "numuCC_background_Emu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta2_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta2_FC_overlay" || ch_name == "numuCC_background_Pmu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta2_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta2_FC_overlay" || ch_name ==  "numuCC_background_nu_theta2_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta2_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta3_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta3_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta3_FC_overlay"     || ch_name == "numuCC_signal_nu_theta3_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta3_FC_overlay" || ch_name == "numuCC_background_Emu_theta3_FC_overlay" || ch_name == "numuCC_background_Pmu_theta3_FC_overlay" || ch_name == "numuCC_background_nu_theta3_FC_overlay"
         || ch_name == "BG_numuCC_theta3_FC_ext"                 || ch_name =="BG_numuCC_theta3_FC_dirt"                 || ch_name == "numuCC_theta3_FC_bnb"
         || ch_name == "BG_numuCC_theta3_FC_ext_v2"              || ch_name =="BG_numuCC_theta3_FC_dirt_v2"              || ch_name == "numuCC_theta3_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) {
      if      (ch_name == "numuCC_signal_Enu_theta3_FC_overlay" || ch_name == "numuCC_background_Enu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta3_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta3_FC_overlay" || ch_name == "numuCC_background_Emu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta3_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta3_FC_overlay" || ch_name == "numuCC_background_Pmu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta3_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta3_FC_overlay" || ch_name ==  "numuCC_background_nu_theta3_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta3_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta4_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta4_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta4_FC_overlay"     || ch_name == "numuCC_signal_nu_theta4_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta4_FC_overlay" || ch_name == "numuCC_background_Emu_theta4_FC_overlay" || ch_name == "numuCC_background_Pmu_theta4_FC_overlay" || ch_name == "numuCC_background_nu_theta4_FC_overlay"
         || ch_name == "BG_numuCC_theta4_FC_ext"                 || ch_name =="BG_numuCC_theta4_FC_dirt"                 || ch_name == "numuCC_theta4_FC_bnb"
         || ch_name == "BG_numuCC_theta4_FC_ext_v2"              || ch_name =="BG_numuCC_theta4_FC_dirt_v2"              || ch_name == "numuCC_theta4_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) {
      if      (ch_name == "numuCC_signal_Enu_theta4_FC_overlay" || ch_name == "numuCC_background_Enu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta4_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta4_FC_overlay" || ch_name == "numuCC_background_Emu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta4_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta4_FC_overlay" || ch_name == "numuCC_background_Pmu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta4_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta4_FC_overlay" || ch_name ==  "numuCC_background_nu_theta4_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta4_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta5_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta5_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta5_FC_overlay"     || ch_name == "numuCC_signal_nu_theta5_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta5_FC_overlay" || ch_name == "numuCC_background_Emu_theta5_FC_overlay" || ch_name == "numuCC_background_Pmu_theta5_FC_overlay" || ch_name == "numuCC_background_nu_theta5_FC_overlay"
         || ch_name == "BG_numuCC_theta5_FC_ext"                 || ch_name =="BG_numuCC_theta5_FC_dirt"                 || ch_name == "numuCC_theta5_FC_bnb"
         || ch_name == "BG_numuCC_theta5_FC_ext_v2"              || ch_name =="BG_numuCC_theta5_FC_dirt_v2"              || ch_name == "numuCC_theta5_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) {
      if      (ch_name == "numuCC_signal_Enu_theta5_FC_overlay" || ch_name == "numuCC_background_Enu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta5_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta5_FC_overlay" || ch_name == "numuCC_background_Emu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta5_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta5_FC_overlay" || ch_name == "numuCC_background_Pmu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta5_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta5_FC_overlay" || ch_name ==  "numuCC_background_nu_theta5_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta5_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta6_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta6_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta6_FC_overlay"     || ch_name == "numuCC_signal_nu_theta6_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta6_FC_overlay" || ch_name == "numuCC_background_Emu_theta6_FC_overlay" || ch_name == "numuCC_background_Pmu_theta6_FC_overlay" || ch_name == "numuCC_background_nu_theta6_FC_overlay"
         || ch_name == "BG_numuCC_theta6_FC_ext"                 || ch_name =="BG_numuCC_theta6_FC_dirt"                 || ch_name == "numuCC_theta6_FC_bnb"
         || ch_name == "BG_numuCC_theta6_FC_ext_v2"              || ch_name =="BG_numuCC_theta6_FC_dirt_v2"              || ch_name == "numuCC_theta6_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) {
      if      (ch_name == "numuCC_signal_Enu_theta6_FC_overlay" || ch_name == "numuCC_background_Enu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta6_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta6_FC_overlay" || ch_name == "numuCC_background_Emu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta6_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta6_FC_overlay" || ch_name == "numuCC_background_Pmu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta6_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta6_FC_overlay" || ch_name ==  "numuCC_background_nu_theta6_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta6_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta7_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta7_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta7_FC_overlay"     || ch_name == "numuCC_signal_nu_theta7_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta7_FC_overlay" || ch_name == "numuCC_background_Emu_theta7_FC_overlay" || ch_name == "numuCC_background_Pmu_theta7_FC_overlay" || ch_name == "numuCC_background_nu_theta7_FC_overlay"
         || ch_name == "BG_numuCC_theta7_FC_ext"                 || ch_name =="BG_numuCC_theta7_FC_dirt"                 || ch_name == "numuCC_theta7_FC_bnb"
         || ch_name == "BG_numuCC_theta7_FC_ext_v2"              || ch_name =="BG_numuCC_theta7_FC_dirt_v2"              || ch_name == "numuCC_theta7_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) {
      if      (ch_name == "numuCC_signal_Enu_theta7_FC_overlay" || ch_name == "numuCC_background_Enu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta7_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta7_FC_overlay" || ch_name == "numuCC_background_Emu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta7_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta7_FC_overlay" || ch_name == "numuCC_background_Pmu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta7_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta7_FC_overlay" || ch_name ==  "numuCC_background_nu_theta7_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta7_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta8_FC_overlay"     || ch_name == "numuCC_signal_Emu_theta8_FC_overlay"     || ch_name == "numuCC_signal_Pmu_theta8_FC_overlay"     || ch_name == "numuCC_signal_nu_theta8_FC_overlay"
         || ch_name == "numuCC_background_Enu_theta8_FC_overlay" || ch_name == "numuCC_background_Emu_theta8_FC_overlay" || ch_name == "numuCC_background_Pmu_theta8_FC_overlay" || ch_name == "numuCC_background_nu_theta8_FC_overlay"
         || ch_name == "BG_numuCC_theta8_FC_ext"                 || ch_name =="BG_numuCC_theta8_FC_dirt"                 || ch_name == "numuCC_theta8_FC_bnb"
         || ch_name == "BG_numuCC_theta8_FC_ext_v2"              || ch_name =="BG_numuCC_theta8_FC_dirt_v2"              || ch_name == "numuCC_theta8_FC_bnb_v2" ){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) {
      if      (ch_name == "numuCC_signal_Enu_theta8_FC_overlay" || ch_name == "numuCC_background_Enu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta8_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta8_FC_overlay" || ch_name == "numuCC_background_Emu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta8_FC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta8_FC_overlay" || ch_name == "numuCC_background_Pmu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta8_FC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta8_FC_overlay" || ch_name ==  "numuCC_background_nu_theta8_FC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta8_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu0_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta0_Pmu_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu0_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu0_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu0_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu0_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu0_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu0_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu0_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu0_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu0_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu0_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu1_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta0_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu1_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu1_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu1_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu1_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu1_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu1_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu1_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu1_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu1_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu1_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu2_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta0_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu2_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu2_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu2_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu2_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu2_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu2_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu2_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu2_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu2_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu2_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta0_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta0_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta0_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu3_theta0_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta0_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta1_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta1_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta1_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu3_theta1_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta1_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta2_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta2_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta2_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu3_theta2_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta2_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta3_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta3_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta3_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu3_theta3_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta3_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta4_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta4_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta4_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu3_theta4_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta4_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta5_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta5_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta5_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu3_theta5_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta5_Pmu_FC_overlay")); }
      else return true;
    } else return false;

  }else if (ch_name == "numuCC_signal_Enu3_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta6_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta6_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta6_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu3_theta6_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta6_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta7_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta7_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta7_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu3_theta7_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta7_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_FC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta8_Pmu_FC_ext"         || ch_name =="BG_numuCC_Enu3_theta8_Pmu_FC_dirt" || ch_name == "numuCC_Enu3_theta8_Pmu_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu3_theta8_Pmu_FC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_FC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta8_Pmu_FC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Pmu0_FC_overlay" || ch_name == "numuCC_background_Pmu0_FC_overlay"
         || ch_name == "BG_numuCC_Pmu0_FC_ext"         || ch_name =="BG_numuCC_Pmu0_FC_dirt" || ch_name == "numuCC_Pmu0_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==0) {
      if      (ch_name == "numuCC_signal_Pmu0_FC_overlay" || ch_name == "numuCC_background_Pmu0_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu0_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu1_FC_overlay" || ch_name == "numuCC_background_Pmu1_FC_overlay"
         || ch_name == "BG_numuCC_Pmu1_FC_ext"         || ch_name =="BG_numuCC_Pmu1_FC_dirt" || ch_name == "numuCC_Pmu1_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==1) {
      if      (ch_name == "numuCC_signal_Pmu1_FC_overlay" || ch_name == "numuCC_background_Pmu1_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu1_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu2_FC_overlay" || ch_name == "numuCC_background_Pmu2_FC_overlay"
         || ch_name == "BG_numuCC_Pmu2_FC_ext"         || ch_name =="BG_numuCC_Pmu2_FC_dirt" || ch_name == "numuCC_Pmu2_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==2) {
      if      (ch_name == "numuCC_signal_Pmu2_FC_overlay" || ch_name == "numuCC_background_Pmu2_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu2_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu3_FC_overlay" || ch_name == "numuCC_background_Pmu3_FC_overlay"
         || ch_name == "BG_numuCC_Pmu3_FC_ext"         || ch_name =="BG_numuCC_Pmu3_FC_dirt" || ch_name == "numuCC_Pmu3_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==3) {
      if      (ch_name == "numuCC_signal_Pmu3_FC_overlay" || ch_name == "numuCC_background_Pmu3_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu3_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu4_FC_overlay" || ch_name == "numuCC_background_Pmu4_FC_overlay"
         || ch_name == "BG_numuCC_Pmu4_FC_ext"         || ch_name =="BG_numuCC_Pmu4_FC_dirt" || ch_name == "numuCC_Pmu4_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==4) {
      if      (ch_name == "numuCC_signal_Pmu4_FC_overlay" || ch_name == "numuCC_background_Pmu4_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu4_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu5_FC_overlay" || ch_name == "numuCC_background_Pmu5_FC_overlay"
         || ch_name == "BG_numuCC_Pmu5_FC_ext"         || ch_name =="BG_numuCC_Pmu5_FC_dirt" || ch_name == "numuCC_Pmu5_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==5) {
      if      (ch_name == "numuCC_signal_Pmu5_FC_overlay" || ch_name == "numuCC_background_Pmu5_FC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu5_FC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_theta_all_0p_FC_overlay" || ch_name == "BG_numuCC_theta_all_0p_FC_ext" || ch_name =="BG_numuCC_theta_all_0p_FC_dirt" || ch_name == "numuCC_theta_all_0p_FC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_Np_FC_overlay" || ch_name == "BG_numuCC_theta_all_Np_FC_ext" || ch_name =="BG_numuCC_theta_all_Np_FC_dirt" || ch_name == "numuCC_theta_all_Np_FC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all_FC_overlay" || ch_name == "BG_numuCC_all_FC_ext" || ch_name =="BG_numuCC_all_FC_dirt" || ch_name == "numuCC_all_FC_bnb") {
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all2_FC_overlay" || ch_name == "BG_numuCC_all2_FC_ext" || ch_name =="BG_numuCC_all2_FC_dirt" || ch_name == "numuCC_all2_FC_bnb") {
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_0p_FC_overlay" || ch_name == "BG_numuCC_0p_FC_ext" || ch_name =="BG_numuCC_0p_FC_dirt" || ch_name == "numuCC_0p_FC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_Np_FC_overlay" || ch_name == "BG_numuCC_Np_FC_ext" || ch_name =="BG_numuCC_Np_FC_dirt" || ch_name == "numuCC_Np_FC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta0_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta0_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta0_PC_overlay"     || ch_name == "numuCC_signal_nu_theta0_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta0_PC_overlay" || ch_name == "numuCC_background_Emu_theta0_PC_overlay" || ch_name == "numuCC_background_Pmu_theta0_PC_overlay" || ch_name == "numuCC_background_nu_theta0_PC_overlay"
         || ch_name == "BG_numuCC_theta0_PC_ext"                 || ch_name =="BG_numuCC_theta0_PC_dirt"                 || ch_name == "numuCC_theta0_PC_bnb"
         || ch_name == "BG_numuCC_theta0_PC_ext_v2"              || ch_name =="BG_numuCC_theta0_PC_dirt_v2"              || ch_name == "numuCC_theta0_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[0] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[1])) {
      if      (ch_name == "numuCC_signal_Enu_theta0_PC_overlay" || ch_name == "numuCC_background_Enu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta0_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta0_PC_overlay" || ch_name == "numuCC_background_Emu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta0_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta0_PC_overlay" || ch_name == "numuCC_background_Pmu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta0_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta0_PC_overlay" || ch_name ==  "numuCC_background_nu_theta0_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta0_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta1_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta1_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta1_PC_overlay"     || ch_name == "numuCC_signal_nu_theta1_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta1_PC_overlay" || ch_name == "numuCC_background_Emu_theta1_PC_overlay" || ch_name == "numuCC_background_Pmu_theta1_PC_overlay" || ch_name == "numuCC_background_nu_theta1_PC_overlay"
         || ch_name == "BG_numuCC_theta1_PC_ext"                 || ch_name =="BG_numuCC_theta1_PC_dirt"                 || ch_name == "numuCC_theta1_PC_bnb"
         || ch_name == "BG_numuCC_theta1_PC_ext_v2"              || ch_name =="BG_numuCC_theta1_PC_dirt_v2"              || ch_name == "numuCC_theta1_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[1] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[2])) {
      if      (ch_name == "numuCC_signal_Enu_theta1_PC_overlay" || ch_name == "numuCC_background_Enu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta1_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta1_PC_overlay" || ch_name == "numuCC_background_Emu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta1_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta1_PC_overlay" || ch_name == "numuCC_background_Pmu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta1_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta1_PC_overlay" || ch_name ==  "numuCC_background_nu_theta1_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta1_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta2_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta2_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta2_PC_overlay"     || ch_name == "numuCC_signal_nu_theta2_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta2_PC_overlay" || ch_name == "numuCC_background_Emu_theta2_PC_overlay" || ch_name == "numuCC_background_Pmu_theta2_PC_overlay" || ch_name == "numuCC_background_nu_theta2_PC_overlay"
         || ch_name == "BG_numuCC_theta2_PC_ext"                 || ch_name =="BG_numuCC_theta2_PC_dirt"                 || ch_name == "numuCC_theta2_PC_bnb"
         || ch_name == "BG_numuCC_theta2_PC_ext_v2"              || ch_name =="BG_numuCC_theta2_PC_dirt_v2"              || ch_name == "numuCC_theta2_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[2] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[3])) {
      if      (ch_name == "numuCC_signal_Enu_theta2_PC_overlay" || ch_name == "numuCC_background_Enu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta2_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta2_PC_overlay" || ch_name == "numuCC_background_Emu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta2_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta2_PC_overlay" || ch_name == "numuCC_background_Pmu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta2_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta2_PC_overlay" || ch_name ==  "numuCC_background_nu_theta2_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta2_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta3_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta3_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta3_PC_overlay"     || ch_name == "numuCC_signal_nu_theta3_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta3_PC_overlay" || ch_name == "numuCC_background_Emu_theta3_PC_overlay" || ch_name == "numuCC_background_Pmu_theta3_PC_overlay" || ch_name == "numuCC_background_nu_theta3_PC_overlay"
         || ch_name == "BG_numuCC_theta3_PC_ext"                 || ch_name =="BG_numuCC_theta3_PC_dirt"                 || ch_name == "numuCC_theta3_PC_bnb"
         || ch_name == "BG_numuCC_theta3_PC_ext_v2"              || ch_name =="BG_numuCC_theta3_PC_dirt_v2"              || ch_name == "numuCC_theta3_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[3] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[4])) {
      if      (ch_name == "numuCC_signal_Enu_theta3_PC_overlay" || ch_name == "numuCC_background_Enu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta3_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta3_PC_overlay" || ch_name == "numuCC_background_Emu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta3_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta3_PC_overlay" || ch_name == "numuCC_background_Pmu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta3_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta3_PC_overlay" || ch_name ==  "numuCC_background_nu_theta3_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta3_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta4_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta4_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta4_PC_overlay"     || ch_name == "numuCC_signal_nu_theta4_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta4_PC_overlay" || ch_name == "numuCC_background_Emu_theta4_PC_overlay" || ch_name == "numuCC_background_Pmu_theta4_PC_overlay" || ch_name == "numuCC_background_nu_theta4_PC_overlay"
         || ch_name == "BG_numuCC_theta4_PC_ext"                 || ch_name =="BG_numuCC_theta4_PC_dirt"                 || ch_name == "numuCC_theta4_PC_bnb"
         || ch_name == "BG_numuCC_theta4_PC_ext_v2"              || ch_name =="BG_numuCC_theta4_PC_dirt_v2"              || ch_name == "numuCC_theta4_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[4] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[5])) {
      if      (ch_name == "numuCC_signal_Enu_theta4_PC_overlay" || ch_name == "numuCC_background_Enu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta4_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta4_PC_overlay" || ch_name == "numuCC_background_Emu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta4_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta4_PC_overlay" || ch_name == "numuCC_background_Pmu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta4_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta4_PC_overlay" || ch_name ==  "numuCC_background_nu_theta4_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta4_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta5_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta5_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta5_PC_overlay"     || ch_name == "numuCC_signal_nu_theta5_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta5_PC_overlay" || ch_name == "numuCC_background_Emu_theta5_PC_overlay" || ch_name == "numuCC_background_Pmu_theta5_PC_overlay" || ch_name == "numuCC_background_nu_theta5_PC_overlay"
         || ch_name == "BG_numuCC_theta5_PC_ext"                 || ch_name =="BG_numuCC_theta5_PC_dirt"                 || ch_name == "numuCC_theta5_PC_bnb"
         || ch_name == "BG_numuCC_theta5_PC_ext_v2"              || ch_name =="BG_numuCC_theta5_PC_dirt_v2"              || ch_name == "numuCC_theta5_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[5] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[6])) {
      if      (ch_name == "numuCC_signal_Enu_theta5_PC_overlay" || ch_name == "numuCC_background_Enu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta5_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta5_PC_overlay" || ch_name == "numuCC_background_Emu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta5_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta5_PC_overlay" || ch_name == "numuCC_background_Pmu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta5_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta5_PC_overlay" || ch_name ==  "numuCC_background_nu_theta5_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta5_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta6_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta6_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta6_PC_overlay"     || ch_name == "numuCC_signal_nu_theta6_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta6_PC_overlay" || ch_name == "numuCC_background_Emu_theta6_PC_overlay" || ch_name == "numuCC_background_Pmu_theta6_PC_overlay" || ch_name == "numuCC_background_nu_theta6_PC_overlay"
         || ch_name == "BG_numuCC_theta6_PC_ext"                 || ch_name =="BG_numuCC_theta6_PC_dirt"                 || ch_name == "numuCC_theta6_PC_bnb"
         || ch_name == "BG_numuCC_theta6_PC_ext_v2"              || ch_name =="BG_numuCC_theta6_PC_dirt_v2"              || ch_name == "numuCC_theta6_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[6] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[7])) {
      if      (ch_name == "numuCC_signal_Enu_theta6_PC_overlay" || ch_name == "numuCC_background_Enu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta6_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta6_PC_overlay" || ch_name == "numuCC_background_Emu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta6_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta6_PC_overlay" || ch_name == "numuCC_background_Pmu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta6_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta6_PC_overlay" || ch_name ==  "numuCC_background_nu_theta6_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta6_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta7_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta7_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta7_PC_overlay"     || ch_name == "numuCC_signal_nu_theta7_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta7_PC_overlay" || ch_name == "numuCC_background_Emu_theta7_PC_overlay" || ch_name == "numuCC_background_Pmu_theta7_PC_overlay" || ch_name == "numuCC_background_nu_theta7_PC_overlay"
         || ch_name == "BG_numuCC_theta7_PC_ext"                 || ch_name =="BG_numuCC_theta7_PC_dirt"                 || ch_name == "numuCC_theta7_PC_bnb"
         || ch_name == "BG_numuCC_theta7_PC_ext_v2"              || ch_name =="BG_numuCC_theta7_PC_dirt_v2"              || ch_name == "numuCC_theta7_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[7] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[8])) {
      if      (ch_name == "numuCC_signal_Enu_theta7_PC_overlay" || ch_name == "numuCC_background_Enu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta7_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta7_PC_overlay" || ch_name == "numuCC_background_Emu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta7_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta7_PC_overlay" || ch_name == "numuCC_background_Pmu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta7_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta7_PC_overlay" || ch_name ==  "numuCC_background_nu_theta7_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta7_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu_theta8_PC_overlay"     || ch_name == "numuCC_signal_Emu_theta8_PC_overlay"     || ch_name == "numuCC_signal_Pmu_theta8_PC_overlay"     || ch_name == "numuCC_signal_nu_theta8_PC_overlay"
         || ch_name == "numuCC_background_Enu_theta8_PC_overlay" || ch_name == "numuCC_background_Emu_theta8_PC_overlay" || ch_name == "numuCC_background_Pmu_theta8_PC_overlay" || ch_name == "numuCC_background_nu_theta8_PC_overlay"
         || ch_name == "BG_numuCC_theta8_PC_ext"                 || ch_name =="BG_numuCC_theta8_PC_dirt"                 || ch_name == "numuCC_theta8_PC_bnb"
         || ch_name == "BG_numuCC_theta8_PC_ext_v2"              || ch_name =="BG_numuCC_theta8_PC_dirt_v2"              || ch_name == "numuCC_theta8_PC_bnb_v2" ){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && (TMath::Cos(muonMomentum.Theta())>=costheta_binning[8] && TMath::Cos(muonMomentum.Theta())<=costheta_binning[9])) {
      if      (ch_name == "numuCC_signal_Enu_theta8_PC_overlay" || ch_name == "numuCC_background_Enu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Enu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu_theta8_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Emu_theta8_PC_overlay" || ch_name == "numuCC_background_Emu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Emu_numuCCinFV"] == (ch_name=="numuCC_signal_Emu_theta8_PC_overlay")); }
      else if (ch_name == "numuCC_signal_Pmu_theta8_PC_overlay" || ch_name == "numuCC_background_Pmu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu_theta8_PC_overlay")); }
      else if (ch_name ==  "numuCC_signal_nu_theta8_PC_overlay" || ch_name ==  "numuCC_background_nu_theta8_PC_overlay") { return (map_cuts_flag["Xs_Ehad_numuCCinFV"]== (ch_name== "numuCC_signal_nu_theta8_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu0_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu0_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu0_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu0_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu0_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu0_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu0_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu0_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu0_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu0_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu0_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu0_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==0 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu0_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu0_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu0_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu1_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu1_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu1_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu1_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu1_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu1_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu1_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu1_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu1_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu1_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu1_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu1_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==1 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu1_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu1_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu1_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu2_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu2_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu2_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu2_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu2_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu2_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu2_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu2_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu2_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu2_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu2_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu2_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==2 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu2_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu2_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu2_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta0_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta0_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta0_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==0) {
      if      (ch_name == "numuCC_signal_Enu3_theta0_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta0_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta0_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta1_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta1_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta1_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==1) {
      if      (ch_name == "numuCC_signal_Enu3_theta1_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta1_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta1_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta2_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta2_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta2_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==2) {
      if      (ch_name == "numuCC_signal_Enu3_theta2_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta2_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta2_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta3_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta3_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta3_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==3) {
      if      (ch_name == "numuCC_signal_Enu3_theta3_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta3_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta3_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta4_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta4_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta4_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==4) {
      if      (ch_name == "numuCC_signal_Enu3_theta4_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta4_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta4_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta5_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta5_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta5_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==5) {
      if      (ch_name == "numuCC_signal_Enu3_theta5_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta5_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta5_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta6_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta6_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta6_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==6) {
      if      (ch_name == "numuCC_signal_Enu3_theta6_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta6_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta6_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta7_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta7_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta7_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==7) {
      if      (ch_name == "numuCC_signal_Enu3_theta7_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta7_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta7_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Enu3_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_PC_overlay"
         || ch_name == "BG_numuCC_Enu3_theta8_Pmu_PC_ext"         || ch_name =="BG_numuCC_Enu3_theta8_Pmu_PC_dirt" || ch_name == "numuCC_Enu3_theta8_Pmu_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Enu_bin==3 && costheta_bin==8) {
      if      (ch_name == "numuCC_signal_Enu3_theta8_Pmu_PC_overlay" || ch_name == "numuCC_background_Enu3_theta8_Pmu_PC_overlay") { return (map_cuts_flag["Xs_Enu_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Enu3_theta8_Pmu_PC_overlay")); }
      else return true;
    } else return false;
  }else if (ch_name == "numuCC_signal_Pmu0_PC_overlay" || ch_name == "numuCC_background_Pmu0_PC_overlay"
         || ch_name == "BG_numuCC_Pmu0_PC_ext"         || ch_name =="BG_numuCC_Pmu0_PC_dirt" || ch_name == "numuCC_Pmu0_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==0) {
      if      (ch_name == "numuCC_signal_Pmu0_PC_overlay" || ch_name == "numuCC_background_Pmu0_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu0_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu1_PC_overlay" || ch_name == "numuCC_background_Pmu1_PC_overlay"
         || ch_name == "BG_numuCC_Pmu1_PC_ext"         || ch_name =="BG_numuCC_Pmu1_PC_dirt" || ch_name == "numuCC_Pmu1_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==1) {
      if      (ch_name == "numuCC_signal_Pmu1_PC_overlay" || ch_name == "numuCC_background_Pmu1_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu1_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu2_PC_overlay" || ch_name == "numuCC_background_Pmu2_PC_overlay"
         || ch_name == "BG_numuCC_Pmu2_PC_ext"         || ch_name =="BG_numuCC_Pmu2_PC_dirt" || ch_name == "numuCC_Pmu2_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==2) {
      if      (ch_name == "numuCC_signal_Pmu2_PC_overlay" || ch_name == "numuCC_background_Pmu2_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu2_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu3_PC_overlay" || ch_name == "numuCC_background_Pmu3_PC_overlay"
         || ch_name == "BG_numuCC_Pmu3_PC_ext"         || ch_name =="BG_numuCC_Pmu3_PC_dirt" || ch_name == "numuCC_Pmu3_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==3) {
      if      (ch_name == "numuCC_signal_Pmu3_PC_overlay" || ch_name == "numuCC_background_Pmu3_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu3_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu4_PC_overlay" || ch_name == "numuCC_background_Pmu4_PC_overlay"
         || ch_name == "BG_numuCC_Pmu4_PC_ext"         || ch_name =="BG_numuCC_Pmu4_PC_dirt" || ch_name == "numuCC_Pmu4_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==4) {
      if      (ch_name == "numuCC_signal_Pmu4_PC_overlay" || ch_name == "numuCC_background_Pmu4_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu4_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_signal_Pmu5_PC_overlay" || ch_name == "numuCC_background_Pmu5_PC_overlay"
         || ch_name == "BG_numuCC_Pmu5_PC_ext"         || ch_name =="BG_numuCC_Pmu5_PC_dirt" || ch_name == "numuCC_Pmu5_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0) && Pmu_bin==5) {
      if      (ch_name == "numuCC_signal_Pmu5_PC_overlay" || ch_name == "numuCC_background_Pmu5_PC_overlay") { return (map_cuts_flag["Xs_Pmu_numuCCinFV"] == (ch_name=="numuCC_signal_Pmu5_PC_overlay")); }
      else return true;
    }  else return false;
  }else if (ch_name == "numuCC_theta_all_PC_overlay" || ch_name == "BG_numuCC_theta_all_PC_ext" || ch_name =="BG_numuCC_theta_all_PC_dirt" || ch_name == "numuCC_theta_all_PC_bnb") {
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_0p_PC_overlay" || ch_name == "BG_numuCC_theta_all_0p_PC_ext" || ch_name =="BG_numuCC_theta_all_0p_PC_dirt" || ch_name == "numuCC_theta_all_0p_PC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_theta_all_Np_PC_overlay" || ch_name == "BG_numuCC_theta_all_Np_PC_ext" || ch_name =="BG_numuCC_theta_all_Np_PC_dirt" || ch_name == "numuCC_theta_all_Np_PC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all_PC_overlay" || ch_name == "BG_numuCC_all_PC_ext" || ch_name =="BG_numuCC_all_PC_dirt" || ch_name == "numuCC_all_PC_bnb") {
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_all2_PC_overlay" || ch_name == "BG_numuCC_all2_PC_ext" || ch_name =="BG_numuCC_all2_PC_dirt" || ch_name == "numuCC_all2_PC_bnb") {
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_0p_PC_overlay" || ch_name == "BG_numuCC_0p_PC_ext" || ch_name =="BG_numuCC_0p_PC_dirt" || ch_name == "numuCC_0p_PC_bnb") {
    if (flag_numuCC && flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_Np_PC_overlay" || ch_name == "BG_numuCC_Np_PC_ext" || ch_name =="BG_numuCC_Np_PC_dirt" || ch_name == "numuCC_Np_PC_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_Np_both_overlay" || ch_name == "BG_numuCC_Np_both_ext" || ch_name =="BG_numuCC_Np_both_dirt" || ch_name == "numuCC_Np_both_bnb") {
    if (flag_numuCC && (!flag_numuCC_1mu0p) && (!flag_nueCC) && (pfeval.reco_muonMomentum[3]>0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_nopi0_nonueCC_FC_bnb" || ch_name == "numuCC_nopi0_nonueCC_FC_numu2nueoverlay"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_FC_overlay_numi" || ch_name == "BG_numuCC_nopi0_nonueCC_FC_ext_numi" || ch_name =="BG_numuCC_nopi0_nonueCC_FC_dirt_numi" || ch_name == "numuCC_nopi0_nonueCC_FC_numi" || ch_name == "numuCC_nopi0_nonueCC_FC_numu2nueoverlay_numi"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_nopi0_nonueCC_PC_bnb" || ch_name == "numuCC_nopi0_nonueCC_PC_numu2nueoverlay"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_nopi0_nonueCC_PC_overlay_numi" || ch_name == "BG_numuCC_nopi0_nonueCC_PC_ext_numi" || ch_name =="BG_numuCC_nopi0_nonueCC_PC_dirt_numi" || ch_name == "numuCC_nopi0_nonueCC_PC_numi" || ch_name == "numuCC_nopi0_nonueCC_PC_numu2nueoverlay_numi"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay" || ch_name =="BG_CCpi0_nonueCC_FC_ext" || ch_name == "BG_CCpi0_nonueCC_FC_dirt" || ch_name == "CCpi0_nonueCC_FC_bnb" || ch_name == "CCpi0_nonueCC_FC_numu2nueoverlay"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_FC_overlay_numi" || ch_name =="BG_CCpi0_nonueCC_FC_ext_numi" || ch_name == "BG_CCpi0_nonueCC_FC_dirt_numi" || ch_name == "CCpi0_nonueCC_FC_numi" || ch_name == "CCpi0_nonueCC_FC_numu2nueoverlay_numi"){
    if (flag_numuCC && flag_FC && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay" || ch_name == "BG_CCpi0_nonueCC_PC_ext" || ch_name == "BG_CCpi0_nonueCC_PC_dirt" || ch_name == "CCpi0_nonueCC_PC_bnb" || ch_name == "CCpi0_nonueCC_PC_numu2nueoverlay"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "CCpi0_nonueCC_PC_overlay_numi" || ch_name == "BG_CCpi0_nonueCC_PC_ext_numi" || ch_name == "BG_CCpi0_nonueCC_PC_dirt_numi" || ch_name == "CCpi0_nonueCC_PC_numi" || ch_name == "CCpi0_nonueCC_PC_numu2nueoverlay_numi"){
    if (flag_numuCC && (!flag_FC) && flag_cc_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_overlay" || ch_name == "BG_NCpi0_nonueCC_ext" || ch_name == "BG_NCpi0_nonueCC_dirt" || ch_name == "NCpi0_nonueCC_bnb" || ch_name == "NCpi0_nonueCC_numu2nueoverlay"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    // if (flag_NC && flag_pi0 && (!flag_nueCC) && flag_FC && (!flag_0p) ) return true; // a test ...
    else return false;
  }else if (ch_name == "NCpi0_nonueCC_overlay_numi" || ch_name == "BG_NCpi0_nonueCC_ext_numi" || ch_name == "BG_NCpi0_nonueCC_dirt_numi" || ch_name == "NCpi0_nonueCC_numi" || ch_name == "NCpi0_nonueCC_numu2nueoverlay_numi"){
    if (flag_NC && flag_pi0 && (!flag_nueCC) ) return true;
    else return false;
  }else if (ch_name == "nueCC_bnb" || ch_name == "nueCC_nueoverlay"){   // side band ...
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay" || ch_name == "nueCC_bnb") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb" || ch_name == "all_but_nueCC_overlay" || ch_name == "all_but_nueCC_ext" || ch_name == "all_but_nueCC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay" || ch_name != "all_but_nueCC_overlay") return true;
    else return false;
  }else if (ch_name == "nueCC_bnb1" || ch_name == "nueCC_nueoverlay1"){
    if (flag_truth_inside &&  ch_name == "nueCC_nueoverlay1" || ch_name == "nueCC_bnb1") return true;
    else return false;
  }else if (ch_name == "all_but_nueCC_bnb1" || ch_name == "all_but_nueCC_overlay1" || ch_name == "all_but_nueCC_ext1" || ch_name == "all_but_nueCC_dirt1"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "all_but_nueCC_overlay1" || ch_name != "all_but_nueCC_overlay1") return true;
    else return false;
  }else if (ch_name == "testA_bnb" || ch_name == "testA_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testA_nueoverlay" || ch_name == "testA_bnb") return true;
    else return false;
  }else if (ch_name == "testA_overlay" || ch_name == "testA_ext" || ch_name == "testA_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testA_overlay" || ch_name != "testA_overlay") return true;
    else return false;
  }else if (ch_name == "testB_bnb" || ch_name == "testB_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testB_nueoverlay" || ch_name == "testB_bnb") return true;
    else return false;
  }else if (ch_name == "testB_overlay" || ch_name == "testB_ext" || ch_name == "testB_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testB_overlay" || ch_name != "testB_overlay") return true;
    else return false;
  }else if (ch_name == "testC_bnb" || ch_name == "testC_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testC_nueoverlay" || ch_name == "testC_bnb") return true;
    else return false;
  }else if (ch_name == "testC_overlay" || ch_name == "testC_ext" || ch_name == "testC_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testC_overlay" || ch_name != "testC_overlay") return true;
    else return false;
  }else if (ch_name == "testD_bnb" || ch_name == "testD_nueoverlay"){
    if (flag_truth_inside &&  ch_name == "testD_nueoverlay" || ch_name == "testD_bnb") return true;
    else return false;
  }else if (ch_name == "testD_overlay" || ch_name == "testD_ext" || ch_name == "testD_dirt"){
    if (!(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && ch_name == "testD_overlay" || ch_name != "testD_overlay") return true;
    else return false;
 // Janet's requests: <600 MeV numuCC PC, FC for three variables = 6 obs channels
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;
  }else if (ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=0 && reco_Enu<600) return true;
    else return false;

  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC && flag_FC && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;
  }else if (ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_600t1500MeV_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_600t1500MeV_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0) && reco_Enu>=600 && reco_Enu<1500) return true;
    else return false;

  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_extra2_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_extra2_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_extra2_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_extra2_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_1mu0p && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_lowEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_lowEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_lowEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_1mu0p && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC2_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC2_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC2_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC3_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC3_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC3_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC4_highEhad_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC4_highEhad_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC4_highEhad_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_tight && (!flag_numuCC_1mu0p) && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // Mike Shaevitz >800 MeV nueCC PC+FC 1 obs channel
  }else if (ch_name == "nueCC_extra_nueoverlay"){
    if (flag_nueCC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext" || ch_name == "BG_nueCC_extra_dirt" || ch_name =="nueCC_extra_bnb"){
    if (flag_nueCC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;

  }else if (ch_name == "nueCC_extra_nueoverlay_fc"){
    if (flag_nueCC && flag_truth_inside && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext_fc" || ch_name == "BG_nueCC_extra_dirt_fc" || ch_name =="nueCC_extra_bnb_fc"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay_fc"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside) && flag_FC) return true;
    else return false;

    // FC and Np
  }else if (ch_name == "nueCC_extra_nueoverlay_fc_np"){
    if (flag_nueCC && flag_truth_inside && flag_FC && (!flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext_fc_np" || ch_name == "BG_nueCC_extra_dirt_fc_np" || ch_name =="nueCC_extra_bnb_fc_np"){
    if (flag_nueCC&& flag_FC && (!flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay_fc_np"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)&& flag_FC && (!flag_0p)) return true;
    else return false;
    // FC and 0p
  }else if (ch_name == "nueCC_extra_nueoverlay_fc_0p"){
    if (flag_nueCC && flag_truth_inside && flag_FC && (flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_ext_fc_0p" || ch_name == "BG_nueCC_extra_dirt_fc_0p" || ch_name =="nueCC_extra_bnb_fc_0p"){
    if (flag_nueCC&& flag_FC && (flag_0p)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC_extra_overlay_fc_0p"){
    if (flag_nueCC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)&& flag_FC && (flag_0p)) return true;
    else return false;

 // cut-based numuCC FC/PC 2 obs channels
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_FC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_FC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_FC_bnb"){
    if (flag_numuCC_cutbased && flag_FC && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;
  }else if (ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_overlay" || ch_name == "BG_numuCC_cutbased_nopi0_nonueCC_PC_ext" || ch_name =="BG_numuCC_cutbased_nopi0_nonueCC_PC_dirt" || ch_name == "numuCC_cutbased_nopi0_nonueCC_PC_bnb"){
    if (flag_numuCC_cutbased && (!flag_FC) && (!flag_nueCC) && (!flag_cc_pi0)) return true;
    else return false;

 // generic selection nu PC+FC 1 obs channel
}else if (ch_name == "generic_nu_overlay" || ch_name == "BG_generic_nu_ext" || ch_name =="BG_generic_nu_dirt" || ch_name == "generic_nu_bnb" ||
          ch_name == "generic_nu_overlay_2" || ch_name == "BG_generic_nu_ext_2" || ch_name =="BG_generic_nu_dirt_2" || ch_name == "generic_nu_bnb_2" ||
          ch_name == "generic_nu_overlay_3" || ch_name == "BG_generic_nu_ext_3" || ch_name =="BG_generic_nu_dirt_3" || ch_name == "generic_nu_bnb_3" ||
          ch_name == "generic_nu_overlay_4" || ch_name == "BG_generic_nu_ext_4" || ch_name =="BG_generic_nu_dirt_4" || ch_name == "generic_nu_bnb_4"){
    if (flag_generic) return true;
    else return false;
 // numuCC selection PC+FC 1 obs channel
  }else if (ch_name == "numuCC_overlay" || ch_name == "BG_numuCC_ext" || ch_name =="BG_numuCC_dirt" || ch_name == "numuCC_bnb"){
    if (flag_numuCC) return true;
    else return false;
  }else if (ch_name == "numuCC_overlay_fc" || ch_name == "BG_numuCC_ext_fc" || ch_name =="BG_numuCC_dirt_fc" || ch_name == "numuCC_bnb_fc"){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_overlay_fc_np" || ch_name == "BG_numuCC_ext_fc_np" || ch_name =="BG_numuCC_dirt_fc_np" || ch_name == "numuCC_bnb_fc_np"){
    if (flag_numuCC && flag_FC && (!flag_0p)) return true;
    else return false;
  }else if (ch_name == "numuCC_overlay_fc_0p" || ch_name == "BG_numuCC_ext_fc_0p" || ch_name =="BG_numuCC_dirt_fc_0p" || ch_name == "numuCC_bnb_fc_0p"){
    if (flag_numuCC && flag_FC && (flag_0p)) return true;
    else return false;
 // cutbased numuCC selection PC+FC 1 obs channel
  }else if (ch_name == "numuCC_cutbased_overlay" || ch_name == "BG_numuCC_cutbased_ext" || ch_name =="BG_numuCC_cutbased_dirt" || ch_name == "numuCC_cutbased_bnb"){
    if (flag_numuCC_cutbased) return true;
    else return false;
 // nueCC 3 variables: n_trakcs, n_showers, gap_n_bad, FC/PC x3 = 6 channels; 4 additional channels
  }else if (ch_name == "nueCC2_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_ext" || ch_name == "BG_nueCC2_FC_dirt" || ch_name =="nueCC2_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC2_PC_nueoverlay" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_ext" || ch_name == "BG_nueCC2_PC_dirt" || ch_name == "nueCC2_PC_bnb"){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC2_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC3_FC_nueoverlay"){
    if (flag_nueCC && flag_FC && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_ext" || ch_name == "BG_nueCC3_FC_dirt" || ch_name =="nueCC3_FC_bnb"){
    if (flag_nueCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_FC_overlay"){
    if (flag_nueCC && flag_FC && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
  }else if (ch_name == "nueCC3_PC_nueoverlay" ){
    if (flag_nueCC && (!flag_FC) && flag_truth_inside) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_ext" || ch_name == "BG_nueCC3_PC_dirt" || ch_name == "nueCC3_PC_bnb"){
    if (flag_nueCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "BG_nueCC3_PC_overlay"){
    if (flag_nueCC && (!flag_FC) && !(eval.truth_isCC==1 && abs(eval.truth_nuPdg)==12 && flag_truth_inside)) return true;
    else return false;
    // add some cuts for Xs related cases ...
  }else if (ch_name == "numuCC_FC_bnb" || ch_name == "BG_numuCC_FC_ext" || ch_name == "BG_numuCC_FC_dirt"
	    || ch_name == "numuCC1_FC_bnb" || ch_name == "BG_numuCC1_FC_ext" || ch_name == "BG_numuCC1_FC_dirt"
	    || ch_name == "numuCC2_FC_bnb" || ch_name == "BG_numuCC2_FC_ext" || ch_name == "BG_numuCC2_FC_dirt"
	    ){
    if (flag_numuCC && flag_FC) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb" || ch_name == "BG_numuCC_PC_ext" || ch_name == "BG_numuCC_PC_dirt"
	    || ch_name == "numuCC1_PC_bnb" || ch_name == "BG_numuCC1_PC_ext" || ch_name == "BG_numuCC1_PC_dirt"
	    || ch_name == "numuCC2_PC_bnb" || ch_name == "BG_numuCC2_PC_ext" || ch_name == "BG_numuCC2_PC_dirt"
	    ){
    if (flag_numuCC && (!flag_FC)) return true;
    else return false;
  }else if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC_background_PC_overlay"
  	    || ch_name == "numuCC1_signal_FC_overlay" || ch_name == "numuCC1_signal_PC_overlay" || ch_name == "numuCC1_background_FC_overlay" || ch_name == "numuCC1_background_PC_overlay"
  	    || ch_name == "numuCC2_signal_FC_overlay" || ch_name == "numuCC2_signal_PC_overlay" || ch_name == "numuCC2_background_FC_overlay" || ch_name == "numuCC2_background_PC_overlay"
  	    ){
    if (ch_name == "numuCC_signal_FC_overlay" || ch_name == "numuCC1_signal_FC_overlay" || ch_name == "numuCC2_signal_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_PC_overlay" || ch_name == "numuCC1_signal_PC_overlay" || ch_name == "numuCC2_signal_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["XsnumuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_FC_overlay" || ch_name == "numuCC1_background_FC_overlay" || ch_name == "numuCC2_background_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_PC_overlay" || ch_name == "numuCC1_background_PC_overlay" || ch_name == "numuCC2_background_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["XsnumuCCinFV"])) return true;
    }
    return false;
  } else if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC_background_Enu_PC_overlay"
	     || ch_name == "numuCC1_signal_Enu_FC_overlay" || ch_name == "numuCC1_signal_Enu_PC_overlay" || ch_name == "numuCC1_background_Enu_FC_overlay" || ch_name == "numuCC1_background_Enu_PC_overlay"
	     || ch_name == "numuCC2_signal_Enu_FC_overlay" || ch_name == "numuCC2_signal_Enu_PC_overlay" || ch_name == "numuCC2_background_Enu_FC_overlay" || ch_name == "numuCC2_background_Enu_PC_overlay"
       || ch_name == "numuCC_signal_Enu_overlay" || ch_name == "numuCC_background_Enu_overlay"
	    ){
    if (ch_name == "numuCC_signal_Enu_FC_overlay" || ch_name == "numuCC1_signal_Enu_FC_overlay" || ch_name == "numuCC2_signal_Enu_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Enu_PC_overlay" || ch_name == "numuCC1_signal_Enu_PC_overlay" || ch_name == "numuCC2_signal_Enu_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Enu_FC_overlay" || ch_name == "numuCC1_background_Enu_FC_overlay" || ch_name == "numuCC2_background_Enu_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Enu_PC_overlay" || ch_name == "numuCC1_background_Enu_PC_overlay" || ch_name == "numuCC2_background_Enu_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_signal_Enu_overlay"){
      if (flag_numuCC && map_cuts_flag["Xs_Enu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Enu_overlay"){
      if (flag_numuCC && (!map_cuts_flag["Xs_Enu_numuCCinFV"])) return true;
    }
    return false;
  // 1D Enu channel with the same inclusive signal definition as the 3D selection
  }else if (ch_name == "numuCC_Enu_mu_FC_bnb" || ch_name == "BG_numuCC_Enu_mu_FC_ext" || ch_name == "BG_numuCC_Enu_mu_FC_dirt"){
    if (flag_numuCC &&   flag_FC  && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8) return true;
    else return false;
  }else if (ch_name == "numuCC_Enu_mu_PC_bnb" || ch_name == "BG_numuCC_Enu_mu_PC_ext" || ch_name == "BG_numuCC_Enu_mu_PC_dirt"){
    if (flag_numuCC && (!flag_FC) && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8) return true;
    else return false;
  } else if (ch_name == "numuCC_signal_Enu_mu_FC_overlay"    || ch_name == "numuCC_signal_Enu_mu_PC_overlay"
         || ch_name == "numuCC_background_Enu_mu_FC_overlay" || ch_name == "numuCC_background_Enu_mu_PC_overlay" ){
    bool pre_cut = flag_numuCC && (!flag_nueCC) && pfeval.reco_muonMomentum[3]>0 && Enu_bin>=0 && Enu_bin<=3 && costheta_bin>=0 && costheta_bin<=8;
    if      (ch_name == "numuCC_signal_Enu_mu_FC_overlay"     && pre_cut &&   flag_FC  &&   map_cuts_flag["Xs_Enu_mu_numuCCinFV"])  { return true; }
    else if (ch_name == "numuCC_signal_Enu_mu_PC_overlay"     && pre_cut && (!flag_FC) &&   map_cuts_flag["Xs_Enu_mu_numuCCinFV"])  { return true; }
    else if (ch_name == "numuCC_background_Enu_mu_FC_overlay" && pre_cut &&   flag_FC  && (!map_cuts_flag["Xs_Enu_mu_numuCCinFV"])) { return true; }
    else if (ch_name == "numuCC_background_Enu_mu_PC_overlay" && pre_cut && (!flag_FC) && (!map_cuts_flag["Xs_Enu_mu_numuCCinFV"])) { return true; }
    return false;
  // ------
  }else if (ch_name == "nueCC_signal_Enu_FC_overlay" || ch_name == "nueCC_signal_Enu_PC_overlay" || ch_name == "nueCC_background_Enu_FC_overlay" || ch_name == "nueCC_background_Enu_PC_overlay"
       || ch_name == "nueCC_signal_Enu_overlay" || ch_name == "nueCC_background_Enu_overlay"
      ){
    if (ch_name == "nueCC_signal_Enu_FC_overlay"){
      if (flag_nueCC && flag_FC && map_cuts_flag["Xs_Enu_nueCCinFV"]) return true;
    }else if (ch_name == "nueCC_signal_Enu_PC_overlay"){
      if (flag_nueCC && (!flag_FC) && map_cuts_flag["Xs_Enu_nueCCinFV"]) return true;
    }else if (ch_name == "nueCC_background_Enu_FC_overlay"){
      if (flag_nueCC && flag_FC && (!map_cuts_flag["Xs_Enu_nueCCinFV"])) return true;
    }else if (ch_name == "nueCC_background_Enu_PC_overlay" ){
      if (flag_nueCC && (!flag_FC) && (!map_cuts_flag["Xs_Enu_nueCCinFV"])) return true;
    }else if (ch_name == "nueCC_signal_Enu_overlay"){
      if (flag_nueCC && map_cuts_flag["Xs_Enu_nueCCinFV"]) return true;
    }else if (ch_name == "nueCC_background_Enu_overlay"){
      if (flag_nueCC && (!map_cuts_flag["Xs_Enu_nueCCinFV"])) return true;
    }
    return false;

  } else if (ch_name == "numuCC_signal_Emu_FC_overlay" || ch_name == "numuCC_signal_Emu_PC_overlay" || ch_name == "numuCC_background_Emu_FC_overlay" || ch_name == "numuCC_background_Emu_PC_overlay"
	     || ch_name == "numuCC1_signal_Emu_FC_overlay" || ch_name == "numuCC1_signal_Emu_PC_overlay" || ch_name == "numuCC1_background_Emu_FC_overlay" || ch_name == "numuCC1_background_Emu_PC_overlay"
	     || ch_name == "numuCC2_signal_Emu_FC_overlay" || ch_name == "numuCC2_signal_Emu_PC_overlay" || ch_name == "numuCC2_background_Emu_FC_overlay" || ch_name == "numuCC2_background_Emu_PC_overlay"
	     ){
    if (ch_name == "numuCC_signal_Emu_FC_overlay" || ch_name == "numuCC1_signal_Emu_FC_overlay" || ch_name == "numuCC2_signal_Emu_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Emu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Emu_PC_overlay" || ch_name == "numuCC1_signal_Emu_PC_overlay" || ch_name == "numuCC2_signal_Emu_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Emu_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Emu_FC_overlay" || ch_name == "numuCC1_background_Emu_FC_overlay" || ch_name == "numuCC2_background_Emu_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Emu_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Emu_PC_overlay" || ch_name == "numuCC1_background_Emu_PC_overlay" || ch_name == "numuCC2_background_Emu_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Emu_numuCCinFV"])) return true;
    }
    return false;

  } else if (ch_name == "numuCC_signal_Ehad_FC_overlay" || ch_name == "numuCC_signal_Ehad_PC_overlay" || ch_name == "numuCC_background_Ehad_FC_overlay" || ch_name == "numuCC_background_Ehad_PC_overlay"
	     || ch_name == "numuCC1_signal_Ehad_FC_overlay" || ch_name == "numuCC1_signal_Ehad_PC_overlay" || ch_name == "numuCC1_background_Ehad_FC_overlay" || ch_name == "numuCC1_background_Ehad_PC_overlay"
	     || ch_name == "numuCC2_signal_Ehad_FC_overlay" || ch_name == "numuCC2_signal_Ehad_PC_overlay" || ch_name == "numuCC2_background_Ehad_FC_overlay" || ch_name == "numuCC2_background_Ehad_PC_overlay"
	    ){
    if (ch_name == "numuCC_signal_Ehad_FC_overlay" || ch_name == "numuCC1_signal_Ehad_FC_overlay" || ch_name == "numuCC2_signal_Ehad_FC_overlay"){
      if (flag_numuCC && flag_FC && map_cuts_flag["Xs_Ehad_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_signal_Ehad_PC_overlay" || ch_name == "numuCC1_signal_Ehad_PC_overlay" || ch_name == "numuCC2_signal_Ehad_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) && map_cuts_flag["Xs_Ehad_numuCCinFV"]) return true;
    }else if (ch_name == "numuCC_background_Ehad_FC_overlay" || ch_name == "numuCC1_background_Ehad_FC_overlay" || ch_name == "numuCC2_background_Ehad_FC_overlay"){
      if (flag_numuCC && flag_FC && (!map_cuts_flag["Xs_Ehad_numuCCinFV"])) return true;
    }else if (ch_name == "numuCC_background_Ehad_PC_overlay" || ch_name == "numuCC1_background_Ehad_PC_overlay" || ch_name == "numuCC2_background_Ehad_PC_overlay"){
      if (flag_numuCC && (!flag_FC) && (!map_cuts_flag["Xs_Ehad_numuCCinFV"])) return true;
    }
    return false;

  }else if (ch_name == "numuCC_FC_bnb_L800MeV" || ch_name == "BG_numuCC_FC_ext_L800MeV" || ch_name == "BG_numuCC_FC_dirt_L800MeV"
	    || ch_name == "numuCC1_FC_bnb_L800MeV" || ch_name == "BG_numuCC1_FC_ext_L800MeV" || ch_name == "BG_numuCC1_FC_dirt_L800MeV"
	    || ch_name == "numuCC2_FC_bnb_L800MeV" || ch_name == "BG_numuCC2_FC_ext_L800MeV" || ch_name == "BG_numuCC2_FC_dirt_L800MeV"
	    ){
    if (flag_numuCC && flag_FC && reco_Enu<800) return true;
    else return false;
  }else if (ch_name == "numuCC_PC_bnb_L800MeV" || ch_name == "BG_numuCC_PC_ext_L800MeV" || ch_name == "BG_numuCC_PC_dirt_L800MeV"
	    || ch_name == "numuCC1_PC_bnb_L800MeV" || ch_name == "BG_numuCC1_PC_ext_L800MeV" || ch_name == "BG_numuCC1_PC_dirt_L800MeV"
	    || ch_name == "numuCC2_PC_bnb_L800MeV" || ch_name == "BG_numuCC2_PC_ext_L800MeV" || ch_name == "BG_numuCC2_PC_dirt_L800MeV"
	    ){
    if (flag_numuCC && (!flag_FC) && reco_Enu<800) return true;
    else return false;
  }else if (ch_name == "numuCC_FC_overlay_L800MeV" || ch_name == "numuCC_PC_overlay_L800MeV"
	    || ch_name == "numuCC1_FC_overlay_L800MeV" || ch_name == "numuCC1_PC_overlay_L800MeV"
	    || ch_name == "numuCC2_FC_overlay_L800MeV" || ch_name == "numuCC2_PC_overlay_L800MeV"   ){
    if (ch_name == "numuCC_FC_overlay_L800MeV" || ch_name == "numuCC1_FC_overlay_L800MeV" || ch_name == "numuCC2_FC_overlay_L800MeV"){
      if (flag_numuCC && flag_FC && reco_Enu<800) return true;
    }else if (ch_name == "numuCC_PC_overlay_L800MeV" || ch_name == "numuCC1_PC_overlay_L800MeV" || ch_name == "numuCC2_PC_overlay_L800MeV" ){
      if (flag_numuCC && (!flag_FC) && reco_Enu<800) return true;
    }
    return false;
  }else if (ch_name == "numuCC_FC_overlay" || ch_name == "numuCC_PC_overlay"
	    || ch_name == "numuCC1_FC_overlay" || ch_name == "numuCC1_PC_overlay"
	    || ch_name == "numuCC2_FC_overlay" || ch_name == "numuCC2_PC_overlay"   ){
    if (ch_name == "numuCC_FC_overlay" || ch_name == "numuCC1_FC_overlay" || ch_name == "numuCC2_FC_overlay"){
      if (flag_numuCC && flag_FC ) return true;
    }else if (ch_name == "numuCC_PC_overlay" || ch_name == "numuCC1_PC_overlay" || ch_name == "numuCC2_PC_overlay" ){
      if (flag_numuCC && (!flag_FC) ) return true;
    }
    return false;
  }else if (ch_name == "nc_pio_energy_FC" || ch_name == "nc_pio_score_FC"
	    || ch_name == "nc_pio_energy_FC_ncpio_overlay" || ch_name == "nc_pio_score_FC_ncpio_overlay"
	    || ch_name == "nc_pio_energy_FC_ncdelta_overlay" || ch_name == "nc_pio_score_FC_ncdelta_overlay"
	    || ch_name == "nc_pio_energy_FC_overlay" || ch_name == "nc_pio_score_FC_overlay"
	    || ch_name == "nc_pio_energy_FC_ext" || ch_name == "nc_pio_score_FC_ext"
	    || ch_name == "nc_pio_energy_FC_dirt" || ch_name == "nc_pio_score_FC_dirt"
	    ){
    if (ch_name == "nc_pio_energy_FC"
	|| ch_name == "nc_pio_energy_FC_ext"
	|| ch_name == "nc_pio_energy_FC_dirt" ){
      if (flag_ncpio_sel && flag_FC) return true;
    }else if (ch_name == "nc_pio_energy_FC_ncpio_overlay" ){
      if (flag_ncpio_sel && flag_FC && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_pio_energy_FC_ncdelta_overlay" ){
      if (flag_ncpio_sel && flag_FC && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_pio_energy_FC_overlay" ){
      if (flag_ncpio_sel && flag_FC && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //if (flag_ncpio_sel && flag_FC) return true;
    }else if (ch_name == "nc_pio_score_FC"
	|| ch_name == "nc_pio_score_FC_ext"
	|| ch_name == "nc_pio_score_FC_dirt"){
      if (flag_FC) return true;
    }else if (ch_name == "nc_pio_score_FC_ncpio_overlay"){
      if (flag_FC && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_pio_score_FC_ncdelta_overlay"){
      if (flag_FC && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_pio_score_FC_overlay"){
      if (flag_FC && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
			&& (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //  if (flag_FC) return true;
    }

    return false;
  // start cuts from Lee's cuts.h
  // NC delta (1g) channels:
  }else if (ch_name == "nc_delta_0p_01" || ch_name == "nc_delta_0p_02" || ch_name == "nc_delta_0p_03" || ch_name == "nc_delta_0p_04"
                    || ch_name == "nc_delta_0p_05" || ch_name == "nc_delta_0p_06" || ch_name == "nc_delta_0p_07" || ch_name == "nc_delta_0p_08"
                    || ch_name == "nc_delta_0p_09" || ch_name == "nc_delta_0p_10" || ch_name == "nc_delta_0p_11" || ch_name == "nc_delta_0p_12"
                    || ch_name == "nc_delta_0p_13" || ch_name == "nc_delta_0p_14" || ch_name == "nc_delta_0p_15" || ch_name == "nc_delta_0p_16"
                    || ch_name == "nc_delta_0p_17" || ch_name == "nc_delta_0p_18" || ch_name == "nc_delta_0p_19" || ch_name == "nc_delta_0p_20"){
            if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
  }else if (ch_name == "nc_delta_Np_01" || ch_name == "nc_delta_Np_02" || ch_name == "nc_delta_Np_03" || ch_name == "nc_delta_Np_04"
            || ch_name == "nc_delta_Np_05" || ch_name == "nc_delta_Np_06" || ch_name == "nc_delta_Np_07" || ch_name == "nc_delta_Np_08"
            || ch_name == "nc_delta_Np_09" || ch_name == "nc_delta_Np_10" || ch_name == "nc_delta_Np_11" || ch_name == "nc_delta_Np_12"
            || ch_name == "nc_delta_Np_13" || ch_name == "nc_delta_Np_14" || ch_name == "nc_delta_Np_15" || ch_name == "nc_delta_Np_16"
            || ch_name == "nc_delta_Np_17" || ch_name == "nc_delta_Np_18" || ch_name == "nc_delta_Np_19" || ch_name == "nc_delta_Np_20"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
  }else if (ch_name == "nc_delta_0p_01_ext" || ch_name == "nc_delta_0p_02_ext" || ch_name == "nc_delta_0p_03_ext" || ch_name == "nc_delta_0p_04_ext"
            || ch_name == "nc_delta_0p_05_ext" || ch_name == "nc_delta_0p_06_ext" || ch_name == "nc_delta_0p_07_ext" || ch_name == "nc_delta_0p_08_ext"
            || ch_name == "nc_delta_0p_09_ext" || ch_name == "nc_delta_0p_10_ext" || ch_name == "nc_delta_0p_11_ext" || ch_name == "nc_delta_0p_12_ext"
            || ch_name == "nc_delta_0p_13_ext" || ch_name == "nc_delta_0p_14_ext" || ch_name == "nc_delta_0p_15_ext" || ch_name == "nc_delta_0p_16_ext"
            || ch_name == "nc_delta_0p_17_ext" || ch_name == "nc_delta_0p_18_ext" || ch_name == "nc_delta_0p_19_ext" || ch_name == "nc_delta_0p_20_ext"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
  }else if (ch_name == "nc_delta_Np_01_ext" || ch_name == "nc_delta_Np_02_ext" || ch_name == "nc_delta_Np_03_ext" || ch_name == "nc_delta_Np_04_ext"
            || ch_name == "nc_delta_Np_05_ext" || ch_name == "nc_delta_Np_06_ext" || ch_name == "nc_delta_Np_07_ext" || ch_name == "nc_delta_Np_08_ext"
            || ch_name == "nc_delta_Np_09_ext" || ch_name == "nc_delta_Np_10_ext" || ch_name == "nc_delta_Np_11_ext" || ch_name == "nc_delta_Np_12_ext"
            || ch_name == "nc_delta_Np_13_ext" || ch_name == "nc_delta_Np_14_ext" || ch_name == "nc_delta_Np_15_ext" || ch_name == "nc_delta_Np_16_ext"
            || ch_name == "nc_delta_Np_17_ext" || ch_name == "nc_delta_Np_18_ext" || ch_name == "nc_delta_Np_19_ext" || ch_name == "nc_delta_Np_20_ext"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
  }else if (ch_name == "nc_delta_0p_01_dirt" || ch_name == "nc_delta_0p_02_dirt" || ch_name == "nc_delta_0p_03_dirt" || ch_name == "nc_delta_0p_04_dirt"
            || ch_name == "nc_delta_0p_05_dirt" || ch_name == "nc_delta_0p_06_dirt" || ch_name == "nc_delta_0p_07_dirt" || ch_name == "nc_delta_0p_08_dirt"
            || ch_name == "nc_delta_0p_09_dirt" || ch_name == "nc_delta_0p_10_dirt" || ch_name == "nc_delta_0p_11_dirt" || ch_name == "nc_delta_0p_12_dirt"
            || ch_name == "nc_delta_0p_13_dirt" || ch_name == "nc_delta_0p_14_dirt" || ch_name == "nc_delta_0p_15_dirt" || ch_name == "nc_delta_0p_16_dirt"
            || ch_name == "nc_delta_0p_17_dirt" || ch_name == "nc_delta_0p_18_dirt" || ch_name == "nc_delta_0p_19_dirt" || ch_name == "nc_delta_0p_20_dirt"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
  }else if (ch_name == "nc_delta_Np_01_dirt" || ch_name == "nc_delta_Np_02_dirt" || ch_name == "nc_delta_Np_03_dirt" || ch_name == "nc_delta_Np_04_dirt"
            || ch_name == "nc_delta_Np_05_dirt" || ch_name == "nc_delta_Np_06_dirt" || ch_name == "nc_delta_Np_07_dirt" || ch_name == "nc_delta_Np_08_dirt"
            || ch_name == "nc_delta_Np_09_dirt" || ch_name == "nc_delta_Np_10_dirt" || ch_name == "nc_delta_Np_11_dirt" || ch_name == "nc_delta_Np_12_dirt"
            || ch_name == "nc_delta_Np_13_dirt" || ch_name == "nc_delta_Np_14_dirt" || ch_name == "nc_delta_Np_15_dirt" || ch_name == "nc_delta_Np_16_dirt"
            || ch_name == "nc_delta_Np_17_dirt" || ch_name == "nc_delta_Np_18_dirt" || ch_name == "nc_delta_Np_19_dirt" || ch_name == "nc_delta_Np_20_dirt"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
  }else if (ch_name == "nc_delta_0p_01_nc_delta_overlay" || ch_name == "nc_delta_0p_01_nc_delta_overlay_add" || ch_name == "nc_delta_0p_02_nc_delta_overlay" || ch_name == "nc_delta_0p_03_nc_delta_overlay" || ch_name == "nc_delta_0p_04_nc_delta_overlay"
            || ch_name == "nc_delta_0p_05_nc_delta_overlay" || ch_name == "nc_delta_0p_06_nc_delta_overlay" || ch_name == "nc_delta_0p_07_nc_delta_overlay" || ch_name == "nc_delta_0p_08_nc_delta_overlay"
            || ch_name == "nc_delta_0p_09_nc_delta_overlay" || ch_name == "nc_delta_0p_10_nc_delta_overlay" || ch_name == "nc_delta_0p_11_nc_delta_overlay" || ch_name == "nc_delta_0p_12_nc_delta_overlay"
            || ch_name == "nc_delta_0p_13_nc_delta_overlay" || ch_name == "nc_delta_0p_14_nc_delta_overlay" || ch_name == "nc_delta_0p_15_nc_delta_overlay" || ch_name == "nc_delta_0p_16_nc_delta_overlay"
            || ch_name == "nc_delta_0p_17_nc_delta_overlay" || ch_name == "nc_delta_0p_18_nc_delta_overlay" || ch_name == "nc_delta_0p_19_nc_delta_overlay" || ch_name == "nc_delta_0p_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_nc_delta_overlay" || ch_name == "nc_delta_Np_01_nc_delta_overlay_add" || ch_name == "nc_delta_Np_02_nc_delta_overlay" || ch_name == "nc_delta_Np_03_nc_delta_overlay" || ch_name == "nc_delta_Np_04_nc_delta_overlay"
            || ch_name == "nc_delta_Np_05_nc_delta_overlay" || ch_name == "nc_delta_Np_06_nc_delta_overlay" || ch_name == "nc_delta_Np_07_nc_delta_overlay" || ch_name == "nc_delta_Np_08_nc_delta_overlay"
            || ch_name == "nc_delta_Np_09_nc_delta_overlay" || ch_name == "nc_delta_Np_10_nc_delta_overlay" || ch_name == "nc_delta_Np_11_nc_delta_overlay" || ch_name == "nc_delta_Np_12_nc_delta_overlay"
            || ch_name == "nc_delta_Np_13_nc_delta_overlay" || ch_name == "nc_delta_Np_14_nc_delta_overlay" || ch_name == "nc_delta_Np_15_nc_delta_overlay" || ch_name == "nc_delta_Np_16_nc_delta_overlay"
            || ch_name == "nc_delta_Np_17_nc_delta_overlay" || ch_name == "nc_delta_Np_18_nc_delta_overlay" || ch_name == "nc_delta_Np_19_nc_delta_overlay" || ch_name == "nc_delta_Np_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_true_Np_nc_delta_overlay" || ch_name == "nc_delta_0p_01_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_delta_0p_01_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_delta_0p_01_true_Np_nc_delta_overlay"));
            }
            return false;

    }else if (ch_name == "nc_delta_Np_01_true_0p_nc_delta_overlay" || ch_name == "nc_delta_Np_01_true_Np_nc_delta_overlay"){
            if (flag_FC && flag_ncdelta_sel && !(flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_delta_Np_01_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_delta_Np_01_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_delta_0p_01_nc_pi0_overlay" || ch_name == "nc_delta_0p_02_nc_pi0_overlay" || ch_name == "nc_delta_0p_03_nc_pi0_overlay" || ch_name == "nc_delta_0p_04_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_05_nc_pi0_overlay" || ch_name == "nc_delta_0p_06_nc_pi0_overlay" || ch_name == "nc_delta_0p_07_nc_pi0_overlay" || ch_name == "nc_delta_0p_08_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_09_nc_pi0_overlay" || ch_name == "nc_delta_0p_10_nc_pi0_overlay" || ch_name == "nc_delta_0p_11_nc_pi0_overlay" || ch_name == "nc_delta_0p_12_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_13_nc_pi0_overlay" || ch_name == "nc_delta_0p_14_nc_pi0_overlay" || ch_name == "nc_delta_0p_15_nc_pi0_overlay" || ch_name == "nc_delta_0p_16_nc_pi0_overlay"
            || ch_name == "nc_delta_0p_17_nc_pi0_overlay" || ch_name == "nc_delta_0p_18_nc_pi0_overlay" || ch_name == "nc_delta_0p_19_nc_pi0_overlay" || ch_name == "nc_delta_0p_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_nc_pi0_overlay" || ch_name == "nc_delta_Np_02_nc_pi0_overlay" || ch_name == "nc_delta_Np_03_nc_pi0_overlay" || ch_name == "nc_delta_Np_04_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_05_nc_pi0_overlay" || ch_name == "nc_delta_Np_06_nc_pi0_overlay" || ch_name == "nc_delta_Np_07_nc_pi0_overlay" || ch_name == "nc_delta_Np_08_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_09_nc_pi0_overlay" || ch_name == "nc_delta_Np_10_nc_pi0_overlay" || ch_name == "nc_delta_Np_11_nc_pi0_overlay" || ch_name == "nc_delta_Np_12_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_13_nc_pi0_overlay" || ch_name == "nc_delta_Np_14_nc_pi0_overlay" || ch_name == "nc_delta_Np_15_nc_pi0_overlay" || ch_name == "nc_delta_Np_16_nc_pi0_overlay"
            || ch_name == "nc_delta_Np_17_nc_pi0_overlay" || ch_name == "nc_delta_Np_18_nc_pi0_overlay" || ch_name == "nc_delta_Np_19_nc_pi0_overlay" || ch_name == "nc_delta_Np_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                        && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
   }else if (ch_name == "nc_delta_0p_01_overlay" || ch_name == "nc_delta_0p_02_overlay" || ch_name == "nc_delta_0p_03_overlay" || ch_name == "nc_delta_0p_04_overlay"
            || ch_name == "nc_delta_0p_05_overlay" || ch_name == "nc_delta_0p_06_overlay" || ch_name == "nc_delta_0p_07_overlay" || ch_name == "nc_delta_0p_08_overlay"
            || ch_name == "nc_delta_0p_09_overlay" || ch_name == "nc_delta_0p_10_overlay" || ch_name == "nc_delta_0p_11_overlay" || ch_name == "nc_delta_0p_12_overlay"
            || ch_name == "nc_delta_0p_13_overlay" || ch_name == "nc_delta_0p_14_overlay" || ch_name == "nc_delta_0p_15_overlay" || ch_name == "nc_delta_0p_16_overlay"
            || ch_name == "nc_delta_0p_17_overlay" || ch_name == "nc_delta_0p_18_overlay" || ch_name == "nc_delta_0p_19_overlay" || ch_name == "nc_delta_0p_20_overlay"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                               && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_overlay" || ch_name == "nc_delta_Np_02_overlay" || ch_name == "nc_delta_Np_03_overlay" || ch_name == "nc_delta_Np_04_overlay"
            || ch_name == "nc_delta_Np_05_overlay" || ch_name == "nc_delta_Np_06_overlay" || ch_name == "nc_delta_Np_07_overlay" || ch_name == "nc_delta_Np_08_overlay"
            || ch_name == "nc_delta_Np_09_overlay" || ch_name == "nc_delta_Np_10_overlay" || ch_name == "nc_delta_Np_11_overlay" || ch_name == "nc_delta_Np_12_overlay"
            || ch_name == "nc_delta_Np_13_overlay" || ch_name == "nc_delta_Np_14_overlay" || ch_name == "nc_delta_Np_15_overlay" || ch_name == "nc_delta_Np_16_overlay"
            || ch_name == "nc_delta_Np_17_overlay" || ch_name == "nc_delta_Np_18_overlay" || ch_name == "nc_delta_Np_19_overlay" || ch_name == "nc_delta_Np_20_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_0p_01_overlay_entire" || ch_name == "nc_delta_0p_02_overlay_entire" || ch_name == "nc_delta_0p_03_overlay_entire" || ch_name == "nc_delta_0p_04_overlay_entire"){
                  if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_delta_Np_01_overlay_entire" || ch_name == "nc_delta_Np_02_overlay_entire" || ch_name == "nc_delta_Np_03_overlay_entire" || ch_name == "nc_delta_Np_04_overlay_entire"){
                  if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01" || ch_name == "nc_delta_Xp_02" || ch_name == "nc_delta_Xp_03" || ch_name == "nc_delta_Xp_04"
            || ch_name == "nc_delta_Xp_05" || ch_name == "nc_delta_Xp_06" || ch_name == "nc_delta_Xp_07" || ch_name == "nc_delta_Xp_08"
            || ch_name == "nc_delta_Xp_09" || ch_name == "nc_delta_Xp_10" || ch_name == "nc_delta_Xp_11" || ch_name == "nc_delta_Xp_12"
            || ch_name == "nc_delta_Xp_13" || ch_name == "nc_delta_Xp_14" || ch_name == "nc_delta_Xp_15" || ch_name == "nc_delta_Xp_16"
            || ch_name == "nc_delta_Xp_17" || ch_name == "nc_delta_Xp_18" || ch_name == "nc_delta_Xp_19" || ch_name == "nc_delta_Xp_20"){
                  if (flag_FC && flag_ncdelta_sel) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_ext" || ch_name == "nc_delta_Xp_02_ext" || ch_name == "nc_delta_Xp_03_ext" || ch_name == "nc_delta_Xp_04_ext"
            || ch_name == "nc_delta_Xp_05_ext" || ch_name == "nc_delta_Xp_06_ext" || ch_name == "nc_delta_Xp_07_ext" || ch_name == "nc_delta_Xp_08_ext"
            || ch_name == "nc_delta_Xp_09_ext" || ch_name == "nc_delta_Xp_10_ext" || ch_name == "nc_delta_Xp_11_ext" || ch_name == "nc_delta_Xp_12_ext"
            || ch_name == "nc_delta_Xp_13_ext" || ch_name == "nc_delta_Xp_14_ext" || ch_name == "nc_delta_Xp_15_ext" || ch_name == "nc_delta_Xp_16_ext"
            || ch_name == "nc_delta_Xp_17_ext" || ch_name == "nc_delta_Xp_18_ext" || ch_name == "nc_delta_Xp_19_ext" || ch_name == "nc_delta_Xp_20_ext"){
                  if (flag_FC && flag_ncdelta_sel) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_dirt" || ch_name == "nc_delta_Xp_02_dirt" || ch_name == "nc_delta_Xp_03_dirt" || ch_name == "nc_delta_Xp_04_dirt"
            || ch_name == "nc_delta_Xp_05_dirt" || ch_name == "nc_delta_Xp_06_dirt" || ch_name == "nc_delta_Xp_07_dirt" || ch_name == "nc_delta_Xp_08_dirt"
            || ch_name == "nc_delta_Xp_09_dirt" || ch_name == "nc_delta_Xp_10_dirt" || ch_name == "nc_delta_Xp_11_dirt" || ch_name == "nc_delta_Xp_12_dirt"
            || ch_name == "nc_delta_Xp_13_dirt" || ch_name == "nc_delta_Xp_14_dirt" || ch_name == "nc_delta_Xp_15_dirt" || ch_name == "nc_delta_Xp_16_dirt"
            || ch_name == "nc_delta_Xp_17_dirt" || ch_name == "nc_delta_Xp_18_dirt" || ch_name == "nc_delta_Xp_19_dirt" || ch_name == "nc_delta_Xp_20_dirt"){
                  if (flag_FC && flag_ncdelta_sel) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_nc_delta_overlay" || ch_name == "nc_delta_Xp_02_nc_delta_overlay" || ch_name == "nc_delta_Xp_03_nc_delta_overlay" || ch_name == "nc_delta_Xp_04_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_05_nc_delta_overlay" || ch_name == "nc_delta_Xp_06_nc_delta_overlay" || ch_name == "nc_delta_Xp_07_nc_delta_overlay" || ch_name == "nc_delta_Xp_08_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_09_nc_delta_overlay" || ch_name == "nc_delta_Xp_10_nc_delta_overlay" || ch_name == "nc_delta_Xp_11_nc_delta_overlay" || ch_name == "nc_delta_Xp_12_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_13_nc_delta_overlay" || ch_name == "nc_delta_Xp_14_nc_delta_overlay" || ch_name == "nc_delta_Xp_15_nc_delta_overlay" || ch_name == "nc_delta_Xp_16_nc_delta_overlay"
            || ch_name == "nc_delta_Xp_17_nc_delta_overlay" || ch_name == "nc_delta_Xp_18_nc_delta_overlay" || ch_name == "nc_delta_Xp_19_nc_delta_overlay" || ch_name == "nc_delta_Xp_20_nc_delta_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_nc_pi0_overlay" || ch_name == "nc_delta_Xp_02_nc_pi0_overlay" || ch_name == "nc_delta_Xp_03_nc_pi0_overlay" || ch_name == "nc_delta_Xp_04_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_05_nc_pi0_overlay" || ch_name == "nc_delta_Xp_06_nc_pi0_overlay" || ch_name == "nc_delta_Xp_07_nc_pi0_overlay" || ch_name == "nc_delta_Xp_08_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_09_nc_pi0_overlay" || ch_name == "nc_delta_Xp_10_nc_pi0_overlay" || ch_name == "nc_delta_Xp_11_nc_pi0_overlay" || ch_name == "nc_delta_Xp_12_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_13_nc_pi0_overlay" || ch_name == "nc_delta_Xp_14_nc_pi0_overlay" || ch_name == "nc_delta_Xp_15_nc_pi0_overlay" || ch_name == "nc_delta_Xp_16_nc_pi0_overlay"
            || ch_name == "nc_delta_Xp_17_nc_pi0_overlay" || ch_name == "nc_delta_Xp_18_nc_pi0_overlay" || ch_name == "nc_delta_Xp_19_nc_pi0_overlay" || ch_name == "nc_delta_Xp_20_nc_pi0_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_delta_Xp_01_overlay" || ch_name == "nc_delta_Xp_02_overlay" || ch_name == "nc_delta_Xp_03_overlay" || ch_name == "nc_delta_Xp_04_overlay"
            || ch_name == "nc_delta_Xp_05_overlay" || ch_name == "nc_delta_Xp_06_overlay" || ch_name == "nc_delta_Xp_07_overlay" || ch_name == "nc_delta_Xp_08_overlay"
            || ch_name == "nc_delta_Xp_09_overlay" || ch_name == "nc_delta_Xp_10_overlay" || ch_name == "nc_delta_Xp_11_overlay" || ch_name == "nc_delta_Xp_12_overlay"
            || ch_name == "nc_delta_Xp_13_overlay" || ch_name == "nc_delta_Xp_14_overlay" || ch_name == "nc_delta_Xp_15_overlay" || ch_name == "nc_delta_Xp_16_overlay"
            || ch_name == "nc_delta_Xp_17_overlay" || ch_name == "nc_delta_Xp_18_overlay" || ch_name == "nc_delta_Xp_19_overlay" || ch_name == "nc_delta_Xp_20_overlay"){
                  if (flag_FC && flag_ncdelta_sel && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // NC Pi0 channels (with NC delta selected events removed):
    }else if (ch_name == "nc_pi0_0p" || ch_name == "nc_pi0_2_0p"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np" || ch_name == "nc_pi0_2_Np"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp" || ch_name == "nc_pi0_2_Xp"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_ext" || ch_name == "nc_pi0_2_0p_ext"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_ext" || ch_name == "nc_pi0_2_Np_ext"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_ext" || ch_name == "nc_pi0_2_Xp_ext"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_dirt" || ch_name == "nc_pi0_2_0p_dirt"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_dirt" || ch_name == "nc_pi0_2_Np_dirt"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_dirt" || ch_name == "nc_pi0_2_Xp_dirt"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_nc_delta_overlay" || ch_name == "nc_pi0_0p_nc_delta_overlay_add" || ch_name == "nc_pi0_2_0p_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_nc_delta_overlay" || ch_name == "nc_pi0_Np_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Np_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_nc_delta_overlay" || ch_name == "nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "nc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_true_Np_nc_delta_overlay" || ch_name == "nc_pi0_Np_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_pi0_Np_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_pi0_Np_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_pi0_0p_true_Np_nc_delta_overlay" || ch_name == "nc_pi0_0p_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="nc_pi0_0p_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="nc_pi0_0p_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "nc_pi0_0p_nc_pi0_overlay" || ch_name == "nc_pi0_2_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_nc_pi0_overlay" || ch_name == "nc_pi0_2_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_nc_pi0_overlay" || ch_name == "nc_pi0_2_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_0p_overlay" || ch_name == "nc_pi0_2_0p_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Np_overlay" || ch_name == "nc_pi0_2_Np_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "nc_pi0_Xp_overlay" || ch_name == "nc_pi0_2_Xp_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // CC pi0 channels (with NC delta and NC pi0 selections removed) (not used for final NC delta fits):
     }else if (ch_name == "cc_pi0_0p" || ch_name == "cc_pi0_2_0p"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np" || ch_name == "cc_pi0_2_Np"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp" || ch_name == "cc_pi0_2_Xp"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_ext" || ch_name == "cc_pi0_2_0p_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_ext" || ch_name == "cc_pi0_2_Np_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_ext" || ch_name == "cc_pi0_2_Xp_ext"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_dirt" || ch_name == "cc_pi0_2_0p_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_dirt" || ch_name == "cc_pi0_2_Np_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_dirt" || ch_name == "cc_pi0_2_Xp_dirt"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_nc_delta_overlay" || ch_name == "cc_pi0_0p_nc_delta_overlay_add" || ch_name == "cc_pi0_2_0p_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_nc_delta_overlay" || ch_name == "cc_pi0_Np_nc_delta_overlay_add" || ch_name == "cc_pi0_2_Np_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_nc_delta_overlay" || ch_name == "cc_pi0_Xp_nc_delta_overlay_add" || ch_name == "cc_pi0_2_Xp_nc_delta_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_nc_pi0_overlay" || ch_name == "cc_pi0_2_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_nc_pi0_overlay" || ch_name == "cc_pi0_2_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_nc_pi0_overlay" || ch_name == "cc_pi0_2_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_0p_overlay" || ch_name == "cc_pi0_2_0p_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Np_overlay" || ch_name == "cc_pi0_2_Np_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "cc_pi0_Xp_overlay" || ch_name == "cc_pi0_2_Xp_overlay"){
                  if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // numuCC channels (with NC delta, NC pi0, and CC pi0 events removed) (used for some constraint tests, not used for final fits):
    }else if (ch_name == "numuCC_noCCpi0_0p"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)  && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_nc_delta_overlay" || ch_name == "numuCC_noCCpi0_0p_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_nc_delta_overlay" || ch_name == "numuCC_noCCpi0_Np_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_nc_delta_overlay" || ch_name == "numuCC_noCCpi0_Xp_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)  && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_0p_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Np_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_noCCpi0_Xp_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    // numuCC channels (with NC delta and NC pi0 events removed):
    }else if (ch_name == "numuCC_0p"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_ext"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_dirt"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel)) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_nc_delta_overlay" || ch_name == "numuCC_0p_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_nc_delta_overlay" || ch_name == "numuCC_Np_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_nc_delta_overlay" || ch_name == "numuCC_Xp_nc_delta_overlay_add"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_true_Np_nc_delta_overlay" || ch_name == "numuCC_Np_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="numuCC_Np_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="numuCC_Np_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "numuCC_0p_true_Np_nc_delta_overlay" || ch_name == "numuCC_0p_true_0p_nc_delta_overlay"){
            if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)){
                return ((is_true_0p(pfeval)==1 && ch_name=="numuCC_0p_true_0p_nc_delta_overlay") || (is_true_0p(pfeval)==0 && ch_name=="numuCC_0p_true_Np_nc_delta_overlay"));
            }
            return false;
    }else if (ch_name == "numuCC_0p_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_nc_pi0_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_overlay"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
                                                       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
          && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
                  return false;
    }else if (ch_name == "numuCC_0p_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && flag_0p) return true;
                  return false;
    }else if (ch_name == "numuCC_Np_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0) && (!flag_0p)) return true;
                  return false;
    }else if (ch_name == "numuCC_Xp_overlay_entire"){
                  if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && (!flag_cc_pi0)) return true;
                  return false;
    // Pi0 1p selections (used for proton-pi0 invariant mass test)
    }else if (ch_name == "nc_pi0_1p" || ch_name == "nc_pi0_1p_ext" || ch_name == "nc_pi0_1p_dirt" || ch_name == "nc_pi0_1p_overlay_entire" || ch_name == "nc_pi0_2_1p" || ch_name == "nc_pi0_2_1p_ext" || ch_name == "nc_pi0_2_1p_dirt" || ch_name == "nc_pi0_2_1p_overlay_entire"){
            if (flag_FC && flag_ncpio_sel && (!flag_ncdelta_sel) && flag_1p) return true;
                  return false;

    }else if (ch_name == "cc_pi0_1p" || ch_name == "cc_pi0_1p_ext" || ch_name == "cc_pi0_1p_dirt" || ch_name == "cc_pi0_1p_overlay_entire" || ch_name == "cc_pi0_2_1p" || ch_name == "cc_pi0_2_1p_ext" || ch_name == "cc_pi0_2_1p_dirt" || ch_name == "cc_pi0_2_1p_overlay_entire"){
            if (flag_FC && flag_cc_pi0 && (!flag_ncpio_sel) && (!flag_ncdelta_sel) && flag_1p) return true;
                  return false;


  // end cuts from Lee's cuts.h
  }else if (ch_name == "nc_delta_energy_FC_0p" || ch_name == "nc_delta_score_FC_0p"
	    || ch_name == "nc_delta_energy_FC_0p_ncpio_overlay" || ch_name == "nc_delta_score_FC_0p_ncpio_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_ncdelta_overlay" || ch_name == "nc_delta_score_FC_0p_ncdelta_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_overlay" || ch_name == "nc_delta_score_FC_0p_overlay"
	    || ch_name == "nc_delta_energy_FC_0p_ext" || ch_name == "nc_delta_score_FC_0p_ext"
	    || ch_name == "nc_delta_energy_FC_0p_dirt" || ch_name == "nc_delta_score_FC_0p_dirt"){

    if (ch_name == "nc_delta_energy_FC_0p" ||  ch_name == "nc_delta_energy_FC_0p_ext" || ch_name == "nc_delta_energy_FC_0p_dirt"){
      if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_ncpio_overlay"){
      if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
						     && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_ncdelta_overlay"){
      if (flag_FC && flag_ncdelta_sel && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_energy_FC_0p_overlay"){
      if (flag_FC && flag_ncdelta_sel && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC && flag_ncdelta_sel && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_FC_0p" ||  ch_name == "nc_delta_score_FC_0p_ext" || ch_name == "nc_delta_score_FC_0p_dirt"){
      if (flag_FC  && flag_0p) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_ncpio_overlay"){
      if (flag_FC  && flag_0p && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
				  && !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_ncdelta_overlay"){
      if (flag_FC  && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_FC_0p_overlay"){
      if (flag_FC  && flag_0p && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC  && flag_0p) return true;
    }

    return false;
  }else if (ch_name == "nc_delta_energy_FC_Np" || ch_name == "nc_delta_score_FC_Np"
	    || ch_name == "nc_delta_energy_FC_Np_ncpio_overlay" || ch_name == "nc_delta_score_FC_Np_ncpio_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_ncdelta_overlay" || ch_name == "nc_delta_score_FC_Np_ncdelta_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_overlay" || ch_name == "nc_delta_score_FC_Np_overlay"
	    || ch_name == "nc_delta_energy_FC_Np_ext" || ch_name == "nc_delta_score_FC_Np_ext"
	    || ch_name == "nc_delta_energy_FC_Np_dirt" || ch_name == "nc_delta_score_FC_Np_dirt"){

    if (ch_name == "nc_delta_energy_FC_Np" ||  ch_name == "nc_delta_energy_FC_Np_ext" || ch_name == "nc_delta_energy_FC_Np_dirt"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_ncpio_overlay"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_ncdelta_overlay"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_energy_FC_Np_overlay"){
      if (flag_FC && flag_ncdelta_sel && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      //      if (flag_FC && flag_ncdelta_sel && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np" ||  ch_name == "nc_delta_score_FC_Np_ext" || ch_name == "nc_delta_score_FC_Np_dirt"){
      if (flag_FC  && (!flag_0p)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_ncpio_overlay"){
      if (flag_FC  && (!flag_0p) && (eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      					&& !(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_ncdelta_overlay"){
      if (flag_FC  && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
    }else if (ch_name == "nc_delta_score_FC_Np_overlay"){
      if (flag_FC  && (!flag_0p) && (!(eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1
      						       && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))))
      	  && (!(eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside))) return true;
      // if (flag_FC  && (!flag_0p)) return true;
    }

    return false;

  //Erin
  }else if (ch_name == "all_bnb" || ch_name == "all_ext" || ch_name == "all_dirt" || ch_name == "all_bnb_LEE"){
    return true;
  }else if (ch_name == "all_bnb_nsbeam"){
    if (flag_nsbeam) return true;
    else return false;
  }else if (ch_name == "nodata_bnb"){
    return false;
  }else if (ch_name == "all_spoverlay" || ch_name == "all_spoverlay_2" || ch_name == "all_spoverlay_3"){
            if ((eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1)) return true;
            return false;
  }else if (ch_name == "all_ncpi0overlay" || ch_name == "all_ncpi0overlay_2" || ch_name == "all_ncpi0overlay_3"){
            if (!(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "all_overlay_sp_ncpi0_BG" || ch_name == "all_overlay_sp_ncpi0_BG_2" || ch_name == "all_overlay_sp_ncpi0_BG_3"){
            if (!(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "all_overlay_sp_BG"){
            if (!(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "generic_nu_ext" || ch_name == "generic_nu_dirt" ||ch_name == "generic_nu_bnb_LEE"){
    if (flag_generic) return true;
    else return false;
  }else if (ch_name == "generic_nu_bnb_nsbeam"){
    if (flag_generic && flag_nsbeam) return true;
    else return false;
  }else if (ch_name == "nodata_bnb"){
    return false;
  }else if (ch_name == "generic_nu_spoverlay" || ch_name == "generic_nu_spoverlay_2" || ch_name == "generic_nu_spoverlay_3"){
            if (flag_generic &&
              (eval.match_completeness_energy>0.1*eval.truth_energyInside && pfeval.truth_single_photon==1)) return true;
            return false;
  }else if (ch_name == "generic_nu_ncpi0overlay" || ch_name == "generic_nu_ncpi0overlay_2" || ch_name == "generic_nu_ncpi0overlay_3"){
            if (flag_generic &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "generic_nu_overlay_sp_ncpi0_BG" || ch_name == "generic_nu_overlay_sp_ncpi0_BG_2" || ch_name == "generic_nu_overlay_sp_ncpi0_BG_3"){
            if (flag_generic &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "generic_nu_overlay_sp_BG"){
            if (flag_generic &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb" || ch_name == "single_photon_ext"
    || ch_name == "single_photon_overlay" || ch_name == "single_photon_dirt" 
    || ch_name == "single_photon_LEE"){
            if (flag_singlephoton_sel) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb" || ch_name == "single_photon_eff_ext"
    || ch_name == "single_photon_eff_overlay" || ch_name == "single_photon_eff_dirt" 
    || ch_name == "single_photon_eff_LEE"){
            if (flag_singlephoton_eff_sel) return true;
            return false;
  }else if (ch_name == "single_shower_bnb" || ch_name == "single_shower_ext"
    || ch_name == "single_shower_overlay" || ch_name == "single_shower_dirt" 
    || ch_name == "single_shower_LEE"){
            if (flag_singleshower_sel) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb" || ch_name == "single_shower_eff_ext"
    || ch_name == "single_shower_eff_overlay" || ch_name == "single_shower_eff_dirt"
    || ch_name == "single_shower_eff_bnb_2" || ch_name == "single_shower_eff_ext_2"
      || ch_name == "single_shower_eff_overlay_2" || ch_name == "single_shower_eff_dirt_2" 
    || ch_name == "single_shower_eff_LEE"){
            if (flag_singleshower_eff_sel) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p" || ch_name == "single_photon_ext_0p"
    || ch_name == "single_photon_overlay_0p" || ch_name == "single_photon_dirt_0p" 
    || ch_name == "single_photon_LEE_0p"){
            if (flag_singlephoton_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p" || ch_name == "single_photon_eff_ext_0p"
    || ch_name == "single_photon_eff_overlay_0p" || ch_name == "single_photon_eff_dirt_0p" 
    || ch_name == "single_photon_eff_LEE_0p"){
            if (flag_singlephoton_eff_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p" || ch_name == "single_shower_ext_0p"
    || ch_name == "single_shower_overlay_0p" || ch_name == "single_shower_dirt_0p" 
    || ch_name == "single_shower_LEE_0p"){
            if (flag_singleshower_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p" || ch_name == "single_shower_eff_ext_0p"
    || ch_name == "single_shower_eff_overlay_0p" || ch_name == "single_shower_eff_dirt_0p" 
    || ch_name == "single_shower_eff_LEE_0p"){
            if (flag_singleshower_eff_sel && flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np" || ch_name == "single_photon_ext_Np"
    || ch_name == "single_photon_overlay_Np" || ch_name == "single_photon_dirt_Np"  
    || ch_name == "single_photon_LEE_Np"){
            if (flag_singlephoton_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np" || ch_name == "single_photon_eff_ext_Np"
    || ch_name == "single_photon_eff_overlay_Np" || ch_name == "single_photon_eff_dirt_Np" 
    || ch_name == "single_photon_eff_LEE_Np"){
            if (flag_singlephoton_eff_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np" || ch_name == "single_shower_ext_Np"
    || ch_name == "single_shower_overlay_Np" || ch_name == "single_shower_dirt_Np" 
    || ch_name == "single_shower_LEE_Np"){
            if (flag_singleshower_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np" || ch_name == "single_shower_eff_ext_Np"
    || ch_name == "single_shower_eff_overlay_Np" || ch_name == "single_shower_eff_dirt_Np" 
    || ch_name == "single_shower_eff_LEE_Np"){
            if (flag_singleshower_eff_sel && !flag_0p) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_BG"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_BG"){
            if (flag_singlephoton_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_BG"){
            if (flag_singleshower_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_BG"){
            if (flag_singleshower_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_0p_BG"){
            if (flag_singlephoton_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_0p_BG"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_0p_BG"){
            if (flag_singleshower_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_0p_BG"){
            if (flag_singleshower_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_Np_BG"){
            if (flag_singlephoton_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_Np_BG"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_Np_BG"){
            if (flag_singleshower_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_Np_BG"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_nodelta"){
            if (flag_singlephoton_sel &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_delta"){
            if (flag_singlephoton_sel &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay" || ch_name == "single_photon_spoverlay_lee"){
            if (flag_singlephoton_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_otpc" || ch_name == "single_photon_spoverlay_otpc_lee"){
            if (flag_singlephoton_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay"){
            if (flag_singlephoton_eff_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay" || ch_name == "single_shower_spoverlay_lee"){
            if (flag_singleshower_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay"){
            if (flag_singleshower_eff_sel &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_0p" || ch_name == "single_photon_spoverlay_lee_0p"){
            if (flag_singlephoton_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_0p"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_0p"){
            if (flag_singleshower_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_0p"){
            if (flag_singleshower_eff_sel && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_Np" || ch_name == "single_photon_spoverlay_lee_Np"){
            if (flag_singlephoton_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_Np"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_Np"){
            if (flag_singleshower_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_Np"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  //nc pi0 overlay
  }else if (ch_name == "single_photon_ncpi0overlay"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_otpc"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              //map_cuts_flag["SPNCPi0Sig"] || 
              map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              //&& !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              //&& pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              //&& pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_BG"){
            if (flag_singlephoton_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay"){
            if (flag_singlephoton_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_BG"){
            if (flag_singlephoton_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay"){
            if (flag_singleshower_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_BG"){
            if (flag_singleshower_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay"){
            if (flag_singleshower_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_BG"){
            if (flag_singleshower_eff_sel &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_0p"){
            if (flag_singlephoton_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_0p_BG"){
            if (flag_singlephoton_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_0p"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_0p_BG"){
            if (flag_singlephoton_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_0p"){
            if (flag_singleshower_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_0p_BG"){
            if (flag_singleshower_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_0p"){
            if (flag_singleshower_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_0p_BG"){
            if (flag_singleshower_eff_sel && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_Np"){
            if (flag_singlephoton_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_Np_BG"){
            if (flag_singlephoton_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_Np"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_Np_BG"){
            if (flag_singlephoton_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_Np"){
            if (flag_singleshower_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_Np_BG"){
            if (flag_singleshower_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_Np"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_Np_BG"){
            if (flag_singleshower_eff_sel && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_nsbeam"){
            if (flag_singlephoton_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_nsbeam"){
            if (flag_singlephoton_eff_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_nsbeam"){
            if (flag_singleshower_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_nsbeam"){
            if (flag_singleshower_eff_sel && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_nsbeam"){
            if (flag_singlephoton_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_nsbeam"){
            if (flag_singlephoton_eff_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_nsbeam"){
            if (flag_singleshower_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_nsbeam"){
            if (flag_singleshower_eff_sel && flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_nsbeam"){
            if (flag_singlephoton_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_nsbeam"){
            if (flag_singlephoton_eff_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_nsbeam"){
            if (flag_singleshower_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_nsbeam"){
            if (flag_singleshower_eff_sel && !flag_0p && flag_nsbeam_photon) return true;
            return false;
  //FC
  }else if (ch_name == "single_photon_bnb_FC" || ch_name == "single_photon_ext_FC"
    || ch_name == "single_photon_overlay_FC" || ch_name == "single_photon_dirt_FC"){
            if (flag_singlephoton_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_FC" || ch_name == "single_photon_eff_ext_FC"
    || ch_name == "single_photon_eff_overlay_FC" || ch_name == "single_photon_eff_dirt_FC"){
            if (flag_singlephoton_eff_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_FC" || ch_name == "single_shower_ext_FC"
    || ch_name == "single_shower_overlay_FC" || ch_name == "single_shower_dirt_FC"){
            if (flag_singleshower_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_FC" || ch_name == "single_shower_eff_ext_FC"
    || ch_name == "single_shower_eff_overlay_FC" || ch_name == "single_shower_eff_dirt_FC"
    || ch_name == "single_shower_eff_bnb_2_FC" || ch_name == "single_shower_eff_ext_2_FC"
      || ch_name == "single_shower_eff_overlay_2_FC" || ch_name == "single_shower_eff_dirt_2_FC"){
            if (flag_singleshower_eff_sel && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_FC" || ch_name == "single_photon_ext_0p_FC"
    || ch_name == "single_photon_overlay_0p_FC" || ch_name == "single_photon_dirt_0p_FC"){
            if (flag_singlephoton_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_FC" || ch_name == "single_photon_eff_ext_0p_FC"
    || ch_name == "single_photon_eff_overlay_0p_FC" || ch_name == "single_photon_eff_dirt_0p_FC"){
            if (flag_singlephoton_eff_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_FC" || ch_name == "single_shower_ext_0p_FC"
    || ch_name == "single_shower_overlay_0p_FC" || ch_name == "single_shower_dirt_0p_FC"){
            if (flag_singleshower_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_FC" || ch_name == "single_shower_eff_ext_0p_FC"
    || ch_name == "single_shower_eff_overlay_0p_FC" || ch_name == "single_shower_eff_dirt_0p_FC"){
            if (flag_singleshower_eff_sel && flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_FC" || ch_name == "single_photon_ext_Np_FC"
    || ch_name == "single_photon_overlay_Np_FC" || ch_name == "single_photon_dirt_Np_FC"){
            if (flag_singlephoton_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_FC" || ch_name == "single_photon_eff_ext_Np_FC"
    || ch_name == "single_photon_eff_overlay_Np_FC" || ch_name == "single_photon_eff_dirt_Np_FC"){
            if (flag_singlephoton_eff_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_FC" || ch_name == "single_shower_ext_Np_FC"
    || ch_name == "single_shower_overlay_Np_FC" || ch_name == "single_shower_dirt_Np_FC"){
            if (flag_singleshower_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_FC" || ch_name == "single_shower_eff_ext_Np_FC"
    || ch_name == "single_shower_eff_overlay_Np_FC" || ch_name == "single_shower_eff_dirt_Np_FC"){
            if (flag_singleshower_eff_sel && !flag_0p && flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_BG_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_BG_FC"){
            if (flag_singleshower_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_0p_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_0p_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_0p_BG_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_0p_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_Np_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_Np_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_Np_BG_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_Np_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_FC"){
            if (flag_singleshower_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_0p_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_0p_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_0p_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_0p_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_Np_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_Np_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_Np_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_Np_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  //nc pi0 overlay FC
  }else if (ch_name == "single_photon_ncpi0overlay_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_BG_FC"){
            if (flag_singlephoton_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_FC"){
            if (flag_singleshower_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_BG_FC"){
            if (flag_singleshower_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_0p_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_0p_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_0p_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_0p_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_0p_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_0p_BG_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_0p_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_0p_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_Np_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_Np_BG_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_Np_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_Np_BG_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_Np_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_Np_BG_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_Np_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_Np_BG_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_nsbeam_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_nsbeam_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_nsbeam_FC"){
            if (flag_singleshower_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_nsbeam_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_nsbeam_FC"){
            if (flag_singlephoton_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_nsbeam_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_nsbeam_FC"){
            if (flag_singleshower_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_nsbeam_FC"){
            if (flag_singleshower_eff_sel && flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_nsbeam_FC"){
            if (flag_singlephoton_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_nsbeam_FC"){
            if (flag_singlephoton_eff_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_nsbeam_FC"){
            if (flag_singleshower_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_nsbeam_FC"){
            if (flag_singleshower_eff_sel && flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  //PC
  }else if (ch_name == "single_photon_bnb_PC" || ch_name == "single_photon_ext_PC"
    || ch_name == "single_photon_overlay_PC" || ch_name == "single_photon_dirt_PC"){
            if (flag_singlephoton_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_PC" || ch_name == "single_photon_eff_ext_PC"
    || ch_name == "single_photon_eff_overlay_PC" || ch_name == "single_photon_eff_dirt_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_PC" || ch_name == "single_shower_ext_PC"
    || ch_name == "single_shower_overlay_PC" || ch_name == "single_shower_dirt_PC"){
            if (flag_singleshower_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_PC" || ch_name == "single_shower_eff_ext_PC"
    || ch_name == "single_shower_eff_overlay_PC" || ch_name == "single_shower_eff_dirt_PC"
    || ch_name == "single_shower_eff_bnb_2_PC" || ch_name == "single_shower_eff_ext_2_PC"
      || ch_name == "single_shower_eff_overlay_2_PC" || ch_name == "single_shower_eff_dirt_2_PC"){
            if (flag_singleshower_eff_sel && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_PC" || ch_name == "single_photon_ext_0p_PC"
    || ch_name == "single_photon_overlay_0p_PC" || ch_name == "single_photon_dirt_0p_PC"){
            if (flag_singlephoton_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_PC" || ch_name == "single_photon_eff_ext_0p_PC"
    || ch_name == "single_photon_eff_overlay_0p_PC" || ch_name == "single_photon_eff_dirt_0p_PC"){
            if (flag_singlephoton_eff_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_PC" || ch_name == "single_shower_ext_0p_PC"
    || ch_name == "single_shower_overlay_0p_PC" || ch_name == "single_shower_dirt_0p_PC"){
            if (flag_singleshower_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_PC" || ch_name == "single_shower_eff_ext_0p_PC"
    || ch_name == "single_shower_eff_overlay_0p_PC" || ch_name == "single_shower_eff_dirt_0p_PC"){
            if (flag_singleshower_eff_sel && flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_PC" || ch_name == "single_photon_ext_Np_PC"
    || ch_name == "single_photon_overlay_Np_PC" || ch_name == "single_photon_dirt_Np_PC"){
            if (flag_singlephoton_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_PC" || ch_name == "single_photon_eff_ext_Np_PC"
    || ch_name == "single_photon_eff_overlay_Np_PC" || ch_name == "single_photon_eff_dirt_Np_PC"){
            if (flag_singlephoton_eff_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_PC" || ch_name == "single_shower_ext_Np_PC"
    || ch_name == "single_shower_overlay_Np_PC" || ch_name == "single_shower_dirt_Np_PC"){
            if (flag_singleshower_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_PC" || ch_name == "single_shower_eff_ext_Np_PC"
    || ch_name == "single_shower_eff_overlay_Np_PC" || ch_name == "single_shower_eff_dirt_Np_PC"){
            if (flag_singleshower_eff_sel && !flag_0p && !flag_FC) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_BG_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_0p_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_0p_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_0p_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_0p_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_Np_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_Np_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_Np_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_Np_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_0p_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_0p_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_0p_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_0p_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_spoverlay_Np_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_photon_eff_spoverlay_Np_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_spoverlay_Np_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "single_shower_eff_spoverlay_Np_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  //nc pi0 overlay PC
  }else if (ch_name == "single_photon_ncpi0overlay_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_BG_PC"){
            if (flag_singleshower_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_0p_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_0p_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_0p_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_0p_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_0p_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_0p_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_0p_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_0p_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_ncpi0overlay_Np_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_overlay_ncpi0_Np_BG_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_ncpi0overlay_Np_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_eff_overlay_ncpi0_Np_BG_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_ncpi0overlay_Np_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_overlay_ncpi0_Np_BG_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_ncpi0overlay_Np_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_shower_eff_overlay_ncpi0_Np_BG_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_nsbeam_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_nsbeam_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_nsbeam_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_nsbeam_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_0p_nsbeam_PC"){
            if (flag_singlephoton_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_0p_nsbeam_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_0p_nsbeam_PC"){
            if (flag_singleshower_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_0p_nsbeam_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_bnb_Np_nsbeam_PC"){
            if (flag_singlephoton_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_photon_eff_bnb_Np_nsbeam_PC"){
            if (flag_singlephoton_eff_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_bnb_Np_nsbeam_PC"){
            if (flag_singleshower_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  }else if (ch_name == "single_shower_eff_bnb_Np_nsbeam_PC"){
            if (flag_singleshower_eff_sel && !flag_FC && !flag_0p && flag_nsbeam) return true;
            return false;
  //nue bdt cut eff
  }else if (ch_name == "single_photon_eff_nue_overlay" ||
            ch_name == "single_photon_eff_nue_dirt" ||
            ch_name == "single_photon_eff_nue_ext" ||
            ch_name == "single_photon_eff_nue_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_nue_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 && flag_singlephoton_eff_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  //ncpi0 bdt cut
  }else if (ch_name == "single_photon_eff_ncpi0_overlay" ||
            ch_name == "single_photon_eff_ncpi0_dirt" ||
            ch_name == "single_photon_eff_ncpi0_ext" ||
            ch_name == "single_photon_eff_ncpi0_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_ncpi0_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
       flag_singlephoton_eff_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_overlay" ||
            ch_name == "single_photon_eff_other_dirt" ||
            ch_name == "single_photon_eff_other_ext" ||
            ch_name == "single_photon_eff_other_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_other_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu && flag_singlephoton_eff_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_overlay" ||
            ch_name == "single_photon_eff_numu_dirt" ||
            ch_name == "single_photon_eff_numu_ext" ||
            ch_name == "single_photon_eff_numu_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_eff_numu_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_eff_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
    }else if (ch_name == "single_photon_eff_pre_overlay" ||
              ch_name == "single_photon_eff_pre_dirt" ||
              ch_name == "single_photon_eff_pre_ext" ||
              ch_name == "single_photon_eff_pre_bnb"){
      if (flag_singlephoton_pre) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_bnb_nsbeam"){
      if (flag_singlephoton_pre && flag_nsbeam_photon) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_spoverlay"){
      if (flag_singlephoton_pre &&
          (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
          map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
          map_cuts_flag["SPNumuCCSig"])) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_overlay_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
           map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
           map_cuts_flag["SPNumuCCSig"])
           && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
           && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
           && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_ncpi0overlay"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;
    }else if (ch_name == "single_photon_eff_pre_overlay_ncpi0_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;

  //nue bdt cut pure
  }else if (ch_name == "single_photon_nue_overlay" ||
            ch_name == "single_photon_nue_dirt" ||
            ch_name == "single_photon_nue_ext" ||
            ch_name == "single_photon_nue_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue) return true;
      return false;
  }else if (ch_name == "single_photon_nue_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_nue_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_nue_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_nue_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_nue_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_singlephoton_nue &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  //ncpi0 bdt cut
  }else if (ch_name == "single_photon_ncpi0_overlay" ||
            ch_name == "single_photon_ncpi0_dirt" ||
            ch_name == "single_photon_ncpi0_ext" ||
            ch_name == "single_photon_ncpi0_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_ncpi0_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
       flag_singlephoton_ncpi0 &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_other_overlay" ||
            ch_name == "single_photon_other_dirt" ||
            ch_name == "single_photon_other_ext" ||
            ch_name == "single_photon_other_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other) return true;
      return false;
  }else if (ch_name == "single_photon_other_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_other_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_other_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_other_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_other_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_singlephoton_other &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_numu_overlay" ||
            ch_name == "single_photon_numu_dirt" ||
            ch_name == "single_photon_numu_ext" ||
            ch_name == "single_photon_numu_bnb"){
    if (flag_singlephoton_pre && flag_singlephoton_numu) return true;
      return false;
  }else if (ch_name == "single_photon_numu_bnb_nsbeam"){
    if (flag_singlephoton_pre && flag_singlephoton_numu && flag_nsbeam_photon) return true;
      return false;
  }else if (ch_name == "single_photon_numu_spoverlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
        (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
        map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
        map_cuts_flag["SPNumuCCSig"])) return true;
      return false;
  }else if (ch_name == "single_photon_numu_overlay_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
         map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
         map_cuts_flag["SPNumuCCSig"])
         && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
         && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
         && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
      return false;
  }else if (ch_name == "single_photon_numu_ncpi0overlay"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
  }else if (ch_name == "single_photon_numu_overlay_ncpi0_BG"){
    if (flag_singlephoton_pre && flag_singlephoton_numu &&
         !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
      return false;
    }else if (ch_name == "single_photon_pre_overlay" ||
              ch_name == "single_photon_pre_dirt" ||
              ch_name == "single_photon_pre_ext" ||
              ch_name == "single_photon_pre_bnb"){
      if (flag_singlephoton_pre) return true;
        return false;
    }else if (ch_name == "single_photon_pre_bnb_nsbeam"){
      if (flag_singlephoton_pre && flag_nsbeam_photon) return true;
        return false;
    }else if (ch_name == "single_photon_pre_spoverlay"){
      if (flag_singlephoton_pre &&
          (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
          map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
          map_cuts_flag["SPNumuCCSig"])) return true;
        return false;
    }else if (ch_name == "single_photon_pre_overlay_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
           map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
           map_cuts_flag["SPNumuCCSig"])
           && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
           && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
           && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
        return false;
    }else if (ch_name == "single_photon_pre_ncpi0overlay"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;
    }else if (ch_name == "single_photon_pre_overlay_ncpi0_BG"){
      if (flag_singlephoton_pre &&
           !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
        return false;

  }else if (ch_name == "single_photon_ncdel_sig_overlay" || ch_name == "single_photon_ncpi0_sig_overlay" ||
    ch_name == "single_photon_ncother_sig_overlay" || ch_name == "single_photon_numucc_sig_overlay" ||
    ch_name == "single_photon_outfv_sig_overlay" || ch_name == "single_photon_bkg_overlay"){
      if (ch_name == "single_photon_ncdel_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNCDeltaSig"]) return true;
      }else if (ch_name == "single_photon_ncpi0_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNCPi0Sig"]) return true;
      }else if (ch_name == "single_photon_ncother_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNCOtherSig"]) return true;
      }else if (ch_name == "single_photon_numucc_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPNumuCCSig"]) return true;
      }else if (ch_name == "single_photon_outfv_sig_overlay"){
        if (flag_singlephoton_sel && map_cuts_flag["SPOutFVSig"]) return true;
      }else if (ch_name == "single_photon_bkg_overlay"){
        if (flag_singlephoton_sel &&
          !map_cuts_flag["SPNCDeltaSig"] && !map_cuts_flag["SPOutFVSig"] &&
          !map_cuts_flag["SPNCPi0Sig"] && !map_cuts_flag["SPNCOtherSig"] &&
          !map_cuts_flag["SPNumuCCSig"]) return true;
      }
      return false;
  // NC Pi0 channels (with single photon selected events removed):
}else if (ch_name == "sp_nc_pi0_0p" || ch_name == "sp_nc_pi0_2_0p" || ch_name == "sp_nc_pi0_3_0p" || ch_name == "sp_nc_pi0_4_0p"
           || ch_name == "sp_nc_pi0_5_0p" || ch_name == "sp_nc_pi0_6_0p" || ch_name == "sp_nc_pi0_7_0p" || ch_name == "sp_nc_pi0_8_0p"
           || ch_name == "sp_nc_pi0_9_0p" || ch_name == "sp_nc_pi0_10_0p" || ch_name == "sp_nc_pi0_11_0p" || ch_name == "sp_nc_pi0_12_0p"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np" || ch_name == "sp_nc_pi0_2_Np" || ch_name == "sp_nc_pi0_3_Np" || ch_name == "sp_nc_pi0_4_Np"
            || ch_name == "sp_nc_pi0_5_Np" || ch_name == "sp_nc_pi0_6_Np" || ch_name == "sp_nc_pi0_7_Np" || ch_name == "sp_nc_pi0_8_Np"
            || ch_name == "sp_nc_pi0_9_Np" || ch_name == "sp_nc_pi0_10_Np" || ch_name == "sp_nc_pi0_11_Np" || ch_name == "sp_nc_pi0_12_Np"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp" || ch_name == "sp_nc_pi0_2_Xp" || ch_name == "sp_nc_pi0_3_Xp" || ch_name == "sp_nc_pi0_4_Xp"
            || ch_name == "sp_nc_pi0_5_Xp" || ch_name == "sp_nc_pi0_6_Xp" || ch_name == "sp_nc_pi0_7_Xp" || ch_name == "sp_nc_pi0_8_Xp"
            || ch_name == "sp_nc_pi0_9_Xp" || ch_name == "sp_nc_pi0_10_Xp" || ch_name == "sp_nc_pi0_11_Xp" || ch_name == "sp_nc_pi0_12_Xp"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_nsbeam" || ch_name == "sp_nc_pi0_2_0p_nsbeam" || ch_name == "sp_nc_pi0_3_0p_nsbeam" || ch_name == "sp_nc_pi0_4_0p_nsbeam"
            || ch_name == "sp_nc_pi0_5_0p_nsbeam" || ch_name == "sp_nc_pi0_6_0p_nsbeam" || ch_name == "sp_nc_pi0_7_0p_nsbeam" || ch_name == "sp_nc_pi0_8_0p_nsbeam"
            || ch_name == "sp_nc_pi0_9_0p_nsbeam" || ch_name == "sp_nc_pi0_10_0p_nsbeam" || ch_name == "sp_nc_pi0_11_0p_nsbeam" || ch_name == "sp_nc_pi0_12_0p_nsbeam"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_nsbeam" || ch_name == "sp_nc_pi0_2_Np_nsbeam" || ch_name == "sp_nc_pi0_3_Np_nsbeam" || ch_name == "sp_nc_pi0_4_Np_nsbeam"
            || ch_name == "sp_nc_pi0_5_Np_nsbeam" || ch_name == "sp_nc_pi0_6_Np_nsbeam" || ch_name == "sp_nc_pi0_7_Np_nsbeam" || ch_name == "sp_nc_pi0_8_Np_nsbeam"
            || ch_name == "sp_nc_pi0_9_Np_nsbeam" || ch_name == "sp_nc_pi0_10_Np_nsbeam" || ch_name == "sp_nc_pi0_11_Np_nsbeam" || ch_name == "sp_nc_pi0_12_Np_nsbeam"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_nsbeam" || ch_name == "sp_nc_pi0_2_Xp_nsbeam" || ch_name == "sp_nc_pi0_3_Xp_nsbeam" || ch_name == "sp_nc_pi0_4_Xp_nsbeam"
            || ch_name == "sp_nc_pi0_5_Xp_nsbeam" || ch_name == "sp_nc_pi0_6_Xp_nsbeam" || ch_name == "sp_nc_pi0_7_Xp_nsbeam" || ch_name == "sp_nc_pi0_8_Xp_nsbeam"
            || ch_name == "sp_nc_pi0_9_Xp_nsbeam" || ch_name == "sp_nc_pi0_10_Xp_nsbeam" || ch_name == "sp_nc_pi0_11_Xp_nsbeam" || ch_name == "sp_nc_pi0_12_Xp_nsbeam"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_ext" || ch_name == "sp_nc_pi0_2_0p_ext" || ch_name == "sp_nc_pi0_3_0p_ext" || ch_name == "sp_nc_pi0_4_0p_ext"
            || ch_name == "sp_nc_pi0_5_0p_ext" || ch_name == "sp_nc_pi0_6_0p_ext" || ch_name == "sp_nc_pi0_7_0p_ext" || ch_name == "sp_nc_pi0_8_0p_ext"
            || ch_name == "sp_nc_pi0_9_0p_ext" || ch_name == "sp_nc_pi0_10_0p_ext" || ch_name == "sp_nc_pi0_11_0p_ext" || ch_name == "sp_nc_pi0_12_0p_ext"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_ext" || ch_name == "sp_nc_pi0_2_Np_ext" || ch_name == "sp_nc_pi0_3_Np_ext" || ch_name == "sp_nc_pi0_4_Np_ext"
            || ch_name == "sp_nc_pi0_5_Np_ext" || ch_name == "sp_nc_pi0_6_Np_ext" || ch_name == "sp_nc_pi0_7_Np_ext" || ch_name == "sp_nc_pi0_8_Np_ext"
            || ch_name == "sp_nc_pi0_9_Np_ext" || ch_name == "sp_nc_pi0_10_Np_ext" || ch_name == "sp_nc_pi0_11_Np_ext" || ch_name == "sp_nc_pi0_12_Np_ext"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_ext" || ch_name == "sp_nc_pi0_2_Xp_ext" || ch_name == "sp_nc_pi0_3_Xp_ext" || ch_name == "sp_nc_pi0_4_Xp_ext"
            || ch_name == "sp_nc_pi0_5_Xp_ext" || ch_name == "sp_nc_pi0_6_Xp_ext" || ch_name == "sp_nc_pi0_7_Xp_ext" || ch_name == "sp_nc_pi0_8_Xp_ext"
            || ch_name == "sp_nc_pi0_9_Xp_ext" || ch_name == "sp_nc_pi0_10_Xp_ext" || ch_name == "sp_nc_pi0_11_Xp_ext" || ch_name == "sp_nc_pi0_12_Xp_ext"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_dirt" || ch_name == "sp_nc_pi0_2_0p_dirt" || ch_name == "sp_nc_pi0_3_0p_dirt" || ch_name == "sp_nc_pi0_4_0p_dirt"
            || ch_name == "sp_nc_pi0_5_0p_dirt" || ch_name == "sp_nc_pi0_6_0p_dirt" || ch_name == "sp_nc_pi0_7_0p_dirt" || ch_name == "sp_nc_pi0_8_0p_dirt"
            || ch_name == "sp_nc_pi0_9_0p_dirt" || ch_name == "sp_nc_pi0_10_0p_dirt" || ch_name == "sp_nc_pi0_11_0p_dirt" || ch_name == "sp_nc_pi0_12_0p_dirt"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_dirt" || ch_name == "sp_nc_pi0_2_Np_dirt" || ch_name == "sp_nc_pi0_3_Np_dirt" || ch_name == "sp_nc_pi0_4_Np_dirt"
            || ch_name == "sp_nc_pi0_5_Np_dirt" || ch_name == "sp_nc_pi0_6_Np_dirt" || ch_name == "sp_nc_pi0_7_Np_dirt" || ch_name == "sp_nc_pi0_8_Np_dirt"
            || ch_name == "sp_nc_pi0_9_Np_dirt" || ch_name == "sp_nc_pi0_10_Np_dirt" || ch_name == "sp_nc_pi0_11_Np_dirt" || ch_name == "sp_nc_pi0_12_Np_dirt"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_dirt" || ch_name == "sp_nc_pi0_2_Xp_dirt" || ch_name == "sp_nc_pi0_3_Xp_dirt" || ch_name == "sp_nc_pi0_4_Xp_dirt"
            || ch_name == "sp_nc_pi0_5_Xp_dirt" || ch_name == "sp_nc_pi0_6_Xp_dirt" || ch_name == "sp_nc_pi0_7_Xp_dirt" || ch_name == "sp_nc_pi0_8_Xp_dirt"
            || ch_name == "sp_nc_pi0_9_Xp_dirt" || ch_name == "sp_nc_pi0_10_Xp_dirt" || ch_name == "sp_nc_pi0_11_Xp_dirt" || ch_name == "sp_nc_pi0_12_Xp_dirt"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_0p_nc_delta_overlay" || ch_name == "sp_nc_pi0_0p_nc_delta_overlay_add" || ch_name == "sp_nc_pi0_2_0p_nc_delta_overlay"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Np_nc_delta_overlay" || ch_name == "sp_nc_pi0_Np_nc_delta_overlay_add" || ch_name == "sp_nc_pi0_2_Np_nc_delta_overlay"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_nc_delta_overlay" || ch_name == "sp_nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "sp_nc_pi0_2_Xp_nc_delta_overlay"){
                if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_nc_pi0_Xp_overlay" || ch_name == "sp_nc_pi0_2_Xp_overlay" || ch_name == "sp_nc_pi0_3_Xp_overlay" || ch_name == "sp_nc_pi0_4_Xp_overlay"
            || ch_name == "sp_nc_pi0_5_Xp_overlay" || ch_name == "sp_nc_pi0_6_Xp_overlay" || ch_name == "sp_nc_pi0_7_Xp_overlay" || ch_name == "sp_nc_pi0_8_Xp_overlay"
            || ch_name == "sp_nc_pi0_9_Xp_overlay" || ch_name == "sp_nc_pi0_10_Xp_overlay" || ch_name == "sp_nc_pi0_11_Xp_overlay" || ch_name == "sp_nc_pi0_12_Xp_overlay"){
                  if (flag_FC && flag_ncpio_sel && (!flag_singlephoton_sel))
                      return true;
                  return false;
  //nc pi0 sp sideband
}else if (ch_name == "sp_bdt_nc_pi0_0p" || ch_name == "sp_bdt_nc_pi0_2_0p" || ch_name == "sp_bdt_nc_pi0_3_0p" || ch_name == "sp_bdt_nc_pi0_4_0p"
          || ch_name == "sp_bdt_nc_pi0_5_0p" || ch_name == "sp_bdt_nc_pi0_6_0p" || ch_name == "sp_bdt_nc_pi0_7_0p" || ch_name == "sp_bdt_nc_pi0_8_0p"
          || ch_name == "sp_bdt_nc_pi0_9_0p" || ch_name == "sp_bdt_nc_pi0_10_0p" || ch_name == "sp_bdt_nc_pi0_11_0p" || ch_name == "sp_bdt_nc_pi0_12_0p"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np" || ch_name == "sp_bdt_nc_pi0_2_Np" || ch_name == "sp_bdt_nc_pi0_3_Np" || ch_name == "sp_bdt_nc_pi0_4_Np"
            || ch_name == "sp_bdt_nc_pi0_5_Np" || ch_name == "sp_bdt_nc_pi0_6_Np" || ch_name == "sp_bdt_nc_pi0_7_Np" || ch_name == "sp_bdt_nc_pi0_8_Np"
            || ch_name == "sp_bdt_nc_pi0_9_Np" || ch_name == "sp_bdt_nc_pi0_10_Np" || ch_name == "sp_bdt_nc_pi0_11_Np" || ch_name == "sp_bdt_nc_pi0_12_Np"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp" || ch_name == "sp_bdt_nc_pi0_2_Xp" || ch_name == "sp_bdt_nc_pi0_3_Xp" || ch_name == "sp_bdt_nc_pi0_4_Xp"
            || ch_name == "sp_bdt_nc_pi0_5_Xp" || ch_name == "sp_bdt_nc_pi0_6_Xp" || ch_name == "sp_bdt_nc_pi0_7_Xp" || ch_name == "sp_bdt_nc_pi0_8_Xp"
            || ch_name == "sp_bdt_nc_pi0_9_Xp" || ch_name == "sp_bdt_nc_pi0_10_Xp" || ch_name == "sp_bdt_nc_pi0_11_Xp" || ch_name == "sp_bdt_nc_pi0_12_Xp"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_2_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_3_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_4_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_6_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_7_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_8_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_10_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_11_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_12_Xp_nsbeam"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_ext" || ch_name == "sp_bdt_nc_pi0_2_0p_ext" || ch_name == "sp_bdt_nc_pi0_3_0p_ext" || ch_name == "sp_bdt_nc_pi0_4_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_5_0p_ext" || ch_name == "sp_bdt_nc_pi0_6_0p_ext" || ch_name == "sp_bdt_nc_pi0_7_0p_ext" || ch_name == "sp_bdt_nc_pi0_8_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_9_0p_ext" || ch_name == "sp_bdt_nc_pi0_10_0p_ext" || ch_name == "sp_bdt_nc_pi0_11_0p_ext" || ch_name == "sp_bdt_nc_pi0_12_0p_ext"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_ext" || ch_name == "sp_bdt_nc_pi0_2_Np_ext" || ch_name == "sp_bdt_nc_pi0_3_Np_ext" || ch_name == "sp_bdt_nc_pi0_4_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_5_Np_ext" || ch_name == "sp_bdt_nc_pi0_6_Np_ext" || ch_name == "sp_bdt_nc_pi0_7_Np_ext" || ch_name == "sp_bdt_nc_pi0_8_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_9_Np_ext" || ch_name == "sp_bdt_nc_pi0_10_Np_ext" || ch_name == "sp_bdt_nc_pi0_11_Np_ext" || ch_name == "sp_bdt_nc_pi0_12_Np_ext"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_ext" || ch_name == "sp_bdt_nc_pi0_2_Xp_ext" || ch_name == "sp_bdt_nc_pi0_3_Xp_ext" || ch_name == "sp_bdt_nc_pi0_4_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_ext" || ch_name == "sp_bdt_nc_pi0_6_Xp_ext" || ch_name == "sp_bdt_nc_pi0_7_Xp_ext" || ch_name == "sp_bdt_nc_pi0_8_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_ext" || ch_name == "sp_bdt_nc_pi0_10_Xp_ext" || ch_name == "sp_bdt_nc_pi0_11_Xp_ext" || ch_name == "sp_bdt_nc_pi0_12_Xp_ext"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_dirt" || ch_name == "sp_bdt_nc_pi0_2_0p_dirt" || ch_name == "sp_bdt_nc_pi0_3_0p_dirt" || ch_name == "sp_bdt_nc_pi0_4_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_5_0p_dirt" || ch_name == "sp_bdt_nc_pi0_6_0p_dirt" || ch_name == "sp_bdt_nc_pi0_7_0p_dirt" || ch_name == "sp_bdt_nc_pi0_8_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_9_0p_dirt" || ch_name == "sp_bdt_nc_pi0_10_0p_dirt" || ch_name == "sp_bdt_nc_pi0_11_0p_dirt" || ch_name == "sp_bdt_nc_pi0_12_0p_dirt"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_dirt" || ch_name == "sp_bdt_nc_pi0_2_Np_dirt" || ch_name == "sp_bdt_nc_pi0_3_Np_dirt" || ch_name == "sp_bdt_nc_pi0_4_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_5_Np_dirt" || ch_name == "sp_bdt_nc_pi0_6_Np_dirt" || ch_name == "sp_bdt_nc_pi0_7_Np_dirt" || ch_name == "sp_bdt_nc_pi0_8_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_9_Np_dirt" || ch_name == "sp_bdt_nc_pi0_10_Np_dirt" || ch_name == "sp_bdt_nc_pi0_11_Np_dirt" || ch_name == "sp_bdt_nc_pi0_12_Np_dirt"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_2_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_3_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_4_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_6_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_7_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_8_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_10_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_11_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_12_Xp_dirt"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_overlay" || ch_name == "sp_bdt_nc_pi0_2_0p_overlay" || ch_name == "sp_bdt_nc_pi0_3_0p_overlay" || ch_name == "sp_bdt_nc_pi0_4_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_5_0p_overlay" || ch_name == "sp_bdt_nc_pi0_6_0p_overlay" || ch_name == "sp_bdt_nc_pi0_7_0p_overlay" || ch_name == "sp_bdt_nc_pi0_8_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_9_0p_overlay" || ch_name == "sp_bdt_nc_pi0_10_0p_overlay" || ch_name == "sp_bdt_nc_pi0_11_0p_overlay" || ch_name == "sp_bdt_nc_pi0_12_0p_overlay"){
                  if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_overlay" || ch_name == "sp_bdt_nc_pi0_2_Np_overlay" || ch_name == "sp_bdt_nc_pi0_3_Np_overlay" || ch_name == "sp_bdt_nc_pi0_4_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_5_Np_overlay" || ch_name == "sp_bdt_nc_pi0_6_Np_overlay" || ch_name == "sp_bdt_nc_pi0_7_Np_overlay" || ch_name == "sp_bdt_nc_pi0_8_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_9_Np_overlay" || ch_name == "sp_bdt_nc_pi0_10_Np_overlay" || ch_name == "sp_bdt_nc_pi0_11_Np_overlay" || ch_name == "sp_bdt_nc_pi0_12_Np_overlay"){
                  if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_2_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_3_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_4_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_5_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_6_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_7_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_8_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_9_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_10_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_11_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_12_Xp_overlay"){
                  if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel))
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_2_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_3_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_4_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_5_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_6_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_7_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_8_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_9_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_10_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_11_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_12_0p_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_2_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_3_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_4_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_6_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_7_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_8_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_10_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_11_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_12_Np_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_2_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_3_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_6_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_7_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_10_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_11_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_12_Xp_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_2_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_3_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_4_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_5_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_6_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_7_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_8_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_9_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_10_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_11_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_12_0p_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_2_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_3_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_4_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_5_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_6_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_7_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_8_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_9_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_10_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_11_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_12_Np_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_2_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_3_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_4_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_6_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_7_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_8_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_10_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_11_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_12_Xp_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_spoverlay_nodelta"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_spoverlay_delta"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_2_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_3_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_4_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_5_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_6_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_7_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_8_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_9_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_10_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_11_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_12_0p_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_4_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_8_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_12_0p_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_2_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_3_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_5_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_6_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_7_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_9_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_10_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_11_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_12_Np_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_12_Np_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_12_Xp_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_12_Xp_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;

  //nc pi0 sp sideband with one shower
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p"
          || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p"
          || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_nsbeam"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_ext"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_ext"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_ext" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_ext"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_dirt"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_dirt"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_dirt"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_overlay"){
                  if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_overlay"){
                  if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_overlay"){
                  if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel))
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_spoverlay_nodelta"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_spoverlay_delta"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_oneshw_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_0p_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Np_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_oneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_oneshw_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_oneshw_12_Xp_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;

  //nc pi0 sp sideband with not one shower
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p"
          || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p"
          || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_nsbeam"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_nsbeam" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_nsbeam"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_ext"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_ext"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_ext"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_ext"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_ext" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_ext"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_dirt"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_dirt"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_dirt"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_dirt" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_dirt"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_overlay"){
                  if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_overlay"){
                  if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p)
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_overlay"){
                  if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel))
                      return true;
                  return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_overlay_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_overlay_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_overlay_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_spoverlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_spoverlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_spoverlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_spoverlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_spoverlay_nodelta"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              (map_cuts_flag["SPOutFVSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_spoverlay_delta"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              map_cuts_flag["SPNCDeltaSig"]) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_ncpi0overlay"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_0p_overlay_ncpi0_BG"
            || ch_name == "sp_bdt_nc_pi0_notoneshw_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_0p_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Np_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) && !flag_0p &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_ncpi0overlay"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;
  }else if (ch_name == "sp_bdt_nc_pi0_notoneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_nc_pi0_notoneshw_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nc_pi0_notoneshw_12_Xp_overlay_ncpi0_BG"){
            if (flag_singlephoton_ncpi0_sel && !flag_singlephoton_oneshw && (!flag_singlephoton_sel) &&
              !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
            return false;

  // numuCC channels (with single photon and NC pi0 events removed):
}else if (ch_name == "sp_numuCC_0p" || ch_name == "sp_numuCC_2_0p" || ch_name == "sp_numuCC_3_0p" || ch_name == "sp_numuCC_4_0p"
          || ch_name == "sp_numuCC_5_0p" || ch_name == "sp_numuCC_6_0p" || ch_name == "sp_numuCC_7_0p" || ch_name == "sp_numuCC_8_0p"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np" || ch_name == "sp_numuCC_2_Np"|| ch_name == "sp_numuCC_3_Np" || ch_name == "sp_numuCC_4_Np"
            || ch_name == "sp_numuCC_5_Np" || ch_name == "sp_numuCC_6_Np" || ch_name == "sp_numuCC_7_Np" || ch_name == "sp_numuCC_8_Np"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp" || ch_name == "sp_numuCC_2_Xp"|| ch_name == "sp_numuCC_3_Xp" || ch_name == "sp_numuCC_4_Xp"
            || ch_name == "sp_numuCC_5_Xp" || ch_name == "sp_numuCC_6_Xp" || ch_name == "sp_numuCC_7_Xp" || ch_name == "sp_numuCC_8_Xp"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_nsbeam" || ch_name == "sp_numuCC_2_0p_nsbeam" || ch_name == "sp_numuCC_3_0p_nsbeam" || ch_name == "sp_numuCC_4_0p_nsbeam"
             || ch_name == "sp_numuCC_5_0p_nsbeam" || ch_name == "sp_numuCC_6_0p_nsbeam" || ch_name == "sp_numuCC_7_0p_nsbeam" || ch_name == "sp_numuCC_8_0p_nsbeam"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_nsbeam" || ch_name == "sp_numuCC_2_Np_nsbeam" || ch_name == "sp_numuCC_3_Np_nsbeam" || ch_name == "sp_numuCC_4_Np_nsbeam"
             || ch_name == "sp_numuCC_5_Np_nsbeam" || ch_name == "sp_numuCC_6_Np_nsbeam" || ch_name == "sp_numuCC_7_Np_nsbeam" || ch_name == "sp_numuCC_8_Np_nsbeam"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_nsbeam" || ch_name == "sp_numuCC_2_Xp_nsbeam" || ch_name == "sp_numuCC_3_Xp_nsbeam" || ch_name == "sp_numuCC_4_Xp_nsbeam"
             || ch_name == "sp_numuCC_5_Xp_nsbeam" || ch_name == "sp_numuCC_6_Xp_nsbeam" || ch_name == "sp_numuCC_7_Xp_nsbeam" || ch_name == "sp_numuCC_8_Xp_nsbeam"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_ext" || ch_name == "sp_numuCC_2_0p_ext" || ch_name == "sp_numuCC_3_0p_ext" || ch_name == "sp_numuCC_4_0p_ext"
             || ch_name == "sp_numuCC_5_0p_ext" || ch_name == "sp_numuCC_6_0p_ext" || ch_name == "sp_numuCC_7_0p_ext" || ch_name == "sp_numuCC_8_0p_ext"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_ext" || ch_name == "sp_numuCC_2_Np_ext" || ch_name == "sp_numuCC_3_Np_ext" || ch_name == "sp_numuCC_4_Np_ext"
             || ch_name == "sp_numuCC_5_Np_ext" || ch_name == "sp_numuCC_6_Np_ext" || ch_name == "sp_numuCC_7_Np_ext" || ch_name == "sp_numuCC_8_Np_ext"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_ext" || ch_name == "sp_numuCC_2_Xp_ext" || ch_name == "sp_numuCC_3_Xp_ext" || ch_name == "sp_numuCC_4_Xp_ext"
             || ch_name == "sp_numuCC_5_Xp_ext" || ch_name == "sp_numuCC_6_Xp_ext" || ch_name == "sp_numuCC_7_Xp_ext" || ch_name == "sp_numuCC_8_Xp_ext"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_dirt" || ch_name == "sp_numuCC_2_0p_dirt" || ch_name == "sp_numuCC_3_0p_dirt" || ch_name == "sp_numuCC_4_0p_dirt"
             || ch_name == "sp_numuCC_5_0p_dirt" || ch_name == "sp_numuCC_6_0p_dirt" || ch_name == "sp_numuCC_7_0p_dirt" || ch_name == "sp_numuCC_8_0p_dirt"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_dirt" || ch_name == "sp_numuCC_2_Np_dirt" || ch_name == "sp_numuCC_3_Np_dirt" || ch_name == "sp_numuCC_4_Np_dirt"
             || ch_name == "sp_numuCC_5_Np_dirt" || ch_name == "sp_numuCC_6_Np_dirt" || ch_name == "sp_numuCC_7_Np_dirt" || ch_name == "sp_numuCC_8_Np_dirt"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_dirt" || ch_name == "sp_numuCC_2_Xp_dirt" || ch_name == "sp_numuCC_3_Xp_dirt" || ch_name == "sp_numuCC_4_Xp_dirt"
             || ch_name == "sp_numuCC_5_Xp_dirt" || ch_name == "sp_numuCC_6_Xp_dirt" || ch_name == "sp_numuCC_7_Xp_dirt" || ch_name == "sp_numuCC_8_Xp_dirt"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_0p_nc_delta_overlay" || ch_name == "sp_numuCC_0p_nc_delta_overlay_add" || ch_name == "sp_numuCC_2_0p_nc_delta_overlay"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Np_nc_delta_overlay" || ch_name == "sp_numuCC_Np_nc_delta_overlay_add" || ch_name == "sp_numuCC_2_Np_nc_delta_overlay"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_nc_delta_overlay" || ch_name == "sp_numuCC_Xp_nc_delta_overlay_add" || ch_name == "sp_numuCC_2_Xp_nc_delta_overlay"){
                if (flag_FC && flag_numuCC && !flag_ncpio_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_numuCC_Xp_overlay" || ch_name == "sp_numuCC_2_Xp_overlay"){
                if (flag_FC && flag_numuCC && (!flag_ncpio_sel) && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p" || ch_name == "sp_bdt_numuCC_2_0p" || ch_name == "sp_bdt_numuCC_3_0p" || ch_name == "sp_bdt_numuCC_4_0p"
             || ch_name == "sp_bdt_numuCC_5_0p" || ch_name == "sp_bdt_numuCC_6_0p" || ch_name == "sp_bdt_numuCC_7_0p" || ch_name == "sp_bdt_numuCC_8_0p"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np" || ch_name == "sp_bdt_numuCC_2_Np" || ch_name == "sp_bdt_numuCC_3_Np" || ch_name == "sp_bdt_numuCC_4_Np"
             || ch_name == "sp_bdt_numuCC_5_Np" || ch_name == "sp_bdt_numuCC_6_Np" || ch_name == "sp_bdt_numuCC_7_Np" || ch_name == "sp_bdt_numuCC_8_Np"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp" || ch_name == "sp_bdt_numuCC_2_Xp" || ch_name == "sp_bdt_numuCC_3_Xp" || ch_name == "sp_bdt_numuCC_4_Xp"
             || ch_name == "sp_bdt_numuCC_5_Xp" || ch_name == "sp_bdt_numuCC_6_Xp" || ch_name == "sp_bdt_numuCC_7_Xp" || ch_name == "sp_bdt_numuCC_8_Xp"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_2_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_3_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_4_Xp_nsbeam"
             || ch_name == "sp_bdt_numuCC_5_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_6_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_7_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_8_Xp_nsbeam"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_ext" || ch_name == "sp_bdt_numuCC_2_0p_ext" || ch_name == "sp_bdt_numuCC_3_0p_ext" || ch_name == "sp_bdt_numuCC_4_0p_ext"
             || ch_name == "sp_bdt_numuCC_5_0p_ext" || ch_name == "sp_bdt_numuCC_6_0p_ext" || ch_name == "sp_bdt_numuCC_7_0p_ext" || ch_name == "sp_bdt_numuCC_8_0p_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_ext" || ch_name == "sp_bdt_numuCC_2_Np_ext" || ch_name == "sp_bdt_numuCC_3_Np_ext" || ch_name == "sp_bdt_numuCC_4_Np_ext"
             || ch_name == "sp_bdt_numuCC_5_Np_ext" || ch_name == "sp_bdt_numuCC_6_Np_ext" || ch_name == "sp_bdt_numuCC_7_Np_ext" || ch_name == "sp_bdt_numuCC_8_Np_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_ext" || ch_name == "sp_bdt_numuCC_2_Xp_ext" || ch_name == "sp_bdt_numuCC_3_Xp_ext" || ch_name == "sp_bdt_numuCC_4_Xp_ext"
             || ch_name == "sp_bdt_numuCC_5_Xp_ext" || ch_name == "sp_bdt_numuCC_6_Xp_ext" || ch_name == "sp_bdt_numuCC_7_Xp_ext" || ch_name == "sp_bdt_numuCC_8_Xp_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_dirt" || ch_name == "sp_bdt_numuCC_2_0p_dirt" || ch_name == "sp_bdt_numuCC_3_0p_dirt" || ch_name == "sp_bdt_numuCC_4_0p_dirt"
             || ch_name == "sp_bdt_numuCC_5_0p_dirt" || ch_name == "sp_bdt_numuCC_6_0p_dirt" || ch_name == "sp_bdt_numuCC_7_0p_dirt" || ch_name == "sp_bdt_numuCC_8_0p_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_dirt" || ch_name == "sp_bdt_numuCC_2_Np_dirt" || ch_name == "sp_bdt_numuCC_3_Np_dirt" || ch_name == "sp_bdt_numuCC_4_Np_dirt"
             || ch_name == "sp_bdt_numuCC_5_Np_dirt" || ch_name == "sp_bdt_numuCC_6_Np_dirt" || ch_name == "sp_bdt_numuCC_7_Np_dirt" || ch_name == "sp_bdt_numuCC_8_Np_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_dirt" || ch_name == "sp_bdt_numuCC_2_Xp_dirt" || ch_name == "sp_bdt_numuCC_3_Xp_dirt" || ch_name == "sp_bdt_numuCC_4_Xp_dirt"
             || ch_name == "sp_bdt_numuCC_5_Xp_dirt" || ch_name == "sp_bdt_numuCC_6_Xp_dirt" || ch_name == "sp_bdt_numuCC_7_Xp_dirt" || ch_name == "sp_bdt_numuCC_8_Xp_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_overlay" || ch_name == "sp_bdt_numuCC_2_Xp_overlay" || ch_name == "sp_bdt_numuCC_3_Xp_overlay" || ch_name == "sp_bdt_numuCC_4_Xp_overlay"
             || ch_name == "sp_bdt_numuCC_5_Xp_overlay" || ch_name == "sp_bdt_numuCC_6_Xp_overlay" || ch_name == "sp_bdt_numuCC_7_Xp_overlay" || ch_name == "sp_bdt_numuCC_8_Xp_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_overlay" || ch_name == "sp_bdt_numuCC_2_0p_overlay" || ch_name == "sp_bdt_numuCC_3_0p_overlay" || ch_name == "sp_bdt_numuCC_4_0p_overlay"
             || ch_name == "sp_bdt_numuCC_5_0p_overlay" || ch_name == "sp_bdt_numuCC_6_0p_overlay" || ch_name == "sp_bdt_numuCC_7_0p_overlay" || ch_name == "sp_bdt_numuCC_8_0p_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_overlay" || ch_name == "sp_bdt_numuCC_2_Np_overlay" || ch_name == "sp_bdt_numuCC_3_Np_overlay" || ch_name == "sp_bdt_numuCC_4_Np_overlay"
             || ch_name == "sp_bdt_numuCC_5_Np_overlay" || ch_name == "sp_bdt_numuCC_6_Np_overlay" || ch_name == "sp_bdt_numuCC_7_Np_overlay" || ch_name == "sp_bdt_numuCC_8_Np_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_2_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_3_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_numuCC_5_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_6_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_7_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_8_Xp_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_2_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_3_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_4_0p_overlay_BG"
             || ch_name == "sp_bdt_numuCC_5_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_6_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_7_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_8_0p_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_2_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_3_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_4_Np_overlay_BG"
             || ch_name == "sp_bdt_numuCC_5_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_6_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_7_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_8_Np_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_spoverlay_nodelta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_spoverlay_delta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  map_cuts_flag["SPNCDeltaSig"]) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_2_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_3_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_4_Xp_spoverlay"
             || ch_name == "sp_bdt_numuCC_5_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_6_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_7_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_8_Xp_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_spoverlay" || ch_name == "sp_bdt_numuCC_2_0p_spoverlay" || ch_name == "sp_bdt_numuCC_3_0p_spoverlay" || ch_name == "sp_bdt_numuCC_4_0p_spoverlay"
             || ch_name == "sp_bdt_numuCC_5_0p_spoverlay" || ch_name == "sp_bdt_numuCC_6_0p_spoverlay" || ch_name == "sp_bdt_numuCC_7_0p_spoverlay" || ch_name == "sp_bdt_numuCC_8_0p_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_spoverlay" || ch_name == "sp_bdt_numuCC_2_Np_spoverlay" || ch_name == "sp_bdt_numuCC_3_Np_spoverlay" || ch_name == "sp_bdt_numuCC_4_Np_spoverlay"
             || ch_name == "sp_bdt_numuCC_5_Np_spoverlay" || ch_name == "sp_bdt_numuCC_6_Np_spoverlay" || ch_name == "sp_bdt_numuCC_7_Np_spoverlay" || ch_name == "sp_bdt_numuCC_8_Np_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_8_Xp_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_8_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_2_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_3_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_5_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_6_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_7_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_8_0p_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_8_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_2_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_3_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_5_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_6_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_7_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_8_Np_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_8_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  //one shower numu sideband
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p" || ch_name == "sp_bdt_numuCC_oneshw_2_0p" || ch_name == "sp_bdt_numuCC_oneshw_3_0p" || ch_name == "sp_bdt_numuCC_oneshw_4_0p"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p" || ch_name == "sp_bdt_numuCC_oneshw_6_0p" || ch_name == "sp_bdt_numuCC_oneshw_7_0p" || ch_name == "sp_bdt_numuCC_oneshw_8_0p"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np" || ch_name == "sp_bdt_numuCC_oneshw_2_Np" || ch_name == "sp_bdt_numuCC_oneshw_3_Np" || ch_name == "sp_bdt_numuCC_oneshw_4_Np"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np" || ch_name == "sp_bdt_numuCC_oneshw_6_Np" || ch_name == "sp_bdt_numuCC_oneshw_7_Np" || ch_name == "sp_bdt_numuCC_oneshw_8_Np"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_nsbeam"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_nsbeam"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_ext"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_ext" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_ext"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_ext"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_ext" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_ext"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_ext"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_ext" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_ext"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_dirt"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_dirt" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_dirt"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_dirt"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_dirt" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_dirt"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_dirt"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_dirt" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_dirt"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_oneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_oneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_oneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_overlay_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_overlay_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_overlay_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_overlay_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_spoverlay_nodelta" ){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_spoverlay_delta" ){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  map_cuts_flag["SPNCDeltaSig"]) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_spoverlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_spoverlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_spoverlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_spoverlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_spoverlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_spoverlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_oneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_oneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_oneshw_8_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;

//not one shower numu sideband
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_nsbeam"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_nsbeam" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_nsbeam"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_ext"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_ext" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_ext"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_ext" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_ext"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_ext" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_ext"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_dirt"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_dirt"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_dirt"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_dirt" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_dirt"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_overlay_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_overlay_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_overlay_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_overlay_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_spoverlay_nodelta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_spoverlay_delta" ){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  map_cuts_flag["SPNCDeltaSig"]) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_spoverlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_spoverlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_spoverlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_spoverlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_spoverlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_ncpi0overlay" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_ncpi0overlay"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_numuCC_notoneshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_numuCC_notoneshw_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_numuCC_notoneshw_8_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_numu_sel && !flag_singlephoton_oneshw && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
              map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
              map_cuts_flag["SPNumuCCSig"])
              && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
              && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
              && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
              && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;

  //single photon other sideband
}else if (ch_name == "sp_bdt_other_0p" || ch_name == "sp_bdt_other_2_0p" || ch_name == "sp_bdt_other_3_0p" || ch_name == "sp_bdt_other_4_0p"
           || ch_name == "sp_bdt_other_5_0p" || ch_name == "sp_bdt_other_6_0p" || ch_name == "sp_bdt_other_7_0p" || ch_name == "sp_bdt_other_8_0p"
          || ch_name == "sp_bdt_other_9_0p" || ch_name == "sp_bdt_other_10_0p" || ch_name == "sp_bdt_other_11_0p" || ch_name == "sp_bdt_other_12_0p"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np" || ch_name == "sp_bdt_other_2_Np" || ch_name == "sp_bdt_other_3_Np" || ch_name == "sp_bdt_other_4_Np"
           || ch_name == "sp_bdt_other_5_Np" || ch_name == "sp_bdt_other_6_Np" || ch_name == "sp_bdt_other_7_Np" || ch_name == "sp_bdt_other_8_Np"
          || ch_name == "sp_bdt_other_9_Np" || ch_name == "sp_bdt_other_10_Np" || ch_name == "sp_bdt_other_11_Np" || ch_name == "sp_bdt_other_12_Np"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp" || ch_name == "sp_bdt_other_2_Xp" || ch_name == "sp_bdt_other_3_Xp" || ch_name == "sp_bdt_other_4_Xp"
             || ch_name == "sp_bdt_other_5_Xp" || ch_name == "sp_bdt_other_6_Xp" || ch_name == "sp_bdt_other_7_Xp" || ch_name == "sp_bdt_other_8_Xp"
            || ch_name == "sp_bdt_other_9_Xp" || ch_name == "sp_bdt_other_10_Xp" || ch_name == "sp_bdt_other_11_Xp" || ch_name == "sp_bdt_other_12_Xp"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_nsbeam" || ch_name == "sp_bdt_other_2_0p_nsbeam" || ch_name == "sp_bdt_other_3_0p_nsbeam" || ch_name == "sp_bdt_other_4_0p_nsbeam"
             || ch_name == "sp_bdt_other_5_0p_nsbeam" || ch_name == "sp_bdt_other_6_0p_nsbeam" || ch_name == "sp_bdt_other_7_0p_nsbeam" || ch_name == "sp_bdt_other_8_0p_nsbeam"
             || ch_name == "sp_bdt_other_9_0p_nsbeam" || ch_name == "sp_bdt_other_10_0p_nsbeam" || ch_name == "sp_bdt_other_11_0p_nsbeam" || ch_name == "sp_bdt_other_12_0p_nsbeam"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_nsbeam" || ch_name == "sp_bdt_other_2_Np_nsbeam" || ch_name == "sp_bdt_other_3_Np_nsbeam" || ch_name == "sp_bdt_other_4_Np_nsbeam"
             || ch_name == "sp_bdt_other_5_Np_nsbeam" || ch_name == "sp_bdt_other_6_Np_nsbeam" || ch_name == "sp_bdt_other_7_Np_nsbeam" || ch_name == "sp_bdt_other_8_Np_nsbeam"
             || ch_name == "sp_bdt_other_9_Np_nsbeam" || ch_name == "sp_bdt_other_10_Np_nsbeam" || ch_name == "sp_bdt_other_11_Np_nsbeam" || ch_name == "sp_bdt_other_12_Np_nsbeam"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_nsbeam" || ch_name == "sp_bdt_other_2_Xp_nsbeam" || ch_name == "sp_bdt_other_3_Xp_nsbeam" || ch_name == "sp_bdt_other_4_Xp_nsbeam"
           || ch_name == "sp_bdt_other_5_Xp_nsbeam" || ch_name == "sp_bdt_other_6_Xp_nsbeam" || ch_name == "sp_bdt_other_7_Xp_nsbeam" || ch_name == "sp_bdt_other_8_Xp_nsbeam"
             || ch_name == "sp_bdt_other_9_Xp_nsbeam" || ch_name == "sp_bdt_other_10_Xp_nsbeam" || ch_name == "sp_bdt_other_11_Xp_nsbeam" || ch_name == "sp_bdt_other_12_Xp_nsbeam"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_ext" || ch_name == "sp_bdt_other_2_0p_ext" || ch_name == "sp_bdt_other_3_0p_ext" || ch_name == "sp_bdt_other_4_0p_ext"
             || ch_name == "sp_bdt_other_5_0p_ext" || ch_name == "sp_bdt_other_6_0p_ext" || ch_name == "sp_bdt_other_7_0p_ext" || ch_name == "sp_bdt_other_8_0p_ext"
               || ch_name == "sp_bdt_other_9_0p_ext" || ch_name == "sp_bdt_other_10_0p_ext" || ch_name == "sp_bdt_other_11_0p_ext" || ch_name == "sp_bdt_other_12_0p_ext"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_ext" || ch_name == "sp_bdt_other_2_Np_ext" || ch_name == "sp_bdt_other_3_Np_ext" || ch_name == "sp_bdt_other_4_Np_ext"
             || ch_name == "sp_bdt_other_5_Np_ext" || ch_name == "sp_bdt_other_6_Np_ext" || ch_name == "sp_bdt_other_7_Np_ext" || ch_name == "sp_bdt_other_8_Np_ext"
             || ch_name == "sp_bdt_other_9_Np_ext" || ch_name == "sp_bdt_other_10_Np_ext" || ch_name == "sp_bdt_other_11_Np_ext" || ch_name == "sp_bdt_other_12_Np_ext"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_ext" || ch_name == "sp_bdt_other_2_Xp_ext" || ch_name == "sp_bdt_other_3_Xp_ext" || ch_name == "sp_bdt_other_4_Xp_ext"
             || ch_name == "sp_bdt_other_5_Xp_ext" || ch_name == "sp_bdt_other_6_Xp_ext" || ch_name == "sp_bdt_other_7_Xp_ext" || ch_name == "sp_bdt_other_8_Xp_ext"
            || ch_name == "sp_bdt_other_9_Xp_ext" || ch_name == "sp_bdt_other_10_Xp_ext" || ch_name == "sp_bdt_other_11_Xp_ext" || ch_name == "sp_bdt_other_12_Xp_ext"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_dirt" || ch_name == "sp_bdt_other_2_0p_dirt" || ch_name == "sp_bdt_other_3_0p_dirt" || ch_name == "sp_bdt_other_4_0p_dirt"
             || ch_name == "sp_bdt_other_5_0p_dirt" || ch_name == "sp_bdt_other_6_0p_dirt" || ch_name == "sp_bdt_other_7_0p_dirt" || ch_name == "sp_bdt_other_8_0p_dirt"
             || ch_name == "sp_bdt_other_9_0p_dirt" || ch_name == "sp_bdt_other_10_0p_dirt" || ch_name == "sp_bdt_other_11_0p_dirt" || ch_name == "sp_bdt_other_12_0p_dirt"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_dirt" || ch_name == "sp_bdt_other_2_Np_dirt" || ch_name == "sp_bdt_other_3_Np_dirt" || ch_name == "sp_bdt_other_4_Np_dirt"
             || ch_name == "sp_bdt_other_5_Np_dirt" || ch_name == "sp_bdt_other_6_Np_dirt" || ch_name == "sp_bdt_other_7_Np_dirt" || ch_name == "sp_bdt_other_8_Np_dirt"
            || ch_name == "sp_bdt_other_9_Np_dirt" || ch_name == "sp_bdt_other_10_Np_dirt" || ch_name == "sp_bdt_other_11_Np_dirt" || ch_name == "sp_bdt_other_12_Np_dirt"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_dirt" || ch_name == "sp_bdt_other_2_Xp_dirt" || ch_name == "sp_bdt_other_3_Xp_dirt" || ch_name == "sp_bdt_other_4_Xp_dirt"
             || ch_name == "sp_bdt_other_5_Xp_dirt" || ch_name == "sp_bdt_other_6_Xp_dirt" || ch_name == "sp_bdt_other_7_Xp_dirt" || ch_name == "sp_bdt_other_8_Xp_dirt"
            || ch_name == "sp_bdt_other_9_Xp_dirt" || ch_name == "sp_bdt_other_10_Xp_dirt" || ch_name == "sp_bdt_other_11_Xp_dirt" || ch_name == "sp_bdt_other_12_Xp_dirt"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_nc_delta_overlay" || ch_name == "sp_bdt_other_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_other_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_nc_delta_overlay" || ch_name == "sp_bdt_other_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_other_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_nc_delta_overlay" || ch_name == "sp_bdt_other_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_other_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_overlay" || ch_name == "sp_bdt_other_2_Xp_overlay" || ch_name == "sp_bdt_other_3_Xp_overlay" || ch_name == "sp_bdt_other_4_Xp_overlay"
             || ch_name == "sp_bdt_other_5_Xp_overlay" || ch_name == "sp_bdt_other_6_Xp_overlay" || ch_name == "sp_bdt_other_7_Xp_overlay" || ch_name == "sp_bdt_other_8_Xp_overlay"
             || ch_name == "sp_bdt_other_9_Xp_overlay" || ch_name == "sp_bdt_other_10_Xp_overlay" || ch_name == "sp_bdt_other_11_Xp_overlay" || ch_name == "sp_bdt_other_2_Xp_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_overlay" || ch_name == "sp_bdt_other_2_0p_overlay" || ch_name == "sp_bdt_other_3_0p_overlay" || ch_name == "sp_bdt_other_4_0p_overlay"
             || ch_name == "sp_bdt_other_5_0p_overlay" || ch_name == "sp_bdt_other_6_0p_overlay" || ch_name == "sp_bdt_other_7_0p_overlay" || ch_name == "sp_bdt_other_8_0p_overlay"
            || ch_name == "sp_bdt_other_9_0p_overlay"|| ch_name == "sp_bdt_other_10_0p_overlay" || ch_name == "sp_bdt_other_11_0p_overlay" || ch_name == "sp_bdt_other_12_0p_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_overlay" || ch_name == "sp_bdt_other_2_Np_overlay" || ch_name == "sp_bdt_other_3_Np_overlay" || ch_name == "sp_bdt_other_4_Np_overlay"
           || ch_name == "sp_bdt_other_5_Np_overlay" || ch_name == "sp_bdt_other_6_Np_overlay" || ch_name == "sp_bdt_other_7_Np_overlay" || ch_name == "sp_bdt_other_8_Np_overlay"
           || ch_name == "sp_bdt_other_9_Np_overlay" || ch_name == "sp_bdt_other_10_Np_overlay" || ch_name == "sp_bdt_other_11_Np_overlay" || ch_name == "sp_bdt_other_12_Np_overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_overlay_BG" || ch_name == "sp_bdt_other_2_Xp_overlay_BG" || ch_name == "sp_bdt_other_3_Xp_overlay_BG" || ch_name == "sp_bdt_other_4_Xp_overlay_BG"
             || ch_name == "sp_bdt_other_5_Xp_overlay_BG" || ch_name == "sp_bdt_other_6_Xp_overlay_BG" || ch_name == "sp_bdt_other_7_Xp_overlay_BG" || ch_name == "sp_bdt_other_8_Xp_overlay_BG"
             || ch_name == "sp_bdt_other_9_Xp_overlay_BG" || ch_name == "sp_bdt_other_10_Xp_overlay_BG" || ch_name == "sp_bdt_other_11_Xp_overlay_BG" || ch_name == "sp_bdt_other_12_Xp_overlay_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_overlay_BG" || ch_name == "sp_bdt_other_2_0p_overlay_BG" || ch_name == "sp_bdt_other_3_0p_overlay_BG" || ch_name == "sp_bdt_other_4_0p_overlay_BG"
             || ch_name == "sp_bdt_other_5_0p_overlay_BG" || ch_name == "sp_bdt_other_6_0p_overlay_BG" || ch_name == "sp_bdt_other_7_0p_overlay_BG" || ch_name == "sp_bdt_other_8_0p_overlay_BG"
             || ch_name == "sp_bdt_other_9_0p_overlay_BG" || ch_name == "sp_bdt_other_10_0p_overlay_BG" || ch_name == "sp_bdt_other_11_0p_overlay_BG" || ch_name == "sp_bdt_other_12_0p_overlay_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_overlay_BG" || ch_name == "sp_bdt_other_2_Np_overlay_BG" || ch_name == "sp_bdt_other_3_Np_overlay_BG" || ch_name == "sp_bdt_other_4_Np_overlay_BG"
             || ch_name == "sp_bdt_other_5_Np_overlay_BG" || ch_name == "sp_bdt_other_6_Np_overlay_BG" || ch_name == "sp_bdt_other_7_Np_overlay_BG" || ch_name == "sp_bdt_other_8_Np_overlay_BG"
             || ch_name == "sp_bdt_other_9_Np_overlay_BG" || ch_name == "sp_bdt_other_10_Np_overlay_BG" || ch_name == "sp_bdt_other_11_Np_overlay_BG" || ch_name == "sp_bdt_other_12_Np_overlay_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_spoverlay" || ch_name == "sp_bdt_other_2_Xp_spoverlay" || ch_name == "sp_bdt_other_3_Xp_spoverlay" || ch_name == "sp_bdt_other_4_Xp_spoverlay"
             || ch_name == "sp_bdt_other_5_Xp_spoverlay" || ch_name == "sp_bdt_other_6_Xp_spoverlay" || ch_name == "sp_bdt_other_7_Xp_spoverlay" || ch_name == "sp_bdt_other_8_Xp_spoverlay"
             || ch_name == "sp_bdt_other_9_Xp_spoverlay" || ch_name == "sp_bdt_other_10_Xp_spoverlay" || ch_name == "sp_bdt_other_11_Xp_spoverlay" || ch_name == "sp_bdt_other_12_Xp_spoverlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_spoverlay" || ch_name == "sp_bdt_other_2_0p_spoverlay" || ch_name == "sp_bdt_other_3_0p_spoverlay" || ch_name == "sp_bdt_other_4_0p_spoverlay"
            || ch_name == "sp_bdt_other_5_0p_spoverlay" || ch_name == "sp_bdt_other_6_0p_spoverlay" || ch_name == "sp_bdt_other_7_0p_spoverlay" || ch_name == "sp_bdt_other_8_0p_spoverlay"
            || ch_name == "sp_bdt_other_9_0p_spoverlay" || ch_name == "sp_bdt_other_10_0p_spoverlay" || ch_name == "sp_bdt_other_11_0p_spoverlay" || ch_name == "sp_bdt_other_12_0p_spoverlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_spoverlay" || ch_name == "sp_bdt_other_2_Np_spoverlay" || ch_name == "sp_bdt_other_3_Np_spoverlay" || ch_name == "sp_bdt_other_4_Np_spoverlay"
            || ch_name == "sp_bdt_other_5_Np_spoverlay" || ch_name == "sp_bdt_other_6_Np_spoverlay" || ch_name == "sp_bdt_other_7_Np_spoverlay" || ch_name == "sp_bdt_other_8_Np_spoverlay"
            || ch_name == "sp_bdt_other_9_Np_spoverlay" || ch_name == "sp_bdt_other_10_Np_spoverlay" || ch_name == "sp_bdt_other_11_Np_spoverlay" || ch_name == "sp_bdt_other_12_Np_spoverlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_2_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_3_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_4_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_other_5_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_6_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_7_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_8_Xp_ncpi0overlay"
             || ch_name == "sp_bdt_other_9_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_10_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_11_Xp_ncpi0overlay" || ch_name == "sp_bdt_other_12_Xp_ncpi0overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                    !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_2_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_3_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_4_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_5_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_6_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_7_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_8_Xp_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_9_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_10_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_11_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_12_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                    !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_ncpi0overlay" || ch_name == "sp_bdt_other_2_0p_ncpi0overlay" || ch_name == "sp_bdt_other_3_0p_ncpi0overlay" || ch_name == "sp_bdt_other_4_0p_ncpi0overlay"
             || ch_name == "sp_bdt_other_5_0p_ncpi0overlay" || ch_name == "sp_bdt_other_6_0p_ncpi0overlay" || ch_name == "sp_bdt_other_7_0p_ncpi0overlay" || ch_name == "sp_bdt_other_8_0p_ncpi0overlay"
             || ch_name == "sp_bdt_other_9_0p_ncpi0overlay" || ch_name == "sp_bdt_other_10_0p_ncpi0overlay" || ch_name == "sp_bdt_other_11_0p_ncpi0overlay" || ch_name == "sp_bdt_other_12_0p_ncpi0overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_2_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_3_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_4_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_5_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_6_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_7_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_8_0p_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_9_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_10_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_11_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_12_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_ncpi0overlay" || ch_name == "sp_bdt_other_2_Np_ncpi0overlay" || ch_name == "sp_bdt_other_3_Np_ncpi0overlay" || ch_name == "sp_bdt_other_4_Np_ncpi0overlay"
             || ch_name == "sp_bdt_other_5_Np_ncpi0overlay" || ch_name == "sp_bdt_other_6_Np_ncpi0overlay" || ch_name == "sp_bdt_other_7_Np_ncpi0overlay" || ch_name == "sp_bdt_other_8_Np_ncpi0overlay"
             || ch_name == "sp_bdt_other_9_Np_ncpi0overlay" || ch_name == "sp_bdt_other_10_Np_ncpi0overlay" || ch_name == "sp_bdt_other_11_Np_ncpi0overlay" || ch_name == "sp_bdt_other_12_Np_ncpi0overlay"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_other_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_2_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_3_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_4_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_5_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_6_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_7_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_8_Np_overlay_ncpi0_BG"
             || ch_name == "sp_bdt_other_9_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_10_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_11_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_other_12_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_other_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  //single photon nue cc sideband
  }else if (ch_name == "sp_bdt_nue_0p" || ch_name == "sp_bdt_nue_2_0p"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np" || ch_name == "sp_bdt_nue_2_Np"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp" || ch_name == "sp_bdt_nue_2_Xp"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_nsbeam" || ch_name == "sp_bdt_nue_2_0p_nsbeam"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_nsbeam" || ch_name == "sp_bdt_nue_2_Np_nsbeam"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_nsbeam" || ch_name == "sp_bdt_nue_2_Xp_nsbeam"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_ext" || ch_name == "sp_bdt_nue_2_0p_ext"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_ext" || ch_name == "sp_bdt_nue_2_Np_ext"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_ext" || ch_name == "sp_bdt_nue_2_Xp_ext"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_dirt" || ch_name == "sp_bdt_nue_2_0p_dirt"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_dirt" || ch_name == "sp_bdt_nue_2_Np_dirt"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_dirt" || ch_name == "sp_bdt_nue_2_Xp_dirt"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_nc_delta_overlay" || ch_name == "sp_bdt_nue_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_nc_delta_overlay" || ch_name == "sp_bdt_nue_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nue_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_overlay" || ch_name == "sp_bdt_nue_2_Xp_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_overlay" || ch_name == "sp_bdt_nue_2_0p_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_overlay" || ch_name == "sp_bdt_nue_2_Np_overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_overlay_BG" || ch_name == "sp_bdt_nue_2_Xp_overlay_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_overlay_BG" || ch_name == "sp_bdt_nue_2_0p_overlay_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_overlay_BG" || ch_name == "sp_bdt_nue_2_Np_overlay_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_spoverlay" || ch_name == "sp_bdt_nue_2_Xp_spoverlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_spoverlay" || ch_name == "sp_bdt_nue_2_0p_spoverlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_spoverlay" || ch_name == "sp_bdt_nue_2_Np_spoverlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;

  }else if (ch_name == "sp_bdt_nue_Xp_ncpi0overlay" || ch_name == "sp_bdt_nue_2_Xp_ncpi0overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_2_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_ncpi0overlay" || ch_name == "sp_bdt_nue_2_0p_ncpi0overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_2_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_ncpi0overlay" || ch_name == "sp_bdt_nue_2_Np_ncpi0overlay"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_2_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  
  //single photon nue cc sideband - all showers
  }else if (ch_name == "sp_bdt_nue_allshw_0p" || ch_name == "sp_bdt_nue_allshw_2_0p"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np" || ch_name == "sp_bdt_nue_allshw_2_Np"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp" || ch_name == "sp_bdt_nue_allshw_2_Xp"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_nsbeam" || ch_name == "sp_bdt_nue_allshw_2_0p_nsbeam"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_nsbeam" || ch_name == "sp_bdt_nue_allshw_2_Np_nsbeam"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_nsbeam" || ch_name == "sp_bdt_nue_allshw_2_Xp_nsbeam"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_nsbeam_photon) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_ext" || ch_name == "sp_bdt_nue_allshw_2_0p_ext"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_ext" || ch_name == "sp_bdt_nue_allshw_2_Np_ext"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_ext" || ch_name == "sp_bdt_nue_allshw_2_Xp_ext"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_dirt" || ch_name == "sp_bdt_nue_allshw_2_0p_dirt"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_dirt" || ch_name == "sp_bdt_nue_allshw_2_Np_dirt"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_dirt" || ch_name == "sp_bdt_nue_allshw_2_Xp_dirt"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_nc_delta_overlay" || ch_name == "sp_bdt_nue_allshw_0p_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_allshw_2_0p_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_nc_delta_overlay" || ch_name == "sp_bdt_nue_allshw_Np_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_allshw_2_Np_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (!flag_0p) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_nc_delta_overlay" || ch_name == "sp_bdt_nue_allshw_Xp_nc_delta_overlay_add" || ch_name == "sp_bdt_nue_allshw_2_Xp_nc_delta_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && (eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && flag_truth_inside)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_overlay" || ch_name == "sp_bdt_nue_allshw_2_Xp_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_overlay" || ch_name == "sp_bdt_nue_allshw_2_0p_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_overlay" || ch_name == "sp_bdt_nue_allshw_2_Np_overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_overlay_BG" || ch_name == "sp_bdt_nue_allshw_2_Xp_overlay_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_overlay_BG" || ch_name == "sp_bdt_nue_allshw_2_0p_overlay_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_overlay_BG" || ch_name == "sp_bdt_nue_allshw_2_Np_overlay_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])
                  && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                  && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                  && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_spoverlay" || ch_name == "sp_bdt_nue_allshw_2_Xp_spoverlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_spoverlay" || ch_name == "sp_bdt_nue_allshw_2_0p_spoverlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_spoverlay" || ch_name == "sp_bdt_nue_allshw_2_Np_spoverlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  (map_cuts_flag["SPNCDeltaSig"] || map_cuts_flag["SPOutFVSig"] ||
                  map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                  map_cuts_flag["SPNumuCCSig"])) return true;
                return false;

  }else if (ch_name == "sp_bdt_nue_allshw_Xp_ncpi0overlay" || ch_name == "sp_bdt_nue_allshw_2_Xp_ncpi0overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Xp_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_allshw_2_Xp_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_ncpi0overlay" || ch_name == "sp_bdt_nue_allshw_2_0p_ncpi0overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_0p_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_allshw_2_0p_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_ncpi0overlay" || ch_name == "sp_bdt_nue_allshw_2_Np_ncpi0overlay"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && (eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;
  }else if (ch_name == "sp_bdt_nue_allshw_Np_overlay_ncpi0_BG" || ch_name == "sp_bdt_nue_allshw_2_Np_overlay_ncpi0_BG"){
                if (flag_singlephoton_nue_sel_allshw && !flag_singlephoton_numu_sel && !flag_singlephoton_ncpi0_sel && (!flag_singlephoton_sel) && !flag_0p &&
                  !(map_cuts_flag["SPNCDeltaSig"] ||
                    map_cuts_flag["SPNCPi0Sig"] || map_cuts_flag["SPNCOtherSig"] ||
                    map_cuts_flag["SPNumuCCSig"])
                    && !(map_cuts_flag["SPOutFVSig"] && pfeval.truth_corr_nuvtxX<260.9 && pfeval.truth_corr_nuvtxX>-0.9
                    && pfeval.truth_corr_nuvtxY<129.0 && pfeval.truth_corr_nuvtxY>-127.1
                    && pfeval.truth_corr_nuvtxZ<1040.9 && pfeval.truth_corr_nuvtxZ>-4.0) 
                    && !(eval.match_completeness_energy>0.1*eval.truth_energyInside
                    && eval.truth_isCC==0 && flag_truth_inside && pfeval.truth_NprimPio==1)) return true;
                return false;

  }else{
    std::cout << "Not sure what cut: " << ch_name << std::endl;
  }

  return false;
}

bool LEEana::get_rw_cut_pass(TString cut, EvalInfo& eval, PFevalInfo& pfeval, TaggerInfo& tagger, KineInfo& kine){
  if(cut == "NCPi0"){
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && !(pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCPi0_Np"){
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && !(is_true_0p(pfeval)) && !(pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCPi0_0p"){
    if (eval.truth_isCC==0 && pfeval.truth_NprimPio==1 && is_true_0p(pfeval) && !(pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCDeltaNp_scale"){
    if(is_NCdelta_sel(tagger, pfeval) && eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && !(is_0p(tagger, kine, pfeval)) ) return true;
    return false;
  }else if(cut == "NCDelta0p_scale"){
    if(is_NCdelta_sel(tagger, pfeval) && eval.truth_isCC==0 && pfeval.truth_NCDelta==1 && is_0p(tagger, kine, pfeval)) return true;
    return false;
  }else if(cut == "NCPi0_NCDelta_Np"){
    if(eval.truth_isCC==0 && !(is_true_0p(pfeval)) && (pfeval.truth_NprimPio==1 || pfeval.truth_NCDelta==1)) return true;
    return false;
  }else if(cut == "NCPi0_NCDelta_0p"){
    if(eval.truth_isCC==0 && is_true_0p(pfeval) && (pfeval.truth_NprimPio==1 || pfeval.truth_NCDelta==1)) return true;
    return false;
  }else{
    std::cout<<"No matching reweighting cut, check reweight configuration file"<<std::endl;
  }
return false;
}

bool LEEana::is_far_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;

  bool flag_numuCC = is_numuCC(tagger);
  bool flag_pi0 = is_pi0(kine, flag_data);
  bool flag_cc_pi0 = is_cc_pi0(kine, flag_data);
  bool flag_NC = is_NC(tagger);

  double reco_Enu = get_reco_Enu_corr(kine, flag_data);

  if ((reco_Enu>=800 && tagger.nue_score >=0) ||
      (tagger.nue_score<=0 && (flag_numuCC || (flag_pi0 && flag_NC) ))) flag = true;
  return flag;
}
bool LEEana::is_near_sideband(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);

  if (reco_Enu < 800 && tagger.nue_score>0 && (reco_Enu>=600 || tagger.nue_score<=7)) flag = true;

  return flag ;
}

bool LEEana::is_LEE_signal(KineInfo& kine, TaggerInfo& tagger, bool flag_data){
  bool flag = false;
  double reco_Enu = get_reco_Enu_corr(kine, flag_data);
  if (reco_Enu < 600 && tagger.nue_score>7) flag = true;
  return flag;
}




bool LEEana::is_truth_nueCC_inside(EvalInfo& eval){
  bool flag = false;

  if (fabs(eval.truth_nuPdg)==12 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;

  return flag;
}

bool LEEana::is_truth_numuCC_inside(EvalInfo& eval){
   bool flag = false;

  if (fabs(eval.truth_nuPdg)==14 && eval.truth_isCC==1 && eval.truth_vtxInside==1)
    flag = true;

  return flag;
}



bool LEEana::is_FC(EvalInfo& eval){
  if (eval.match_isFC){
    return true;
  }else{
    return false;
  }
}

bool LEEana::is_cc_pi0(KineInfo& kine, bool flag_data){

  bool flag = false;

  if (flag_data){
    if (kine.kine_pio_mass>0){
      //     TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      // TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      //TLorentzVector pio = p1 + p2;
      //pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;

      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300)
	flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 ) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
      flag = true;
  }

  return flag;
}


bool LEEana::is_pi0(KineInfo& kine, bool flag_data){
  bool flag = false;

  if (flag_data){
    if (kine.kine_pio_mass>0){
      //      TLorentzVector p1(kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Sin(kine.kine_pio_theta_1/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_1/180.*3.1415926), kine.kine_pio_energy_1*TMath::Cos(kine.kine_pio_theta_1/180.*3.1415926), kine.kine_pio_energy_1);
      //TLorentzVector p2(kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Cos(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Sin(kine.kine_pio_theta_2/180.*3.1415926)*TMath::Sin(kine.kine_pio_phi_2/180.*3.1415926), kine.kine_pio_energy_2*TMath::Cos(kine.kine_pio_theta_2/180.*3.1415926), kine.kine_pio_energy_2);
      // TLorentzVector pio = p1 + p2;
      // pio *= em_charge_scale;
      double pio_mass = kine.kine_pio_mass * em_charge_scale;

      if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1* em_charge_scale > 40 && kine.kine_pio_energy_2* em_charge_scale > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && pio_mass > 22 && pio_mass < 300)
	flag = true;
    }
  }else{
    if ((kine.kine_pio_flag==1 && kine.kine_pio_vtx_dis < 9 || kine.kine_pio_flag==2) && kine.kine_pio_energy_1 > 40 && kine.kine_pio_energy_2 > 25 && kine.kine_pio_dis_1 < 110 && kine.kine_pio_dis_2 < 120 && kine.kine_pio_angle > 0 && kine.kine_pio_angle < 174  && kine.kine_pio_mass > 22 && kine.kine_pio_mass < 300)
      flag = true;
  }

  return flag;
}


bool LEEana::is_NCpio_sel(TaggerInfo& tagger_info, KineInfo& kine){ // includes all cuts except FC
  bool flag = false;
  if (tagger_info.nc_pio_score > 1.816 && tagger_info.numu_cc_flag >=0 && kine.kine_pio_energy_1 > 0. && kine.kine_pio_energy_2 > 0.) flag = true;
  return flag;
}
bool LEEana::is_NCdelta_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){ // includes all cuts except FC
  bool flag = false;
  if (tagger_info.nc_delta_score > 2.61 && tagger_info.numu_cc_flag >=0 && pfeval.reco_showerKE > 0) flag = true;
  return flag;
}

//Erin
bool LEEana::is_singlephoton_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.4 && tagger_info.single_photon_other_score > 0.2 &&
    tagger_info.single_photon_ncpi0_score > -0.05 && tagger_info.single_photon_nue_score > -1.0 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 && tagger_info.single_photon_nue_score > -3.0 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
bool LEEana::is_singleshower_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.4 && tagger_info.single_photon_other_score > 0.2 &&
    tagger_info.single_photon_ncpi0_score > -0.05 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
bool LEEana::is_singleshower_eff_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 &&
    tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_numu_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    //tagger_info.single_photon_numu_score < 0.4 && //pure
    tagger_info.single_photon_numu_score < 0.1 && //eff
    //tagger_info.single_photon_numu_score < -3.0 && //numu pure
    tagger_info.single_photon_numu_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_other_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score < -0.4 &&
    tagger_info.single_photon_other_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_ncpi0_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    //tagger_info.single_photon_numu_score > 0.4 && tagger_info.single_photon_other_score > 0.2 &&
    //tagger_info.single_photon_ncpi0_score < -0.05 && //pure
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score < -0.4 && //eff
    //tagger_info.single_photon_numu_score > -1.0 && tagger_info.single_photon_other_score > -0.4 &&
    //tagger_info.single_photon_ncpi0_score < -0.4 && //nc pi0 eff
    tagger_info.single_photon_ncpi0_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_nue_sel(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 && tagger_info.single_photon_nue_score < -3.0 &&
    tagger_info.shw_sp_n_20br1_showers==1 &&
    tagger_info.single_photon_nue_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_singlephoton_nue_sel_allshw(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0 &&
    tagger_info.single_photon_numu_score > 0.1 && tagger_info.single_photon_other_score > -0.4 &&
    tagger_info.single_photon_ncpi0_score > -0.4 && tagger_info.single_photon_nue_score < -3.0 &&
    tagger_info.single_photon_nue_score > -20.0) {flag = true;}
  return flag;
}

bool LEEana::is_nsbeam(PFevalInfo& pfeval, EvalInfo& eval){
  bool flag = false;

  double delta_time_calc = -9999.;
  //Merge Peaks
  double gap=18.936;
  double Shift=0;
  if (pfeval.run >= 17380){ Shift=2916.0; }
  else if (pfeval.run >= 13697){ Shift = 3147.3;}//3166.1;}
  else if (pfeval.run >= 10812){ Shift = 3568.5; }
  else if (pfeval.run >= 8321){ Shift = 3610.7;}
  else if (pfeval.run >= 5800){ Shift = 3164.4;}
  else if (pfeval.run >= 0){ Shift = 3168.9;}
  //else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}
  //if(run>8000 && run<10812){Shift=3610.7; }
  //if(run>=10812 && run <12500){Shift=3568.5; }
  double TThelp=pfeval.evtTimeNS-Shift+gap*0.5;
  double TT_merged = -9999.;

  //merge peaks
  if(TThelp>=0 && TThelp<gap*81.0){
    TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
  }

  delta_time_calc = TT_merged;

  if (abs(delta_time_calc) < 5.0) {flag = true;}
  return flag;
}

bool LEEana::is_nsbeam_photon(PFevalInfo& pfeval, EvalInfo& eval){
  bool flag = false;

  double delta_time_calc = -9999.;
  //Merge Peaks
  double gap=18.936;
  double Shift=0;
  if (pfeval.run >= 17380){ Shift=2916.0; }
  else if (pfeval.run >= 13697){ Shift = 3147.3;}//3166.1;}
  else if (pfeval.run >= 10812){ Shift = 3568.5; }
  else if (pfeval.run >= 8321){ Shift = 3610.7;}
  else if (pfeval.run >= 5800){ Shift = 3164.4;}
  else if (pfeval.run >= 0){ Shift = 3168.9;}
  //else if (pfeval.run > 0 ){ Shift = 3166.0;}//3168.9;}
  //if(run>8000 && run<10812){Shift=3610.7; }
  //if(run>=10812 && run <12500){Shift=3568.5; }
  double TThelp=pfeval.evtTimeNS-Shift+gap*0.5;
  double TT_merged = -9999.;

  //merge peaks
  if(TThelp>=0 && TThelp<gap*81.0){
    TT_merged=(TThelp-(int((TThelp)/gap))*gap)-gap*0.5;
  }

  delta_time_calc = TT_merged;

  if (delta_time_calc > -6.6 && delta_time_calc < 3.4) {flag = true;}
  return flag;
}

//break down selection
bool LEEana::is_singlephoton_pre(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20mev_showers > 0 &&
    pfeval.reco_nuvtxX>5.0 && pfeval.reco_nuvtxX<250.0) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_numu_score > 0.4) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_other(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_other_score > 0.2) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_ncpi0_score > -0.05 ) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_nue_score > -1.0) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_numu(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_numu_score > 0.1) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_other(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_other_score > -0.4) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_ncpi0(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_ncpi0_score > -0.4) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_eff_nue(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.single_photon_nue_score > -3.0) {flag = true;}
  return flag;
}
bool LEEana::is_singlephoton_oneshw(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;
  if (tagger_info.shw_sp_n_20br1_showers==1) {flag = true;}
  return flag;
}
//


bool LEEana::is_NC(TaggerInfo& tagger_info){
  bool flag = false;
  if ((!tagger_info.cosmict_flag) && tagger_info.numu_score < 0)
    flag = true;

  return flag;
}


bool LEEana::is_numuCC(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9)
    flag = true;

  return flag;
}

bool LEEana::is_numuCC_tight(TaggerInfo& tagger_info, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0)
    flag = true;

  return flag;
}

bool LEEana::is_0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
          //if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>0) Nproton++; //Erin: CHANGE, "actual 0p" aka no 35 MeV threshold
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
      if(Nproton==0) flag = true;
  }

  return flag;
}

bool LEEana::is_1p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
    // 1 lepton <=1 proton 0 charged pion
    // 1 lepton guaranteed by numu cc flag
    // using pi0 flag to remove pi0 component in channel definition
    int Nproton = 0;
    int Npion = 0;
    for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
        int pdgcode = kine.kine_particle_type->at(i);
        if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
        if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
    if(Nproton==1) flag = true;
  }

  return flag;
}


bool LEEana::is_0pi(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0){
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
      if(Npion==0) flag = true;
  }

  return flag;
}

bool LEEana::is_numuCC_1mu0p(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval){
  bool flag = false;

  if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){
      // 1 lepton <=1 proton 0 charged pion
      // 1 lepton guaranteed by numu cc flag
      // using pi0 flag to remove pi0 component in channel definition
      int Nproton = 0;
      int Npion = 0;
      for(size_t i=0; i<kine.kine_energy_particle->size(); i++)
      {
          int pdgcode = kine.kine_particle_type->at(i);
          if(abs(pdgcode)==2212 && kine.kine_energy_particle->at(i)>35) Nproton++; // KE threshold: 50 MeV, 1.5 cm?
          if(abs(pdgcode)==211 && kine.kine_energy_particle->at(i)>10) Npion++; // KE threshold: 10 MeV
      }
      if(Nproton==0) flag = true;
  }

  return flag;
}


bool LEEana::is_numuCC_lowEhad(TaggerInfo& tagger_info, KineInfo& kine, PFevalInfo& pfeval, bool flag_data){
    bool flag = false;

    if (tagger_info.numu_cc_flag>=0 && tagger_info.numu_score > 0.9 && pfeval.reco_muonMomentum[3]>0){

      double reco_Enu = get_reco_Enu_corr(kine, flag_data);

      Float_t Ehadron = reco_Enu - pfeval.reco_muonMomentum[3]*1000.;
      if(Ehadron<200) // MeV
        {
	  flag = true;
        }
    }
    return flag;
}

bool LEEana::is_numuCC_cutbased(TaggerInfo& tagger_info){
  bool flag = false;

  if (tagger_info.numu_cc_flag==1 && tagger_info.cosmict_flag==0)
    flag = true;

  return flag;
}


bool LEEana::is_nueCC(TaggerInfo& tagger_info){
  bool flag = false;
  // default 7.0
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 7.0)
    //  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score <= 7.0 && tagger_info.nue_score > 0)
    flag = true;

  return flag;
}

bool LEEana::is_loosenueCC(TaggerInfo& tagger_info){
  bool flag = false;
  if (tagger_info.numu_cc_flag >=0 && tagger_info.nue_score > 4.0)
    flag = true;

  return flag;
}

bool LEEana::is_generic(EvalInfo& eval){
  // not very useful for the main analysis
  bool flag = is_preselection(eval);

  flag = flag && (eval.stm_clusterlength > 15);
  return flag;
}

bool LEEana::is_preselection(EvalInfo& eval){ // == T_BDTvars.numu_cc_flag >= 0
  bool flag = false;

  // match code ...
  int tmp_match_found = eval.match_found;
  if (eval.is_match_found_int){
    tmp_match_found = eval.match_found_asInt;
  }

  if (tmp_match_found == 1 && eval.stm_eventtype != 0 && eval.stm_lowenergy ==0 && eval.stm_LM ==0 && eval.stm_TGM ==0 && eval.stm_STM==0 && eval.stm_FullDead == 0 && eval.stm_clusterlength >0) flag = true;


  return flag;
}

int LEEana::mcc8_pmuon_costheta_bin(float pmuon, float costh){

    if (costh>=-1 and costh<-0.5) {
      if (pmuon>=0 and pmuon<0.18) return 1;
      else if (pmuon>=0.18 and pmuon<0.30) return 2;
      else if (pmuon>=0.30 and pmuon<0.45) return 3;
      else if (pmuon>=0.45 and pmuon<0.77) return 4;
      else if (pmuon>=0.77 and pmuon<2.5) return 5;
      else return -10000;
    }
    else if (costh>=-0.5 and costh<0){
      if (pmuon>=0 and pmuon<0.18) return 6;
      else if (pmuon>=0.18 and pmuon<0.30) return 7;
      else if (pmuon>=0.30 and pmuon<0.45) return 8;
      else if (pmuon>=0.45 and pmuon<0.77) return 9;
      else if (pmuon>=0.77 and pmuon<2.5) return 10;
      else return -10000;
    }
    else if (costh>0 and costh<0.27){
      if (pmuon>=0 and pmuon<0.18) return 11;
      else if (pmuon>=0.18 and pmuon<0.30) return 12;
      else if (pmuon>=0.30 and pmuon<0.45) return 13;
      else if (pmuon>=0.45 and pmuon<0.77) return 14;
      else if (pmuon>=0.77 and pmuon<2.5) return 15;
      else return -10000;
    }
    else if (costh>=0.27 and costh<0.45){
      if (pmuon>=0 and pmuon<0.30) return 16;
      else if (pmuon>=0.30 and pmuon<0.45) return 17;
      else if (pmuon>=0.45 and pmuon<0.77) return 18;
      else if (pmuon>=0.77 and pmuon<2.5) return 19;
      else return -10000;
    }
    else if (costh>=0.45 and costh<0.62){
      if (pmuon>=0 and pmuon<0.30) return 20;
      else if (pmuon>=0.30 and pmuon<0.45) return 21;
      else if (pmuon>=0.45 and pmuon<0.77) return 22;
      else if (pmuon>=0.77 and pmuon<2.5) return 23;
      else return -10000;
    }
    else if (costh>=0.62 and costh<0.76){
      if (pmuon>=0 and pmuon<0.30) return 24;
      else if (pmuon>=0.30 and pmuon<0.45) return 25;
      else if (pmuon>=0.45 and pmuon<0.77) return 26;
      else if (pmuon>=0.77 and pmuon<2.5) return 27;
      else return -10000;
    }
    else if (costh>=0.76 and costh<0.86){
      if (pmuon>=0 and pmuon<0.30) return 28;
      else if (pmuon>=0.30 and pmuon<0.45) return 29;
      else if (pmuon>=0.45 and pmuon<0.77) return 30;
      else if (pmuon>=0.77 and pmuon<1.28) return 31;
      else if (pmuon>=1.28 and pmuon<2.5) return 32;
      else return -10000;
    }
    else if (costh>=0.86 and costh<0.94){
      if (pmuon>=0 and pmuon<0.30) return 33;
      else if (pmuon>=0.30 and pmuon<0.45) return 34;
      else if (pmuon>=0.45 and pmuon<0.77) return 35;
      else if (pmuon>=0.77 and pmuon<1.28) return 36;
      else if (pmuon>=1.28 and pmuon<2.5) return 37;
      else return -10000;
    }
    else if (costh>=0.94 and costh<1.00){
      if (pmuon>=0 and pmuon<0.30) return 38;
      else if (pmuon>=0.30 and pmuon<0.45) return 39;
      else if (pmuon>=0.45 and pmuon<0.77) return 40;
      else if (pmuon>=0.77 and pmuon<1.28) return 41;
      else if (pmuon>=1.28 and pmuon<2.5) return 42;
      else return -10000;
    }
    else return -10000;
}

// return bin number of a variable in alternative binning choice
// bin edges are defined in config (default: ./configurations/alt_var_xbins.txt)
int LEEana::alt_var_index(std::string var1, float val1, std::string var2, float val2, std::string config){

  if (map_var_hist.empty()) {
    // load all variables from configuration
    std::vector<float> xbins;
    float val;
    std::string varname, line;
    std::ifstream in(config);
    while ( std::getline(in, line) ) {
      std::istringstream ss(line);
      ss >> varname;
      xbins.clear();
      while (ss >> val) {
        xbins.push_back(val);
      }
      map_var_hist[varname] = TH1F(varname.c_str(), varname.c_str(), xbins.size()-1, xbins.data());
    }
    in.close();
  }

  auto varhist1 = map_var_hist[var1];
  auto varhist2 = map_var_hist[var2];
  int bin1 = varhist1.FindBin(val1);
  int bin2 = varhist2.FindBin(val2);
  int nBins1 = varhist1.GetNbinsX();
  int nBins2 = varhist2.GetNbinsX();
  if( bin1>0 and bin1<=nBins1 and bin2>0 and bin2<=nBins2){
    return bin1 + (bin2-1)*nBins1;
  }

  return -1;
}

#endif
