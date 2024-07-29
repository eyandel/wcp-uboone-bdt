#!/bin/bash

# Make box around text @climagic
function box() { t="$1xxxx";c=${2:-=}; echo ${t//?/$c}; echo "$c $1 $c"; echo ${t//?/$c}; }

usage() { echo "Usage: $0 [-c <path to pruned cv file>] [-s <path to folder containing pruned syst files>] [-o <path to output folder>] [-p \"prefix name for the output CV file\"]" 1>&2; exit 1; }

data=false
while getopts ":c:s:o:p:h" opt; do
  case "$opt" in
    c)
      inputcv="$OPTARG"
      ;;
    s)
      inputsystfolder="$OPTARG"
      ;;
    o)
      outfolder="$OPTARG"
      ;;
    p)
      outcvname="$OPTARG"
      ;;
    h | *)
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$inputsystfolder" ]]; then
  echo "No input syst folder provided. Assuming data!"
  data=true
fi

if [[ -z "$outfolder" || -z "$outcvname" ]]; then
  echo "Either no output folder (with -o) or no output prefix provided (with -p). Exiting!"
  usage
  exit 1
fi

if [[ ! -e "$inputcv" ]]; then
  echo "Provide valid input cv file with -c - Exiting!"
  usage
  exit 1
fi

if [[ ! -d "$inputsystfolder" ]]; then
  echo "Could not find input syst folder provided by -s - Exiting!"
  usage
  exit 1
fi

if [[ ! -d "$outfolder" ]]; then
  echo "Making output directory"
  mkdir -p "$outfolder"
fi

if [[ ! -d $outfolder"/systs"  && ! "$data" ]]; then
  echo "Making output directory for systs"
  mkdir -p $outfolder"/systs"
fi


systs=("UBGenieFluxSmallUni" \
       "expskin_FluxUnisim" "horncurrent_FluxUnisim" "kminus_PrimaryHadronNormalization"  \
       "kplus_PrimaryHadronFeynmanScaling" "kzero_PrimaryHadronSanfordWang" "nucleoninexsec_FluxUnisim" "nucleonqexsec_FluxUnisim" "nucleontotxsec_FluxUnisim" \
       "piminus_PrimaryHadronSWCentralSplineVariation" "pioninexsec_FluxUnisim" "pionqexsec_FluxUnisim" "piontotxsec_FluxUnisim" "piplus_PrimaryHadronSWCentralSplineVariation" \
       "reinteractions_piminus_Geant4" "reinteractions_piplus_Geant4" "reinteractions_proton_Geant4")

# process CV first
box "Processing CV"

box "BDT Convert"
bdt_convert $inputcv cv_bdt.tmp.root

cvfile=$outfolder"/"$outcvname".root"

box "Converting CV and writing to "$cvfile
convert_cv_spec "cv_bdt.tmp.root" $cvfile
rm cv_bdt.tmp.root

# process systs
if [[ !"$data" ]]; then
  box "Processing for systematics"

  for syst in "${systs[@]}"; do
    box $syst
    merge_xf $cvfile $inputsystfolder"/"$syst".root" $outfolder"/systs/"$syst".root" $syst &
    sleep 5
  done
  wait
fi
