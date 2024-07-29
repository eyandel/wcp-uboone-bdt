#!/bin/bash

# Make box around text @climagic
function box() { t="$1xxxx";c=${2:-=}; echo ${t//?/$c}; echo "$c $1 $c"; echo ${t//?/$c}; }

usage() { echo "Usage: $0 [-i <path to input txt files>] [-o <path to output directory>] [-n] (optional : for numi)" 1>&2; exit 1; }

numi=false
while getopts ":i:o:nh" opt; do
  case "$opt" in
    i)
      inputtxt="$OPTARG"
      ;;
    o)
      outdir="$OPTARG"
      ;;
    n)
      numi=true
      ;;
    h | *)
      usage
      exit 1
      ;;
  esac
done

if [[ ! -e "$inputtxt" || -z "$outdir" ]]; then
  echo "Either an invalid input text file provided with -i or an invalid output directory with -o - Exiting!"
  usage
  exit 1
fi

if [[ ! -d "$outdir" ]]; then
  echo "Making output directory"
  mkdir -p "$outdir"
fi

if [[ "$numi" ]]; then
  box "Processing for NuMI"
fi

box "Pruning CV to "$outdir"/cv.root"
prune_checkout_trees $inputtxt $outdir"/cv.root"

# process systs
systs=("UBGenieFluxSmallUni" \
       "expskin_FluxUnisim" "horncurrent_FluxUnisim" "kminus_PrimaryHadronNormalization"  \
       "kplus_PrimaryHadronFeynmanScaling" "kzero_PrimaryHadronSanfordWang" "nucleoninexsec_FluxUnisim" "nucleonqexsec_FluxUnisim" "nucleontotxsec_FluxUnisim" \
       "piminus_PrimaryHadronSWCentralSplineVariation" "pioninexsec_FluxUnisim" "pionqexsec_FluxUnisim" "piontotxsec_FluxUnisim" "piplus_PrimaryHadronSWCentralSplineVariation" \
\
       "reinteractions_piminus_Geant4" "reinteractions_piplus_Geant4" "reinteractions_proton_Geant4")

box "Launching threads for pruning the systematic files"
for syst in "${systs[@]}"; do
  echo $syst
  if [[ "$numi" ]]; then
    prune_weightsep24_trees $inputtxt $outdir"/"$syst".root" $syst &
  else
    prune_weightsep24_trees_numi $inputtxt $outdir"/"$syst".root" $syst &
  fi
  sleep 10
done
wait
