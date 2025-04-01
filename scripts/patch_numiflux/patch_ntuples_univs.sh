#!/bin/bash

indir="/home/nitish/data/uboone/inputs/flugg_studies/patch/existing/"
outdir=`echo $indir | sed -e 's/existing/new/g'`
mkdir -p $outdir

mkdir -p tempmacro
cp utils.h tempmacro/.

for fin in `fd -g "kminus*.root" $indir`; do
  fout=`echo $fin | sed -e 's/existing/new/g'`
  foutdir=`dirname "$fout"`
  foutname=`basename "$foutdir"`
  mkdir -p "$foutdir"

  box "$fin"
  isrhc="true"
  isfullosc="false"
  if [[ "$fin" =~ "fhc" ]]; then
    isrhc="false"
  fi
  if [[ "$fin" =~ "fullosc" ]]; then
    isfullosc="true"
  fi
  cp patch_ntuples_univs.C "tempmacro/patch_ntuples_univs_$foutname.C"
  sed -i 's/patch_ntuples_univs/patch_ntuples_univs_'"$foutname"'/g' "tempmacro/patch_ntuples_univs_$foutname.C"
  time root -b -q 'tempmacro/patch_ntuples_univs_'"$foutname"'.C+("'"$fin"'", "'"$fout"'", '"$isrhc"', '"$isfullosc"')' &
done

wait
rm -rf tempmacro
box "Finished!"
