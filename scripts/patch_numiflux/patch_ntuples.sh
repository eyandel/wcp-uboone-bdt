#!/bin/bash

indir="/home/nitish/data/uboone/inputs/flugg_studies/patch/existing/"
outdir=`echo $indir | sed -e 's/existing/new/g'`
mkdir -p $outdir

mkdir -p tempmacro
cp utils.h tempmacro/.

for fin in `fd -g "checkout*.root" $indir`; do
  fout=`echo $fin | sed -e 's/existing/new/g'`
  foutdir=`dirname "$fout"`
  foutname=`basename "$fout" | sed -e 's/checkout.*_\(run.*\).root/\1/g'`
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
  mkdir -p tempmacro
  cp patch_ntuples.C "tempmacro/patch_ntuples_$foutname.C"
  sed -i 's/patch_ntuples/patch_ntuples_'"$foutname"'/g' "tempmacro/patch_ntuples_$foutname.C"
  time root -b -q 'tempmacro/patch_ntuples_'"$foutname"'.C+("'"$fin"'", "'"$fout"'", '"$isrhc"', '"$isfullosc"')' &
done

wait
rm -rf tempmacro
box "Finished!"
