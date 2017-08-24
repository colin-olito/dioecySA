#!/bin/bash

export LANG=C # to get back ABCabc instead of AaBbCc

cd analyzeH/
echo analyzeH
for dir in [A-Z]*/; do
  echo $dir
  if [ -d "$dir" ]; then
    cd $dir
    Rscript ../../scripts/run_mcmc.R Herm
    Rscript ../../scripts/run_maps.R Herm
    cd ../
  fi 
done
cd ../

cd analyzeD/
echo analyzeD
for dir in [A-Z]*/; do
  echo $dir
  if [ -d "$dir" ]; then
    cd $dir
    Rscript ../../scripts/run_mcmc.R Dio
    Rscript ../../scripts/run_maps.R Dio
    cd ../
  fi 
done
cd ../

cd analyze4/
echo analyze4
for dir in [A-Z]*/; do
  echo $dir
  if [ -d "$dir" ]; then
    cd $dir
    Rscript ../../scripts/run_mcmc.R 4
    Rscript ../../scripts/run_maps.R 4
    cd ../
  fi 
done
cd summaries/
Rscript tip_root.R
cd ../
cd ../
