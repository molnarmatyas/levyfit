#!/bin/bash

make clean
make exe/onedim_EbE_or_Eavg_fit.exe

nevt_avgs=(1 50 100 200 500 1000 5000 10000)

echo "Fitting for energy 3p0"
for avg in "${nevt_avgs[@]}"; do
  echo "Fitting for nevt_avg: ${avg}"
  exe/onedim_EbE_or_Eavg_fit.exe 11 "3p0" 1 10000 ${avg} 1 &> fit_log_3p0_nevtavg${avg}.log # avg. by nevt events
done

echo "Fitting for energy 27"
for avg in "${nevt_avgs[@]}"; do
  echo "Fitting for nevt_avg: ${avg}"
  exe/onedim_EbE_or_Eavg_fit.exe 11 "27" 1 10000 ${avg} 1 &> fit_log_3p0_nevtavg${avg}.log # avg. by nevt events
done

echo "All done."