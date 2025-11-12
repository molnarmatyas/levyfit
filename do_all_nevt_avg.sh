#!/bin/bash

make clean
make exe/onedim_EbE_or_Eavg_fit.exe

nevt_avgs=(1 10 50 100 200 500 1000 5000 10000)
energies=("3p0" "7p7" "27")

for energy in "${energies[@]}"; do
  echo "Fitting for energy ${energy}"
  for avg in "${nevt_avgs[@]}"; do
    echo "Fitting for nevt_avg: ${avg}"
    exe/onedim_EbE_or_Eavg_fit.exe 11 "${energy}" 1 10000 ${avg} 1 &> fit_log_${energy}_nevtavg${avg}.log # avg. by nevt events
  done
done

echo "All done."