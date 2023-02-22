#!/bin/zsh

# Clock time
SECONDS=0
# Exit if any command fails
set -e
# Check for R dependencies
R/packages.R
# Check for files in data directory
FNUM=$(find data/marx_traced -name '*marx.RData' -print | wc -l)
if [[ $FNUM -lt 64 ]]; then
  # Download data from osf
  # https://osf.io/3emwf/files/osfstorage
  R/download-data.R
fi
# Preprocess data
R/preprocess.R
# Run study and visualize results
FNUM=$(find data/marx_traced -name '*marx.RData' -print | wc -l)
if [[ ! $FNUM -lt 64 ]]; then
  echo 'Found previous marker traces:'
  ls -1q data/marx_traced/*marx.RData | head -n 10
  echo '\nClassify markers?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      unset itr
      echo '\nMax iterations for jackknife resampling?'
      echo 'Note: computation cost is high! (recommended < 50 iterations)'
      read 'itr?Number: '
      while [[ -z ${itr} ]]; do
        read 'itr?Number: '
      done
      unset prop
      echo '\nProportion of markers to keep during jackknife resampling?'
      echo 'Note: enter proportion between (0-1; 0.9=default)'
      read 'prop?Proportion: '
      while [[ -z ${prop} ]]; do
        read 'prop?Proportion: '
      done
      unset kcluster
      echo '\nNumber of k clusters for Gaussian Mixture Modelling?'
      echo 'Note: recommended k = 14 clusters'
      read 'kcluster?Clusters: '
      while [[ -z ${kcluster} ]]; do
        read 'kcluster?Clusters: '
      done
      unset thresh
      echo '\nThreshold gradient for marker classification?'
      echo 'Note: markers with cluster centriods below this gradient'
      echo '  will be classified as subducted (recommended 3-5 C/km)'
      read 'thresh?Threshold: '
      while [[ -z ${thresh} ]]; do
        read 'thresh?Threshold: '
      done
      unset ptpathfilter
      echo '\nPT path filter to use for marker classification?'
      echo 'Note: options are "maxP" (default) "maxT" and "maxPT"'
      read 'ptpathfilter?Filter: '
      while [[ -z ${ptpathfilter} ]]; do
        read 'ptpathfilter?Filter: '
      done
      unset ncores
      echo '\nAvailable cores on this machine:'
      nproc --all
      echo 'Please enter number of cores for parallel computing'
      read 'ncores?Cores: '
      while [[ -z ${ncores} ]]; do
        read 'ncores?Cores: '
      done
      R/classify.R $itr $prop $kcluster $thresh $ptpathfilter $ncores
      echo '\nMarker classification successful!'
      # Print clock time
      t=$SECONDS
      printf '\nTime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      break
    elif [[ $p == 'no' ]]; then
      # Print clock time
      t=$SECONDS
      printf '\nTime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      break
    else
      read 'p?yes/no: '
    fi
  done
  echo '\nVisualize results?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      unset ncores
      echo '\nAvailable cores on this machine:'
      nproc --all
      echo 'Please enter number of cores for parallel computing'
      read 'ncores?Cores: '
      while [[ -z ${ncores} ]]; do
        read 'ncores?Cores: '
      done
      unset thresh
      echo '\nThreshold gradient for marker classification?'
      echo 'Note: markers with cluster centriods below this gradient'
      echo '  will be classified as subducted (recommended 3-5 C/km)'
      read 'thresh?Threshold: '
      while [[ -z ${thresh} ]]; do
        read 'thresh?Threshold: '
      done
      unset ptpathfilter
      echo '\nPT path filter to use for marker classification?'
      echo 'Note: options are "maxP" (default) "maxT" and "maxPT"'
      read 'ptpathfilter?Filter: '
      while [[ -z ${ptpathfilter} ]]; do
        read 'ptpathfilter?Filter: '
      done
      R/marx-vis.R $ncores $thresh $ptpathfilter
      R/summary.R $ptpathfilter
      echo '\nCopying data and figures to draft folder for manuscript'
      cp data/*.csv draft/assets/r/
      cp data/*.RData draft/assets/r/
      cp -r figs/* draft/assets/figs/
      echo '\nFinished!'
      # Print clock time
      t=$SECONDS
      printf '\nTime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    elif [[ $p == 'no' ]]; then
      echo 'okay bye'
      # Print clock time
      t=$SECONDS
      printf '\nTime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
fi