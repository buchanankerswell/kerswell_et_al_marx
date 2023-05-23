#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)
# Set download timeout to 10min
cat(rep('~', 80), sep='')
cat('\nSetting download timeout to 10 minutes ...')
options(timeout = 1200)
# Download data from osf
# https://osf.io/3emwf/files/osfstorage
data.url <-
  'https://files.osf.io/v1/resources/3emwf/providers/osfstorage/63753f6db07665090b617b8f/?zip='
# Download .zip file
cat('\nDownloading data from osf ...')
cat('\nurl: https://osf.io/3emwf/files/osfstorage')
cat('\nThis may take > 5 minutes ...')
download.file(data.url, 'data.zip', quiet = T)
# Extract
unzip('data.zip', exdir = 'data')
# Remove .zip file
file.remove('data.zip')
# Write log
cat('\ndownload-data.R complete!\n')
sink()