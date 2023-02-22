#!/usr/bin/env Rscript

# Make dir for log
dir.create('data', showWarnings = F)

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Quiet loading
sshhh <- function(p) {
  suppressWarnings(
    suppressPackageStartupMessages({
      cat('\nChecking required package:', p)
      require(p, quietly = T, character.only = TRUE)
    })
  )
}

# Install dependencies method
using <- function(...) {
  # Try loading required packages
  pkgs <- unlist(list(...))
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  # Try installing required packages
  if(length(need) > 0) {
    install.packages(need)
    lapply(need, sshhh)
  }
  # Check installs
  req <- unlist(lapply(pkgs, sshhh))
  need <- pkgs[req == FALSE]
  if(length(need) > 0) {
    cat('\n', rep('~', 80), sep = '')
    cat('\nFailed to install packages:\n')
    writeLines(need)
    cat('\nFor troubleshooting tips see:')
    cat('\nhttps://github.com/buchanankerswell/kerswell_kohn_backarc')
    cat('\n', rep('~', 80), '\n', sep = '')
    stop()
  }
}

# Package list
package.list <- c(
  'reshape2',
  'plot3D',
  'ggplotify',
  'Hmisc',
  'mclust',
  'magrittr',
  'tidyr',
  'readr',
  'stringr',
  'tibble',
  'dplyr',
  'purrr',
  'furrr',
  'future',
  'ggforce',
  'ggnewscale',
  'ggrepel',
  'patchwork',
  'gridExtra',
  'progress',
  'viridis',
  'metR',
  'colorspace'
)

cat(rep('~', 80), sep='')
cat('\nChecking for required R packages ...')
using(package.list)

# Write log
cat('\nAll required packages installed and available!')
cat('\npackages.R complete!\n')

# Print session info
cat(rep('~', 80), '\n', sep='')
sessionInfo()

sink()