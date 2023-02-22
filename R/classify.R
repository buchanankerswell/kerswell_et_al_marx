#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Load functions and libraries
cat(rep('~', 80), sep='')
cat('\nLoading packages and functions')
source('R/functions.R')

# Test arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  cat('\nNo arguments passed to R/classify.R')
  n <- 0
  p <- 0.9
  k <- 14
  gradient.threshold <- 3
  pt.path.filter <- 'maxP'
  cores <- availableCores()-2
  cat('\nUsing defaults')
  cat('\nJackknife resampling iterations  :', n)
  cat('\nJackknife resampling proportion  :', p)
  cat('\nNumber of cluster components     :', k)
  cat('\nThermal gradient threshold       :', gradient.threshold)
  cat('\nPT path filter                   :', pt.path.filter)
} else if (length(args) != 0) {
  n <- suppressWarnings(as.integer(args[1]))
  if(n < 0 | n > 100) {
    n <- 30
    cat('\nJackknife resampling iterations must be between [0-100]')
    cat('\nDefaulting to', n)
  }
  p <- suppressWarnings(as.numeric(args[2]))
  if(p < 0 | p > 1) {
    p <- 0.9
    cat('\nJackknife resampling proportion must be between [0-1]')
    cat('\nDefaulting to', p)
  }
  k <- suppressWarnings(as.integer(args[3]))
  if(k < 1 | k > 20) {
    k <- 10
    cat('\nNumber of cluster components must be between [0-20]')
    cat('\nDefaulting to', k)
  }
  gradient.threshold <- suppressWarnings(as.numeric(args[4]))
  if(gradient.threshold < 0 | gradient.threshold > 5) {
    gradient.threshold <- 3
    cat('\nThreshold thermal gradient must be between [0-5] C/km')
    cat('\nDefaulting to', gradient.threshold)
  }
  pt.path.filter <- suppressWarnings(args[5])
  if(!(pt.path.filter %in% c('maxT', 'maxP', 'maxPT'))) {
    pt.path.filter <- 'maxP'
    cat('\nUnrecognized PT path filter. Please use "maxT", "maxP", or "maxPT"')
    cat('\nDefaulting to', pt.path.filter)
  }
  cores <- suppressWarnings(as.integer(args[6]))
  if(cores < 1 | cores > availableCores()) {
    cores <- availableCores()-2
    cat('\nUsing default number of cores: ', cores, sep = '')
  }
}

# Get filepaths
paths <- list.files('data/marx_traced', pattern = '-marx.RData', full.names = T)
models <- str_extract(paths, 'cd.[0-9]+')
cat('\nFound traced marx for ', length(models), ' numerical experiments:\n', sep = '')
cat('\nClassifying markers')
cat('\nDate                             :', as.character(Sys.Date()))
cat('\nJackknife resampling iterations  :', n)
cat('\nJackknife resampling proportion  :', p)
cat('\nNumber of cluster components     :', k)
cat('\nThermal gradient threshold       :', gradient.threshold)
cat('\nPT path filter                   :', pt.path.filter)
cat('\n', rep('~', 80), '\n', sep='')

# Function to parallelize
fun <- function(model, path, n, p, gradient.threshold, pt.path.filter, k) {
  cat('\nClassifying markers [', model, ']', sep = '')
  # Load markers
  marx <- load_marx(path)
  # Classify markers
  marx.class <- marx_classify(marx, model, gradient.threshold, pt.path.filter, k)
  # Jackknife sampling
  jk <- jknife(marx, model, gradient.threshold, pt.path.filter, n, p, k)
  # Save
  assign(paste0(model, '.marx.classified'), list('jk' = jk, 'marx.class' = marx.class))
  save(
    list = paste0(model, '.marx.classified'),
    file = paste0(
      'data/marx_classified', '/', model, '-marx-classified.RData'
    )
  )
}

# Set parallel plan
cat('\nParallel computing with ', cores, ' cores', sep = '')
plan(multisession, workers = cores)

# Classify markers
future_walk2(
  models,
  paths,
  ~fun(.x, .y, n, p, gradient.threshold, pt.path.filter, k),
  .progress = T,
  .options = furrr_options(seed = T)
)

# Write log
cat('\nclassify.R complete!\n')
sink()
