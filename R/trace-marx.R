#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Load functions and libraries
cat(rep('~', 80), sep='')
cat('\nLoading packages and functions ...')
source('R/functions.R')

# Test arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  stop('\nNo arguments passed to R/trace-marx.R', call.=FALSE)
} else if (length(args) != 0) {
  prns <- list.files(args[1], pattern = '.prn', full.names = T)
  if(length(prns) == 0) {
    stop(paste0('\nNo .prn files found in filepath: ', args[1]), call.=FALSE)
  }
  cores <- suppressWarnings(as.integer(args[2]))
  if(cores < 1 | cores > availableCores()) {
    cores <- availableCores()-2
    cat('\nUsing default number of cores: ', cores, sep = '')
  }
}

# Compile filepaths
files <-
  tibble(path = prns) %>%
  mutate(
    model = str_extract(path, 'cd[a-z][0-9]+'),
    tstep = as.integer(str_sub(str_extract(path, '_[0-9]+'), 2))/10
  )
files.ordered <- files[order(files$tstep),]

# Models to trace
mods <- str_extract(files.ordered$path, 'cd[a-z]{1}[0-9]{2,3}')
models <- unique(mods)
cat('\nFound .prn files for ', length(models), ' numerical experiments', sep = '')
print(files.ordered, n = 10)

# Trace markers
fun_marx <-
  function(model) {
    cat('\nTracing marx for [', model, ']', sep = '')
    fpaths <- files[files$model == model,]$path
    read_prn(
      prn.paths = fpaths,
      marx.est = 2.5e4,
      area = c(5e5, 1.26e6, 1.75e4, 2.85e4),
      markers = T,
      grid = F
    )
  }

# Set parallel plan
cat('\nParallel computing with ', cores, ' cores', sep = '')
plan(multisession, workers = cores)

# Parallel computing
future_walk(models, fun_marx, .progress = T, .options = furrr_options(seed = T))

# Write log
cat('\ntrace-marx.R complete!\n')
sink()
