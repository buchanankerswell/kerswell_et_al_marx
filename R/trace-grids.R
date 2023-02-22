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

# Get model names
models <- prns %>% str_extract('cd.[0-9]+') %>% unique()
cat('\nFound .prn files for ', length(models), ' numerical experiments', sep = '')
writeLines(models)

# Trace grids
fun_grids <- function(model) {
  cat('\nTracing grids for [', model, ']', sep = '')
  # File paths (prn binaries)
  prns <- list.files(args[1], pattern = model, full.names = T)
  # Compile filepaths
  files <-
    tibble(path = prns) %>%
    mutate(
      model = str_extract(path, 'cd[a-z][0-9]+'),
      tstep = as.integer(str_sub(str_extract(path, '_[0-9]+'), 2))/10
    )
  files.ordered <- files[order(files$tstep),]
  # Get timesteps
  marx <- load_marx(paste0('data/marx_traced/', model, '-marx.RData'))
  times <-
    tibble(tstep = marx$tstep-1, time = marx$time/1e6) %>%
    slice(
      c(
        which.min(abs(time - 5)),
        which.min(abs(time - 10)),
        which.min(abs(time - 15)),
        which.min(abs(time - 20))
      )
    )
  fpaths <- files.ordered %>% filter(tstep %in% times$tstep)
  read_prn(
    prn.paths = fpaths$path,
    marx.est = 2.5e4,
    area = c(5e5, 1.26e6, 1.75e4, 2.85e4),
    markers = F,
    grid = T
  )
}

# Set parallel plan
cat('\nParallel computing with ', cores, ' cores', sep = '')
plan(multisession, workers = cores)

# Parallel computing
future_walk(models, fun_grids, .progress = T, .options = furrr_options(seed = T))

# Write log
cat('\ntrace-grids.R complete!\n')
sink()
