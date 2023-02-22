#!/usr/bin/env Rscript

# Helper function
library(magrittr)
drive_upload_folder <- function(folder, drive_path) {
  # Make figs dir in drive
  googledrive::drive_mkdir('figs', drive_path, overwrite = T)
  # Get contents of local folder
  contents <- fs::dir_info(folder, type = c("file", "dir"))
  contents.files <- dplyr::filter(contents, type == 'file')$path
  # Upload files
  purrr::walk(contents.files,
    ~googledrive::drive_upload(x, path = paste0(drive_path, '/figs/'))
  )
  # Upload files recursively within folders
  contents.dirs <- dplyr::filter(contents, type == 'directory')$path
  purrr::walk(contents.dirs, function(x) {
    googledrive::drive_mkdir(
      paste0(stringr::str_remove(x, paste0(folder, '/'))),
      paste0(drive_path, '/figs/'),
      overwrite = T
    )
    contents <- fs::dir_info(x, type = c("file", "dir"))
    contents.files <- dplyr::filter(contents, type == 'file')$path
    purrr::walk(contents.files, function(y) {
      googledrive::drive_upload(
        y,
        path = paste0(drive_path, '/figs/', stringr::str_remove(y, paste0(folder, '/')))
      )
    })
  })
}

options(gargle_oauth_email = 'buck.kerswell@gmail.com')

# Test arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
path_output <- args[2]
path_figs <- args[3]
if(args[4] == 'upload') {
  trackdown::upload_file(
    file = file,
    hide_code = T,
    path_output = path_output,
    rich_text_par = list(rgb_color = list(red = 217/255, green = 234/255, blue = 211/255))
  )
  drive_upload_folder(args[3], 'trackdown')
} else if(args[4] == 'update') {
  trackdown::update_file(
    file = file,
    hide_code = T,
    path_output = path_output,
    rich_text_par = list(rgb_color = list(red = 217/255, green = 234/255, blue = 211/255))
  )
} else if (args[4] == 'download') {
  trackdown::download_file(file = file)
} else if(is.null(args[4])) {
  stop('please specify "upload", "update", or "download" as the second Rscript argument')
} else if(!is.null(args[4]) && !(args[4] %in% c('upload', 'update'))) {
  stop('please specify "upload", "update", or "download" as the second Rscript argument')
}