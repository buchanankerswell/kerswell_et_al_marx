#!/usr/bin/env Rscript

# Load packages
# Quiet loading
sshhh <- function(p){
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, character.only=TRUE)
    )
  )}

# Package list
p.list <-
  c(
    'Hmisc',
    'scales',
    'plot3D',
    'mclust',
    'magrittr',
    'tidyr',
    'readr',
    'stringr',
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

# auto-load quietly
sapply(p.list, sshhh)

# Read binary (.prn) files and trace markers
read_prn <-
  function(
    prn.paths,
    marx.est = 2.5e4,
    area = c(500000, 1260000, 17500, 28500),
    markers = TRUE,
    grid = TRUE
  ){
  if(grid == FALSE & markers == FALSE) {
    stop('grid and markers cannot both be FALSE')
  }
  # List prn files
  fpaths <- prn.paths
  forder <-
    order(
      fpaths %>%
      map_int(
        ~.x %>%
        str_extract('[0-9]+.prn+') %>%
        str_extract('[0-9]+') %>%
        as.integer()
      )
    )
  fpaths.ordered <- fpaths[forder]
  fnames <- fpaths.ordered %>% str_extract('cd[a-z][0-9]+_[0-9]+')
  f.mod <- fpaths.ordered[1] %>% str_extract('cd[a-z][0-9]+')
  # Coordinates of the sampling area, [m]
  xmin <- area[1]
  xmax <- area[2]
  zmin <- area[3]
  zmax <- area[4]
  if(markers) {
    # Create marker arrays (n.markers x n.tsteps)
    mmm <- matrix(NA, marx.est) # marker global index
    mty <- matrix(NA, marx.est, length(fpaths)) # type
    mti <- matrix(NA, marx.est, length(fpaths)) # time, [yr]
    mxx <- matrix(NA, marx.est, length(fpaths)) # x, [m]
    mzz <- matrix(NA, marx.est, length(fpaths)) # y, [m]
    mtk <- matrix(NA, marx.est, length(fpaths)) # T, [K]
  }
  # Create grid array
  grids <- vector('list', length(fpaths))
  # Counters
  mfind <- 0 # Marker counter
  tstep <- 1 # .prn (tstep) counter
  # Main loop
  for(path in fpaths.ordered) {
    # Filename
    f <- path
    # Model name
    f.name <- str_extract(f, 'cd.[0-9]+.[0-9]+')
    f.mod <- str_extract(f, 'cd.[1-9]+')
    # Open connection
    f.prn <- file(f, 'rb')
    # Read sizes of variables
    readBin(f.prn, 'integer', 4, 1, signed = F)
    # Read model parameters
    # Grid resolution
    xnumx <- readBin(f.prn, 'integer', 1, 8)
    znumz <- readBin(f.prn, 'integer', 1, 8)
    # Markers per cell
    mnumx <- readBin(f.prn, 'integer', 1, 8)
    mnumz <- readBin(f.prn, 'integer', 1, 8)
    # Number of markers
    marknum <- readBin(f.prn, 'integer', 1, 8)
    # Model sizes
    xsize <- readBin(f.prn, 'numeric', 1, 8)
    zsize <- readBin(f.prn, 'numeric', 1, 8)
    # Pressure value
    pinit <- readBin(f.prn, 'numeric', 5, 8)
    # Gravity
    gx <- readBin(f.prn, 'numeric', 1, 8)
    gz <- readBin(f.prn, 'numeric', 1, 8)
    # Number of rocks
    rocknum <- readBin(f.prn, 'integer', 1, 4)
    # Number of Boundary conditions
    boundnum <- readBin(f.prn, 'integer', 1, 8)
    # Stage, time
    stg <- readBin(f.prn, 'integer', 1, 4)
    timesum <- readBin(f.prn, 'numeric', 1, 8)
    # Skip rock properties
    curpos <- 4+2*4+16*8+rocknum*(8*24+4)
    seek(f.prn, curpos, 'start')
    if(grid == TRUE) {
      # Initialize Matrices
      pr <- matrix(NA, znumz, xnumx) # Pressure [Pa]
      vx <- matrix(NA, znumz, xnumx) # Velocity [m/s]
      vz <- matrix(NA, znumz, xnumx) # Velocity [m/s]
      exx <- matrix(NA, znumz, xnumx) # Strain rate [1/s]
      ezz <- matrix(NA, znumz, xnumx) # Strain rate [1/s]
      exz <- matrix(NA, znumz, xnumx) # Strain rate [1/s]
      sxx <- matrix(NA, znumz, xnumx) # Stress [Pa]
      szz <- matrix(NA, znumz, xnumx) # Stress [Pa]
      sxz <- matrix(NA, znumz, xnumx) # Stress [Pa]
      ro <- matrix(NA, znumz, xnumx) # Density [kg/m^3]
      nu <- matrix(NA, znumz, xnumx) # Viscosity [Pa s]
      nd <- matrix(NA, znumz, xnumx) # ?
      mu <- matrix(NA, znumz, xnumx) # Standard viscosity for node [Pa s]
      ep <- matrix(NA, znumz, xnumx) # Surface trace
      et <- matrix(NA, znumz, xnumx) # Free array
      pr0 <- matrix(NA, znumz, xnumx) # Last cycle pressure [Pa]
      prb <- matrix(NA, znumz, xnumx) # ?
      dv <- matrix(NA, znumz, xnumx) # ?
      tk <- matrix(NA, znumz, xnumx) # Temperature [K]
      cp <- matrix(NA, znumz, xnumx) # Heat capacity [J/kg]
      kt <- matrix(NA, znumz, xnumx) # Thermal conductivity [Wt/m K]
      ht <- matrix(NA, znumz, xnumx) # Heat sources [W/m^3]
      eii <- matrix(1, znumz, xnumx)*1e-16 # Strain rate tensor [1/s]
      sii <- matrix(1, znumz, xnumx)*1e+4 # Stress tensor [Pa]
      # Progress bar
      pb.nodes <-
        progress_bar$new(
          format = paste0(
            'Reading Nodes [',
            f.name,
            '] [:bar] :percent in: :elapsed'),
          total = xnumx*znumz,
          clear = F,
          width = 100
        )
      # Read nodes information
      for(i in seq_len(xnumx)) {
        for(j in seq_len(znumz)) {
          vbuf <- readBin(f.prn, 'numeric', 3, 4)
          pr[j,i] <- vbuf[1]
          vx[j,i] <- vbuf[2]
          vz[j,i] <- vbuf[3]
          vbuf1 <- readBin(f.prn, 'integer', 3, 8)
          vbuf2 <- readBin(f.prn, 'numeric', 16, 4)
          exx[j,i] <- vbuf2[1]
          ezz[j,i] <- vbuf2[2]
          exz[j,i] <- vbuf2[3]
          sxx[j,i] <- vbuf2[4]
          szz[j,i] <- vbuf2[5]
          sxz[j,i] <- vbuf2[6]
          ro[j,i] <- vbuf2[7]
          nu[j,i] <- vbuf2[8]
          nd[j,i] <- vbuf2[9]
          mu[j,i] <- vbuf2[10]
          ep[j,i] <- vbuf2[11]
          et[j,i] <- vbuf2[12]
          pr0[j,i] <- vbuf2[13]
          prb[j,i] <- vbuf2[14]
          dv[j,i] <- vbuf2[15]
          tk[j,i] <- vbuf2[16]
          vbuf3 <- readBin(f.prn, 'integer', 1, 8)
          vbuf4 <- readBin(f.prn, 'numeric', 3, 4)
          cp[j,i] <- vbuf4[1]
          kt[j,i] <- vbuf4[2]
          ht[j,i] <- vbuf4[3]
          pb.nodes$tick()
        }
      }
    }
    # Skip all nodes
    curpos2 <- curpos+(4*22+8*4)*xnumx*znumz
    seek(f.prn, curpos2, 'start')
    # Read gridline positions
    gx <- readBin(f.prn, 'numeric', xnumx, 4)
    gz <- readBin(f.prn, 'numeric', znumz, 4)
    if(grid == TRUE) {
      # Progress bar
      pb.eii <-
        progress_bar$new(
          format =
            paste0('Stress & Strain [',
                   f.name,
                   '] [:bar] :percent in: :elapsed'),
          total = (xnumx*znumz)-2440,
          clear = FALSE,
          width = 100
        )
      for(i in seq_len(xnumx-2)) {
        for(j in seq_len(znumz-2)) {
          eii[j+1,i+1] <-
            (exz[j+1,i+1]^2+((exx[j+1,i+1]+exx[j+2,i+1]+exx[j+1,i+2]+exx[j+2,i+2])/4)^2)^0.5;
          sii[j+1,i+1] <-
            (sxz[j+1,i+1]^2+((sxx[j+1,i+1]+sxx[j+2,i+1]+sxx[j+1,i+2]+sxx[j+2,i+2])/4)^2)^0.5;
          pb.eii$tick()
        }
      }
      # Save grid
      grd.vars <-
        purrr::map(list(nu, tk), ~{
          colnames(.x) <- gx
          rownames(.x) <- gz
          .x
        })
      var.names <- c('nu', 'tk')
      # Make grid object
      assign(paste0('grid.', f.name), set_names(grd.vars, var.names))
      # Save
      grids[[tstep]] <- get(paste0('grid.', f.name))
      attr(grids[[tstep]], 'time') <- timesum
    }
    if(markers == TRUE) {
      # Skip all boundary conditions
      curpos3 <- curpos2+(xnumx+znumz)*4+(4*4+8*3)*(boundnum-1)
      # Progress bar
      pb.marx <-
        progress_bar$new(
          format =
            paste0(
              'Parsing ',
              marknum,
              ' Markers [',
              f.name,
              '] [:bar] :percent in: :elapsed'
            ),
          total = marknum,
          clear = FALSE,
          width = 100
        )
      # Find Markers
      if(mfind == 0){
        seek(f.prn, curpos3, 'start')
        for(m in seq_len(marknum)){
          mbuf <- readBin(f.prn, 'numeric', 9, 4) %>% replace(is.na(.), Inf)
          mt <- readBin(f.prn, 'integer', 1, 1, signed = F)
          mx <- mbuf[1]
          mz <- mbuf[2]
          mk <- mbuf[3]
          # Save markers from the interest area
          if(mx>=xmin & mx<=xmax & mz>=zmin & mz<=zmax & mt>1){
            mfind <- mfind+1
            mmm[mfind] <- m # global index
            mty[mfind,tstep] <- mt # type
            mti[mfind,tstep] <- timesum # time [yr]
            mxx[mfind,tstep] <- mx # x [m]
            mzz[mfind,tstep] <- mz # z [m]
            mtk[mfind,tstep] <- mk # T [K]
          }
          pb.marx$tick()
        }
        if(mfind == 0){
          break
        }
      } else {
        cat('\nUpdating PT paths of ', mfind, ' markers [', f.name, ']', sep = '')
        # Update PT-paths of selected markers
        for(mi in seq_len(mfind)) {
          m <- mmm[mi]
          curpos4 <- curpos3+(m-1)*(9*4+1)
          seek(f.prn, curpos4, 'start')
          mbuf <- readBin(f.prn, 'numeric', 9, 4) %>% replace(is.na(.), Inf)
          mt <- readBin(f.prn, 'integer', 1, 1, signed = F)
          mx <- mbuf[1]
          mz <- mbuf[2]
          mk <- mbuf[3]
          # Save marker data
          mty[mi,tstep] <- mt # type
          mti[mi,tstep] <- timesum # time [yr]
          mxx[mi,tstep] <- mx # x [m]
          mzz[mi,tstep] <- mz # z [m]
          mtk[mi,tstep] <- mk # T [K]
        }
      }
      if(mfind == 0) {
        break
      }
    }
    # Close connection
    close(f.prn)
    # .prn (tstep) counter
    tstep <- tstep + 1
  }
  # Save markers
  if(markers) {
    assign(
      paste0('marx.', f.mod),
      # Combine into one tibble
      map2(
        list(mti, mxx, mzz, mtk, mty),
        c('time', 'x', 'z', 'T', 'type'),
        ~{
          .x %>%
          as_tibble(rownames = 'id') %>%
          pivot_longer(-id, names_to = 'tstep', values_to = .y) %>%
          mutate(
            'tstep' = as.integer(str_extract(tstep, '[0-9]+')),
            'id' = as.integer(id)
          )
      }) %>%
      reduce(left_join, by = c('id', 'tstep')) %>%
      group_by(id) %>%
      drop_na()
    )
  }
  # Write objects to .RData files
  if(markers == TRUE & grid == TRUE) {
    assign(f.mod, list(grid = set_names(grids, fnames), marx = get(paste0('marx.', f.mod))))
    save(list = f.mod, file = paste0('data/marx_traced/', f.mod, '-grids-marx.RData'))
  } else if(markers == TRUE & grid == FALSE) {
    assign(f.mod, list(marx = get(paste0('marx.', f.mod))))
    save(list = f.mod, file = paste0('data/marx_traced/', f.mod, '-marx.RData'))
  } else if(markers == FALSE & grid == TRUE) {
    assign(f.mod, list(grid = set_names(grids, fnames)))
    save(list = f.mod, file = paste0('data/marx_traced/', f.mod, '-grids.RData'))
  }
}

# Load .RData file and process markers
load_marx <- function(path) {
  # Model name
  mod <- str_extract(path, 'cd.[0-9]+')
  cat('\nLoading markers [', mod, ']', sep = '')
  # Load .RData file
  load(path)
  # Transform units into Ma, km, C, and GPa
  # Compute thermal gradients in C/km
  # Compute depth from surface by subtracting 18km from marker z position
  get(mod)$marx %>%
  filter(type > 1) %>%
  group_by(id) %>%
  slice_min(type) %>%
  ungroup() %>%
  mutate(
    time = time/1e6,
    x = x/1e3,
    z = (z/1e3)-18,
    T = T-273,
    P = z*1/35,
    grad = T/z
  )
}

# Gaussian mixture modelling (Scrucca et al., 2016) to classify recovered rocks
marx_classify <-
  function(
    marx,
    model,
    gradient.threshold = 3,
    pt.path.filter = 'maxP',
    k = 14
  ) {
  # Find tstep when toe of accretionary wedge is within
  # 500 km of convergence region
  marx.500km.from.convergence.region <-
    marx %>%
    filter(tstep > 1) %>%
    group_by(id) %>%
    slice_min(z) %>%
    ungroup() %>%
    filter(z < 0 & x <= 1e3)
  # If the toe of the accretionary wedge does not reach the convergence region
  # then use the maximum timestep for the classifier
  if(nrow(marx.500km.from.convergence.region) > 0) {
    tstep.cutoff <- min(marx.500km.from.convergence.region$tstep)
  } else {
    tstep.cutoff <- max(marx$tstep)
  }
  # If the accretionary wedge reaches the convergence region before 20Ma
  # then use 20Ma as the cutoff for the classifier
  if(marx$time[marx$tstep == tstep.cutoff][1] <= 20) {
    tstep.cutoff.classify <- marx$tstep[which.min(abs(marx$time - 20))]
  } else {
    tstep.cutoff.classify <- tstep.cutoff
  }
  cat(
    '\nmodel:          ', model,
    '\ntstep cutoff:   ', tstep.cutoff,
    '\ntime cutoff:    ', round(marx$time[marx$tstep == tstep.cutoff][1], 1), ' Ma',
    '\ntstep classify: ', tstep.cutoff.classify,
    '\ntime classify:  ', round(marx$time[marx$tstep == tstep.cutoff.classify][1], 1),
    ' Ma',
    sep = ''
  )
  # Filter markers
  marx.tstep.cutoff.classify <- filter(marx, tstep <= tstep.cutoff.classify)
  # Compute features
  fts.cutoff.classify <- marx_ft(marx.tstep.cutoff.classify)
  # Save marker ids
  ids <- filter(fts.cutoff.classify, pt.path.position == pt.path.filter)$id
  # Fit Eigenvalue decomposition models and select best using BIC
  X <- select(filter(fts.cutoff.classify, pt.path.position == pt.path.filter), T, z)
  # Try clustering
  cat('\nGaussian mixture model (GMM) clustering with the mclust package ...')
  cat('\nFitting eigenvalue decomposition models and selecting best using BIC ...')
  mcl.bic <- try(mclustBIC(X, G = k, verbose = F))
  # If clustering doesn't converge or throws error
  if(class(mcl.bic) == 'try-error') {
    # Keep trying i times
    i <- 1 # Counter
    imax <- 10
    while(class(mcl.bic) == 'try-error' & i <= imax) {
      if(i < imax) {
        cat('\nPrevious decomposition failed to converge ...')
        cat('\nFitting Eigenvalue decomposition again ...')
        mcl.bic <- try(mclustBIC(X, G = k, verbose = F))
        i <- i + 1
      } else {
        stop('Clustering could not converge ...')
      }
    }
  }
  cat('\nFitting GMM ...')
  # GMM clustering using model picked by BIC
  mcl <- Mclust(X, x = mcl.bic)
  # Summarise classes
  gmm.class <- tibble(id = ids, gmm.class = mcl$classification)
  # Get parameter centroids and calculate threshold PT line
  cents <-
    t(mcl$parameters$mean) %>%
    as_tibble() %>%
    mutate(
      gmm.class = 1:n(),
      aboveGrad = T >= gradient.threshold*z,
      below1300C = T <= 1300,
      below120km = z <= 120,
      .before = 'T'
    )
  # Classify group as recovered
  recov <- filter(cents, aboveGrad & below1300C & below120km)
  # Summarise and add recovered class
  cat('\nSummarising marker features ...')
  if(tstep.cutoff == tstep.cutoff.classify) {
    fts.summary.cutoff.classify <-
      fts.cutoff.classify %>%
      left_join(gmm.class, by = 'id') %>%
      mutate(recovered = ifelse(gmm.class %in% recov$gmm.class, TRUE, FALSE)) %>%
      ungroup()
    attr(fts.summary.cutoff.classify, 'tstep.cutoff') <- tstep.cutoff.classify
    attr(fts.summary.cutoff.classify, 'time.cutoff') <-
      marx$time[marx$tstep == tstep.cutoff.classify][1]
    fts.summary.cutoff <- fts.summary.cutoff.classify
    attr(fts.summary.cutoff, 'tstep.cutoff') <- tstep.cutoff
    attr(fts.summary.cutoff, 'time.cutoff') <- marx$time[marx$tstep == tstep.cutoff][1]
    cat('\nReturning classified markers ...')
    return(fts.summary.cutoff.classify)
  } else if (tstep.cutoff != tstep.cutoff.classify) {
    fts.summary.cutoff.classify <-
      fts.cutoff.classify %>%
      left_join(gmm.class, by = 'id') %>%
      mutate(recovered = ifelse(gmm.class %in% recov$gmm.class, TRUE, FALSE)) %>%
      ungroup()
    attr(fts.summary.cutoff.classify, 'tstep.cutoff') <-
      tstep.cutoff.classify
    attr(fts.summary.cutoff.classify, 'time.cutoff') <-
      marx$time[marx$tstep == tstep.cutoff.classify][1]
    fts.summary.cutoff <-
      marx_ft(filter(marx, tstep <= tstep.cutoff)) %>%
      left_join(gmm.class, by = 'id') %>%
      mutate(recovered = ifelse(gmm.class %in% recov$gmm.class, TRUE, FALSE)) %>%
      ungroup()
    attr(fts.summary.cutoff, 'tstep.cutoff') <- tstep.cutoff
    attr(fts.summary.cutoff, 'time.cutoff') <- marx$time[marx$tstep == tstep.cutoff][1]
    cat('\nReturning classified markers ...')
    return(fts.summary.cutoff)
  }
}

# Summarizes marker features at different positions along PT path
marx_ft <- function(marx) {
  # Summarize marker features at maxT maxP and maxPT
  bind_rows(
    marx %>%
    group_by(id) %>%
    slice_max(z) %>%
    ungroup() %>%
    mutate(pt.path.position = 'maxP') %>%
    distinct(id, z, .keep_all = T),
    marx %>%
    group_by(id) %>%
    slice_max(T) %>%
    ungroup() %>%
    mutate(pt.path.position = 'maxT') %>%
    distinct(id, T, .keep_all = T),
    marx %>%
    group_by(id) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      tstep = NA,
      time = NA,
      x = NA,
      z = distinct(slice_max(group_by(marx, id), z), id, z)$z,
      T = distinct(slice_max(group_by(marx, id), T), id, T)$T,
      P = z*1/35,
      grad = T/z,
      pt.path.position = 'maxPT'
    )
  )
}

# Find peaks (modes) in density by computing second derivative
find_modes <-
  function(x, n = 10, peak.proportion.threshold = 65) {
    d <- density(x)
    modes <-
      as_tibble(data.frame(d[c('x', 'y')])[c(F, diff(diff(d$y)>=0)<0),]) %>%
      arrange(desc(y)) %>%
      mutate(mode = 1:n(), .before = x)
    modes <- filter(modes, y > max(modes$y)/peak.proportion.threshold) %>% slice(1:n)
    return(modes)
  }

# Calculate ratios, statistics, and P-T-x-z of markers
marx_stats <- function(marx.class) {
  if(nrow(marx.class[marx.class$recovered,]) == 0) {
    return(NULL)
  }
  # Compute density to find peaks in PT
  modes.P <- find_modes(marx.class$P[marx.class$recovered])
  modes.T <- find_modes(marx.class$T[marx.class$recovered])
  modes.Grad <- find_modes(marx.class$grad[marx.class$recovered])
  s <-
    marx.class %>%
    summarise(n = n(), recovered = sum(recovered, na.rm = T)) %>%
    mutate(
      rec.ratio = recovered/n,
      mode1.P = modes.P$x[which.max(modes.P$y)],
      mode2.P = modes.P$x[which.max(modes.P$x)],
      mode1.T = modes.T$x[which.max(modes.T$y)],
      mode2.T = modes.T$x[which.max(modes.T$x)],
      mode1.grad = modes.Grad$x[which.max(modes.Grad$y)],
      mode2.grad = modes.Grad$x[which.max(modes.Grad$x)]
    )
  return(s)
}

# Jacknife sampling of marx_classify()
jknife <-
  function(
    marx,
    model,
    gradient.threshold = 5,
    pt.path.filter = 'maxP',
    n = 20,
    p = 0.90,
    k = 14
  ) {
  # Jacknife sampling
  smpl <-
    purrr::map(seq_len(n),
      ~{
        cat('\nJknife sample ', .x, ' [', model, ']', sep = '')
        # Subsample
        smp.ids <- sample(unique(marx$id), length(unique(marx$id))*p)
        marx.smp <- marx[marx$id %in% smp.ids,]
        # Classify markers
        marx.class.smp <- marx_classify(marx.smp, model, gradient.threshold, pt.path.filter, k)
        # Compute summary statistics and cumulative probability curves
        marx_stats(filter(marx.class.smp, pt.path.position == pt.path.filter))
      }
    ) %>%
    compact()
  if(length(smpl) == 0) {
    return(NULL)
  }
  stats <-
    map_df(smpl, ~.x, .id = 'run') %>%
    summarise(
      n.est = mean(n, na.rm = T),
      n.sigma = sd(n, na.rm = T),
      recovered.est = mean(recovered, na.rm = T),
      recovered.sigma = sd(recovered, na.rm = T),
      rec.ratio.est = mean(rec.ratio, na.rm = T),
      rec.ratio.sigma = sd(rec.ratio, na.rm = T),
      mode1.P.est = mean(mode1.P),
      mode1.P.sigma = sd(mode1.P),
      mode2.P.est = mean(mode2.P),
      mode2.P.sigma = sd(mode2.P),
      mode1.T.est = mean(mode1.T),
      mode1.T.sigma = sd(mode1.T),
      mode2.T.est = mean(mode2.T),
      mode2.T.sigma = sd(mode2.T),
      mode1.grad.est = mean(mode1.grad),
      mode1.grad.sigma = sd(mode1.grad),
      mode2.grad.est = mean(mode2.grad),
      mode2.grad.sigma = sd(mode2.grad)
    )
  return(stats)
}

# Compute density in 2d
compute_dens_2d <- function(x, n = 100) {
  k.dens.fts <- MASS::kde2d(x$T, x$P, n = n)
  dens.summary <-
    expand.grid(k.dens.fts$x, k.dens.fts$y) %>%
    as_tibble() %>%
    rename(T = Var1, P = Var2) %>%
    mutate(
      k.dens = as.vector(k.dens.fts$z),
      cnt = nrow(x) / sum(k.dens) * k.dens
    )
  return(dens.summary)
}

# Use to flatten output of Hmisc::rcorr()
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  tibble(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut],
    p = pmat[ut]
  )
}

# Summarise density
get_dens_summary <- function(dset) {
  modes.T <- find_modes(get(dset)$T)
  modes.P <- find_modes(get(dset)$P)
  modes.grad <- find_modes(get(dset)$grad)
  tibble(
    rec.ratio.est = NA,
    mode1.P.est = modes.P$x[which.max(modes.P$y)],
    mode2.P.est = modes.P$x[which.max(modes.P$x)],
    mode1.T.est = modes.T$x[which.max(modes.T$y)],
    mode2.T.est = modes.T$x[which.max(modes.T$x)],
    mode1.grad.est = modes.grad$x[which.max(modes.grad$y)],
    mode2.grad.est = modes.grad$x[which.max(modes.grad$x)]
  ) %>%
  pivot_longer(everything(), names_to = 'var', values_to = paste0(dset, '.thresh'))
}

# Draw grids
draw_grid <-
  function(
    node,
    model = NULL,
    marx = NULL,
    time = NULL,
    box = c(up = -18, down = 300, left = 0, right = 2000),
    base.size = 12,
    p.type = 'viscosity',
    mk.p.type = 'recovered',
    cls.col = 'batlow',
    bk.alpha = 0.6,
    bk.col = rgb(1, 1, 1, 1),
    mk.size = 0.2,
    iso.alpha = 1,
    iso.col = 'black',
    iso.text.col = 'white',
    iso.skip = 2,
    iso.size = 3,
    rec.col = 'gold',
    sub.col = 'darkorange3',
    guide.pal = 'grays',
    guide.rev = TRUE,
    guide.dir = 'horizontal',
    guide.just = 'left',
    guide.pos = 'bottom',
    guide.width = 2
  ) {
  # Get time cutoff
  if(p.type == 'blank') {
    p <-
      node$tk %>%
      reshape2::melt() %>%
      rename(x = Var2, z = Var1, tk = value) %>%
      ggplot() +
      geom_contour(
        aes(x = x/1000, y = (z/1000), z = tk - 273),
        color = iso.col,
        breaks = c(0, seq(100, 1300, 200)),
        size = 0.3,
        alpha = iso.alpha
      ) +
      geom_text_contour(
        aes(x = x/1000, y = (z/1000), z = tk - 273),
        size = iso.size,
        skip = iso.skip,
        color = iso.text.col,
        breaks = c(0, seq(100, 1300, 200)),
        label.placer = label_placer_fraction(c(0.2))
      ) +
      labs(
        title = paste0('markers (', model, ') ', round(time/1e6), ' Ma'),
        x = 'distance (km)',
        y = 'distance (km)'
      ) +
      coord_equal(
        expand = F,
        xlim = c(box[3], box[4]),
        ylim = c(box[2], box[1])
      ) +
      scale_x_continuous(breaks = seq(600, 1800, 400)) +
      theme_minimal(base_size = base.size) +
      theme(
        plot.margin = margin(1, 1, 1, 1),
        legend.box.margin = margin(),
        legend.margin = margin(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'grey50'),
        legend.position = guide.pos
      )
  } else if(p.type == 'temperature') {
    p <-
      node$tk %>%
      reshape2::melt() %>%
      rename(x = Var2, z = Var1, tk = value) %>%
      ggplot() +
      geom_contour_fill(aes(x = x/1000, y = z/1000, z = tk - 273), alpha = bk.alpha) +
      geom_contour(
        aes(x = x/1000, y = (z/1000), z = tk - 273),
        color = iso.col,
        breaks = c(0, seq(100, 1300, 200)),
        size = 0.3,
        alpha = iso.alpha
      ) +
      geom_text_contour(
        aes(x = x/1000, y = (z/1000), z = tk - 273),
        size = iso.size,
        skip = iso.skip,
        color = iso.text.col,
        breaks = c(0, seq(100, 1300, 200)),
        label.placer = label_placer_fraction(c(0.2))
      ) +
      labs(
        title = paste0('temperature (', model, ') ', round(time/1e6), ' Ma'),
        x = 'distance (km)',
        y = 'distance (km)'
      ) +
      scale_fill_continuous_sequential(
        guide.pal,
        name = 'temperature (˚C)',
        rev = guide.rev,
        breaks = MakeBreaks(bins = 6),
        na.value = 'grey50',
        guide =
          guide_colorstrip(
            inside = T,
            title.position = ifelse(guide.dir == 'horizontal', 'top', 'top'),
            title.vjust = ifelse(guide.dir == 'horizontal', 1, 0),
            barwidth = unit(guide.width, 'in')
          )
      ) +
      coord_equal(
        expand = F,
        xlim = c(box[3], box[4]),
        ylim = c(box[2], box[1])
      ) +
      theme_minimal(base_size = base.size) +
      theme(
        plot.margin = margin(1, 1, 1, 1),
        legend.box.margin = margin(),
        legend.margin = margin(),
        panel.grid = element_blank(),
        legend.position = guide.pos
      )
  } else if(p.type == 'viscosity') {
    p <-
      node$tk %>%
      reshape2::melt() %>%
      rename(x = Var2, z = Var1, tk = value) %>%
      ggplot() +
      geom_contour_fill(
        data = rename(reshape2::melt(node$nu), x = Var2, z = Var1, nu = value),
        aes(x = x/1000, y = (z/1000), z = log10(nu)),
        alpha = bk.alpha
      ) +
      geom_contour(
        aes(x = x/1000, y = (z/1000), z = tk - 273),
        color = iso.col,
        breaks = c(0, seq(100, 1300, 200)),
        size = 0.3,
        alpha = iso.alpha
      ) +
      geom_text_contour(
        aes(x = x/1000, y = (z/1000), z = tk - 273),
        size = iso.size,
        skip = iso.skip,
        color = iso.text.col,
        breaks = c(0, seq(100, 1300, 200)),
        label.placer = label_placer_fraction(c(0.2))
      ) +
      labs(
        title = paste0('log viscosity (', model, ') ', round(time), ' Ma'),
        x = 'distance (km)',
        y = 'distance (km)'
      ) +
      scale_fill_continuous_sequential(
        guide.pal,
        name = 'log10 viscosity (Pa s)',
        rev = guide.rev,
        breaks = MakeBreaks(bins = 6),
        na.value = 'grey50',
        guide =
          guide_colorstrip(
            inside = T,
            title.position = ifelse(guide.dir == 'horizontal', 'top', 'top'),
            title.vjust = ifelse(guide.dir == 'horizontal', 1, 0),
            barwidth = unit(guide.width, 'in')
          )
      ) +
      coord_equal(
        expand = F,
        xlim = c(box[3], box[4]),
        ylim = c(box[2], box[1])
      ) +
      scale_x_continuous(breaks = seq(600, 1800, 400)) +
      theme_minimal(base_size = base.size) +
      theme(
        plot.margin = margin(1, 1, 1, 1),
        legend.box.margin = margin(),
        legend.margin = margin(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = 'grey50'),
        legend.position = guide.pos
      )
  }
  if(!is.null(marx)) {
    if(mk.p.type == 'recovered') {
      if(nrow(marx[marx$time == time & !marx$recovered,]) > 0) {
        p <-
          p +
          geom_point(
            data = marx[marx$time == time & !marx$recovered,],
            aes(x = x, y = z+18, color = 'no'),
            shape = 15,
            size = mk.size
          )
      }
      if(nrow(marx[marx$time == time & marx$recovered,]) > 0) {
        p <-
          p +
          geom_point(
            data = marx[marx$time == time & marx$recovered,],
            aes(x = x, y = z+18, color = 'yes'),
            shape = 15,
            size = mk.size
          ) +
          geom_contour(
            aes(x = x/1000, y = (z/1000), z = tk - 273),
            color = iso.col,
            breaks = c(0, seq(100, 1300, 200)),
            size = 0.3,
            alpha = iso.alpha
          ) +
          geom_text_contour(
            aes(x = x/1000, y = (z/1000), z = tk - 273),
            size = iso.size,
            skip = iso.skip,
            color = iso.text.col,
            breaks = c(0, seq(100, 1300, 200)),
            label.placer = label_placer_fraction(c(0.2))
          ) +
          scale_color_manual(
            name = 'recovered',
            values = c('yes' = rec.col, 'no' = sub.col),
            guide =
              guide_legend(
                direction = guide.dir,
                title.position = ifelse(guide.dir == 'horizontal', 'top', 'top'),
                label.position = ifelse(guide.dir == 'horizontal', 'right', 'right'),
                title.vjust = ifelse(guide.dir == 'horizontal', 0, 0),
                override.aes = list(alpha = 1, size = 5)
              )
          )
      }
    } else if(mk.p.type == 'class') {
      if(nrow(marx[marx$time == time & !marx$recovered,]) > 0) {
        p <-
          p +
          geom_point(
            data = marx[marx$time == time & !marx$recovered,],
            aes(x = x, y = z+18, color = gmm.class),
            shape = 15,
            size = mk.size
          )
      }
      if(nrow(marx[marx$time == time & marx$recovered,]) > 0) {
        p <-
          p +
          geom_point(
            data = marx[marx$time == time & marx$recovered,],
            aes(x = x, y = z+18, color = gmm.class),
            shape = 15,
            size = mk.size
          ) +
          geom_contour(
            aes(x = x/1000, y = (z/1000), z = tk - 273),
            color = iso.col,
            breaks = c(0, seq(100, 1300, 200)),
            size = 0.3,
            alpha = iso.alpha
          ) +
          geom_text_contour(
            aes(x = x/1000, y = (z/1000), z = tk - 273),
            size = iso.size,
            skip = iso.skip,
            color = iso.text.col,
            breaks = c(0, seq(100, 1300, 200)),
            label.placer = label_placer_fraction(c(0.2))
          ) +
          scale_color_continuous_sequential(
            cls.col,
            rev = F,
            na.value = 'transparent',
            name = 'classifier cluster groups',
            guide =
              guide_colorstrip(
                inside = T,
                title.position = 'top',
                title.vjust = 1,
                barwidth = unit(3, 'in')
              )
          )
      }
    } else {
        p <- p
    }
  }
  return(p)
}
