#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Load functions and libraries
cat(rep('~', 80), sep='')
cat('\nLoading packages and functions ...\n\n')
source('R/functions.R')
load('data/preprocessed.RData')

# Test arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  cat('\nNo arguments passed to R/marx-vis.R')
  cores <- availableCores()-2
  gradient.threshold <- 3
  pt.path.filter <- 'maxP'
  cat('\nUsing defaults')
  cat('\nCores for parallel computing     :', cores)
  cat('\nThermal gradient threshold       :', gradient.threshold)
  cat('\nPT path filter                   :', pt.path.filter, '\n')
} else if (length(args) != 0) {
  cores <- suppressWarnings(as.integer(args[1]))
  if(cores < 1 | cores > availableCores()) {
    cores <- availableCores()-2
    cat('\nCores for parallel computing must be between [0-', availableCores(), ']', sep = '')
    cat('\nDefaulting to', cores, '\n')
  }
  gradient.threshold <- suppressWarnings(as.numeric(args[2]))
  if(gradient.threshold < 0 | gradient.threshold > 5) {
    gradient.threshold <- 3
    cat('\nThreshold thermal gradient must be between [0-5] C/km')
    cat('\nDefaulting to', gradient.threshold, '\n')
  }
  pt.path.filter <- suppressWarnings(args[3])
  if(!(pt.path.filter %in% c('maxT', 'maxP', 'maxPT'))) {
    pt.path.filter <- 'maxP'
    cat('\nUnrecognized PT path filter. Please use "maxT", "maxP", or "maxPT"')
    cat('\nDefaulting to', pt.path.filter, '\n')
  }
}

# Load marker and grid data
paths_marx <- list.files('data/marx_classified', pattern = 'classified.RData', full.names = T)
paths_grids <- list.files('data/marx_traced', pattern = 'grids.RData', full.names = T)
models <- str_extract(paths_marx, 'cd.[0-9]+')

# Test input files
if(length(paths_marx) != length(paths_grids)) {
  stop('Number of classified marx and grid files are different!')
} else if(length(paths_marx) < 1) {
  stop(
    paste0('No classified marx files found in data/marx_classified !!')
  )
} else if(length(paths_grids) < 1) {
  stop('No grids files found in data/marx_traced/!')
}

# Main plotting function to parallelize
fun <- function(model, path_marx, path_grid) {
  cat('\nPlotting figures for model [', model, ']', sep = '')
  library(patchwork)
  # Load classified markers data
  # Filter out rock types after dehydration/melting
  load(path_marx)
  time.cutoff <- attr(get(paste0(model, '.marx.classified'))$marx.class, 'time.cutoff')
  marx.class <-
    get(paste0(model, '.marx.classified'))$marx.class %>%
    mutate(time = ifelse(is.na(time), time.cutoff, time)) %>%
    filter(time <= time.cutoff)
  marx.class.filtered.pt.path <- filter(marx.class, pt.path.position == pt.path.filter)
  # Summarise marker cdfP
  cdfP.marx <-
    marx.class.filtered.pt.path %>%
    filter(recovered) %>%
    select(P) %>%
    arrange(P) %>%
    mutate(prob = (row_number()-1)/n())
  # Load marx data
  marx <-
    load_marx(paste0('data/marx_traced/', model, '-marx.RData')) %>%
    left_join(select(marx.class.filtered.pt.path, id, recovered, gmm.class), by = 'id')
  sample.size <- 200
  if (
    sample.size > nrow(marx.class[marx.class$recovered,]) |
    sample.size > nrow(marx.class[!marx.class$recovered,])
  ) {
    sample.size <-
      min(nrow(marx.class[marx.class$recovered,]), nrow(marx.class[!marx.class$recovered,]))
  }
  marx.class.sample <-
    bind_rows(
      marx.class %>%
      filter(id %in% sample(marx.class$id[marx.class$recovered], sample.size)) %>%
      mutate(recovered = ifelse(recovered, 'yes', 'no')),
      marx.class %>%
      filter(id %in% sample(marx.class$id[!marx.class$recovered], sample.size)) %>%
      mutate(recovered = ifelse(recovered, 'yes', 'no'))
    ) %>%
    arrange(id, tstep)
  # Compute 2d marker density
  if(nrow(marx.class.filtered.pt.path[marx.class.filtered.pt.path$recovered,]) > 0){
    k.dens.fts <-
      MASS::kde2d(
        marx.class.filtered.pt.path$T[marx.class.filtered.pt.path$recovered],
        marx.class.filtered.pt.path$P[marx.class.filtered.pt.path$recovered],
        n = 30
      )
    fts.dens <-
      expand.grid(k.dens.fts$x, k.dens.fts$y) %>%
      as_tibble() %>%
      rename(T = Var1, P = Var2) %>%
      mutate(
        k.dens = as.vector(k.dens.fts$z),
        cnt =
          nrow(marx.class.filtered.pt.path[marx.class.filtered.pt.path$recovered,]) /
          sum(k.dens) * k.dens
      )
  } else {
    fts.dens <- NULL
  }
  # Load nodes data (grids)
  load(path_grid)
  # Nodes data
  grid1 <- get(paste0(model))$grid[[1]]
  grid2 <- get(paste0(model))$grid[[2]]
  grid3 <- get(paste0(model))$grid[[3]]
  grid4 <- get(paste0(model))$grid[[4]]
  tsteps <- names(get(paste0(model))$grid)
  tstep1 <- as.integer(str_sub(str_extract(tsteps[1], '_[0-9]+'), 2))
  tstep2 <- as.integer(str_sub(str_extract(tsteps[2], '_[0-9]+'), 2))
  tstep3 <- as.integer(str_sub(str_extract(tsteps[3], '_[0-9]+'), 2))
  tstep4 <- as.integer(str_sub(str_extract(tsteps[4], '_[0-9]+'), 2))
  time1 <- attr(grid1, 'time')/1e6
  time2 <- attr(grid2, 'time')/1e6
  time3 <- attr(grid3, 'time')/1e6
  time4 <- attr(grid4, 'time')/1e6
  # Contour limits
  break.points.marx <- MakeBreaks(bins = 6)
  break.points.ag18 <- MakeBreaks(bins = 5)
  break.points.pd15 <- MakeBreaks(binwidth = 5)
  break.lims.marx <- c(round(max(fts.dens$cnt)/10), NA)
  break.lims.ag18 <- c(round(max(ag18.dens$cnt)/6), NA)
  break.lims.pd15 <- c(round(max(pd15.dens$cnt)/12), NA)
  if(break.lims.marx[1] < 1) {break.lims.marx[1] <- 1}
  # Colors
  ag18.col <- 'grays'
  pd15.col <- 'burg'
  marx.col <- 'emrld'
  cls.col <- 'viridis'
  rec.col <- 'gold'
  sub.col <- 'gold4'
  guide.pal <- 'grays'
  base.size <- 16
  # Max PT conditions plot
  n <- 10
  thresh <- 65
  d.marx.P <- density(marx.class.filtered.pt.path$P[marx.class.filtered.pt.path$recovered])
  d.pd15.P <- density(pd15$P)
  d.ag18.P <- density(ag18$P)
  d.marx.T <- density(marx.class.filtered.pt.path$T[marx.class.filtered.pt.path$recovered])
  d.pd15.T <- density(pd15$T)
  d.ag18.T <- density(ag18$T)
  p1 <-
    ggplot() +
    geom_contour_fill(
      data = ag18.dens,
      aes(T, P, z = cnt),
      breaks = break.points.ag18
    ) +
    scale_fill_continuous_sequential(
      ag18.col,
      breaks = break.points.ag18,
      limits = break.lims.ag18,
      rev = F,
      na.value = 'transparent',
      name = 'ag18',
      guide =
        guide_colorstrip(
          inside = T,
          label = F,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(0.5, 'in'),
          order = 3
        )
    ) +
    geom_contour2(
      data = pd15.dens,
      aes(T, P, z = cnt, color = ..level..),
      size = 1,
      breaks = break.points.pd15,
    ) +
    scale_color_continuous_sequential(
      pd15.col,
      breaks = break.points.pd15,
      limits = break.lims.pd15,
      rev = F,
      na.value = 'transparent',
      name = 'pd15',
      guide =
        guide_colorstrip(
          inside = T,
          label = F,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(0.5, 'in'),
          order = 2
        )
    ) +
    new_scale_fill() +
    geom_point(
      data = ungroup(filter(marx.class.filtered.pt.path, recovered)),
      aes(T, P),
      size = 1
    ) +
    geom_contour_fill(
      data = fts.dens,
      aes(T, P, z = cnt)
    ) +
    geom_contour_tanaka(
      data = fts.dens,
      aes(T, P, z = cnt),
      smooth = 1,
      sun.angle = 30
    ) +
    scale_fill_continuous_sequential(
      marx.col,
      breaks = break.points.marx,
      limits = break.lims.marx,
      na.value = 'transparent',
      rev = F,
      name =
        paste0('marker count n=', nrow(filter(marx.class.filtered.pt.path, recovered))),
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(2.8, 'in'),
          order = 1
        )
    ) +
    geom_abline(slope = 1/700, intercept = 0, linewidth = 0.25) +
    geom_abline(slope = 1/175, intercept = 0, linewidth = 0.25) +
    geom_path(data = pl, aes(T, P), size = 1, color = 'white') +
    geom_path(data = atg, aes(T, P), size = 1, color = 'white') +
    annotate(
      'text',
      x = head(atg,1)$T,
      y = head(atg,1)$P,
      label = 'antigorite',
      angle = 50,
      vjust = -0.3,
      hjust = -0.1,
      size = 5,
      color = 'white'
    ) +
    annotate(
      'text',
      x = tail(pl,1)$T,
      y = tail(pl,1)$P,
      label = 'eclogite',
      angle = 24,
      vjust = -0.25,
      hjust = 1.3,
      size = 5,
      color = 'white'
    ) +
    annotate(
      'label',
      x = 1000,
      y = 1000/700,
      label = '20',
      size = 5,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    annotate(
      'label',
      x = 175*4,
      y = 4,
      label = '5',
      size = 5,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    annotate(
      'label',
      label = 'a',
      x = Inf,
      y = Inf,
      size = 9,
      hjust = 1,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 4)) +
    scale_x_continuous(breaks = seq(0, 1000, 200)) +
    labs(x = 'temperature (˚C)', y = 'pressure (GPa)') +
    theme_dark(base_size = base.size) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      legend.spacing.x = unit(0.06, 'in'),
      panel.grid = element_blank()
    )
  p2 <-
    ggplot() +
    geom_line(
      data = arrange(pd15, P) %>% mutate(cdfP = (row_number()-1)/n()),
      aes(P, cdfP, color = 'pd15')
    ) +
    geom_line(
      data = arrange(ag18, P) %>% mutate(cdfP = (row_number()-1)/n()),
      aes(P, cdfP, color = 'ag18')
    ) +
    geom_line(
      data = drop_na(cdfP.marx),
      aes(P, prob, color = 'markers')
    ) +
    labs(y = 'CDF', x = 'P (GPa)') +
    guides(color = 'none') +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 4)) +
    scale_color_manual(
      values = c(
        'pd15' = sequential_hcl(8, palette = pd15.col)[8],
        'ag18' = sequential_hcl(8, palette = ag18.col)[1],
        'markers' = sequential_hcl(8, palette = marx.col)[8]
      )
    ) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), position = 'right') +
    scale_x_continuous(breaks = seq(0, 4, 1)) +
    theme_dark(base_size = 9) +
    theme(
      axis.text = element_text(color = 'white'),
      axis.ticks = element_line(color = 'white'),
      axis.title = element_text(color = 'white'),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  p3 <-
    ggplot() +
    geom_path(
      data = data.frame(d.pd15.P[c('x', 'y')]),
      aes(x, y),
      color = sequential_hcl(8, palette = pd15.col)[8]
    ) +
    geom_path(
      data = data.frame(d.ag18.P[c('x', 'y')]),
      aes(x, y),
      color = sequential_hcl(8, palette = ag18.col)[1]
    ) +
    geom_path(
      data = data.frame(d.marx.P[c('x', 'y')]),
      aes(x, y),
      color = sequential_hcl(8, palette = marx.col)[8]
    ) +
    labs(x = 'P (GPa)', y = NULL) +
    coord_cartesian(xlim = c(0, 4)) +
    scale_x_continuous(breaks = seq(0, 4, 1)) +
    scale_y_continuous(position = 'right') +
    theme_dark(base_size = 9) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(color = 'white'),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  p3.b <-
    ggplot() +
    geom_path(
      data = data.frame(d.pd15.T[c('x', 'y')]),
      aes(x, y),
      color = sequential_hcl(8, palette = pd15.col)[8]
    ) +
    geom_path(
      data = data.frame(d.ag18.T[c('x', 'y')]),
      aes(x, y),
      color = sequential_hcl(8, palette = ag18.col)[1]
    ) +
    geom_path(
      data = data.frame(d.marx.T[c('x', 'y')]),
      aes(x, y),
      color = sequential_hcl(8, palette = marx.col)[8]
    ) +
    labs(x = 'T (˚C)', y = 'PDF') +
    coord_cartesian(xlim = c(0, 1000)) +
    scale_x_continuous(breaks = seq(0, 1000, 250)) +
    scale_y_continuous(position = 'right') +
    theme_dark(base_size = 10) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = 'white'),
      axis.ticks.x = element_line(color = 'white'),
      axis.title.x = element_text(color = 'white'),
      axis.title.y = element_text(color = 'white'),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  pp1 <-
    p1 +
    inset_element(p2, 0.11, 0.55, 0.394, 0.805, align_to = 'full') +
    inset_element(p3, 0.11, 0.795, 0.363, 0.995, align_to = 'full') +
    inset_element(p3.b, 0.354, 0.734, 0.63, 0.995, align_to = 'full') +
    plot_layout(guides = 'collect') &
    theme(
      plot.margin = margin(2, 2, 2, 2),
      legend.box.margin = margin(),
      legend.margin = margin(),
      legend.justification = 'top',
      legend.box.just = 'top',
      legend.direction = 'horizontal',
      legend.position = 'bottom',
      legend.spacing.x = unit(0.06, 'in')
    )
  p4 <-
    draw_grid(
      node = grid1,
      model = model,
      marx = marx,
      time = time1,
      box = c(0, 300, 500, 2000),
      p.type = 'viscosity',
      iso.text.col = 'white',
      rec.col = rec.col,
      sub.col = sub.col,
      guide.pal = guide.pal,
      base.size = base.size
    ) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
    ) +
    annotate(
      'label',
      x = Inf,
      y = Inf,
      size = 5,
      hjust = 1,
      vjust = 0,
      label = paste0(round(time1, 1), ' Ma'),
      fill = ifelse(time.cutoff >= time1, 'white', 'red4'),
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    )
  p5 <-
    draw_grid(
      node = grid2,
      model = model,
      marx = marx,
      time = time2,
      box = c(0, 300, 500, 2000),
      p.type = 'viscosity',
      iso.text.col = 'white',
      rec.col = rec.col,
      sub.col = sub.col,
      guide.pal = guide.pal,
      base.size = base.size
    ) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
    ) +
    annotate(
      'label',
      x = Inf,
      y = Inf,
      size = 5,
      hjust = 1,
      vjust = 0,
      label = paste0(round(time2, 1), ' Ma'),
      fill = ifelse(time.cutoff >= time2, 'white', 'red4'),
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    )
  p6 <-
    draw_grid(
      node = grid3,
      model = model,
      marx = marx,
      time = time3,
      box = c(0, 300, 500, 2000),
      p.type = 'viscosity',
      iso.text.col = 'white',
      rec.col = rec.col,
      sub.col = sub.col,
      guide.pal = guide.pal,
      base.size = base.size
    ) +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
    ) +
    annotate(
      'label',
      x = Inf,
      y = Inf,
      size = 5,
      hjust = 1,
      vjust = 0,
      label = paste0(round(time3, 1), ' Ma'),
      fill = ifelse(time.cutoff >= time3, 'white', 'red4'),
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    )
  p7 <-
    draw_grid(
      node = grid4,
      model = model,
      marx = marx,
      time = time4,
      box = c(0, 300, 500, 2000),
      p.type = 'viscosity',
      rec.col = rec.col,
      sub.col = sub.col,
      guide.pal = guide.pal,
      base.size = base.size
    ) +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_blank()
    ) +
    annotate(
      'label',
      x = Inf,
      y = Inf,
      size = 5,
      hjust = 1,
      vjust = 0,
      label = paste0(round(time4, 1), ' Ma'),
      fill = ifelse(time.cutoff >= time4, 'white', 'red4'),
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    )
  # Composition
  pp2 <-
    (p4 +
      scale_y_continuous(breaks = seq(100, 300, 100)) +
      theme(axis.text.x = element_blank()) +
      annotate(
        'label',
        label = 'b',
        x = Inf,
        y = -Inf,
        size = 9,
        hjust = 1,
        vjust = 1,
        fill = 'grey90',
        label.padding = unit(0.02, 'in'),
        label.r = unit(0, 'in')
      ) +
      guides(color = 'none', fill = 'none')
    ) /
    (p5 +
      guides(color = 'none', fill = 'none') +
      scale_y_continuous(breaks = seq(100, 300, 100)) +
      theme(axis.text.x = element_blank())
    ) /
    (p6 +
      guides(color = 'none', fill = 'none') +
      scale_y_continuous(breaks = seq(100, 300, 100)) +
      theme(axis.text.x = element_blank())
    ) /
    (p7 +
      scale_y_continuous(breaks = seq(100, 300, 100))
    ) &
    theme(
      plot.margin = margin(1, 1, 1, 1),
      legend.margin = margin(),
      legend.box.margin = margin(),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      legend.key = element_rect(color = NA, fill = 'grey50')
    )
  pp3 <- 
    (pp1 | pp2) &
    theme(
      legend.margin = margin(),
      legend.box.margin = margin(),
      legend.box.just = 'top',
      legend.justification = 'top',
      plot.margin = margin(2, 2, 2, 2),
      legend.position = 'bottom'
    )
  suppressWarnings(
    ggsave(
      paste0('figs/marx_comp', '/', model, '-marx-comp.png'),
      plot = pp3,
      device = 'png',
      width = 9,
      height = 5
    )
  )
  p8 <-
    ggplot() +
    geom_path(
      data = 
        marx[
          marx$id %in% filter(marx.class.sample, recovered == 'no')$id &
          marx$time <= time.cutoff
        ,],
      aes(T, P, color = 'no', group = id)
    ) +
    geom_path(
      data =
        marx[
          marx$id %in% filter(marx.class.sample, recovered == 'yes')$id &
          marx$time <= time.cutoff
        ,],
      aes(T, P, color = 'yes', group = id)
    ) +
    scale_color_manual(
      name = 'recovered',
      values = c('yes' = rec.col, 'no' = sub.col),
      guide =
        guide_legend(
          override.aes = list(linewidth = 5),
          title.position = 'top',
          label.position = 'right',
          title.vjust = 0,
          order = 2
        )
    ) +
    new_scale_color() +
    geom_point(
      data = filter(marx.class.sample, recovered == 'yes' & pt.path.position != 'maxPT'),
      aes(T, P, color = pt.path.position),
      size = 0.8
    ) +
    scale_color_manual(
      name = 'PT path position',
      values = c('maxT' = 'white', 'maxP' = 'black', 'maxPT' = 'deeppink'),
      guide =
        guide_legend(
          override.aes = list(size = 5),
          title.position = 'top',
          label.position = 'right',
          title.vjust = 0,
          order = 1
        )
    ) +
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 7)) +
    geom_abline(slope = 1/700, intercept = 0, linewidth = 0.25) +
    geom_abline(slope = 1/175, intercept = 0, linewidth = 0.25) +
    geom_path(data = pl, aes(T, P), size = 1, color = 'white') +
    geom_path(data = atg, aes(T, P), size = 1, color = 'white') +
    annotate(
      'text',
      x = head(atg,1)$T,
      y = head(atg,1)$P,
      label = 'antigorite',
      angle = 35,
      vjust = -0.2,
      hjust = -0.1,
      size = 5,
      color = 'white'
    ) +
    annotate(
      'text',
      x = tail(pl,1)$T,
      y = tail(pl,1)$P,
      label = 'eclogite',
      angle = 14,
      vjust = -0.4,
      hjust = 1.3,
      size = 5,
      color = 'white'
    ) +
    annotate(
      'label',
      x = 1000,
      y = 1000/700,
      label = '20',
      size = 5,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    annotate(
      'label',
      x = 175*4,
      y = 4,
      label = '5',
      size = 5,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    annotate(
      'label',
      label = 'b',
      x = Inf,
      y = Inf,
      size = 9,
      hjust = 1,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    labs(x = 'temperature (˚C)', y = 'pressure (GPa)') +
    theme_dark(base_size = base.size) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      legend.spacing.x = unit(0.06, 'in'),
      panel.grid = element_blank()
    )
  p9 <-
    ggplot() +
    geom_density(
      data = filter(marx.class.sample, recovered == 'yes' & pt.path.position != 'maxPT'),
      aes(P, color = pt.path.position, group = pt.path.position),
      fill = NA,
      position = 'identity',
      alpha = 0.6
    ) +
    scale_x_continuous(
      position = 'top',
      labels = function(x){sprintf('%.0f', x)},
      limits = c(0, 3),
      breaks = seq(0, 3, 1)
    ) +
    labs(x = 'P (GPa)', y = 'PDF') +
    scale_color_manual(
      name = 'PT path position',
      values = c('maxT' = 'white', 'maxP' = 'black'),
      guide = 'none'
    ) +
    coord_flip() +
    theme_dark(base_size = 9) +
    theme(
      axis.text = element_text(color = 'white'),
      axis.title = element_text(color = 'white'),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  p10 <-
    ggplot() +
    geom_density(
      data = filter(marx.class.sample, recovered == 'yes' & pt.path.position != 'maxPT'),
      aes(T, color = pt.path.position, group = pt.path.position),
      fill = NA,
      position = 'identity',
      alpha = 0.6
    ) +
    scale_x_continuous(
      position = 'top',
      labels = function(x){sprintf('%.0f', x)},
      limits = c(0, 800),
      breaks = seq(0, 800, 200)
    ) +
    labs(x = 'T (˚C)', y = 'PDF') +
    scale_color_manual(
      name = 'PT path position',
      values = c('maxT' = 'white', 'maxP' = 'black'),
      guide = 'none'
    ) +
    coord_flip() +
    theme_dark(base_size = 9) +
    theme(
      axis.text = element_text(color = 'white'),
      axis.ticks = element_blank(),
      axis.title = element_text(color = 'white'),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  pp4 <-
    p8 +
    ggtitle('classification') +
    inset_element(p9, 0.11, 0.525, 0.404, 0.755, align_to = 'full') +
    inset_element(p10, 0.11, 0.745, 0.428, 0.945, align_to = 'full') +
    plot_layout(guides = 'collect') &
    theme(
      plot.margin = margin(2, 2, 2, 2),
      plot.title = element_text(margin = margin()),
      legend.box.margin = margin(),
      legend.margin = margin(),
      legend.justification = 'top',
      legend.box.just = 'top',
      legend.direction = 'horizontal',
      legend.position = 'bottom',
      legend.spacing.x = unit(0.06, 'in')
    )
  p11 <-
    ggplot() +
    geom_path(
      data =
        marx[
          marx$id %in% filter(marx.class.sample, recovered == 'no')$id &
          marx$time <= time.cutoff
        ,],
      aes(T, P, color = gmm.class, group = id)
    ) +
    geom_path(
      data =
        marx[
          marx$id %in% filter(marx.class.sample, recovered == 'yes')$id &
          marx$time <= time.cutoff
        ,],
      aes(T, P, color = gmm.class, group = id)
    ) +
    scale_color_continuous_sequential(
      cls.col,
      breaks = MakeBreaks(bins = length(unique(marx.class.sample$gmm.class))),
      rev = F,
      na.value = 'transparent',
      name = 'classifier clustering groups',
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(3, 'in'),
          order = 1
        )
    ) +
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 7)) +
    geom_abline(slope = 1/700, intercept = 0, linewidth = 0.25) +
    geom_abline(slope = 1/175, intercept = 0, linewidth = 0.25) +
    geom_path(data = pl, aes(T, P), size = 1, color = 'white') +
    geom_path(data = atg, aes(T, P), size = 1, color = 'white') +
    annotate(
      'text',
      x = head(atg,1)$T,
      y = head(atg,1)$P,
      label = 'antigorite',
      angle = 35,
      vjust = -0.2,
      hjust = -0.1,
      size = 5,
      color = 'white'
    ) +
    annotate(
      'text',
      x = tail(pl,1)$T,
      y = tail(pl,1)$P,
      label = 'eclogite',
      angle = 14,
      vjust = -0.4,
      hjust = 1.3,
      size = 5,
      color = 'white'
    ) +
    annotate(
      'label',
      x = 1000,
      y = 1000/700,
      label = '20',
      size = 5,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    annotate(
      'label',
      x = 175*4,
      y = 4,
      label = '5',
      size = 5,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    annotate(
      'label',
      label = 'a',
      x = Inf,
      y = Inf,
      size = 9,
      hjust = 1,
      vjust = 1,
      fill = 'grey90',
      label.padding = unit(0.02, 'in'),
      label.r = unit(0, 'in')
    ) +
    labs(x = 'temperature (˚C)', y = 'pressure (GPa)') +
    theme_dark(base_size = base.size) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      legend.spacing.x = unit(0.06, 'in'),
      panel.grid = element_blank()
    )
  p12 <-
    ggplot() +
    geom_hline(yintercept = gradient.threshold, color = 'grey90') +
    annotate(
      'text',
      x = -Inf,
      y = gradient.threshold,
      vjust = -0.2,
      hjust = 0,
      size = 2.5,
      color = 'grey90',
      label = 'classifier threshold'
    ) +
    geom_boxplot(
      data = filter(marx.class.sample, pt.path.position == pt.path.filter),
      aes(gmm.class, grad, fill = gmm.class, color = gmm.class, group = gmm.class),
      outlier.shape = NA,
      size = 0.5
    ) +
    scale_y_continuous(
      position = 'right',
      labels = function(x){sprintf('%.0f', x)},
      limits = c(0, 15),
      breaks = seq(0, 15, 5)
    ) +
    labs(x = 'cluster', y = 'Grad (˚C/km)') +
    guides(color = 'none', fill = 'none') +
    scale_color_continuous_sequential(
      cls.col,
      breaks = MakeBreaks(bins = length(unique(marx.class.sample$gmm.class))),
      rev = F,
      na.value = 'transparent',
      name = 'classifier clustering groups',
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(3, 'in'),
          order = 1
        )
    ) +
    scale_fill_continuous_sequential(
      cls.col,
      breaks = MakeBreaks(bins = length(unique(marx.class.sample$gmm.class))),
      rev = F,
      na.value = 'transparent',
      name = 'classifier clustering groups',
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(3, 'in'),
          order = 1
        )
    ) +
    theme_dark(base_size = 9) +
    theme(
      axis.text = element_text(color = 'white'),
      axis.title = element_text(color = 'white'),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  p13 <-
    ggplot() +
    geom_hline(yintercept = 120, color = 'grey90') +
    annotate(
      'text',
      x = -Inf,
      y = 120,
      vjust = -0.2,
      hjust = 0,
      size = 2.5,
      color = 'grey90',
      label = 'classifier threshold'
    ) +
    geom_boxplot(
      data = filter(marx.class.sample, pt.path.position == pt.path.filter),
      aes(gmm.class, z, fill = gmm.class, color = gmm.class, group = gmm.class),
      outlier.shape = NA,
      size = 0.5
    ) +
    labs(x = 'cluster', y = 'Z (km)') +
    guides(color = 'none', fill = 'none') +
    scale_color_continuous_sequential(
      cls.col,
      breaks = MakeBreaks(bins = length(unique(marx.class.sample$gmm.class))),
      rev = F,
      na.value = 'transparent',
      name = 'classifier clustering groups',
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(3, 'in'),
          order = 1
        )
    ) +
    scale_fill_continuous_sequential(
      cls.col,
      breaks = MakeBreaks(bins = length(unique(marx.class.sample$gmm.class))),
      rev = F,
      na.value = 'transparent',
      name = 'classifier clustering groups',
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(3, 'in'),
          order = 1
        )
    ) +
    scale_y_continuous(
      position = 'right',
      labels = function(x){sprintf('%.0f', x)},
      limits = c(0, 250),
      breaks = seq(0, 250, 50)
    ) +
    theme_dark(base_size = 9) +
    theme(
      axis.text = element_text(color = 'white'),
      axis.title = element_text(color = 'white'),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'grey90', fill = NA),
      plot.background = element_blank()
    )
  pp5 <-
    p11 +
    ggtitle('GMM clustering') +
    inset_element(p12, 0.11, 0.525, 0.415, 0.755, align_to = 'full') +
    inset_element(p13, 0.11, 0.745, 0.428, 0.945, align_to = 'full') +
    plot_layout(guides = 'collect') &
    theme(
      plot.margin = margin(2, 2, 2, 2),
      plot.title = element_text(margin = margin()),
      legend.box.margin = margin(),
      legend.margin = margin(),
      legend.justification = 'top',
      legend.box.just = 'top',
      legend.direction = 'horizontal',
      legend.position = 'bottom',
      legend.spacing.x = unit(0.06, 'in')
    )
  pp6 <- pp5 + pp4 +
    plot_layout(guides = 'keep') &
    theme(plot.margin = margin(2, 2, 2, 2))
  suppressWarnings(
    ggsave(
      paste0('figs/class_comp', '/', model, '-class-comp.png'),
      plot = pp6,
      device = 'png',
      width = 9,
      height = 5
    )
  )
}

# Set parallel plan
plan(multisession, workers = cores)

# Parallel computing
future_pwalk(
  list(models, paths_marx, paths_grids),
  ~fun(..1, ..2, ..3),
  .progress = T,
  .options = furrr_options(seed = TRUE)
)

# Write log
cat('\nmarx-vis.R complete!\n')
sink()
