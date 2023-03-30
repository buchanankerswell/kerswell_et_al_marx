#!/usr/bin/env Rscript

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Load functions and libraries
cat(rep('~', 80), sep='')
cat('\nLoading packages and functions ...')
source('R/functions.R')
load('data/preprocessed.RData')

# Test arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  cat('\nNo arguments passed to R/summary.R')
  pt.path.filter <- 'maxP'
  cat('\nPT path filter:', pt.path.filter)
} else if (length(args) != 0) {
  pt.path.filter <- suppressWarnings(args[1])
  if(!(pt.path.filter %in% c('maxT', 'maxP', 'maxPT'))) {
    pt.path.filter <- 'maxP'
    cat('\nUnrecognized PT path filter. Please use "maxT", "maxP", or "maxPT"')
    cat('\nDefaulting to', pt.path.filter)
  }
}

# Get marker and grid data paths
cat('\nLoading data ...')
paths_marx <- list.files('data/marx_classified', pattern = 'classified.RData', full.names = T)
paths_grids <- list.files('data/marx_traced', pattern = 'grids.RData', full.names = T)
models <- str_extract(paths_marx, 'cd.[0-9]+')
numerical.model.parameters.kerswell.etal.2021 <-
  read_csv('data/numerical-model-parameters.csv', show_col_types = F)

# Test input files
if (length(paths_marx) != length(paths_grids)) {
  stop('Number of classified marx and grid files are different!')
} else if (length(paths_marx) < 1) {
  stop(
    paste0('No classified marx files found in data/marx_classified !!')
  )
} else if (length(paths_grids) < 1) {
  stop('No grids files found in data/marx_traced/!')
}

# Summarise classification results
mods.summary <-
  select(
    numerical.model.parameters.kerswell.etal.2021,
    model,
    thermal.parameter,
    coupling.depth,
    upper.plate.thickness,
    oceanic.plate.age,
    convergence.velocity
  )

cat('\nCompiling markers data from all models ...')
# Save as list
mlist <-
  map(paths_marx, ~{
    load(.x)
    get(paste0(str_extract(.x, 'cd[a-z][0-9]+'), '.marx.classified'))
  }) %>%
  set_names(models)

# Summarize markers results
marx.class.summary.filtered.pt.path <-
  map_df(mlist, ~.x$marx.class, .id = 'model') %>%
  filter(pt.path.position == pt.path.filter)

# Summarise jackknife samples
if (!is.null(reduce(map(mlist, ~.x$jk), ~.x))) {
  stats.summary <-
    map_df(mlist, ~.x$jk, .id = 'model') %>%
    left_join(mods.summary, by = 'model')
  write_csv(stats.summary, 'data/classification-stats-summary.csv')
} else {
  stats.summary <- NULL
}

# Clean up environment
rm(mlist)

# Write marker and stats summaries to csv
if(
  length(list.files('data/marx_classified', pattern = '.csv')) < 64
  ) {
  cat('\nWriting markers results to data/marx_classified')
  walk(
    models,
    ~write_csv(
      filter(marx.class.summary.filtered.pt.path, model == .x),
      paste0('data/marx_classified', '/', .x, '-marx-classified.csv')
    )
  )
}

# Save markers data used for plotting and analysis
save(marx.class.summary.filtered.pt.path, file = 'data/marx-class-summary.RData')

# Summarise density peaks (modes) and other threshold values
pd15.dens.summary <- get_dens_summary('pd15')
ag18.dens.summary <- get_dens_summary('ag18')

# Compute 2d density for recovered markers
dens.summary <-
  compute_dens_2d(filter(marx.class.summary.filtered.pt.path, recovered), n = 80)

# Summarise marker cdfP
cdfP.summary <-
  marx.class.summary.filtered.pt.path %>%
  filter(recovered) %>%
  select(P) %>%
  arrange(P) %>%
  mutate(prob = (row_number()-1)/n())

# Contour limits
break.points.marx <- MakeBreaks(bins = 6)
break.points.ag18 <- MakeBreaks(bins = 5)
break.points.pd15 <- MakeBreaks(binwidth = 5)
break.lims.marx <- c(round(max(dens.summary$cnt)/8, -1), NA)
break.lims.ag18 <- c(round(max(ag18.dens$cnt)/6), NA)
break.lims.pd15 <- c(round(max(pd15.dens$cnt)/12), NA)
if(break.lims.marx[1] < 1) {break.lims.marx[1] <- 1}

# Colors
ag18.col <- 'grays'
pd15.col <- 'burg'
marx.col <- 'emrld'
rec.col <- 'gold'
sub.col <- 'gold4'
guide.pal <- 'grays'

# Summarise marker cdfs
cat('\nDrawing plots ...')

# Plot marker summary results
# Compute density modes
n <- 10
thresh <- 65
d.marx.P <- density(
  marx.class.summary.filtered.pt.path$P[marx.class.summary.filtered.pt.path$recovered]
)
d.pd15.P <- density(pd15$P)
d.ag18.P <- density(ag18$P)
d.marx.T <- density(
  marx.class.summary.filtered.pt.path$T[marx.class.summary.filtered.pt.path$recovered]
)
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
    aes(T, P, z = cnt, color = after_stat(level)),
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
    data = filter(marx.class.summary.filtered.pt.path, recovered),
    aes(T, P),
    size = 1,
    alpha = 0.1
  ) +
  geom_contour_fill(
    data = dens.summary,
    aes(T, P, z = cnt)
  ) +
  geom_contour_tanaka(
    data = dens.summary,
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
      paste0(
        'marker count n=',
        nrow(filter(marx.class.summary.filtered.pt.path, recovered))
      ),
    guide =
      guide_colorstrip(
        inside = T,
        title.position = 'top',
        title.vjust = 1,
        barwidth = unit(2.8, 'in'),
        order = 1
      )
  ) +
  geom_abline(slope = 1/700, intercept = 0, size = 0.25) +
  geom_abline(slope = 1/175, intercept = 0, size = 0.25) +
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
    vjust = -0.5,
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
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 4)) +
  scale_x_continuous(breaks = seq(0, 1000, 200)) +
  labs(x = 'temperature (˚C)', y = 'pressure (GPa)') +
  theme_dark(base_size = 16) +
  theme(
    legend.title = element_text(margin = margin(0, 0, -5, 0)),
    legend.spacing.x = unit(0, 'in'),
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
    data = drop_na(cdfP.summary),
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
  theme_dark(base_size = 10) +
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
  theme_dark(base_size = 10) +
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
  ggtitle('global markers') +
  inset_element(p2, 0.11, 0.495, 0.39, 0.755, align_to = 'full') +
  inset_element(p3, 0.11, 0.745, 0.353, 0.945, align_to = 'full') +
  inset_element(p3.b, 0.343, 0.685, 0.63, 0.945, align_to = 'full') +
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
cat('\nSaving plot to figs/summary/marx-summary-comp.png')
suppressWarnings(
  ggsave(
    plot = pp1,
    file = 'figs/summary/marx-summary-comp.png',
    device = 'png',
    width = 4.5,
    height = 5,
    units = 'in'
  )
)

# Plot rock record
# Compute density modes
p4 <-
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
    aes(T, P, z = cnt, color = after_stat(level)),
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
  geom_point(data = ag18, aes(T, P, shape = 'ag18'), size = 2) +
  geom_point(data = pd15, aes(T, P, shape = 'pd15'), size = 1.5) +
  scale_shape_manual(
    name = 'samples',
    values = c('ag18' = 0, 'pd15' = 16),
    breaks = c('pd15', 'ag18'),
    guide =
      guide_legend(
        override.aes = list(size = 5),
        title.position = 'top',
        label.position = 'right',
        title.vjust = 0,
        order = 1
      )
  ) +
  geom_abline(slope = 1/700, intercept = 0, size = 0.25) +
  geom_abline(slope = 1/175, intercept = 0, size = 0.25) +
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
    vjust = -0.5,
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
  coord_cartesian(xlim = c(0, 1000), ylim = c(0, 4)) +
  scale_x_continuous(breaks = seq(0, 1000, 200)) +
  labs(x = 'temperature (˚C)', y = 'pressure (GPa)') +
  theme_dark(base_size = 16) +
  theme(
    legend.title = element_text(margin = margin(-3.8, 0, -5, 0)),
    legend.spacing.x = unit(0, 'in'),
    panel.grid = element_blank()
  )
p5 <-
  ggplot() +
  geom_line(
    data = arrange(pd15, P) %>% mutate(cdfP = (row_number()-1)/n()),
    aes(P, cdfP, color = 'pd15')
  ) +
  geom_line(
    data = arrange(ag18, P) %>% mutate(cdfP = (row_number()-1)/n()),
    aes(P, cdfP, color = 'ag18')
  ) +
  labs(y = 'CDF', x = 'pressure (GPa)') +
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
  theme_dark(base_size = 10) +
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
p6 <-
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
  labs(x = 'pressure (GPa)', y = NULL) +
  coord_cartesian(xlim = c(0, 4)) +
  scale_x_continuous(breaks = seq(0, 4, 1)) +
  scale_y_continuous(position = 'right') +
  theme_dark(base_size = 10) +
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
p6.b <-
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
pp2 <-
  p4 +
  ggtitle('rock record') +
  inset_element(p5, 0.11, 0.495, 0.39, 0.755, align_to = 'full') +
  inset_element(p6, 0.11, 0.745, 0.353, 0.945, align_to = 'full') +
  inset_element(p6.b, 0.343, 0.685, 0.63, 0.945, align_to = 'full') +
  plot_layout(guides = 'collect') &
  theme(
    plot.margin = margin(2, 2, 2, 2),
    plot.title = element_text(margin = margin()),
    legend.box.margin = margin(0, 0, 12.1, 0),
    legend.margin = margin(),
    legend.justification = 'top',
    legend.box.just = 'top',
    legend.direction = 'horizontal',
    legend.position = 'bottom',
    legend.spacing.x = unit(0.06, 'in')
  )
cat('\nSaving plot to figs/summary/rox-comp.png')
# Save
suppressWarnings(
  ggsave(
    'figs/summary/rox-comp.png',
    plot = pp2,
    device = 'png',
    width = 4.5,
    height = 5
 )
)

# Rearrange dataframe for correlations and transform units into
# celcius, celcius per kilometer, percent, gigapascal
if (!is.null(stats.summary)) {
  d.corr <-
    stats.summary %>%
    select(
      upper.plate.thickness,
      oceanic.plate.age,
      convergence.velocity,
      rec.ratio.est,
      mode1.P.est,
      mode2.P.est,
      mode1.T.est,
      mode2.T.est,
      mode1.grad.est,
      mode2.grad.est
    ) %>%
    mutate(
      # Transform into typical units (celcius, celcius per kilometer, gigapascal, percent)
      rec.ratio.est = rec.ratio.est*100,
      mode1.P.est = mode1.P.est,
      mode2.P.est = mode2.P.est,
      mode1.T.est = mode1.T.est,
      mode2.T.est = mode2.T.est,
      mode1.grad.est = mode1.grad.est,
      mode2.grad.est = mode2.grad.est
    ) %>%
    rename(
      'OP age' = oceanic.plate.age,
      'Velocity' = convergence.velocity,
      'UP thickness' = upper.plate.thickness,
      'P mode1' = mode1.P.est,
      'P mode2' = mode2.P.est,
      'T mode1' = mode1.T.est,
      'T mode2' = mode2.T.est,
      'Grad mode1' = mode1.grad.est,
      'Grad mode2' = mode2.grad.est,
      'Rec rate' = rec.ratio.est
    )
  d.modes <-
    stats.summary %>%
    select(
      model,
      upper.plate.thickness,
      oceanic.plate.age,
      convergence.velocity,
      rec.ratio.est,
      mode1.P.est,
      mode2.P.est,
      mode1.T.est,
      mode2.T.est,
      mode1.grad.est,
      mode2.grad.est
    ) %>%
    mutate(
      # Transform into typical units (celcius, celcius per kilometer, gigapascal, percent)
      rec.ratio.est = rec.ratio.est*100,
      mode1.P.est = mode1.P.est,
      mode2.P.est = mode2.P.est,
      mode1.T.est = mode1.T.est,
      mode2.T.est = mode2.T.est,
      mode1.grad.est = mode1.grad.est,
      mode2.grad.est = mode2.grad.est
    ) %>%
    rename(
      'OP age' = oceanic.plate.age,
      'Velocity' = convergence.velocity,
      'UP thickness' = upper.plate.thickness,
      'P mode1' = mode1.P.est,
      'P mode2' = mode2.P.est,
      'T mode1' = mode1.T.est,
      'T mode2' = mode2.T.est,
      'Grad mode1' = mode1.grad.est,
      'Grad mode2' = mode2.grad.est,
      'Rec rate' = rec.ratio.est
    ) %>%
    pivot_longer(
      -c(model, `UP thickness`, `OP age`, Velocity),
      names_to = 'var',
      values_to = 'varval'
    ) %>%
    pivot_longer(-c(var, varval, model), names_to = 'setting', values_to = 'settingval') %>%
    mutate(
      setting = factor(setting, levels = c('OP age', 'Velocity', 'UP thickness')),
      var =
        factor(var, levels =
          c(
            'P mode1',
            'P mode2',
            'T mode1',
            'T mode2',
            'Grad mode1',
            'Grad mode2',
            'Rec rate'
          )
        )
    )
  # Calculate correlations with boundary conditions
  # use spearmans to test nonlinear correlations (specifically: monotonically-increasing)
  # https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
  corr <- d.corr %>% as.matrix() %>% Hmisc::rcorr(type = 'spearman')
  corr$P[is.na(corr$P)] <- 1
  # Flatten correlation matrix into dataframe
  corr.df <-
    flattenCorrMatrix(corr$r, corr$P) %>%
    filter(row %in% c('UP thickness', 'OP age', 'Velocity')) %>%
    filter(!(column %in% c('UP thickness', 'OP age', 'Velocity'))) %>%
    mutate(
      row = factor(row, levels = c('UP thickness', 'Velocity', 'OP age')),
      column =
        factor(column, levels =
          c(
            'P mode1',
            'P mode2',
            'T mode1',
            'T mode2',
            'Grad mode1',
            'Grad mode2',
            'Rec rate'
          )
        ),
      lvl =
        map_chr(.$p, ~{
          if (.x <= 0.001) {'***'}
          else if (.x <= 0.01) {'**'}
          else if (.x <= 0.05) {'*'}
          else {'-'}
      })
    )
  p7 <-
    ggplot() +
    ggtitle('monotonic correlation test (always increasing or decreasing)') +
    geom_tile(
      data = mutate(corr.df, cor = ifelse(abs(p) < 0.05, cor, NA)),
      aes(column, row, fill = cor),
      color = 'grey20',
      size = 0.2
    ) +
    geom_text(
      data = corr.df,
      aes(column, row, label = paste0(lvl, '\n', ifelse(abs(p) < 0.05, round(cor, 2), ''))),
      size = 5
    ) +
    coord_cartesian(expand = F) +
    scale_fill_continuous_diverging(
      'blue-red3',
      rev = T,
      limits = c(-0.7, 0.7),
      breaks = seq(-0.7, 0.7, 0.35),
      na.value = 'grey50',
      name = bquote("Spearman's correlation coefficient"~rho),
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(3.5, 'in')
        )
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(
      plot.margin = margin(2, 2, 2, 2),
      panel.background = element_rect(color = 'black'),
      legend.box.margin = margin(),
      legend.margin = margin(),
      legend.justification = 'top',
      legend.box.just = 'top',
      legend.position = 'bottom',
      axis.text.y = element_text(angle = 90, hjust = 0.5, vjust = -8),
      axis.ticks = element_blank()
    )
  p7.b <-
    ggplot() +
    ggtitle('distributions of PT modes and recovery rates') +
    geom_boxplot(
      data = d.modes,
      aes(varval, settingval, group = settingval),
      width = 10,
      outlier.size = 0.5
    ) +
    facet_grid(vars(setting), vars(var), scales = 'free', switch = 'y') +
    scale_y_continuous(position = 'right') +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 16) +
    theme(
      plot.margin = margin(2, 2, 2, 2),
      panel.background = element_rect(color = 'black', fill = 'grey50'),
      panel.grid = element_blank(),
      legend.box.margin = margin(),
      legend.margin = margin(),
      legend.justification = 'top',
      legend.box.just = 'top',
      legend.position = 'bottom',
      axis.text.x = element_text(angle = -60, hjust = 0.5),
      axis.ticks = element_blank()
    )
  pp.corr <-
    p7.b / p7
  # Save
  suppressWarnings(
    ggsave(
      'figs/summary/thermo-kinematic-correlations-summary.png',
      plot = pp.corr,
      device = 'png',
      width = 9,
      height = 10
    )
  )
}

# Compute density for subsets of numerical models
plots.dens <-
  map(
    c(
      'young-slow',
      'young-fast',
      'young-thick',
      'young-thin',
      'old-slow',
      'old-fast',
      'old-thick',
      'old-thin'
    ),
    ~{
    m <-
      marx.class.summary.filtered.pt.path %>%
      left_join(
        select(
          numerical.model.parameters.kerswell.etal.2021,
          model,
          oceanic.plate.age,
          convergence.velocity,
          upper.plate.thickness
        ),
        by = 'model'
      )
      if (.x == 'young-slow') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(32.6, 55.0), convergence.velocity %in% c(40, 66))
      } else if (.x == 'young-fast') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(32.6, 55.0), convergence.velocity %in% c(80, 100))
      } else if (.x == 'young-thick') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(32.6, 55.0), upper.plate.thickness %in% c(78, 94))
      } else if (.x == 'young-thin') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(32.6, 55.0), upper.plate.thickness %in% c(46, 62))
      } else if (.x == 'old-slow') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(85.0, 110.0), convergence.velocity %in% c(40, 66))
      } else if (.x == 'old-fast') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(85.0, 110.0), convergence.velocity %in% c(80, 100))
      } else if (.x == 'old-thick') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(85.0, 110.0), upper.plate.thickness %in% c(78, 94))
      } else if (.x == 'old-thin') {
        m <- 
          m %>%
          filter(oceanic.plate.age %in% c(85.0, 110.0), upper.plate.thickness %in% c(46, 62))
      }
    d.pd15.P <- density(pd15$P)
    d.ag18.P <- density(ag18$P)
    d.marx.P <- density(m$P[m$recovered])
    d.pd15.T <- density(pd15$T)
    d.ag18.T <- density(ag18$T)
    d.marx.T <- density(m$T[m$recovered])
    dens.summary <- compute_dens_2d(filter(m, recovered), n = 80)
    cdfP.summary <-
      m %>%
      filter(recovered) %>%
      select(P) %>%
      arrange(P) %>%
      mutate(prob = (row_number()-1)/n())
    break.lims.marx <- c(round(max(dens.summary$cnt)/10, -1), NA)
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
        aes(T, P, z = cnt, color = after_stat(level)),
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
        data = filter(m, recovered),
        aes(T, P),
        size = 1,
        alpha = 0.1
      ) +
      geom_contour_fill(
        data = dens.summary,
        aes(T, P, z = cnt)
      ) +
      geom_contour_tanaka(
        data = dens.summary,
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
      name = paste0('marker count n=', nrow(filter(m, recovered))),
      guide =
        guide_colorstrip(
          inside = T,
          title.position = 'top',
          title.vjust = 1,
          barwidth = unit(2.8, 'in'),
          order = 1
        )
    ) +
    geom_abline(slope = 1/700, intercept = 0, size = 0.25) +
    geom_abline(slope = 1/175, intercept = 0, size = 0.25) +
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
      vjust = -0.5,
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
    coord_cartesian(xlim = c(0, 1000), ylim = c(0, 4)) +
    scale_x_continuous(breaks = seq(0, 1000, 200)) +
    labs(x = 'temperature (˚C)', y = 'pressure (GPa)') +
    theme_dark(base_size = 16) +
    theme(
      legend.title = element_text(margin = margin(0, 0, -5, 0)),
      legend.spacing.x = unit(0, 'in'),
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
        data = drop_na(cdfP.summary),
        aes(P, prob, color = 'markers')
      ) +
      labs(y = 'CDF', x = 'pressure (GPa)') +
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
      theme_dark(base_size = 10) +
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
      labs(x = 'pressure (GPa)', y = NULL) +
      coord_cartesian(xlim = c(0, 4)) +
      scale_x_continuous(breaks = seq(0, 4, 1)) +
      scale_y_continuous(position = 'right') +
      theme_dark(base_size = 10) +
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
      ggtitle(paste0(str_replace(.x, '-', ' '), ' systems')) +
      inset_element(p2, 0.11, 0.495, 0.39, 0.755, align_to = 'full') +
      inset_element(p3, 0.11, 0.745, 0.353, 0.945, align_to = 'full') +
      inset_element(p3.b, 0.343, 0.685, 0.63, 0.945, align_to = 'full') +
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
    cat('\nSaving plot to figs/summary/marx-summary-comp-', .x, '.png', sep = '')
    suppressWarnings(
      ggsave(
        plot = pp1,
        file = paste0('figs/summary/marx-summary-comp-', .x, '.png'),
        device = 'png',
        width = 4.5,
        height = 5,
        units = 'in'
      )
    )
    pp1 <-
      p1 +
      ggtitle(paste0(str_replace(.x, '-', ' '), ' systems')) +
      inset_element(p2, 0.11, 0.495, 0.39, 0.755, align_to = 'full') +
      inset_element(p3, 0.11, 0.745, 0.353, 0.945, align_to = 'full') +
      inset_element(p3.b, 0.343, 0.685, 0.63, 0.945, align_to = 'full') +
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
    pp1
    }
  )

# Composition plot
p <-
(plots.dens[[1]] + plots.dens[[2]] + plot_layout(guides = 'keep')) /
(plots.dens[[3]] + plots.dens[[4]] + plot_layout(guides = 'keep')) +
plot_layout(guides = 'keep')
cat('\nSaving plot to figs/summary/rox-marx-model-comp-young.png')
suppressWarnings(
  ggsave(
    'figs/summary/rox-marx-model-comp-young.png',
    plot = p,
    device = 'png',
    width = 9,
    height = 10
  )
)

# Composition plot
p <-
(plots.dens[[5]] + plots.dens[[6]] + plot_layout(guides = 'keep')) /
(plots.dens[[7]] + plots.dens[[8]] + plot_layout(guides = 'keep')) +
plot_layout(guides = 'keep')
cat('\nSaving plot to figs/summary/rox-marx-model-comp-old.png')
suppressWarnings(
  ggsave(
    'figs/summary/rox-marx-model-comp-old.png',
    plot = p,
    device = 'png',
    width = 9,
    height = 10
  )
)

# Write log
cat('\nsummary.R complete!\n')
sink()
