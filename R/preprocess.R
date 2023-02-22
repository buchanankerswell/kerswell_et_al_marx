#!/usr/bin/env Rscript

# Create directories
dir.create('figs', showWarnings = F)
dir.create('draft/assets/figs', showWarnings = F)
dir.create('figs/marx_comp', showWarnings = F)
dir.create('figs/class_comp', showWarnings = F)
dir.create('figs/summary', showWarnings = F)
dir.create('data/marx_classified', showWarnings = F)
dir.create('data/marx_traced', showWarnings = F)

# Capture output
sink(file = paste0('data/log-', Sys.Date()), append = T, type = 'output', split = T)

# Load functions and libraries
cat(rep('~', 80), sep='')
cat('\nLoading packages and functions ...')
source('R/functions.R')

# Read Penniston-Dorland et al., 2015 dataset
cat('\nProcessing data/pd15.csv')
pd15 <- read_csv('data/pd15.csv', show_col_types = F)
# compute pd15 density in 2d
k.dens.pd15 <- MASS::kde2d(pd15$T, pd15$P, n = 10, lims = c(c(0, 1000), c(0, 4)))
pd15.dens <-
  expand.grid(k.dens.pd15$x, k.dens.pd15$y) %>%
  as_tibble() %>%
  rename(T = Var1, P = Var2) %>%
  mutate(
    k.dens = as.vector(k.dens.pd15$z),
    cnt = nrow(pd15) / sum(k.dens) * k.dens
  )

# Read Agard et al., 2018 dataset
cat('\nProcessing data/ag18.csv')
ag18 <- read_csv('data/ag18.csv', show_col_types = F)
# compute ag18 density in 2d
k.dens.ag18 <- MASS::kde2d(ag18$T, ag18$P, n = 10, lims = c(c(0, 1000), c(0, 4)))
ag18.dens <-
  expand.grid(k.dens.ag18$x, k.dens.ag18$y) %>%
  as_tibble() %>%
  rename(T = Var1, P = Var2) %>%
  mutate(
    k.dens = as.vector(k.dens.ag18$z),
    cnt = nrow(ag18) / sum(k.dens) * k.dens
  )

# Define metamorphic reactions and facies
# Antigorite reaction Schmidt & Poli (1998)
cat('\nDefining metamorphic reactions')
atg <-
  tibble(P = seq(-0.1, 4.1, length.out = 50)) %>%
  mutate(
    z = P*35,
    T =
      ifelse(
        z <= 63,
        751.5 + (6.008*z) - (3.469e-2*z^2),
        1013.2 - ((6.039e-5)*z) - ((4.289e-3)*(z^2))
      )
  ) %>%
  mutate(T = T-273)
# Eclogitization reaction Ito & Kennedy (1971)
pl <-
  tibble(T = seq(0, 1373, length.out = 50)) %>%
  mutate(
    P = ((T * 20.0)-1460.0)/1e4,
    z = P * 35
  ) %>%
  mutate(P = ifelse(P <= 1.2, 1.2, P)) %>%
  filter(T >= 473) %>%
  mutate(T = T-273)

# Clean up environment
cat('\nCleaning up environment')
rm(list = lsf.str())
rm(p.list)

# Save
cat('\nSaving data to: data/preprocessed.RData')
save.image('data/preprocessed.RData')

# Write log
cat('\npreprocess.R complete!\n')
sink()