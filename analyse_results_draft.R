# Setup -------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F, scipen = 19)

# Introduction ------------------------------------------------------------

# This is a draft script for generating plots from a NUMTdumper find 
# run. It currently works only for non-multiplicative terms, but if you 
# have any multiplicative terms don't worry, they'll just be ignored.
# As long as you have the below packages installed, it should just run.
# The plots are exactly the same format as figs 2 and 3 in the paper, see
# https://www.biorxiv.org/content/10.1101/2020.06.17.157347v1.full

# Load libraries ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# Load variables ----------------------------------------------------------

termdescriptions = list(
  'library_n' = 'Minimum absolute ASV abundance by library',
  'library_p' = 'Minimum relative ASV abundance by library',
  'library|clade_n' = 'Minimum absolute ASV abundance by library and clade',
  'library|clade_p' = 'Minimum relative ASV abundance by library and clade',
  'library|clade_n' = 'Minimum absolute ASV abundance by library and taxon',
  'library|taxon_p' = 'Minimum relative ASV abundance by library and taxon',
  'total_n' = 'Minimum absolute ASV abundance',
  'total_p' = 'Minimum relative ASV abundance',
  'total|clade_n' = 'Minimum absolute ASV abundance by clade',
  'total|clade_p' = 'Minimum relative ASV abundance by clade',
  'total|taxon_n' = 'Minimum absolute ASV abundance by taxon',
  'total|taxon_p' = 'Minimum relative ASV abundance by taxon'
)

# Load data ---------------------------------------------------------------

### REPLACE WITH THE PATH TO YOUR RESULTS FILE HERE ###
results <- read.csv("_results.csv",
                    check.names = F)

# Organise data -----------------------------------------------------------

# Reshape the data matrix and separate out the variable details
idvars <- colnames(results)[c(1:which(colnames(results) == 'asvs_total'), ncol(results))]
results <- pivot_longer(results, -all_of(idvars), 
                        names_to = c('data', 'part', 'stat'), 
                        names_sep = '_', 
                        values_to = 'value')

# Extract non-multiplicative terms and generate the two plot sets ---------
  # Extract
oneterms <- filter(results, !grepl('\\*', term)) %>%
  pivot_longer(contains('_threshold'), 
               names_to = c('category', 'metric', 't'), 
               names_sep = '_', 
               values_to = 'threshold') %>%
  filter(paste(category, metric, sep = '_') == term) %>%
  select(-c(category, metric, t)) %>%
  mutate(termdescription = unlist(termdescriptions[term]))

# To extract the passing ASVs at a certain threshold value, use the 
# number shown on the plot as the resultindex for a NUMTdumper dump run

  # Plots 1
unite(oneterms, 'part', data, part) %>%
  filter(part %in% c('targets_rejected', 'nontargets_retained'), stat == 'p') %>%
  ggplot(aes(x = threshold, y = value, col = part)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = resultindex), nudge_y = .1, size = 2.5) +  
  scale_colour_discrete(name = '', 
                        labels = c('Verified Non-Authentic ASVs\nRetained (false positives)', 
                                   'Verified Authentic ASVs\nRejected (false negatives)')) + 
  labs(y = '% Retained/Rejected') + 
  facet_wrap( ~ termdescription, ncol = 1, scales = 'free_x', strip.position = 'bottom') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), legend.position = 'bottom',
        strip.placement = 'outside', strip.background = element_blank())

  # Plots 2
filter(oneterms, stat == 'estimate') %>%
  pivot_wider(names_from = part, values_from = value) %>%
  ggplot(aes(x = threshold, y = retained)) +
  geom_ribbon(aes(ymin = 0, ymax = total, fill = data), alpha = .5) +
  geom_line(aes(col = data)) +
  geom_point(aes(col = data)) + 
  geom_text(aes(y = retained, label = resultindex), nudge_y = 10, size = 2.5) + 
  scale_colour_discrete(name = '', 
                        labels = c('Estimated Retained\nNon-Authentic ASVs', 
                                   'Estimated Retained\nAuthentic ASVs')) +
  scale_fill_discrete(name = '', 
                      labels = c('Estimated Initial\nNon-Authentic ASVs', 
                                 'Estimated Initial\nAuthentic ASVs')) + 
  labs(y = 'Number of ASVs') + 
  facet_wrap( ~ termdescription, ncol = 1, scales = 'free_x', strip.position = 'bottom') +
  theme_bw() + 
  theme(axis.title.x = element_blank(), legend.position = 'bottom',
        strip.placement = 'outside', strip.background = element_blank())


