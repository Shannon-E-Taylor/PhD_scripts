

library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')
library('patchwork')
library('ggpubr')
library('stringr')
source('config.R')

data = read.csv('../data/all_overviews.csv')
# data <- as.data.frame(apply(data,2, str_remove_all, " "))

data$Text <- str_remove_all(data$Text, " ")

data <- data %>% 
  mutate(Text = recode(Text, psm_length = 'PSM_length')) %>% 
  mutate(Text = recode(Text, PSM_lnegth = 'PSM_length')) %>% 
  mutate(Text = recode(Text, PSM_lenght = 'PSM_length')) %>% 
  mutate(Text = recode(Text, psm = 'PSM_length')) %>% 
  mutate(Text = recode(Text, Embryo_length = 'embryo_length'))
  
# Extract lengths
lengths <- data[data$Text %in% c(
  'embryo_length', 'PSM_length', 'PSM_width', 'somite_length',
  'species', 'somite_stage', 
  'somite_width', 'paraxial-mesoderm', 'head', 'hw'),
  c('image_id', 'Text', 'length', 'species', 'somite_count')]

lengths$length <- as.numeric(lengths$length)
lengths$somite_count <- as.numeric(lengths$somite_count)


plt_params <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'none'), 
  scale_x_discrete(labels=c("AC","RC"))
)


overview_measurements <- as.data.frame(pivot_wider(lengths, 
                                     names_from = 'Text', 
                                     values_from = 'length', 
                                     values_fn = mean))
                                     

# overview_measurements$species <- as.character(overview_measurements$species)
# na.omit(overview_measurements, cols = c('species'))
overview_measurements$species <- factor(overview_measurements$species,
                                        levels = c('A. calliptera',
                                                   'R. chillingali',
                                                   'M. zebra', ''))

overview_measurements <- overview_measurements[
  overview_measurements$species %in% c("A. calliptera", "R. chillingali"),
  ]

overview_measurements <- overview_measurements[overview_measurements$somite_count > 0, ]

somite_thresh <- 7
pos_const = 10

mean_chilli_EL = mean(overview_measurements[overview_measurements['somite_count'] < somite_thresh & overview_measurements['species'] == 'R. chillingali', ]$embryo_length, na.rm = TRUE)
mean_calli_EL = mean(overview_measurements[overview_measurements['somite_count'] < somite_thresh & overview_measurements['species'] == 'A. calliptera', ]$embryo_length, na.rm = TRUE)

mean_chilli_psm = mean(overview_measurements[overview_measurements['somite_count'] < somite_thresh & overview_measurements['species'] == 'R. chillingali', ]$PSM_length, na.rm = TRUE)
mean_calli_psm = mean(overview_measurements[overview_measurements['somite_count'] < somite_thresh & overview_measurements['species'] == 'A. calliptera', ]$PSM_length, na.rm = TRUE)


EL <- ggplot(
  overview_measurements[overview_measurements$somite_count < somite_thresh, ], 
             aes(x = species, y = embryo_length, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE) + geom_jitter(width = 0.3, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 2600/pos_const, size = 3,
          label = eval(paste(c(round(mean_chilli_EL / mean_calli_EL, 2)), 'fold\nchange'))) +
  ylab('Embryo length (\u00B5m)') + 
  ylim(0, 2600) + 
  ggtitle('A. Embryo length')

prop.psm <- ggplot(
  overview_measurements[overview_measurements$somite_count < somite_thresh, ], 
  aes(x = species, y = PSM_length / embryo_length * 100, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE) + geom_jitter(width = 0.3, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 50/pos_const, size = 3,
           label = eval(paste(c(round(mean_chilli_EL / mean_calli_EL, 2)), 'fold\nchange'))) +
  ylab('PSM length (% embryo length)') + 
  ylim(0, 50) + 
  ggtitle('C. PSM length (% embryo)')

psm <- ggplot(overview_measurements[overview_measurements$somite_count < somite_thresh, ], 
              aes(x = species, y = PSM_length, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE) + 
  geom_jitter(width = 0.3, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 900/pos_const, size = 3,
           label = eval(paste(c(round(mean_chilli_psm / mean_calli_psm, 2)), 'fold\nchange'))) +
  ylab('PSM length (\u00B5m)') +
  ylim(0, 900) + 
  ggtitle('B. PSM length') 

EL + psm + prop.psm & stat_compare_means(method = "t.test")

EL <- ggplot(
  overview_measurements[overview_measurements$somite_count < somite_thresh, ], 
  aes(x = somite_count, y = embryo_length, color = species, fill= species)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = cichlid_colours) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + ylab('Embryo length (\u00B5m)') + 
  ylim(0, 2600) + 
  ggtitle('A. Embryo length')

prop.psm <- ggplot(
  overview_measurements[overview_measurements$somite_count < somite_thresh, ], 
  aes(x = somite_count, y = PSM_length / embryo_length * 100, color = species, fill= species)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = cichlid_colours) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
   ylab('PSM length (% embryo length)') + 
  ylim(0, 50) + 
  ggtitle('C. PSM length (% embryo)')

psm <- ggplot(overview_measurements[overview_measurements$somite_count < somite_thresh, ], 
              aes(x = somite_count, y = PSM_length, color = species, fill= species)) + 
  geom_point(size = 1) + 
  ylab('PSM length (\u00B5m)') +
  ylim(0, 900) + 
  scale_color_manual(values = cichlid_colours) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  ggtitle('B. PSM length') 

EL + psm + prop.psm
