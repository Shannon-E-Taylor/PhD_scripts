
library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')
library('patchwork')

source('config.R')
# contrast colour #7A9E9F

# Read in, and proccess, the data
volumes <- read.csv('../data/volumes.csv')
volumes$somite_count <- as.integer(volumes$somite_count)

# Remove data lacking somite count and incorrect species 
volumes <- volumes[volumes$somite_count > 0 ,]#
volumes <- volumes[volumes$species %in% c('A. calliptera', 'R. chillingali'), ]

# somite volume 0 means I didn't measure it 
# (likely because somite is not complete in image)
volumes[volumes$mean_volume_of_somites == 0, ]$mean_volume_of_somites = NA

# and get rid of mis-annotated data
to_remove = c('Merian_xx', # doesn't actually have psm - mislabelled 
              'Rickshaw_1', 'Rickshaw_2') # staging of these is wrong. 
volumes = volumes[volumes$incomplete != TRUE,]
volumes = volumes[!volumes$embryo_id %in% c(to_remove), ]

volumes <- assign_class(
  volumes, 
  AC_midseg_threshold = AC_midseg_threshold, 
  RC_midseg_threshold = RC_midseg_threshold, 
  AC_lateseg_threshold = AC_lateseg_threshold, 
  RC_lateseg_threshold = RC_lateseg_threshold
)

data = read.csv('../data/all_overviews.csv')

# Extract lengths
lengths <- data[data$Text %in% c(
  'embryo_length', 'PSM_length', 'PSM_width', 'somite_length',
  'species', 'somite_stage', 
  'somite_width', 'paraxial-mesoderm', 'head', 'hw'),
  c('image_id', 'Text', 'length', 'species', 'somite_count')]

overview_measurements <- pivot_wider(lengths, 
                                     names_from = 'Text', 
                                     values_from = 'length', values_fn = mean)

overview_measurements$species <- factor(overview_measurements$species, 
                                        levels = c('A. calliptera', 'R. chillingali', 'M. zebra'))

overview_measurements <- overview_measurements[overview_measurements$species %in% c('A. calliptera', 'R. chillingali'), ]

overview_measurements <- overview_measurements[overview_measurements$somite_count > 0, ]



# Read in, and process, the cell counts data 
cell_counts <- read.csv('../data/cell_numbers.csv')
cell_counts$somite_count <- as.integer(cell_counts$somite_count)

# Remove data lacking somite count and incorrect species 
cell_counts <- cell_counts[cell_counts$somite_count > 0 ,]#
cell_counts <- cell_counts[cell_counts$species %in% c('A. calliptera', 'R. chillingali'), ]



mean_chilli_vol = mean(volumes[volumes['somite_count'] < 15 & volumes['species'] == 'R. chillingali', ]$psm_volume, na.rm = TRUE)
mean_calli_vol = mean(volumes[volumes['somite_count'] < 15 & volumes['species'] == 'A. calliptera', ]$psm_volume, na.rm = TRUE)

mean_chilli_EL = mean(overview_measurements[overview_measurements['somite_count'] < 15 & overview_measurements['species'] == 'R. chillingali', ]$embryo_length, na.rm = TRUE)
mean_calli_EL = mean(overview_measurements[overview_measurements['somite_count'] < 15 & overview_measurements['species'] == 'A. calliptera', ]$embryo_length, na.rm = TRUE)

mean_chilli_psm = mean(overview_measurements[overview_measurements['somite_count'] < 15 & overview_measurements['species'] == 'R. chillingali', ]$PSM_length, na.rm = TRUE)
mean_calli_psm = mean(overview_measurements[overview_measurements['somite_count'] < 15 & overview_measurements['species'] == 'A. calliptera', ]$PSM_length, na.rm = TRUE)



mean_chilli_count = mean(cell_counts[cell_counts$somite_count < 15 & cell_counts$species == 'R. chillingali', ]$cells_in_psm, na.rm = TRUE)
mean_calli_count = mean(cell_counts[cell_counts$somite_count < 15 & cell_counts$species == 'A. calliptera', ]$cells_in_psm, na.rm = TRUE)



plt_params <- list(
  scale_color_manual(values = cichlid_colours), 
  scale_fill_manual(values = cichlid_colours), 
  theme_bw(),
  theme(legend.position = 'none'), 
  scale_x_discrete(labels=c("AC","RC"))# , 
  # stat_compare_means(method = "t.test", vjust = 0.4)
)


pos_const = 10

vol <- ggplot(volumes[volumes$somite_count < 15, ], 
              aes(x = species, y = psm_volume/1000000, 
                  color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE, notch = TRUE) + 
  geom_jitter(width = 0.1, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 9/pos_const, size = 3, 
           label = eval(paste(c(round(mean_chilli_vol / mean_calli_vol, 2)), 'fold\nchange'))) + 
  ylab('PSM volume (10\u2076 \u00B5m\u00b3)') +
  scale_y_continuous(limits = c(0, 9)) + 
  ggtitle('D. PSM volume')

vol

cells <- ggplot(cell_counts[cell_counts$somite_count < 15, ], 
                aes(x = species, y = cells_in_psm, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE) + geom_jitter(width = 0.1, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 10000/pos_const, size = 3, 
           label = eval(paste(c(round(mean_chilli_count / mean_calli_count, 2)), 'fold\nchange'))) + 
  ylab('Cells in PSM') + 
  ylim(0, 10000) + 
  ggtitle('G. Cells in PSM')


EL <- ggplot(overview_measurements[overview_measurements$somite_count < 15, ], 
             aes(x = species, y = embryo_length, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE) + geom_jitter(width = 0.3, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 2600/pos_const, size = 3,
           label = eval(paste(c(round(mean_chilli_EL / mean_calli_EL, 2)), 'fold\nchange'))) +
  ylab('Embryo length (\u00B5m)') + 
  ylim(0, 2600) + 
  ggtitle('Embryo length')




psm <- ggplot(overview_measurements[overview_measurements$somite_count < 15, ], 
              aes(x = species, y = PSM_length, color = species, fill= species)) + 
  geom_boxplot(alpha = 0.3, outliers = FALSE) + 
  geom_jitter(width = 0.3, height = 0, size = 1) + 
  plt_params + 
  stat_summary(fun.y=mean, geom="point", fill = 'red', color = 'red', shape = 3) + 
  annotate("text", x = 1.5, y = 900/pos_const, size = 3,
           label = eval(paste(c(round(mean_chilli_psm / mean_calli_psm, 2)), 'fold\nchange'))) +
  ylab('PSM length (\u00B5m)') +
  ylim(0, 900) + 
  ggtitle('A. PSM length') 

library('ggpubr')

psm_len <- ggplot(overview_measurements, 
                  aes(x = somite_count, y = PSM_length, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylim(0, 900) + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ylab('PSM length (\u00B5m)') + 
  xlab('Somite count') + 
  ggtitle('B. PSM length') + 
  stat_poly_line(
    data = overview_measurements[
      overview_measurements$somite_count < 15,], 
    se = FALSE) + 
  stat_poly_line(data = overview_measurements[overview_measurements$somite_count > 15,], se = FALSE)

som_len <- ggplot(overview_measurements, 
                  aes(x = somite_count, y = somite_length, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylab('Nascent somite\nlength (\u00B5m)') + 
  xlab('Somite count') + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ggtitle('C. Nascent somite length') + 
  stat_poly_line(data = overview_measurements[overview_measurements$somite_count < 40,], se = FALSE)

psm + psm_len + som_len + 
  plot_layout(design = 'ABBCC', guides = 'collect') & 
  theme(legend.position = 'none') 


a_guess = 100000
r_guess = 0

get_exponential_coefficents <- function(df, variable){
  
  exponential_model = nls(
    as.formula(paste(variable, "~ a*exp(r*somite_count)")), 
    data = df, 
    start = list(a = a_guess, r = r_guess) 
  )
  # extract parameters. 
  exponential_coef = coef(exponential_model)
  a = exponential_coef[1]
  r = exponential_coef[2]
  
  return(c(a, r))
  
}

out <- get_exponential_coefficents(
  volumes[volumes$species == 'R. chillingali' & 
            volumes$class == 'early', ], 
  variable = 'psm_volume'
)
a_RC_early_psm = out[[1]]
r_RC_early_psm = out[[2]]

out = get_exponential_coefficents(
  volumes[volumes$species == 'R. chillingali' & 
            volumes$class != 'early', ], 
  variable = 'psm_volume'
)
a_RC_midlate_psm = out[[1]]
r_RC_midlate_psm = out[[2]]


out = get_exponential_coefficents(
  volumes[volumes$species == 'A. calliptera' & 
            volumes$class == 'early', ], 
  variable = 'psm_volume'
)
a_AC_early_psm = out[[1]]
r_AC_early_psm = out[[2]]

out = get_exponential_coefficents(
  volumes[volumes$species == 'A. calliptera' & 
            volumes$class != 'early', ], 
  variable = 'psm_volume'
)


a_AC_midlate_psm = out[[1]]
r_AC_midlate_psm = out[[2]]

exp_function <- function(x, a, r) {
  a * exp(r * x) / 1000000
}

psm_vol <- ggplot(volumes, 
                  aes(x = somite_count, y = psm_volume/1000000, 
                      color = species), 
) + 
  geom_point(size = 1) + 
  theme_bw() + 
  scale_color_manual(values = cichlid_colours) +
  geom_function(fun = exp_function,
                args = list(a = a_RC_midlate_psm, r = r_RC_midlate_psm),
                color = cichlid_colours[[2]], 
                linewidth = 1,
                xlim = c(17, 40)) +
  geom_function(fun = exp_function,
                args = list(a = a_AC_midlate_psm, r = r_AC_midlate_psm),
                color = cichlid_colours[[1]], 
                linewidth = 1,
                xlim = c(17, 32)) +
  stat_poly_line(data= volumes[volumes$somite_count < 15, ], se = FALSE) + 
  ylab('PSM volume (10\u2076 \u00B5m\u00b3)') +
  xlab('Somite count') + 
  ggtitle('E. PSM volume') + 
  scale_y_continuous(limits = c(0, 9)) + 
  xlim(0, 40)+ 
  theme(legend.position = 'none')

psm_vol

out = get_exponential_coefficents(
  volumes[#volumes$species == 'R. chillingali' & 
    volumes$class == 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_RC_som_early = out[[1]]
r_RC_som_early = out[[2]]

out = get_exponential_coefficents(
  volumes[#volumes$species == 'A. calliptera' &
    volumes$class == 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_AC_som_early = out[[1]]
r_AC_som_early = out[[2]]

out = get_exponential_coefficents(
  volumes[#volumes$species == 'R. chillingali' & 
    volumes$class != 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_RC_som_midlate = out[[1]]
r_RC_som_midlate = out[[2]]

out = get_exponential_coefficents(
  volumes[#volumes$species == 'A. calliptera' &
    volumes$class != 'early', ], 
  variable = 'mean_volume_of_somites'
)

a_AC_som_midlate = out[[1]]
r_AC_som_midlate = out[[2]]

som_vol <- ggplot(volumes, 
                  aes(x = somite_count, 
                      y = mean_volume_of_somites/1000000, 
                      color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  scale_color_manual(values = cichlid_colours) +
  geom_function(fun = exp_function, 
                args = list(a = a_AC_som_midlate, r = r_AC_som_midlate), 
                color = cichlid_colours[[1]],
                linewidth = 1,
                xlim = c(16, 30)
  ) + 
  geom_function(fun = exp_function, 
                args = list(a = a_RC_som_midlate, r = r_RC_som_midlate), 
                color = cichlid_colours[[2]],
                linewidth = 1,
                xlim = c(16, 40)) + 
  ylab('Nascent somite \nvolume (10\u2076 \u00B5m\u00b3)') + 
  xlab('Somite count') + 
  xlim(0, 40) + 
  ylim(0, .25)+ 
  ggtitle('F. Somite volume') + 
  theme(legend.position = 'none')

som_vol


ymax = 1000





psm_cells <- ggplot(cell_counts, 
                    aes(x = somite_count, y = cells_in_psm, color = species), 
) + 
  geom_point(size = 1) + 
  theme_bw() + 
  scale_color_manual(values = cichlid_colours) +
  ylab('Cells in PSM') + 
  xlab('Somite count') + 
  ggtitle('H. Cells in PSM') + 
  xlim(0, 40) + 
  ylim(0, 10000) + 
  theme(legend.position = 'bottom') + 
  stat_poly_line(data = cell_counts[cell_counts$somite_count > 15,], se = FALSE)


psm_cells

design = '
ABBCC
DEEFF
GHHII
'

psm + psm_len +  som_len + 
  vol + psm_vol + som_vol + 
  cells + psm_cells + 
  plot_layout(design = design, # guides = 'collect', 
              axes = 'collect_y')
  
ggsave('../graphs/PSM_size.png', width = 10, height = 8, dpi =600)

library('effectsize')

ggplot(overview_measurements, 
       aes(x = somite_count, y = somite_length, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylab('Nascent somite\nlength (\u00B5m)') + 
  xlab('Somite count') + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ggtitle('C. Nascent somite length') + 
  stat_poly_line(se = FALSE) + stat_poly_eq(use_label("eq", "R2"))


stat.lm <- lm(
  somite_length~somite_count*species, 
  data = overview_measurements)
anova(stat.lm)
eta_squared(stat.lm)

stat.lm <- lm(
  psm_volume~somite_count+species, 
  data = volumes[volumes$somite_count < 15,])
anova(stat.lm)
eta_squared(stat.lm)

stat.lm <- lm(
  log(mean_volume_of_somites)~somite_count+species, 
  data = volumes[volumes$somite_count < 40,])
anova(stat.lm)
eta_squared(stat.lm)


stat.lm <- lm(
  PSM_length~somite_count*species, 
  data = overview_measurements[overview_measurements$somite_count>15,])
anova(stat.lm)
eta_squared(stat.lm)

stat.lm <- lm(
  cells_in_psm~somite_count*species, 
  data = cell_counts[cell_counts$somite_count>15,])
anova(stat.lm)
eta_squared(stat.lm)

stat.lm <- lm(
  PSM_length~somite_count*species, 
  data = overview_measurements[overview_measurements$somite_count>15,])
anova(stat.lm)
eta_squared(stat.lm)
summary(stat.lm)

ggplot(overview_measurements, 
       aes(x = somite_count, y = PSM_length, color = species)) + 
  geom_point(size = 1) + 
  theme_bw() + 
  ylim(0, 900) + 
  scale_color_manual(values = cichlid_colours) + 
  theme(legend.position = 'none') + 
  ylab('PSM length (\u00B5m)') + 
  xlab('Somite count') + 
  ggtitle('B. PSM length') + 
  stat_poly_eq(
    data = overview_measurements[
      overview_measurements$somite_count > 15,], 
    use_label("eq", "R2"))
