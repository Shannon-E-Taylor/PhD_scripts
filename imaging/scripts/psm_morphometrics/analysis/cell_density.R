
library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')
library('stringr')
library('patchwork')
library('broom')
library('latex2exp')
library('ggpmisc')

library('ggpubr')


# Import variables from config file 
source("config.R") 



density_data = read.csv('../data/Density_across_psm.csv')

density_data <- density_data[density_data$species %in% c('R. chillingali', 
                                                         'A. calliptera'), ]
density_data$density = density_data$total_cell_number / density_data$psm_volume


density_data$somite_count <- as.integer(density_data$somite_count)
density_data = density_data[density_data$somite_count > 1, ]
density_data$stage <- density_data$somite_count


density_data$max_stage = 32
density_data$max_stage[density_data$species == 'R. chillingali'] = 38

density_data$propsom = density_data$somite_count / density_data$max_stage

density_data <- density_data[!is.na(density_data$somite_count), ]

tokeep = read.csv('../data/cell_counts_to_keep.csv')


to_exclude = tokeep[!tokeep$KEEP,]$psm_id

density_data <- density_data[!density_data$psm_id%in% to_exclude,]

density_data <- density_data[!is.na(density_data$somite_count), ]

density_data <- assign_class(density_data, 
                             AC_midseg_threshold, 
                             RC_midseg_threshold, 
                             AC_lateseg_threshold, 
                             RC_lateseg_threshold
)

dmin = 0.001
dmax = 0.004

density_data <- density_data[density_data$species %in% c(#'R. chillingali', 
  'A. calliptera'), ]


cols_to_keep <- c('psm_id', 'overview_id', 'embryo_id', 
                  'species', 'somite_count', 
                  'psm_volume', 'total_cell_number', 
                  'class', 
                  # 'density', 
                  'max_stage', 'propsom'
)
# Select columns starting with 'region' and those in cols_to_keep
density_data_sub <- density_data %>%
  select(all_of(cols_to_keep), starts_with('Region'))


# Melt the dataframe
melted_df <- density_data_sub %>%
  pivot_longer(
    cols = starts_with("Region"), # Columns starting with 'region'
    names_to = c("Region_number", ".value"), # Split column names into parts
    names_pattern = "Region_(\\d+)_(.*)"    # Regex to extract parts
  ) %>%
  mutate(
    Region_number = case_when(
      Region_number == 1 ~ 1,
      Region_number == 2 ~ 2,
      Region_number == 6 ~ 3,
      Region_number == 24 ~ 4,
      Region_number == 120 ~ 5,
      TRUE ~ as.numeric(Region_number) # Keep other values as-is, if any
    )
  )

model <- aov(density ~ propsom +class, 
             data = density_data)
summary(model)



cell_density <- ggplot(density_data, 
                       aes(x = propsom*100, 
                           y = density, 
                           color = species)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = c('black')) + 
  stat_poly_line(se = FALSE) + 
  stat_poly_eq(use_label("R2")) + 
  xlim(0, 100) + 
  ylim(dmin, dmax)+ 
  xlab('% somitogenesis') + 
  ylab('Density (cells per \U00B5m\U00B3)') + 
  ggtitle('E. Density over time')+ 
  theme(legend.position = 'none')


cell_density

# melted_df <- melted_df[
#   melted_df$species %in% c('R. chillingali', 'A. calliptera'), 
#   ]


summary(aov(density~propsom+Region_number, 
            data = melted_df))






density_across_stage <- ggplot(
  melted_df[melted_df$class=='early',], 
  aes(x = Region_number, 
      y = density, colour = somite_count, 
      shape = species, 
      group = psm_id)) + 
  # geom_jitter(width = 0.1, height = 0) + 
  geom_point(size = 1)+ 
  geom_line() + 
  scale_colour_viridis_c(name="Somite count") + theme_bw() + 
  # facet_grid(~species) +  
  ylab('Density (cells per \U00B5m\U00B3)') + 
  xlab('Region number\n(1: most anterior; 5: most posterior)')  + 
  ylim(dmin, dmax)


density_across_stage

###
model <- aov(density~propsom+Region_number, 
             data = melted_df[melted_df$class == 'early', ])
# anova(model)
anova.model <- anova(model)
summary(model)
# residuals are normal
hist(residuals(model),
     col="darkgray")



stat.lm <- lm(density~propsom+Region_number,
              data = melted_df[melted_df$class == 'early', ])
anova(stat.lm)
# eta_squared(stat.lm)

stat.lm <- lm(density~Region_number,
              data = melted_df[melted_df$class=='early',])
anova(stat.lm)
# eta_squared(stat.lm)


# residulas vs fitted data looks OK
plot(fitted(model),
     residuals(model))

density_plt <- cell_density +  density_across_stage + 
  plot_layout(widths = c(1,1), 
              axes = 'collect', 
              axis_titles = "collect") + 
  plot_annotation(
    title = 'PSM density increases over time and differs at early segmentation'
  )


density_plt

w = 10
h = 3

ggsave('../graphs/Cell_density.png', density_plt, width = 5, height = h, dpi = 300)

# Compute the ANOVA for each class



density_by_class <- ggplot(
  melted_df, 
  aes(x = Region_number, y = density, 
      # colour = 'black', 
      group = psm_id)) + 
  # geom_jitter(width = 0.1, height = 0) + 
  geom_point(color = 'black')+ 
  geom_line(color = 'black') + 
  # scale_colour_viridis_c(name="somite count") + 
  theme_bw() + 
  facet_grid(~class) +  ylab('Density (cells per \U00B5m\U00B3)') + 
  xlab('Region number (1: most anterior; 5: most posterior)')  + 
  ylim(dmin, dmax) + 
  ggtitle('F. Density gradient at early segmentation') #+ 

density_by_class

density_plt <- cell_density +  density_by_class + 
  plot_layout(widths = c(1.5, 3), 
              axes = 'collect', 
              axis_titles = "collect") 

density_plt


ggsave('../graphs/density_by_class.png', density_plt, 
       width = 8, height = 3, dpi = 300)

mean(melted_df[(melted_df$class=='early') & 
                 (melted_df$Region_number=='1'), 
               ]$density)


mean(melted_df[(melted_df$class=='early') & 
                 (melted_df$Region_number=='5'), 
]$density)

mean(melted_df[(melted_df$class=='early') & 
                 (melted_df$Region_number=='1'), 
]$density) / mean(melted_df[(melted_df$class=='early') & 
                              (melted_df$Region_number=='5'), 
]$density)

