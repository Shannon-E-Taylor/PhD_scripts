


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


#############
# MORPHOMETRICS # 
#################



psm_data = read.csv('../data/computed_psm_dimensions.csv')

psm_data$somite_count = as.integer(psm_data$somite_count)

psm_data = psm_data[
  (psm_data$species %in% c('A. calliptera', 'R. chillingali')) & (psm_data$somite_count > 0), 
]

drops = c('somite_count_y', 'species_y', 'somite_count_x', 'species_x', 
          # 'overview_id', 
          'X', 'xsize', 'zsize'#, #'psm_id'
)
psm_data = psm_data[, !names(psm_data) %in% drops]

psm_data = psm_data[!is.na(psm_data$species), ]

# psm_data <- na.omit(psm_data, cols= 'species')

psm_data$max_stage = 32
psm_data[psm_data$species == 'R. chillingali', ]$max_stage = 38

psm_data$propsom = psm_data$somite_count / psm_data$max_stage


# exclude outlier. 
# psm_data <- psm_data[!(psm_data$mpz > 200 & psm_data$somite_count > 20), ]

# method and params in config 
psm_data = assign_class(psm_data, 
                        AC_midseg_threshold, 
                        RC_midseg_threshold, 
                        AC_lateseg_threshold, 
                        RC_lateseg_threshold
)



cichlid_colours = c("A. calliptera"="#C3BE46FF", "R. chillingali"="#4D4D4DFF")



plot_absolute_stage <- function(data, y, species_diff){
  
  plt <- ggplot(data, 
                aes(x = somite_count, y = get(y), 
                    colour = species)
                ) +
    theme_bw() + 
    xlim(0, 40) + 
    xlab('Somite Stage') + 
    theme(legend.position = "bottom") + 
    scale_colour_manual(values = cichlid_colours) + 
    scale_y_continuous(expand = expansion(mult = 0.1)) 
  
  if (species_diff){
    plt <- plt + 
      stat_poly_line(se = FALSE) + 
      stat_poly_eq(use_label(c('eq', 'R2')), 
                   label.x = "right", 
                   label.y = c(0.95, 0.85))  + 
      geom_point()
  }
  
  else {
    plt <- plt + 
      stat_poly_line(se = FALSE, color = 'black', show.legend=FALSE) + 
      stat_poly_eq(color = 'black', 
                   use_label(c('eq', 'R2'))) + 
      geom_point(show.legend=FALSE)
  }
  
  
  return(plt)
}

psm <- plot_absolute_stage(psm_data, 'psm', species_diff = TRUE) + 
  ylab('PSM length (\u00B5m)')  

noto <- plot_absolute_stage(psm_data, 'noto', species_diff = TRUE) + 
  ylab('notocord length (\u00B5m)') 

mpz <- plot_absolute_stage(psm_data, 'mpz', species_diff = TRUE) + 
  ylab('MPZ length (\u00B5m)')

width <- plot_absolute_stage(psm_data, 'psm_width', species_diff = FALSE) + 
  ylab('PSM width (\u00B5m)') + theme(legend.position='none')

depth <- plot_absolute_stage(psm_data, 'psm_depth', species_diff = FALSE) + 
  ylab('PSM depth (\u00B5m)') + theme(legend.position='none')

som <- plot_absolute_stage(psm_data, 'som', species_diff = TRUE) + 
  ylab('Somite length (\u00B5m)') 

len.width <- plot_absolute_stage(psm_data, 'len_width', species_diff = FALSE) + 
  ylab('PSM length : width') 

len.depth <- plot_absolute_stage(psm_data, 'len_depth', species_diff = FALSE) + 
  ylab('PSM length : depth')

width.depth <- plot_absolute_stage(psm_data, 'width_depth', species_diff = FALSE) + 
  ylab('PSM width : depth') 

n <- plot_absolute_stage(psm_data, 'n', species_diff = FALSE) + 
  ylab('Notocord width (\u00B5m)')

som.psm <- plot_absolute_stage(psm_data, 'psm_som', TRUE) + 
  ylab('PSM : somite length') 

psm + noto + mpz + 
  width + depth + som +  
  len.width + len.depth + width.depth + 
  n + som.psm + 
  plot_layout(ncol = 3, axes = "collect", guides = 'collect') & 
  theme(legend.position = 'bottom')


ggsave('../graphs/PSM_mophometrics_supp_data.png', width = 10, height = 10, dpi = 300)

#########
# PCA 
#########



psm_data_for_pca = psm_data
psm_data_for_pca = na.omit(psm_data_for_pca)

psm_data_names =  names(psm_data_for_pca)
keeps = which(!psm_data_names %in% c('species', 'somite_count', 'class', 'propsom', 'overview_id', 'psm_id', 'max_stage'))



numerical_data = psm_data_for_pca[keeps]

res.pca <- prcomp(numerical_data, scale = TRUE)

library('factoextra')

fviz_eig(res.pca, addlabels = TRUE)


pca_values = data.frame(predict(res.pca))

pca_values$species = psm_data_for_pca$species
pca_values$somite_count = psm_data_for_pca$somite_count

pca_values = assign_class(pca_values, 
                          AC_midseg_threshold, 
                          RC_midseg_threshold, 
                          AC_lateseg_threshold, 
                          RC_lateseg_threshold
)



# Calculate the hulls for each group
hull_cyl <- pca_values %>%
  group_by(species, class) %>%
  slice(chull(PC1, PC2))


# extract explained variance for each PC
PC1_perc = round(get_eigenvalue(res.pca)[1,2], digits=0)
PC2_perc = round(get_eigenvalue(res.pca)[2,2], digits=0)
PC3_perc = round(get_eigenvalue(res.pca)[3,2], digits=0)


pca <- ggplot() + 
  # Plot light grey points for all data
  geom_point(data = select(pca_values, c('PC1', 'PC2', 'species')), 
             aes(x = PC1, y = PC2, shape = species), 
             color = 'darkgrey', 
             alpha = 1, size = 1) + 
  # Draw polygons around groups (hull) with desired aesthetics
  geom_polygon(data = hull_cyl, 
               aes(x = PC1, y = PC2, linetype = species), 
               color = 'darkgrey', alpha = 0) + 
  # Add the colored points on top, with somite_count as color
  geom_point(data = pca_values, 
             aes(x = PC1, y = PC2, 
                 color = species, 
                 shape = species)) + 
  # Color scale for somite_count
  # scale_color_viridis_c() + 
  scale_colour_manual(values = cichlid_colours) + 
  theme_bw() + 
  facet_grid( ~ class) + 
  xlab(sprintf("PC1: %.f%% explained variance", PC1_perc)) + 
  ylab(sprintf("PC2: %.f%% explained variance", PC2_perc)) + 
  theme(legend.position = 'bottom') + 
  ggtitle('E. PCA separates embryos by stage and species')

pca

ggsave('../graphs/PCA.png', width = 7, height = 4, dpi = 300)



library('gridGraphics')

design <- '
ABC
DDE'

width + depth + mpz + pca + guide_area() + 
  plot_layout(design = design, axes = "collect", guides = 'collect') & 
  theme(legend.position = 'bottom')


ggsave('../graphs/dimensions_and_pca.png', 
       width = w, height = h*2, 
       dpi = 300)



design <- '
ABCD
EEEF'

width <- width + ggtitle('A. PSM width decreases')+ 
  theme(legend.position = 'none')
depth <- depth + ggtitle('B. PSM depth increases')+ 
  theme(legend.position = 'none')
mpz <- mpz + ggtitle('C. MPZ length differs')+ 
  theme(legend.position = 'none')
som.psm <- som.psm + ggtitle('D. Somite ratio differs')+ 
  theme(legend.position = 'none')


width + depth + mpz + som.psm + 
  pca + plot_spacer() + # cell_density + density_across_stage + 
  plot_layout(design = design, 
              axes = "collect_y", 
              guides = 'collect'
              ) & 
  theme(legend.position = 'bottom')

  
ggsave('../graphs/Figure3_PSM_dimensions_full.png',  width = 12, height = 6, dpi = 600)




ggplot(psm_data, 
       aes(y = mpz/psm, x = somite_count, color = species)) + 
         geom_point()

stat.lm <- lm(
  psm_som ~somite_count+species, 
  data = psm_data)
anova(stat.lm)
eta_squared(stat.lm)

