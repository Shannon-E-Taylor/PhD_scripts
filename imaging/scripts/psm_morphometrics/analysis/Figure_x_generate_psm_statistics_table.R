

options(knitr.kable.NA = '')

library('kableExtra')

library('ggplot2')
library('data.table')
library('tidyr')
library('dplyr')


path_to_paper = '../../../../../Apps/Overleaf/Axial development in cichilds/tables/'

list.files(path_to_paper)



# # set working directory to source file location. 
# this.dir <- dirname(parent.frame(2)$ofile)
# setwd(this.dir)


#############
# FUNCTIONS #
#############

compute_statistics <- function(data, variable){
  
  # variable~stage  
  least_means_regression_1 <- lm(data[[variable]] ~ somite_count, data = data)
  R2_1 = summary(least_means_regression_1)$r.squared
  slope_1 = summary(least_means_regression_1)$coefficients[2,1]
  std_error_slope_1 = summary(least_means_regression_1)$coefficients[2,2]
  CI = compute_CI(slope_1, std_error_slope_1)[[1]]
  # single variable anova. 
  anova_1 <- anova(least_means_regression_1)
  p_value <- anova_1$`Pr(>F)`[[1]]
  F_value <- anova_1$`F value`[[1]]
  
  # variable~stage  
  
  least_means_regression_2 <- lm(data[[variable]] ~ somite_count + species, data = data)
  R2_2 = summary(least_means_regression_2)$r.squared
  slope_2 = least_means_regression_2$coefficients[[2]]    
  std_error_slope_2 = summary(least_means_regression_2)$coefficients[2,2]
  
  sp_difference = least_means_regression_2$coefficients[[3]]
  sp_difference_error = summary(least_means_regression_2)$coefficients[3,2]
  
  p_species = anova(least_means_regression_2)$`Pr(>F)`[2]
  p_count = anova(least_means_regression_2)$`Pr(>F)`[1]
  
  
  # choose better model..
  av <- anova(least_means_regression_1, least_means_regression_2)
  
  stats_out <- c(
    # variable~stage
    'R2' = R2_1, 
    'slope' = slope_1, 
    'slope_error' = std_error_slope_1,
    'F' = F_value, 
    'p' = p_value, 
    # variable~stage
    'RSS_1' = av$RSS[[1]], 
    'RSS_2' = av$RSS[[2]], 
    'p_model2_better' = av$`Pr(>F)`[[2]],
    'R2_2' = R2_2, 
    'slope_2' = slope_2, 
    'slope_error_2' = std_error_slope_2, 
    'species_difference' = sp_difference, 
    'species_difference_error' = sp_difference_error, 
    'p_count' = p_count, 
    'p_species' = p_species
  )
  
  # test for interaction variable if species is a significant stratifier 
  if (p_species < 0.05){
    least_means_regression_3 <- lm(
      data[[variable]] ~ somite_count + species + species*somite_count, 
      data = data)
    R2_3 = summary(least_means_regression_3)$r.squared
    p_species = anova(least_means_regression_3)$`Pr(>F)`[2]
    p_count = anova(least_means_regression_3)$`Pr(>F)`[1]
    p_interaction = anova(least_means_regression_3)$`Pr(>F)`[3]
    
    stats_to_append <- c(
      'R2_3' = R2_3, 
      'p_species' = p_species, 
      'p_count' = p_count, 
      'p_interaction' = p_interaction
    )
    
    stats_out = append(stats_out, stats_to_append)
  }
  
  table <- as.data.frame(t(stats_out))
  
  # Set appropriate column names using the vector names
  colnames(table) <- names(stats_out)
  rownames(table) <- variable
  return(table)
}

compute_CI <- function(mean, sd){
  lb = round(mean - (1.96*sd), digits = 2)
  ub = round(mean + (1.96*sd), digits = 2)
  ub_ = max(ub, lb)
  lb_ = min(ub, lb)
  return (paste(lb_, ' \u2014 ', ub_, sep='') )
}

compute_statistics_display <- function(p_value){
  if (is.na(p_value)){
    val = NA
  }
  else if (p_value > 0.05){
    val = 'n.s.'
  }
  else if (p_value < 0.001){
    val = formatC(as.numeric(p_value), format = "e", digits = 1)
  }
  else{
    val = round(p_value, digits = 2)
  }
  return(val)
}

###############
# IMPORT DATA #
###############

psm_data = read.csv('../data/computed_psm_dimensions.csv')

psm_data$somite_count = as.integer(psm_data$somite_count)

psm_data = psm_data[
  (psm_data$species %in% c('A. calliptera', 'R. chillingali')) & (psm_data$somite_count > 0), 
]

drops = c('somite_count_y', 'species_y', 'somite_count_x', 'species_x', 
          'overview_id', 'X', 'xsize', 'zsize', 'psm_id')
psm_data = psm_data[, !names(psm_data) %in% drops]

psm_data$mpz = psm_data$psm - psm_data$noto
psm_data$len_width = psm_data$psm / psm_data$psm_width
psm_data$len_depth = psm_data$psm / psm_data$psm_depth
psm_data$width_depth = psm_data$psm_width / psm_data$psm_depth

psm_data <- na.omit(psm_data, cols= 'species')

psm_data$max_stage = 32
psm_data[psm_data$species == 'R. chillingali', ]$max_stage = 38

psm_data$propsom = psm_data$somite_count / psm_data$max_stage


# exclude outlier. 
psm_data <- psm_data[!(psm_data$mpz > 200 & psm_data$somite_count > 20), ]

AC_midseg_threshold=15
RC_midseg_threshold=17

AC_lateseg_threshold=25
RC_lateseg_threshold=30


# set Class (stage) by somite stage. 
psm_data$class = 'early'

psm_data[
  (psm_data$species == 'A. calliptera') & 
    (psm_data$somite_count > AC_midseg_threshold)
  , ]$class = 'mid'

psm_data[
  (psm_data$species == 'R. chillingali') & 
    (psm_data$somite_count > RC_midseg_threshold)
  , ]$class = 'mid'


psm_data[
  (psm_data$species == 'A. calliptera') & 
    (psm_data$somite_count > AC_lateseg_threshold)
  , ]$class = 'late'

psm_data[
  (psm_data$species == 'R. chillingali') & 
    (psm_data$somite_count > RC_lateseg_threshold)
  , ]$class = 'late'

psm_data$class <- factor(psm_data$class, levels = c('early', 'mid', 'late'))

stats_out = data.frame()

variable_names = c(
  'psm', 'noto', 'psm_width', 'psm_depth', 
  'mpz', 'som', 'n', 
  'len_width', 'len_depth', 'width_depth')

for (variable in variable_names){
  df_out = compute_statistics(psm_data, variable)
  stats_out <- dplyr::bind_rows(stats_out, data.frame(df_out))
  
}

p_thresh = 0.05


p_values_1 = stats_out$p < p_thresh
stats_out$p = sapply(stats_out$p, compute_statistics_display)


p_values_2 = stats_out$p_species < p_thresh
stats_out$p_species = sapply(stats_out$p_species, compute_statistics_display)

# this is the same as p3, so delete.... 
stats_out$p_model2_better = NULL


p_values_3 = (!is.na(stats_out$p_interaction) & stats_out$p_interaction < p_thresh)
stats_out$p_interaction = sapply(stats_out$p_interaction, compute_statistics_display)



# other stats too... 
stats_out$p_count =  sapply(stats_out$p_count, compute_statistics_display)
stats_out$p_species.1 =  sapply(stats_out$p_species.1, compute_statistics_display)
stats_out$p_count.1 =  sapply(stats_out$p_count.1, compute_statistics_display)


# one 
stats_out$CI <- mapply(compute_CI, stats_out$slope, stats_out$slope_error)
stats_out <- stats_out %>% relocate(CI, .after=slope)
stats_out$slope_error = NULL

stats_out$CI2 <- mapply(compute_CI, stats_out$slope_2, stats_out$slope_error_2)
stats_out <- stats_out %>% relocate(CI2, .after=slope_2)
stats_out$slope_error_2 = NULL

stats_out$CI3 <- mapply(compute_CI, stats_out$species_difference, stats_out$species_difference_error)
stats_out <- stats_out %>% relocate(CI3,.after=species_difference )
stats_out$species_difference_error = NULL

stats_out$RSS_1 <- round(stats_out$RSS_1, digits = 0)
stats_out$RSS_2 <- round(stats_out$RSS_2, digits = 0)




smaller_residual = stats_out$RSS_1 > stats_out$RSS_2



# get column indicies last
p1 <- which(names(stats_out) == 'p') + 1
p2 <- which(names(stats_out) == 'p_species') + 1
p3 <- which(names(stats_out) == 'p_interaction') + 1


where_RSS1 <- which(names(stats_out) == 'RSS_1') + 1
where_RSS2<- which(names(stats_out) == 'RSS_2') + 1







stat_name_header = c(
  # variable ~ stage
  'R\u00B2', 'slope', '95 CI',  'F statistic', 'p-value',
  'RSS model 1', 'RSS model 2', # 'P model 2 better', 
  # variable ~ stage + species
  'R\u00B2', 'slope', '95 CI',
  'species difference', '95 CI species diff.',
  'p count', 'p species', 
  # variable ~ stage:species
  'R\u00B2', 'p species', 'p count', 
  'p interaction' 
)

row_names <- c(
  'PSM length', 'notocord length', 'PSM width', 'PSM depth', 
  'MPZ length', 'Somite length', 'notocord width', 
  'length:width', 'length:depth', 'depth:width'
)

rownames(stats_out) <- row_names


model_header = c('', 
                 'variable~stage' = 5, 
                 'comparison' = 2, 
                 'variable~stage + species'=7,
                 'variable~stage:species'=4
)

#format(round(a, digits=2), nsmall = 2) 


stats_out %>%
  kbl(
    caption = 'stats',
    digits = 2,
    booktabs = TRUE,
    escape = FALSE,
    col.names = stat_name_header, 
    row.names = TRUE, 
    align = "c", 
    format = 'latex'
  ) %>%
  add_header_above(model_header,
                   # background = "#D3D3D3", 
                   line = TRUE, bold = TRUE
  ) %>%
  column_spec(p1, background = ifelse(p_values_1, "#CAFF70", "white")) %>%
  column_spec(p2, background = ifelse(p_values_2, "#CAFF70", "white")) %>%
  column_spec(p3, background = ifelse(p_values_3, "#CAFF70", "white")) %>%
  
  column_spec(where_RSS1,
              background = ifelse(!smaller_residual, "#9fc5e8", "white")
  ) %>%
  
  column_spec(where_RSS2,
              background = ifelse(smaller_residual, "#9fc5e8", "white")
  ) %>%
  
  kable_classic(full_width = FALSE, 
                latex_options = c("hold_position", 
                                  'repeat_header', 'scale_down'
                )) %>%
  landscape() %>% 
  row_spec(0, bold = FALSE, color = "black") %>%
  column_spec(1, bold = TRUE, color = "black")%>%
  
  save_kable(paste0(path_to_paper, 'PSM_statistics_all.tex'))



stats_out[, 0:14] %>%
  kbl(
    caption = 'stats',
    digits = 2,
    booktabs = TRUE,
    escape = FALSE,
    col.names = stat_name_header[0:14], 
    row.names = TRUE, 
    align = "c", 
    format = 'latex'
  ) %>%
  add_header_above(model_header[0:4],
                   # background = "#D3D3D3", 
                   line = TRUE, bold = TRUE
  ) %>%
  column_spec(p1, background = ifelse(p_values_1, "#CAFF70", "white")) %>%
  column_spec(p2, background = ifelse(p_values_2, "#CAFF70", "white")) %>%
  # column_spec(p3, background = ifelse(p_values_3, "#CAFF70", "white")) %>%
  
  column_spec(where_RSS1,
              background = ifelse(!smaller_residual, "#9fc5e8", "white")
  ) %>%
  
  column_spec(where_RSS2,
              background = ifelse(smaller_residual, "#9fc5e8", "white")
  ) %>%
  
  kable_classic(full_width = FALSE, 
                latex_options = c("hold_position", 
                                  'repeat_header', 'scale_down'
                )) %>%
  landscape() %>% 
  row_spec(0, bold = FALSE, color = "black") %>%
  column_spec(1, bold = TRUE, color = "black")%>%
  save_kable(paste0(path_to_paper, 'PSM_statistics.tex'))
  
caption = 'Results from statistical tests'

stats_out[, 14:18] %>%
  kbl(
    caption = 'stats',
    digits = 2,
    booktabs = TRUE,
    escape = FALSE,
    col.names = stat_name_header[14:18], 
    row.names = TRUE, 
    align = "c", 
    format = 'latex'
  ) %>%
  add_header_above(model_header[4:5],
                   # background = "#D3D3D3", 
                   line = TRUE, bold = TRUE
  ) %>%
  column_spec(p1, background = ifelse(p_values_1, "#CAFF70", "white")) %>%
  column_spec(p2, background = ifelse(p_values_2, "#CAFF70", "white")) %>%
  # column_spec(p3, background = ifelse(p_values_3, "#CAFF70", "white")) %>%
  
  column_spec(where_RSS1,
              background = ifelse(!smaller_residual, "#9fc5e8", "white")
  ) %>%
  
  column_spec(where_RSS2,
              background = ifelse(smaller_residual, "#9fc5e8", "white")
  ) %>%
  
  kable_classic(full_width = FALSE, 
                latex_options = c("hold_position", 
                                  'repeat_header', 'scale_down'
                )) %>%
  landscape() %>% 
  row_spec(0, bold = FALSE, color = "black") %>%
  column_spec(1, bold = TRUE, color = "black")%>%
  save_kable(paste0(path_to_paper, 'species_interactions.tex'))


