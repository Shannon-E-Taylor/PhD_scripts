library('ggplot2')
library('ggpmisc')
library('effectsize')
library('patchwork')

data <- read.csv('../data/somite_counts.csv')
cichlid_colours = c('#C3BE46FF', '#4D4D4DFF',  '#00CCFF') # calliptera, Ramphos, zebras


datesfemales_to_avoid = c("13/03/2023", "24/03/2024")

data <- data[!data$date %in% datesfemales_to_avoid, ]

data_sub <- data[
  !(data$species=='Astatotilapia calliptera' & data$min_age >65), 
  ]

plt <- ggplot(data, 
       aes(x = min_age, y = somites, 
           color = species)) + 
  geom_hline(yintercept=39, 
             color = cichlid_colours[2], 
             linetype = 'longdash', 
             linewidth = 1, 
             alpha = 0.5) + 
  geom_hline(yintercept=33, 
             color = cichlid_colours[1], 
             linetype = 'longdash', 
             linewidth = 1, 
             alpha = 0.5) + 
  geom_point(size = 2, pch = 20) + 
  stat_poly_line(data = data_sub, 
                 se = FALSE, 
                 linewidth = 1, method = 'lm') + 
  scale_colour_manual(values = cichlid_colours) + 
  stat_poly_eq(data = data_sub, 
               use_label("eq", "R2"), 
               label.x = "right", label.y = "bottom", 
               method = 'lm') + 
  theme_bw() + 
  ylab('Somite number') + 
  xlab('Time (hpf)') + 
  ggtitle('B. Somitogenesis rate') + 
  theme(legend.position = 'bottom', 
        legend.margin=margin(-10, 0, 0, 0))

plt 

# ggsave('../graphs/somitogenesis_rate.png', 
#        plt, 
#        dpi = 300, width = 1500, height = 900, 
#        units = 'px')


data <- read.csv('../data/pharyngula_stage_somite_counts.csv')
data <- data[data$somite_count > 0, ]


metadata <- read.csv('../data/pharyngula_stage_scans_metadata.csv')
metadata$species <- NULL

data <- merge(metadata, data, by.x = 'image_id', by.y = 'id')

data<- data[data$species %in% c('A. calliptera', 'R. chillingali'), ]


ac_mean_som <- mean(data[data$species == 'A. calliptera',]$somite_count)
ac_sd_som <- sd(data[data$species == 'A. calliptera',]$somite_count)

rc_mean_som <- mean(data[data$species == "R. chillingali",]$somite_count)
rc_sd_som <- sd(data[data$species == "R. chillingali",]$somite_count)

ac <- paste("A. calliptera: ", round(ac_mean_som, 1), "", 
            "\u00B1", " ", round(ac_sd_som, 1))

rc <- paste("R. chillingali: ", round(rc_mean_som, 1), "", 
            "\u00B1", " ", round(rc_sd_som, 1)) 

counts <- ggplot(data, aes(x = somite_count, fill = species)) + 
  geom_histogram(binwidth = 1, color = "white", position = "stack",) + 
  scale_fill_manual(values = cichlid_colours) + 
  theme_bw() + 
  theme(legend.position = 'bottom') + 
  annotate('text', x=-Inf, y=Inf, vjust=1.2, hjust=-.1, 
           label = paste(ac, rc, sep = '\n'), 
           # color = cichlid_colours[1]
           ) + 
  # annotate('text', x=Inf, y=Inf, vjust=2, hjust=1.1, 
  #          label = paste(rc), 
  #          # color = cichlid_colours[1]
  # ) + 
  ggtitle("A. Somite counts") + 
  labs(x = "Somite number", y = "Frequency") + 
  theme(legend.position = 'bottom') 

counts + plt + 
  plot_layout(design = 'ABB') & 
  theme(legend.position = 'bottom') 

h = 3.3

ggsave('../graphs/somitogenesis_rate.png', 
       width = 10, height = 4, 
       dpi = 600)
# 
# stat.lm <- lm(somites ~ min_age*species, 
#               data = data_sub)
# anova(stat.lm)
# eta_squared(stat.lm)
# 
# # supp figure 
# ggplot(data, 
#        aes(x = min_age, y = somites, 
#            color = species, shape = date)) + 
#   geom_point() + 
#   scale_shape_manual(values=1:21) +
#   stat_poly_line(data = data_sub, se = FALSE, linewidth = 0.1) + 
#   scale_colour_manual(values = cichlid_colours) + 
#   # stat_poly_eq(data = data_sub, use_label("eq", "R2")) + 
#   theme_bw() 
