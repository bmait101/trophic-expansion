
# This script tests for effects of the longitudinal gradient (PC1)
# on fish community responces. Specifically:

# 1) fish abundance
# 2) fish density
# 3) fish taxa richness
# 4) fish density as a function of fish taxa richness

#===============================================================================
# Load libraries and custom functions
#===============================================================================
library(tidyverse)
library(here)
library(mgcv)
library(cowplot)

source("code/fx_theme_pub.R", chdir = TRUE)
theme_set(theme_Publication())
#===============================================================================
# Load data
#===============================================================================

data_all <- read_csv(here("data-derived", "fish_structure_summary.csv"))
gradient <- read_csv(here("data-derived", "data_PCA_results.csv")) 

# Summarize data by site (across years) for analysis and plotting
data <- data_all %>% 
  filter(sample_year %in% c(2016, 2017, 2018)) %>% 
  select(-stream_name) %>% 
  unite(site_yr, sample_year, site_id, remove = FALSE)%>%
  left_join(gradient, by = "site_id") %>% 
  mutate_if(is.character, as.factor)
data

# Summarize data by site (across years) for analysis and plotting
data_summ <- data %>% 
  select(-PC1, -PC2, -site_group) %>% 
  group_by(site_id) %>% 
  summarise(n = n(),
            mean_abund = mean(abund, na.rm = TRUE),
            sd_abund = sd(abund, na.rm = TRUE),
            se_abund = sd_abund/sqrt(n),
            mean_density = mean(density, na.rm = TRUE),
            sd_density = sd(density, na.rm = TRUE),
            se_density = sd_density/n, 
            mean_richness = mean(richness, na.rm = TRUE),
            sd_richness = sd(richness, na.rm = TRUE), 
            se_richness = sd_richness/n, 
            mean_shannon = mean(shannon, na.rm = TRUE),
            sd_shannon = sd(shannon, na.rm = TRUE), 
            se_shannon = sd_shannon/n, 
            mean_simpson = mean(simpson, na.rm = TRUE),
            sd_simpson = sd(simpson, na.rm = TRUE), 
            se_simpson = sd_simpson/n) %>%
  left_join(gradient, by = "site_id") %>% 
  mutate_if(is.character, as.factor)
data_summ

#===============================================================================
# Fish taxa richness
#===============================================================================
# Assumptions:
plot(data$richness ~ data $PC1)
plot(log(data$richness) ~ data $PC1)

# No zero inflation. 

hist(data$richness)
hist(log(data$richness))
shapiro.test(data$richness)
shapiro.test(log10(data$richness))

# Fit models
# Fit linear model
lm_rich <- lm(data = data, richness ~ PC1)
summary(lm_rich)

lm_rich_log <- lm(data = data, log(richness) ~ PC1)
summary(lm_rich_log)

glm_rich <- glm(richness ~ PC1, family = Gamma(link = "log"), data = data)
summary(glm_rich)

gam_rich <- data %>% gam(richness ~ s(PC1, k = 5), data = ., method = "REML")
gam_richre <- data %>% gam(richness ~ s(PC1, k = 5) + s(site_id, bs="re"), data = ., method = "REML")
summary(gam_rich)
plot(gam_rich, residuals = TRUE, pch = 1, cex = 1)
gratia::appraise(gam_rich)


summary(lm_rich_log)

p.richness <- 
  data %>%
  ggplot(aes(x = PC1, y = richness)) +    
  geom_smooth(method = "lm", color = "black") + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_log10(limits = c(1,50), breaks=c(1,10,50)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", y = "log Fish taxa richness", 
       fill = "Stream System") + 
  annotate(geom = "text", x = 0.5, y = 50, parse = TRUE, size = 4, hjust = 0,
           label = as.character(expression(paste(F['1,38']==117.61,", ",italic(P),' < ',0.001)))) +
  annotate(geom = "text", x = 0.5, y = 30, parse = TRUE, size = 4,hjust = 0,
           label = as.character(expression(paste(R['adj']^2==0.75)))) +
  theme(legend.position = c(0.8, 0.2))
p.richness

# ggsave(filename = "figs/fish_community_richness.pdf", 
#        plot = p.richness, device = cairo_pdf, 
#        units = "in", width = 5.5, height = 5)

#===============================================================================
# Fish absolute abundance
#===============================================================================
# Assumptions:
plot(data$abund ~ data $PC1)
plot(log(data$abund) ~ data $PC1)

hist(data$abund)
hist(log10(data$abund))
shapiro.test(data$abund)
shapiro.test(log(data$abund))

# Some zero in flation, but the log transformation seems to have worked out

# Fit models
# Fit linear model
lm_abund <- lm(data = data, abund ~ PC1)
summary(lm_abund)

lm_abund_log <- lm(data = data, log(abund) ~ PC1)
summary(lm_abund_log)

# Fit GAM
gam_abund <- data %>% gam(abund ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_abund)
gratia::draw(gam_abund)

# # Get model fit and 95% CIs and bind to raw data
# gam_abund_pred <- lm_abund_log %>% 
#   predict(newdata = data, se.fit = TRUE) %>% 
#   as_tibble() %>% 
#   mutate(lwr = fit - 2 * se.fit,
#          upr = fit + 2 * se.fit) %>% 
#   select(-se.fit)

# Custom plot using GAM predictions
summary(lm_abund_log)

p.abundance <- 
  data %>%
  ggplot(aes(x = PC1, y = abund)) +    
  geom_smooth(method = "lm", color = "black") + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_log10(limits = c(10,1000), breaks=c(10,100,1000)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", y = "log Fish abundance (n0. indiv.)",
       fill = "Stream System") +
  annotate(geom = "text", x = 0.5, y = 1000, parse = TRUE, size = 4, hjust = 0, vjust = 1,
           label = as.character(expression(paste(F['1,38']==211.7,", ",italic(P),' < ',0.001)))) +
  annotate(geom = "text", x = 0.5, y = 700, parse = TRUE, size = 4,hjust = 0, vjust = 1,
           label = as.character(expression(paste(R['adj']^2==0.85)))) 
p.abundance

# ggsave(filename = "figs/fish_community_abund.pdf", 
#        plot = p.abundance, device = cairo_pdf, 
#        units = "in", width = 5.5, height = 5)

#===============================================================================
# Fish density
#===============================================================================

# Assumptions:
plot(data$density ~ data$PC1)
plot(log10(data$density) ~ data$PC1)

hist(data$density)
hist(log(data$density))
shapiro.test(data$density)
shapiro.test(log10(data$density))

# Fit models
# Fit linear model
lm_density <- lm(data = data, density ~ PC1)
summary(lm_density)

lm_density_log <- lm(data = data, log(density) ~ PC1)
summary(lm_density_log)

# Fit GAM
gam_density <- data %>% gam(density ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_density)

# Get model fit and 95% CIs and bind to raw data
density_pred <- lm_density_log %>% 
  predict(newdata = data, se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit)

# Custom plot using GAM predictions
summary(lm_density_log)

p.desnity <- 
  data %>%  
  ggplot(aes(x = PC1, y = density)) +    
  geom_smooth(method = "lm", color = "black") + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_log10(limits = c(0.01,0.15), breaks=c(0.01,0.07, 0.15)) +
  labs(title = "", x = "Long. gradient (PC1)", y = bquote('log Fish density ('*'no' ~ m^-2*')'), 
       fill = "Stream System") +
  annotate(geom = "text", x = 0.5, y = 0.15, parse = TRUE, size = 4, hjust = 0, vjust = 1,
           label = as.character(expression(paste(F['1,38']==4.28,", ",italic(P)==0.05)))) +
  annotate(geom = "text", x = 0.5, y = 0.12, parse = TRUE, size = 4,hjust = 0, vjust = 1,
           label = as.character(expression(paste(R['adj']^2==0.08)))) 
p.desnity

# ggsave(filename = "figs/fish_community_density.pdf", 
#        plot = p.desnity, device = cairo_pdf, 
#        units = "in", width = 5.5, height = 5)


#===============================================================================
# Fish diveristy
#===============================================================================

hist(data$shannon)
hist(log(data$shannon))
shapiro.test(data$shannon)
shapiro.test(sqrt(data$shannon))

# Fit models
# Fit linear model
lm_shannon <- lm(data = data, shannon ~ PC1)
summary(lm_shannon)

lm_shannon_log <- lm(data = data, log10(shannon) ~ PC1)
summary(lm_shannon_log)

# Fit GAM
gam_shannon <- data %>% gam(shannon ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_shannon)
plot(gam_shannon)


p.shannon <- 
  data %>% 
  ggplot(aes(x = PC1, y = shannon)) +    
  geom_smooth(method = "lm", color = "black") + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") +   
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3), breaks = seq(0,3,1)) +
  labs(title = "", x = "Long. gradient (PC1)", y = "Shannon diversity", 
       fill = "Stream System") +
  annotate(geom = "text", x = 0.5, y = 3, parse = TRUE, size = 4, hjust = 0, vjust = 1,
           label = as.character(expression(paste(F['1,38']==66.65,", ",italic(P),' <',0.001)))) +
  annotate(geom = "text", x = 0.5, y = 3*0.9, parse = TRUE, size = 4,hjust = 0, vjust = 1,
           label = as.character(expression(paste(R['adj']^2==0.63)))) 
p.shannon



#===============================================================================
# Panel Plot
#===============================================================================
# ggpubr::ggarrange(p.richness, p.abundance, p.desnity,
#                   nrow = 1, ncol = 3,
#                   common.legend = TRUE, legend = "top")

panel <- plot_grid(
  p.richness + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.desnity + theme(plot.margin=unit(c(1,1,1,1),"mm"),
                    legend.position = "none"), 
  p.shannon + theme(plot.margin=unit(c(1,1,1,1),"mm"),
                    legend.position = "none"),
  labels = c("A", "B", "C"),
  #hjust = -1,
  nrow = 1, ncol = 3, 
  align = "vh")
panel

# ggsave(filename = "figs/fish-community_panel.pdf", plot = last_plot(), 
#       device = cairo_pdf,  units = "in", width = 10, height = 3)
#===============================================================================
