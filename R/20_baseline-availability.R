# Baseline availability analyses

## prep 

# packages
library(tidyverse) 
library(here)
library(cowplot)

# library(magrittr)  
library(mgcv)  
library(gratia)

source(here("R", "fx_theme_pub.R"))
theme_set(theme_Publication())

## Data -------------

gradient <- read_csv(here("out", "data_PCA_results.csv"))
chla <- read_csv(here("data", "chla.csv"))
afdm <- read_csv(here("data", "afdm.csv"))


# Wrangle  ----------------------------

# chla
chla <- chla %>% 
  # Remove Horse Creek Sites
  filter(!site_id %in% c("HC00","HC01","HC02","HC-DOS","HC-SLB")) %>% 
  # Remove the spring 2017 sample event
  filter(sample_hitch != "H0") %>% 
  # Get rid of missing data
  drop_na(BeforeAcid_ul) %>% 
  # Merge the PC1
  left_join(gradient, by = "site_id") %>% 
  # Wrnagle and calculate response metrics
  mutate(sample_hitch = fct_recode(sample_hitch, 
                                   June = "H1", Jul = "H2", Aug = "H3"),
         sample_hitch = fct_relevel(sample_hitch, 
                                    levels = c("June","Jul","Aug")), 
         sample_year = as.factor(sample_year), 
         # Response metrics
         dilution_fct = mlSamp / (mlSamp + mlEtoh), 
         chla_ugmL = 2.05 * (BeforeAcid_ul-AfterAcid_ul) * 
           (ExtractVol_L/filtered_vol_L) * 
           (1/dilution_fct) / 1000, 
         chla_mgmL = chla_ugmL * (1/1000), 
         chla_ug_cm2 = 2.05 * (BeforeAcid_ul-AfterAcid_ul) * 
           (ExtractVol_L/area_cm2) * 
           (1/dilution_fct)) %>% 
  select(filter_id, sample_type, sample_year, sample_hitch, PC1, 
         chla_ugmL, 
         chla_mgmL, 
         chla_ug_cm2)

# AFDM
afdm <- afdm %>%
  filter(!site_id %in% c("HC00","HC01","HC02","HC-DOS","HC-SLB")) %>% 
  filter(sample_hitch != "H0") %>% 
  drop_na(dry_mass_g) %>% 
  left_join(gradient, by = "site_id") %>% 
  mutate(sample_hitch = fct_recode(sample_hitch, 
                                   June = "H1", Jul = "H2", Aug = "H3"),
         sample_hitch = fct_relevel(sample_hitch, 
                                    levels = c("June","Jul","Aug")), 
         sample_year = as.factor(sample_year), 
         # Reponse mmetrics
         AFDM_g = dry_mass_g - ash_mass_g, 
         AFDM_mg_L = ((AFDM_g / vol_filtered_ml) * 1000) * 1000, 
         AFDM_ug_L = ((AFDM_g / vol_filtered_ml) * 1000), 
         AFDM_mg_cm2 = (AFDM_g / area_cm2) * 1000) %>% 
  select(sample_id, sample_type, sample_year, sample_hitch, PC1, 
         AFDM_mg_L, 
         AFDM_ug_L, 
         AFDM_mg_cm2)

# Model fitting: chla seston -----------------------------

hist(chla$chla_ugmL)
m1 <- gam(chla_ugmL ~ 
            sample_hitch +
            s(PC1, by = sample_hitch, k = 3), 
  data = chla %>% filter(sample_type == "seston") , 
  method = "REML")

summary(m1)
appraise(m1)
draw(m1, residuals = TRUE)

# Get model fit and 95% CIs and bind to raw data
fit_gam <- m1 %>% 
  predict(newdata = chla %>% filter(sample_type == "seston"), se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit) %>% 
  bind_cols(chla %>% filter(sample_type == "seston"))

# Plot predictions with raw data
p.chla.ses <- 
  fit_gam %>%
  ggplot(aes(x = PC1, y = chla_ugmL, 
             fill = sample_hitch, color = sample_hitch)) +
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(aes(PC1, fit, linetype = sample_hitch), size = 1) +
  geom_point(aes(shape = sample_year), 
             color = "black", size = 2, alpha = 0.5, shape = 21) + 
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_viridis_d(option = "C") +
  scale_fill_viridis_d(option = "C") +
  # scale_x_continuous(expand = c(0, 0), limits = c(0,9), breaks = seq(0,9,1)) + 
  # scale_y_continuous(limits = c(-0.0008, 0.01), 
  #                    labels = scales::number_format(accuracy = 0.001)) +
  labs(x = "PC1", color = "Month", fill = 'Month', shape = "Year",
       y = expression(paste("Seston chl ", alpha, " (", mu, "g mL"^-1, ")")))  +
  guides(linetype = 'none', 
         color = 'none', 
         fill = guide_legend(override.aes = list(shape = 21))) + 
  annotate(geom = "text", label = "'Deviance expl.' == '74.9%'",
           x = min(data$PC1), y = 0.0095, hjust = 0, parse = TRUE, size = 4) 
p.chla.ses

# Model fitting: chla fbom ----

m2 <- gam(chla_ugmL ~ 
            sample_hitch + 
            s(PC1, by = sample_hitch, k = 3), 
  data = chla %>% filter(sample_type == "FBOM"), 
  method = "REML")

summary(m2)
appraise(m2)
draw(m2, residuals = TRUE)

# Get model fir and 95% CIs and bind to raw data
fit_gam <- m2 %>% 
  predict(newdata = chla %>% filter(sample_type == "FBOM"), se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit) %>% 
  bind_cols(chla %>% filter(sample_type == "FBOM"))

# Plot predictions with raw data
p.chla.fbom <- 
  fit_gam %>%
  ggplot(aes(x = PC1, y = chla_ugmL,  
             fill = sample_hitch, color = sample_hitch)) +
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(aes(PC1, fit, linetype = sample_hitch), size = 1) +
  geom_point(aes(shape  = sample_year), color = "black", size = 2, alpha = 0.5, shape = 21) + 
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_viridis_d(option = "C", guide = 'none') +
  scale_fill_viridis_d(option = "C") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,9), breaks = seq(0,9,1)) + 
  scale_y_continuous(limits = c(0,1.5)) + 
  labs(x = "PC1", color = "Month", fill = 'Month',
       y = expression(paste("FBOM chl ", alpha, " (", mu, "g mL"^-1, ")")))  +
  guides(linetype = 'none') + 
  annotate(geom = "text", label = "'Deviance expl.' == '79.5%'",
           x = min(data$PC1), y = 1.4, hjust = 0, parse = TRUE, size = 4) 
p.chla.fbom

# Model fitting: chla biofilm ----

# GAM and predictions
m3 <- gam(chla_ug_cm2 ~ 
            sample_hitch + 
            s(PC1, by = sample_hitch, k = 3), 
          data = chla %>% filter(sample_type == "biofilm") , 
          method = "REML")

summary(m3)

fit_gam <- m3 %>% 
  predict(newdata = chla %>% filter(sample_type == "biofilm"),se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit) %>% 
  bind_cols(chla %>% filter(sample_type == "biofilm") )

# Plot predictions with raw data
p.chla.biofilm <- 
  fit_gam %>%
  ggplot(aes(x = PC1, y = chla_ug_cm2, 
             fill = sample_hitch, color = sample_hitch)) +
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(aes(PC1, fit, linetype = sample_hitch), size = 1) +
  geom_point(aes(shape  = sample_year), color = "black", size = 2, alpha = 0.5, shape = 21) + 
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_viridis_d(option = "C", guide = 'none') +
  scale_fill_viridis_d(option = "C") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,9), breaks = seq(0,9,1)) + 
  scale_y_continuous(limits = c(0,3.5), 
                     labels = scales::number_format(accuracy = 0.1)) + 
  labs(x = "PC1", color = "Month", fill = 'Month',
       y = expression(paste("Biofilm chl ", alpha, " (", mu, "g cm"^-2, ")")))  +
  guides(linetype = 'none') + 
  annotate(geom = "text", label = "'Deviance expl.' == '83.7%'",
           x = 0, y = 0.01, hjust = 0, parse = TRUE, size = 4) 
p.chla.biofilm

# Model fitting: AFDM seston ----
# Lets check the 3 lm assumptions (Gaussian distribution, no interactions, linear)

m4 <- gam(AFDM_mg_L ~ 
               sample_hitch + 
               s(PC1, by = sample_hitch, k = 3), 
             data =  afdm %>% filter(sample_type == "seston") , 
             method = "REML")

summary(m4)

fit_gam <- m4 %>% 
  predict(newdata = afdm %>% filter(sample_type == "seston") , type = "link", se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit) %>% 
  bind_cols(afdm %>% filter(sample_type == "seston") )

p.afdm.ses <- 
  fit_gam %>%
  ggplot(aes(x = PC1, y = AFDM_mg_L, 
             fill = sample_hitch, color = sample_hitch)) +
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(aes(PC1, fit, linetype = sample_hitch), size = 1) +
  geom_point(aes(shape = sample_year), color = "black", size = 2, alpha = 0.5, shape = 21) + 
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_viridis_d(option = "C", guide = 'none') +
  scale_fill_viridis_d(option = "C") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,9), breaks = seq(0,9,1)) + 
  scale_y_continuous(limits = c(-0.2,2.25)) + 
  labs(x = "PC1", color = "Month", fill = 'Month',
       y = expression(paste("Seston AFDM (mg L"^-1, ")"))) +
  guides(linetype = 'none') + 
  annotate(geom = "text", label = "'Deviance expl.' == '79.4%'",
           x = min(data$PC1), y = 2.2, hjust = 0, parse = TRUE, size = 4) 
p.afdm.ses

# Model fitting: AFDM fbom ----

m5 <- gam(
  AFDM_ug_L ~ sample_hitch + s(PC1, by = sample_hitch, k = 3), 
  data = afdm %>% filter(sample_type == "fbom") , method = "REML"
)

summary(m5)

fit_gam <- m5 %>% 
  predict(newdata = afdm %>% filter(sample_type == "fbom"), se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit) %>% 
  bind_cols(afdm %>% filter(sample_type == "fbom"))

# Plot predictions with raw data
p.afdm.fbom <- 
  fit_gam %>%
  ggplot(aes(x = PC1, y = AFDM_ug_L, 
             fill = sample_hitch, color = sample_hitch)) +
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(aes(PC1, fit, linetype = sample_hitch), size = 1) +
  geom_point(aes(shape = sample_year), color = "black", size = 2, alpha = 0.5, shape = 21) + 
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_viridis_d(option = "C", guide = 'none') +
  scale_fill_viridis_d(option = "C") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,9), breaks = seq(0,9,1)) +
  scale_y_continuous(limits = c(-0.1,1.5)) +
  labs(x = "PC1", color = "Month", fill = 'Month',
       y = expression(paste("FBOM AFDM (mg L"^-1, ")"))) +
  guides(linetype = 'none') + 
  annotate(geom = "text", label = "'Deviance expl.' == '73.3%'",
           x = min(data$PC1), y = 1.5, hjust = 0, parse = TRUE, size = 4) 
p.afdm.fbom

# Model fitting: AFDM biofilm ----

m6 <- gam(
  AFDM_mg_cm2 ~ sample_hitch + s(PC1, by = sample_hitch, k = 3), 
  data =  afdm %>% filter(sample_type == "biofilm") , method = "REML"
)

summary(m6)

fit_gam <- m6 %>% 
  predict(newdata = afdm %>% filter(sample_type == "biofilm"),se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit) %>% 
  bind_cols(afdm %>% filter(sample_type == "biofilm"))

# Plot predictions with raw data
p.afdm.bio <- 
  fit_gam %>%
  ggplot(aes(x = PC1, y = AFDM_mg_cm2, 
             fill = sample_hitch, color = sample_hitch)) +
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.3, color = NA) +
  geom_line(aes(PC1, fit, linetype = sample_hitch), size = 1) +
  geom_point(aes(shape = sample_year), color = "black", size = 2, alpha = 0.5, shape = 21) + 
  #scale_shape_manual(values = c(21, 24)) +
  scale_color_viridis_d(option = "C", guide = 'none') +
  scale_fill_viridis_d(option = "C") +
  scale_x_continuous(expand = c(0, 0), limits = c(0,9), breaks = seq(0,9,1)) + 
  scale_y_continuous(limits = c(-0.1,1.7)) + 
  labs(x = "PC1", color = "Month", fill = 'Month',
       y = expression(paste("Biofilm AFDM (mg cm"^-2, ")"))) +
  guides(linetype = 'none') + 
  annotate(geom = "text", label = "'Deviance expl.' == '87.4%'",
           x = min(data$PC1), y = 1.6, hjust = 0, parse = TRUE, size = 4) 
p.afdm.bio


# Panel Plot ------------------------------

p.panel <- plot_grid(
  p.afdm.ses+ theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"mm")), 
  p.afdm.fbom+ theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"mm")), 
  p.afdm.bio+ theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"mm")), 
  p.chla.ses + theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"mm")), 
  p.chla.fbom + theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"mm")), 
  p.chla.biofilm + theme(legend.position="none", plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A", "B", "C", "D", "E", "F"),
  #hjust = -1,
  ncol = 3, nrow = 2, 
  align = "vh"
)

# extract the legend from one of the plots
legend <- get_legend( p.chla.ses )

# Plot
panel <- plot_grid(p.panel, legend, rel_widths = c(3, .3))
panel

# Save it
ggsave(filename = here("out", "resource_availability.pdf"), 
       plot = panel, device = cairo_pdf,
       units = "in", width = 10, height = 5.5)

