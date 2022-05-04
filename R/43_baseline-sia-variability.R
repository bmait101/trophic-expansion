
# Baseline resource isotopic variability

# Libraries ----
library(tidyverse)  
library(here)
library(mgcv)
library(gratia)
library(cowplot)
library(broom)

# source("code/fx_theme_pub.R", chdir = TRUE)
theme_set(theme_bw()) 

## Data -------------------------------------

source(here("R", "40_prep-sia-data.R"))
lg <- read_csv(here("out", "data_PCA_results.csv"))

# Tidy up data
baseline <- sia_baselines %>% 
  #filter(taxon_code %in% c("biofilm", "detritus")) %>%
  mutate(sample_hitch = fct_recode(sample_hitch, 
                                   June = "H1", Jul = "H2", Aug = "H3"),
         sample_hitch = fct_relevel(sample_hitch, 
                                    levels = c("June","Jul","Aug")), 
         sample_year = as.factor(sample_year)) %>% 
  select(-stream_name) %>% 
  mutate(yr_site = paste(sample_year, site_id, sep = "_")) %>% 
  filter(yr_site != "2017_LR05")  

baseline <- baseline %>%
  left_join(lg, by = "site_id")

baseline_var <- 
  baseline %>% 
  group_by(sample_year, site_id) %>%
  summarise(d13C_range = max(d13C)-min(d13C),
            d13C_CV = sd(d13C)/mean(d13C),
            d15N_range = max(d15N)-min(d15N),
            d15N_CV = sd(d15N)/mean(d15N)) %>%
  left_join(lg, by = "site_id")
baseline_var

# Cor tests  -----------------------------------------
tidy(cor.test(baseline$d13C, baseline$PC1))
tidy(cor.test(baseline_var$d13C_range, baseline_var$PC1))
tidy(cor.test(baseline_var$d13C_CV, baseline_var$PC1))
tidy(cor.test(baseline$d15N, baseline$PC1))
tidy(cor.test(baseline_var$d15N_range, baseline_var$PC1))
tidy(cor.test(baseline_var$d15N_CV, baseline_var$PC1))


# Pooled d13C mean  ------------------------------------
# Assumptions:
# plot(baseline$d13C ~ baseline$PC1)
# hist(baseline$d13C)

# Fit models
lm_rich <- lm(data = baseline, d13C ~ PC1)
lm_rich %>% broom::glance()

# Fit GAM
gam_rich <- baseline %>% gam(d13C ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_rich)
draw(gam_rich, residuals = TRUE)

# Get model fit and 95% CIs and bind to raw data
# rich_pred <- gam_rich %>% 
#   predict(newdata = baseline, se.fit = TRUE) %>% 
#   as_tibble() %>% 
#   mutate(lwr = fit - 2 * se.fit,
#          upr = fit + 2 * se.fit) %>% 
#   select(-se.fit)

p.bl.c.mean <- baseline %>%
  #bind_cols(rich_pred) %>% 
  ggplot(aes(PC1, d13C)) +
  # geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.5, fill = "grey") +
  # geom_line(aes(PC1, fit), size = 1, color = "black") +
  ggpubr::stat_cor(method = "pearson", label.x = 0.5, label.y = -12, vjust = 1) + 
  geom_point(aes(fill = resource, shape = resource), color = "black", size = 2) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(-36,-12), breaks = seq(-36,-12,6)) +
  scale_shape_manual(values = c(21,22,23,24,25), guide = 'none')  +
  scale_color_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") + 
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression(~{delta}^13*C~'(\u2030)'),
       fill = "Resource") + 
  guides(fill = guide_legend(override.aes = list(shape = c(21,22,23,24,25))))
p.bl.c.mean

# Pooled d13C range ---- -------------------------------------------------
# plot(baseline_var$d13C_range ~ baseline_var$PC1)
# hist(baseline_var$d13C_range)

# Fit models
# Fit linear model
lm_rich <- lm(data = baseline_var, d13C_range ~ PC1)
lm_rich %>% broom::glance()

# Fit GAM
# gam_rich <- baseline_var %>% gam(d13C_range ~ s(PC1, k = 3), data = ., method = "REML")
# summary(gam_rich)
# 
# BIC(lm_rich, gam_rich)
# 
# # Get model fit and 95% CIs and bind to raw data
# rich_pred <- gam_rich %>% 
#   predict(newdata = baseline_var, se.fit = TRUE) %>% 
#   as_tibble() %>% 
#   mutate(lwr = fit - 2 * se.fit,
#          upr = fit + 2 * se.fit) %>% 
#   select(-se.fit)

p.bl.c.range <- 
  baseline_var %>%
  #bind_cols(rich_pred) %>% 
  ggplot(aes(PC1, d13C_range)) +
  ggpubr::stat_cor(method = "pearson", label.x = 0.5, label.y = 16, vjust = 1) + 
  #geom_smooth(se = TRUE, method = "lm", color = "black") + 
  # geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.5, fill = "grey") +
  # geom_line(aes(PC1, fit), size = 1, color = "black") +
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(4,16), breaks = seq(4,16,2)) +
  # ggrepel::geom_text_repel(aes(label = site_id)) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression(~{delta}^13*C~'Range (\u2030)'),
       fill = "Stream System")
p.bl.c.range


# Pooled d13C CV ------------------------------------------------------
# plot(baseline_var$d13C_CV ~ baseline_var$PC1)
# hist(baseline_var$d13C_CV)
# shapiro.test(baseline_var$d13C_CV)
# cor.test(baseline_var$d13C_CV, baseline_var$PC1)

# Fit models
# Fit linear model
lm_rich <- lm(data = baseline_var, d13C_CV ~ PC1)
lm_rich %>% broom::glance()

# Fit GAM
gam_rich <- baseline_var %>% gam(d13C_CV ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_rich)

BIC(lm_rich, gam_rich)

# Get model fit and 95% CIs and bind to raw data
rich_pred <- gam_rich %>% 
  predict(newdata = baseline_var, se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit)

p.bl.c.cv <- 
  baseline_var %>%
  #bind_cols(rich_pred) %>% 
  ggplot(aes(PC1, d13C_CV)) +
  ggpubr::stat_cor(method = "pearson", label.x = 0.5, label.y = -0.04, vjust = 1) + 
  # geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.5, fill = "grey") +
  # geom_line(aes(PC1, fit), size = 1, color = "black") +
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(-0.18,-0.04), breaks = seq(-0.18,-0.04,0.02)) +
  # ggrepel::geom_text_repel(aes(label = site_id)) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression(~{delta}^13*C~'CV'),
       fill = "Stream System")

p.bl.c.cv

# Pooled d15N mean -------------------------------------------------------------

# plot(baseline$d15N ~ baseline$PC1)
# hist(baseline$d15N)
# shapiro.test(baseline$d15N)

# Fit models
# Fit linear model
lm_NR <- lm(data = baseline, d15N ~ PC1)
summary(lm_NR)
tidy(lm_NR)

# Fit GAM
gam_rich <- baseline %>% gam(d15N ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_rich)

AIC(lm_NR, gam_rich)

lm_NR_pred <- lm_NR %>% 
  predict(newdata = baseline, se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit)

p.bl.n.mean <- 
  baseline %>%
  bind_cols(lm_NR_pred) %>% 
  ggplot(aes(PC1, d15N)) +
  # ggpubr::stat_cor(method = "pearson", label.x = 0.5, label.y = 9, vjust = 1) + 
  geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.5, fill = "grey") +
  geom_line(aes(PC1, fit), size = 1, color = "black") +
  geom_point(aes(fill = resource, shape = resource), color = "black", size = 2) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(-1,9), breaks = seq(0,9,3)) +
  scale_shape_manual(values = c(21,22,23,24,25), guide = 'none')  +
  scale_color_brewer(palette = "Set1") + 
  scale_fill_brewer(palette = "Set1") + 
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression(~{delta}^15*N~'(\u2030)'),
       fill = "Resource") + 
  annotate(geom = "text", x = 0.5, y = 9, 
           parse = TRUE, size = 4, hjust = 0, vjust = 1,
           label = as.character(expression(
             paste(F['1,384']==138.9,", ",italic(P),' < ',0.001)))) +
  annotate(geom = "text", x = 0.5, y = 9*0.9, 
           parse = TRUE, size = 4,hjust = 0, vjust = 1,
           label = as.character(expression(
             paste(R['adj']^2==0.26)))) + 
  guides(fill = guide_legend(override.aes = list(shape = c(21,22,23,24,25))))
p.bl.n.mean


# Pooled d15N range ----

plot(baseline_var$d15N_range ~ baseline_var$PC1)
# hist(baseline_var$d15N_range)
# shapiro.test(baseline_var$d15N_range)
# cor.test(baseline_var$d15N_range, baseline_var$PC1)

# Fit linear model
lm_rich <- lm(data = baseline_var, d15N_range ~ PC1)
summary(lm_rich)


# Get model fit and 95% CIs and bind to raw data
rich_pred <- lm_rich %>%
  predict(newdata = baseline_var, se.fit = TRUE) %>%
  as_tibble() %>%
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>%
  select(-se.fit)

p.bl.n.range <- 
  baseline_var %>%
  # bind_cols(rich_pred) %>%
  ggplot(aes(PC1, d15N_range)) +
  # ggpubr::stat_cor(method = "pearson", label.x = 0.5, label.y = 6, vjust = 1) + 
  # geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.5, fill = "grey") +
  # geom_line(aes(PC1, fit), size = 1, color = "black") +
  #geom_smooth(color = "black", method = "lm")  +
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(2,6), breaks = seq(2,6,1)) +
  # ggrepel::geom_text_repel(aes(label = site_id)) + 
  # ggrepel::geom_text_repel(aes(label = sample_year), hjust= 1) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression({delta}^15*N~'Range (\u2030)'),
       fill = "Stream System") 
# annotate(geom = "text", x = 0.5, y = 6*.98, parse = TRUE, size = 3, hjust = 0,
#          label = as.character(expression(paste(F['1,25']==4.676,", ",italic(P)==0.04)))) +
# annotate(geom = "text", x = 0.5, y = 6*0.9, parse = TRUE, size = 3, hjust = 0,
#          label = as.character(expression(paste(R['adj']^2==0.36))))
p.bl.n.range

# Pooled d15N CV ----

# plot(baseline_var$d15N_range ~ baseline_var$PC1)
# hist(baseline_var$d15N_range)
# shapiro.test(baseline_var$d15N_range)
# cor.test(baseline_var$d15N_range, baseline_var$PC1)

# Fit linear model
#baseline_var <- baseline_var %>% filter(d15N_CV <1)
lm_rich <- lm(data = baseline_var, d15N_CV ~ PC1)
summary(lm_rich)

# Fit GAM
gam_rich <- baseline_var %>% gam(d15N_CV ~ s(PC1, k = 3), data = ., method = "REML")
summary(gam_rich)

BIC(lm_rich, gam_rich)

# Get model fit and 95% CIs and bind to raw data
rich_pred <- gam_rich %>% 
  predict(newdata = baseline_var, se.fit = TRUE) %>% 
  as_tibble() %>% 
  mutate(lwr = fit - 2 * se.fit,
         upr = fit + 2 * se.fit) %>% 
  select(-se.fit)

summary(gam_rich)
model_label <- c("s(PC1, 1.9)",
                 "'Deviance expl.' == '42.2%'")

p.bl.n.cv <- 
  baseline_var %>%
  # bind_cols(rich_pred) %>% 
  #filter(!site_id == "MB01") %>% 
  ggplot(aes(PC1, d15N_CV)) +
  ggpubr::stat_cor(method = "pearson", label.x = 0.5, label.y = 1, vjust = 1) + 
  # geom_ribbon(aes(x = PC1, ymin = lwr, ymax = upr), alpha = 0.5, fill = "grey") +
  # geom_line(aes(PC1, fit), size = 1, color = "black") +
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks = seq(0,1,1)) +
  # ggrepel::geom_text_repel(aes(label = sample_year)) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression({delta}^15*N~'CV'),
       fill = "Stream System") 
  # annotate(geom = "text", 
  #          x = 0.5, y = c(2, 2*0.9), vjust = 1,
  #          hjust = 0, label = model_label, parse = TRUE, size = 3)
p.bl.n.cv

# Panel Plot ---- 
# Build panels
p.panel <- cowplot::plot_grid(
  p.bl.c.mean + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.bl.c.range + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.bl.c.cv + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.bl.n.mean + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.bl.n.range + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.bl.n.cv + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A", "B", "C", "D", "E", "F"),
  ncol = 3, nrow = 2, 
  align = "vh"
)
#p.panel

# extract the legend from one of the plots
legend1 <- cowplot::get_legend( p.bl.c.mean )
legend2 <- cowplot::get_legend( p.bl.c.range )

legend <- cowplot::plot_grid(legend1, legend2, nrow = 2, ncol = 1)

# Plot
p.panel.full <- cowplot::plot_grid(p.panel, legend, rel_widths = c(3, .4))
p.panel.full

# Save it
ggsave(filename = here("out", "resource_sia_var.pdf"), 
       plot = p.panel.full, device = cairo_pdf,
       units = "in", width = 10, height = 5.5)


