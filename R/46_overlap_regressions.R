library(tidyverse)
library(mgcv)
library(gratia)
source("code/fx_theme_pub.R", chdir = TRUE)
theme_set(theme_Publication()) 


spp_ovlp4 <- spp_ovlp3 %>% 
  mutate(species = factor(species), 
         site_id = factor(site_id))


# overlap ------------------------------------------------------------------

hist(spp_ovlp4$mean_overlap)
hist(log(spp_ovlp4$mean_overlap)+2)
hist(spp_ovlp4$overlap)
hist(log(spp_ovlp4$overlap))


tmp <- lm(overlap ~ PC1 + species, data = spp_ovlp4)
summary(tmp)


k <- 3

m1 <- gam(log(mean_overlap) ~ 
            s(PC1, k = k) + 
            s(site_id, bs = "re") + 
            s(species, bs = "re"), 
          # family = Gamma(link = "log"),
          method = "REML", 
          data = spp_ovlp4)

k.check(m1)
summary(m1)
appraise(m1)
draw(m1, residuals = TRUE)


m2 <- gam(log(mean_overlap) ~ 
            s(PC1, k = k) + 
            s(richness, k = k) + 
            s(site_id, bs = "re") + 
            s(species, bs = "re"), 
          # family = Gamma(link = "log"),
          method = "REML", 
          data = spp_ovlp4)

k.check(m2)
summary(m2)
appraise(m2)
draw(m2, residuals = TRUE)

AIC(m1, m2) %>% arrange(AIC)




# Predict from fitted model
new_data <- with(spp_ovlp4, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  richness = 1, 
  site_id = levels(spp_ovlp4$site_id)[[1]],
  species = levels(spp_ovlp4$species)[[1]]
  
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m2, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m2))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions
model_label <- c("s(PC1, 1.54)")
p.pc1 <- new_data %>% 
  ggplot(aes(x=PC1,y=exp(fit))) + 
  geom_ribbon(aes(ymin=exp(lower), ymax=exp(upper), x = PC1), alpha=.5, fill="grey")+ 
  geom_line() + 
  geom_point(data=spp_ovlp4, aes(x=PC1,y=mean_overlap, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(expand=c(0,0), limits=c(0,9), breaks = seq(0,9,1)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,1.1), breaks = seq(0,1,.25)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = "Niche overlap", 
       fill = "Stream system") 



# Predict from fitted model
new_data <- with(spp_ovlp4, tibble(
  richness = seq(min(richness), max(richness), length.out = 200), 
  PC1 = 1, 
  site_id = levels(spp_ovlp4$site_id)[[1]],
  species = levels(spp_ovlp4$species)[[1]]
  
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m2, 
      newdata = new_data, 
      terms="s(richness)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m2))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions
model_label <- c("s(richness, 1.92)")
p.rich <- new_data %>% 
  ggplot(aes(x=richness,y=exp(fit))) + 
  geom_ribbon(aes(ymin=exp(lower), ymax=exp(upper), x = richness), alpha=.5, fill="grey")+ 
  geom_line() +
  geom_point(data=spp_ovlp4, aes(x=richness, y=mean_overlap, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  # scale_x_continuous(expand=c(0,0), limits=c(0,9), breaks = seq(0,9,1)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,1.2), breaks = seq(0,1,.25)) +
  labs(title = "", x = "Fish species richness", 
       y = "Niche overlap", 
       fill = "Stream system")





p.panel.ovp <- cowplot::plot_grid(
  p.pc1 + theme(legend.position = c(.7, .8), plot.margin=unit(c(1,1,1,1),"mm")), 
  p.rich + theme(legend.position = 'none', plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A","B"),
  ncol = 2, nrow = 1, 
  align = "vh"
)
p.panel.ovp

# # Save it
ggsave(here("figs1", "niche_overlap.pdf"), plot = p.panel.ovp,
       units = "in", width = 9, height = 4)


# 3d

m3 <- gam(log(mean_overlap) ~ 
            s(PC1, k = k) + 
            s(richness, k = k) + 
            ti(PC1, richness) + 
            s(site_id, bs = "re") + 
            s(species, bs = "re"), 
          # family = Gamma(link = "log"),
          method = "REML", 
          data = spp_ovlp4)
m4 <- gam(mean_overlap ~ 
            s(PC1, k = k) + 
            s(richness, k = k) + 
            ti(PC1, richness) + 
            s(site_id, bs = "re") + 
            s(species, bs = "re"), 
          family = Gamma(link = "log"),
          method = "REML", 
          data = spp_ovlp4)

AIC(m3, m4)
k.check(m3)
summary(m3)
appraise(m3)
draw(m3, residuals = TRUE)

AIC(m1, m2, m3) %>% arrange(AIC)

model_2_p <- predict_gam(m3, exclude_terms = "s(site_id,species)")
model_2_p <- predict_gam(m4)
model_2_p

model_2_p %>%
  ggplot(aes(PC1, richness, z = fit)) +
  geom_raster(aes(fill = fit)) +
  geom_contour(colour = "white") + 
  viridis::scale_fill_viridis()
