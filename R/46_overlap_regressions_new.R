library(tidyverse)
library(mgcv)
library(tidymv)
library(mgcViz)
library(cowplot)
library(visreg)
source("code/fx_theme_pub.R", chdir = TRUE)
theme_set(theme_Publication()) 


spp_ovlp <- spp_ovlp %>% 
  mutate(species = factor(species), 
         site_id = factor(site_id))


# overlap ------------------------------------------------------------------

hist(spp_ovlp$mean_overlap)
hist(log(spp_ovlp$mean_overlap)+2)


k <- 3


m1 <- gam(mean_overlap ~ 
            s(PC1, k = k) + 
            s(site_id, bs = "re") + 
            s(species, bs = "re"), 
          family = Gamma(link = "log"),
          method = "REML", 
          data = spp_ovlp)

k.check(m1)
summary(m1)
appraise(m1)
draw(m1, residuals = TRUE)


m2 <- gam(mean_overlap ~ 
            s(PC1, k = k) + 
            s(richness, k = k) + 
            s(site_id, bs = "re") + 
            s(species, bs = "re"), 
          family = Gamma(link = "log"),
          method = "REML", 
          data = spp_ovlp)

k.check(m2)
summary(m2)
appraise(m2)
draw(m2, residuals = TRUE)

AIC(m1, m2) %>% arrange(AIC)








# Panel Plot ------------------------------------------------------------------
# Build panels
p.panel.ovp <- cowplot::plot_grid(
  a + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  b + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A", "B"),
  ncol = 2, nrow = 1, 
  align = "vh"
)
p.panel.ovp

# # Save it
ggsave(filename = "figs/fig_4_fish_ovlp.pdf", plot = p.panel.ovp, device = cairo_pdf,
       units = "in", width = 9, height = 4)










spp_ovlp_native <- spp_ovlp %>% 
  left_join(fish_meta, by = "taxon_code") %>% 
  filter(origin == "Native")
spp_ovlp_nonnative <- spp_ovlp %>% 
  left_join(fish_meta, by = "taxon_code") %>% 
  filter(origin == "NonNative")

#  ------------------------------------------------------------------
# Natives only
#  ------------------------------------------------------------------

# overlap ------------------------------------------------------------------
# Assumptions:
plot(spp_ovlp_native$mean_overlap ~ spp_ovlp_native$PC1)
plot(spp_ovlp_native$mean_overlap ~ spp_ovlp_native$richness)

# Models ------------------------------------------------------------------
m1 <- glm(data = spp_ovlp_native, 
          mean_overlap ~ PC1, 
          family = Gamma(link = "log"))
summary(m1)

m2 <- glm(data = spp_ovlp_native, 
          mean_overlap ~ PC1 + richness,
          family = Gamma(link = "log"))
summary(m2)

m3 <- spp_ovlp_native %>% 
  gam(mean_overlap ~ s(PC1, bs = "tp", k = 3),
      data = ., method = "REML", family = Gamma(link = "log"))
summary(m3)
plot(m3, pages=1, residuals = TRUE, pch = 21)

m4 <- spp_ovlp_native %>% 
  gam(mean_overlap ~ s(PC1, bs = "tp", k = 3) + s(richness, bs = "tp", k = 3),
      data = ., method = "REML", family = Gamma(link = "log"))
summary(m4)
plot(m4, pages=1, residuals = TRUE, pch = 21)
visreg(m4, "PC1", gg=TRUE, data = spp_ovlp_native)
visreg(m4, "richness", gg=TRUE, data = spp_ovlp_native)
visreg2d(m4, "PC1", "richness", plot.type = "rgl", data = spp_ovlp_native)


# Compare models ---------------------------------------------------------------
AIC(m1, m2, m3, m4) %>% arrange(AIC)
BIC(m1, m2, m3, m4) %>% arrange(BIC)

# Build plots ------------------------------------------------------------------

c <- 
  visreg(m4, "PC1", 
         gg=TRUE, 
         data = spp_ovlp_native, 
         line = list(col = "black"), 
         points = list(col = "black", alpha = 0.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,0.5), breaks = seq(-2,0.5,0.5)) +
  labs(title = "", x = "Long. gradient (PC1)") 
c

d <- 
  visreg(m4, "richness", 
         gg=TRUE, 
         data = spp_ovlp_native, 
         line = list(col = "black"), 
         points = list(col = "black", alpha = 0.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(1.5,14.5), breaks = seq(2,14,2)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,0.5), breaks = seq(-2,0.5,0.5)) +
  labs(title = "", x = "Fish taxa richness") 
d

# Panel Plot ------------------------------------------------------------------
# Build panels
p.panel.ovp <- cowplot::plot_grid(
  c + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  d + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("C", "D"),
  ncol = 2, nrow = 1, 
  align = "vh"
)
p.panel.ovp

# # Save it
ggsave(filename = "figs/fig_4b_fish_ovlp_natived.pdf", plot = p.panel.ovp, device = cairo_pdf,
       units = "in", width = 9, height = 4)


#  ------------------------------------------------------------------
# Non-natived only
#  ------------------------------------------------------------------

# overlap ------------------------------------------------------------------
# Assumptions:
plot(spp_ovlp_nonnative$mean_overlap ~ spp_ovlp_nonnative$PC1)
plot(spp_ovlp_nonnative$mean_overlap ~ spp_ovlp_nonnative$richness)

# Models ------------------------------------------------------------------
m1 <- glm(data = spp_ovlp_nonnative, 
          mean_overlap ~ PC1, 
          family = Gamma(link = "log"))
summary(m1)

m2 <- glm(data = spp_ovlp_nonnative, 
          mean_overlap ~ PC1 + richness,
          family = Gamma(link = "log"))
summary(m2)

m3 <- spp_ovlp_nonnative %>% 
  gam(mean_overlap ~ s(PC1, bs = "tp", k = 3),
      data = ., method = "REML", family = Gamma(link = "log"))
summary(m3)
plot(m3, pages=1, residuals = TRUE, pch = 21)

m4 <- spp_ovlp_nonnative %>% 
  gam(mean_overlap ~ s(PC1, bs = "tp", k = 3) + s(richness, bs = "tp", k = 3),
      data = ., method = "REML", family = Gamma(link = "log"))
summary(m4)
plot(m4, pages=1, residuals = TRUE, pch = 21)
visreg(m4, "PC1", gg=TRUE, data = spp_ovlp_nonnative)
visreg(m4, "richness", gg=TRUE, data = spp_ovlp_nonnative)
visreg2d(m4, "PC1", "richness", plot.type = "rgl", data = spp_ovlp_nonnative)


# Compare models ---------------------------------------------------------------
AIC(m1, m2, m3, m4) %>% arrange(AIC)
BIC(m1, m2, m3, m4) %>% arrange(BIC)

# Build plots ------------------------------------------------------------------

e <- 
  visreg(m4, "PC1", 
         gg=TRUE, 
         data = spp_ovlp_nonnative, 
         line = list(col = "black"), 
         points = list(col = "black", alpha = 0.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,9), breaks = seq(0,9,1)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,0.5), breaks = seq(-2,0.5,0.5)) +
  labs(title = "", x = "Long. gradient (PC1)") 
e

f <- 
  visreg(m4, "richness", 
         gg=TRUE, 
         data = spp_ovlp_nonnative, 
         line = list(col = "black"), 
         points = list(col = "black", alpha = 0.5)) + 
  scale_x_continuous(expand=c(0,0), limits=c(1.5,14.5), breaks = seq(2,14,2)) +
  scale_y_continuous(expand=c(0,0), limits=c(-2,0.5), breaks = seq(-2,0.5,0.5)) +
  labs(title = "", x = "Fish taxa richness") 
f

# Panel Plot ------------------------------------------------------------------
# Build panels
p.panel.ovp <- cowplot::plot_grid(
  a + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  b + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  c + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  d + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  e + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  f + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A","B","C","D", "E", "F"),
  ncol = 2, nrow = 3, 
  align = "vh"
)
p.panel.ovp

# # Save it
ggsave(filename = "figs/fig_4c_fish_ovlp_panel.pdf", plot = p.panel.ovp, device = cairo_pdf,
       units = "in", width = 9, height = 12)
