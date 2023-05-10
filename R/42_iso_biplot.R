
# Isotope bi-plot across all sites (ms figure 4)

source(here("R", "40_prep-sia-data.R"))


dat <- sia_tidy %>% 
  filter(sample_year %in% c(2016, 2017)) %>% 
  select(-length_mm, -sia_sample_id) %>% 
  filter(yr_site != "2017_LR05")  # lost most of the samples so removed 
dat

levels(dat$compartment)
levels(dat$resource)

dat_base_means <- dat %>% 
  filter(compartment == "baseline") %>% 
  group_by(site_id, resource) %>% 
  summarise(mean_c = mean(d13C), 
            sd_c   = sd(d13C),
            mean_n = mean(d15N), 
            sd_n   = sd(d15N)
  ) 
dat_base_means

dat %>% 
  filter(compartment %in% c("invert")) %>% 
  mutate(taxon_code = if_else(taxon_code=="Amphipoda","Gammaride", as.character(taxon_code))) %>% 
  mutate(taxon_code = if_else(taxon_code=="Gammarus","Gammaride", as.character(taxon_code))) %>% 
  left_join(meta_invert) %>% 
  group_by(taxon_code, ffg) %>% tally %>% print(n = Inf)

dat_invert_means <- 
  dat %>% 
  mutate(taxon_code = if_else(taxon_code=="Amphipoda","Gammaride", as.character(taxon_code))) %>% 
  mutate(taxon_code = if_else(taxon_code=="Gammarus","Gammaride", as.character(taxon_code))) %>% 
  filter(compartment %in% c("invert")) %>% 
  left_join(meta_invert) %>% 
  group_by(site_id, resource, ffg) %>% 
  summarise(mean_c = mean(d13Cc), 
            sd_c   = sd(d13Cc),
            mean_n = mean(d15N), 
            sd_n   = sd(d15N)
            ) %>% 
  rename(taxon_code = ffg) %>% 
  mutate(resource = taxon_code) %>% 
  select(-taxon_code)
  
dat_invert_means

dat_cons_means <- 
  dat %>% 
  filter(compartment %in% c("fish")) %>% 
  group_by(site_id, resource, taxon_code) %>% 
  summarise(mean_c = mean(d13Cc), 
            sd_c   = sd(d13Cc),
            mean_n = mean(d15N), 
            sd_n   = sd(d15N)
  ) %>% 
  select(-taxon_code) %>% 
  bind_rows(dat_invert_means) %>% 
  bind_rows(dat_base_means) %>% 
  mutate(resource = as.factor(resource)) %>% 
  mutate(resource = fct_recode(resource, 
                               Fish      = "fish", 
                               Collector = "collector", 
                               Filterer  = "filterer",
                               Grazer    = "grazer",
                               Predator  = "predator",
                               Shredder  = "shredder",
                               `Epilithic algae`   = "biofilm",
                               Detritus  = "detritus",
                               `Aquatic phototroph` = "photo",
                               `Terrestrial C3`     = "terrC3",
                               `Terrestrial C4`     = "terrC4"
                               )) %>% 
  mutate(resource = fct_relevel(resource, 
                                levels = c("Fish",
                                           "Predator",
                                           "Collector",
                                           "Grazer",
                                           "Shredder",
                                           "Filterer",
                                           "Epilithic algae",
                                           "Detritus",
                                           "Aquatic phototroph",
                                           "Terrestrial C3",
                                           "Terrestrial C4"
                                           ))
         ) %>% 
  mutate(site_id = fct_relevel(site_id, levels = c(
    "MB01","SW01","LR01","MB02","SW02","LR02","LR03","MB03",
    "LR04","LR05","SW03","LR06","MB04","LR07","SW04","LR08"))
  ) %>% 
  rename(Compartment = resource) 
dat_cons_means

max(dat_cons_means$mean_n)

dat_cons_means %>% 
  #filter(! Compartment %in% c("Predator","Collector","Grazer","Shredder","Filterer")) %>% 
  ggplot(aes(x = mean_c, y = mean_n, shape = Compartment, fill = Compartment)) +
  geom_errorbar(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n)) +
  geom_errorbarh(aes(xmin = mean_c - sd_c, xmax = mean_c + sd_c)) +
  geom_point(size = 1.5) + 
  scale_shape_manual(values = c(21,
                                21,23,22,24,25,
                                23,21,22,24,25
                                )) + 
  scale_fill_manual(values = c("white",
                                "black","black","black","black","black",
                                "olivedrab2","tan4","darkgreen","brown","brown"
                               )) + 
  scale_x_continuous(limits = c(-35, -15), breaks = seq(-35, -15, 5)) + 
  scale_y_continuous(limits = c(-2, 18), breaks = seq(0, 18, 6)) + 
  labs(x = expression(~{delta}^13*C~'(\u2030)'), 
       y = expression(~{delta}^15*N~'(\u2030)')) + 
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1)) + 
  lemon::facet_rep_wrap(vars(factor(site_id)), ncol = 4) + 
  theme(axis.line = element_blank())


# Save it
path <- here::here("out", "r1_sia_biplot_v2")
ggsave(glue::glue("{path}.pdf"), plot = last_plot(), 
       width = 8.5, height = 7.5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)

