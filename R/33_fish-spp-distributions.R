
# Plot fish species occurance and rel. abundance along gradient

## Prep --------------------

# libraries
library(tidyverse)

# some labels
fam_levels <- c("Catostomidae","Centrarchidae","Clupeidae",
                "Cyprinidae","Fundulidae","Gasterosteridae",
                "Ictaluridae","Percidae","Salmonidae","Sciaenidae")
origin_levels <- c("Native","NonNative")

## Data ------------------

fish_spp_abundnace <- read_csv(here("out", "fish_spp_counts.csv"))

den <- fish_spp_abundnace %>% 
  filter(stream_name != "Horse") %>% 
  group_by(site_id, taxon_code) %>%
  summarise(density = mean(density))

data <- data_field %>% 
  group_by(site_id, taxon_code) %>%
  summarise(count = sum(count)) %>% 
  left_join(gradient, by = "site_id") %>% 
  left_join(fish, by = "taxon_code") %>% 
  mutate(site_group = factor(site_group, levels = c("Upstream", "Midstream", "Downstream")), 
         family = factor(family, levels = fam_levels), 
         origin = factor(origin, levels = origin_levels)) %>% 
  left_join(den, by = c("site_id", "taxon_code"))
data

# Plotting ----------------------------------------------------------------

# Plot
spp_dist <- 
  data %>% 
  group_by(common_name) %>% 
  ggplot(aes(x = fct_reorder(common_name, PC1, .desc = TRUE), y = PC1)) + 
  coord_flip(clip = "off") +
  geom_point(aes(fill = trophic_group, shape = origin, size = density), color = "black") +
  scale_shape_manual(values = c(21, 24), name = "Origin") +
  labs(y = "Long. Gradient (PC1)",
       x = "Fish species", size = "Density",
       fill = "Trophic Group") +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  theme(legend.key.size = unit(0.5, "cm"), 
        legend.title = element_text(face = "plain"), 
        panel.grid.major = element_line(color = "grey", size = 0.25))
spp_dist


# ggsave(filename = "figs/spp-distribution.png", plot = last_plot(), 
#        units = "in", width = 11, height = 6)


#===============================================================================
# Panel Plot with other data
#===============================================================================


top_row <- plot_grid(
  p.richness + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.desnity + theme(plot.margin=unit(c(1,1,1,1),"mm"),
                    legend.position = "none"), 
  p.shannon + theme(plot.margin=unit(c(1,1,1,1),"mm"),
                    legend.position = "none"),
  labels = c("A", "B", "C"), ncol = 3, nrow = 1,
  align = "vh"
  )

panel <- plot_grid(top_row, spp_dist, labels = c('', 'D'), ncol = 1, rel_heights = c(0.7, 1))
panel

# ggsave(filename = here("out", "fish_comm_panel.pdf"), 
#        plot = panel, device = cairo_pdf, 
#        units = "in", width = 11, height = 8)


