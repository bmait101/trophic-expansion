library(tidyverse)
library(here)
library(vegan)
library(mgcv)
library(cowplot)
source(here("code", "fx_theme_pub.R"))
theme_set(theme_Publication())

# Data ------------------------------------------------------
bugs <- read_csv(here("data-raw", "sample_log_bugs.csv")) 
bugs_meta <- read_csv(here("data-raw", "metadata_bugs.csv")) %>% 
  select(taxon_code, ffg, ffg_shd, ffg_col, ffg_grz, ffg_prd)
gradient <- read_csv(here("data-derived", "data_PCA_results.csv"))

data_sites <- read_csv("data-derived/data_PCA_results.csv") %>% 
  filter(site_id != "LR01") 

# Wrangle
bugs <- bugs %>% 
  filter(stream_name != "Mono") %>%  
  filter(sample_year == "2016") %>%  
  filter(site_id != "LR01")%>%        
  filter(! site_id %in% c("SW-DG","SW-SRR","SW-WR")) %>% 
  filter(gear_type=="dnet") %>% 
  mutate(taxon_code = if_else(taxon_code == "Glossiphoniidae", "Hirudinea", taxon_code)) %>% 
  mutate(taxon_code = if_else(taxon_code == "Gammarus", "Gammaride", taxon_code)) %>% 
  mutate(taxon_code = if_else(taxon_code == "Amphipoda", "Gammaride", taxon_code)) %>% 
  select(sample_year, stream_name, site_id, sample_hitch, gear_type, taxon_code, count_individuals) %>% 
  group_by(site_id, taxon_code) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(taxon_code) %>% 
  drop_na(taxon_code) %>% 
  filter(taxon_code != "UNKN") %>% 
  #filter(n() > 1) %>% 
  left_join(bugs_meta, by = "taxon_code") %>% 
  filter(! ffg %in% c("terrestrial-collector","terrestrial-herbivore","terrestrial-predator"))
bugs

# Counts 
# bugs %>%  
#   group_by(taxon_code) %>% 
#   tally() %>% print(n = Inf) 
# 
# bugs %>%
#   group_by(ffg, taxon_code) %>% 
#   tally() %>% print(n = Inf) 
# 
# bugs %>%
#   group_by(site_id) %>% 
#   tally() %>% print(n = Inf) 

# Prep --------------------------------

# Cal richness and abundance and combine
site_richness <-
    enframe(vegan::specnumber(data_field_wide_matrix)) %>%
    rename(richness = value) %>% 
    select(-name)

bug_summary <- bugs %>%
  group_by(site_id) %>%
  summarise(.groups = "drop",
            abund = sum(n)) %>% 
  bind_cols(site_richness) %>% 
  left_join(gradient, by = "site_id")
bug_summary


# Make a sites by species matirx (wide data)
data_bug_wide <- bugs %>% 
  group_by(site_id, taxon_code) %>%
  summarise(.groups = "drop",
            abund = sum(n)) %>% 
  # spread obs to make wide table
  spread(taxon_code, abund, fill = 0) 
data_bug_wide

# Set the abundance matrix to its own df
data_field_wide_matrix <- as.matrix(data_bug_wide[2:ncol(data_bug_wide)])

# All values >c0 are retuned as 1
data_field_wide_matrix[data_field_wide_matrix>0] <- 1

# Assiagn to df and name with species names
sites_present <- 
  enframe(colSums(data_field_wide_matrix)) %>% 
  rename(taxon_code = name, occupancy = value)
#sites_present$taxon_code <- rownames(sites_present)

# Check it
sites_present <- sites_present %>% 
  left_join(bugs_meta, by = "taxon_code") %>% 
  select(ffg, taxon_code, occupancy) %>% 
  arrange(desc(occupancy)) %>% 
  arrange(ffg, taxon_code) %>% 
  rename(FFG = ffg, Species = taxon_code, Occupancy=occupancy) %>% 
  write_csv(here("results", "bug_species_occupancy.csv"))

# Rarefaction --------------------------------

# Species richness increases with sample size, and differences in richness
# actually may be caused by differences in sample size.
# Thus, we rarefy species richness to the same number of individuals. 
# But differences seem to be small vs. non-rarefied estimates

# Vector of site_ids
sitesLabels <- data_bug_wide$site_id  
# Set the abundance matrix to its own df
BCI <- as.matrix(data_bug_wide[2:ncol(data_bug_wide)])
rownames(BCI) <- data_bug_wide$site_id 
BCI

rarecurve(BCI, col = "black")
rare_rich <- rarefy(BCI, min(rowSums(BCI)))
rare_rich_tbl <- rare_rich %>% 
  enframe() %>%
  rename(site_id = name, rare_rich = value)

bug_summary <- bug_summary %>% 
  left_join(rare_rich_tbl, by = "site_id")
bug_summary

lm_rich <- lm(data = bug_summary, rare_rich ~ PC1)
summary(lm_rich)

p.richness <- 
  bug_summary %>%
  ggplot(aes(x = PC1, y = rare_rich)) +    
  ggpubr::stat_cor(method = "pearson", label.x = 0.5) + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  # scale_y_log10(limits = c(1,50), breaks=c(1,10,50)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", y = "Invertebrte taxa richness", 
       fill = "Stream System")
p.richness


# FFGs ----------------------------------------------------
ffgs <- bugs %>% 
  group_by(site_id) %>% 
  summarise(mean_shd = mean(ffg_shd, na.rm = TRUE), 
            mean_col = mean(ffg_col, na.rm = TRUE), 
            mean_grz = mean(ffg_grz, na.rm = TRUE), 
            mean_prd = mean(ffg_prd, na.rm = TRUE)
            ) %>% 
  left_join(gradient, by = "site_id") 
ffgs



# Fit GAM
gam_KUD95 <- ffgs %>% 
  gam(mean_shd ~ s(PC1, bs = "tp", k = 3), 
      data = ., method = "REML", family = gaussian)
summary(gam_KUD95)
plot(gam_KUD95, pages=1, residuals = TRUE, pch = 21)

fit_gam <- gam_KUD95 %>%  
  predict(newdata = ffgs, type = "link", se.fit = TRUE) %>% 
  as_tibble() %>% 
  rename(fit_gam = fit) %>% 
  mutate(lwr_gam = fit_gam - 2 * se.fit,
         upr_gam = fit_gam + 2 * se.fit) %>% 
  select(-se.fit)
summary(gam_KUD95) # gam style summary of fitted model

model_label <- c("s(PC1, 1.64)",
                 "'Deviance expl.' == '69.5%'")

p.shd <- 
  ffgs %>% 
  bind_cols(fit_gam) %>% 
  ggplot(aes(x = PC1, y = mean_shd)) +    
  geom_ribbon(aes(x = PC1, ymin = lwr_gam, ymax = upr_gam), alpha = 0.5, fill = "grey") +
  geom_line(aes(PC1, fit_gam), size = 1, color = "black") +
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  #scale_y_continuous(breaks=seq(1.8,2.1,0.1)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  annotate(geom = "text", x = 3, y = c(3, 3*0.9), 
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 4) + 
  labs(title = "", x = "Long. gradient (PC1)", y = "Shredder FFG trait affinity", 
       fill = "Stream System")
p.shd

# Fit linear model
lm_KUD95 <- ffgs %>% lm(mean_col ~ PC1, data = .)
summary(lm_KUD95)
p.col <- 
  ffgs %>% 
  ggplot(aes(x = PC1, y = mean_col)) +    
  geom_smooth(method = "lm", color = "black") + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(breaks=seq(1.8,2.1,0.1)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", y = "Collector FFG trait affinity", 
       fill = "Stream System") +
  annotate(geom = "text", x = 0.5, y = 2.1, parse = TRUE, size = 4, hjust = 0, vjust = 1,
           label = as.character(expression(paste(F['1,14']==5.6,", ",italic(P)==0.3)))) +
  annotate(geom = "text", x = 0.5, y = 2.1*0.98, parse = TRUE, size = 4,hjust = 0, vjust = 1,
           label = as.character(expression(paste(R['adj']^2==0.23))))
p.col

# Fit linear model
lm_KUD95 <- ffgs %>% lm(mean_grz ~ PC1, data = .)
summary(lm_KUD95)
p.grz <- 
  ffgs %>% 
  ggplot(aes(x = PC1, y = mean_grz)) +    
  geom_smooth(method = "lm", color = "black") + 
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  # scale_y_continuous(breaks=seq(1.8,2.1,0.1)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", y = "Grazer FFG trait affinity", 
       fill = "Stream System") +
  annotate(geom = "text", x = 4, y = 1.9, parse = TRUE, size = 4, hjust = 0, vjust = 1,
           label = as.character(expression(paste(F['1,14']==4.7,", ",italic(P)==0.05)))) +
  annotate(geom = "text", x = 4, y = 1.9*0.96, parse = TRUE, size = 4,hjust = 0, vjust = 1,
           label = as.character(expression(paste(R['adj']^2==0.20))))
p.grz

p.prd <- 
  ffgs %>% 
  ggplot(aes(x = PC1, y = mean_prd)) +    
  geom_point(aes(fill = stream_name), color = "black", size = 3, shape = 21) +
  ggpubr::stat_cor(method = "pearson", label.x = 0.5) + 
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  scale_y_continuous(limits = c(2.5,2.9), breaks=seq(2.5,2.9,0.1)) +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  labs(title = "", x = "Long. gradient (PC1)", y = "Predator FFG trait affinity", 
       fill = "Stream System")
p.prd



middle_row <- plot_grid(p.shd + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
                        p.col + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
                        labels = c('B', 'C'))

bottom_row <- plot_grid(p.grz + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
                        p.prd + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
                        labels = c('D', 'E'))

panel <- plot_grid(p.richness, middle_row, bottom_row, labels = c('A', '', ''), ncol = 1)

ggsave(filename = here("figs1", "invert_ffgs.pdf"), 
       plot = panel, device = cairo_pdf, 
       units = "in", width = 8, height = 10)




# NMDS ---------------------------------------------------------------------

# Set dataframs
mat_fish_occup <- as.matrix((BCI > 0) + 0)
data_sites <- data_sites %>% 
  arrange(site_id) %>% 
  select(-site_id, -PC1, -PC2) %>% 
  as.data.frame()
rownames(data_sites) <- rownames(mat_fish_occup)

# run nmds
nmds_fish_occup <- metaMDS(mat_fish_occup,
                           distance = "bray", trace = FALSE, 
                           k = 3,maxit = 999, trymax = 500,wascores = TRUE)

# Check stress
nmds_fish_occup$stress

# plot the results
plot(nmds_fish_occup, type = "t")

# Shepards test/goodness of fit
stressplot(nmds_fish_occup) # Produces a Shepards diagram

# Plotting NMDS
fish.envfit <- envfit(nmds_fish_occup, data_sites, permutations = 999) # this fits environmental vectors
fish.spp.fit <- envfit(nmds_fish_occup, mat_fish_occup, permutations = 999) # this fits species vectors

site.scrs <- as.data.frame(scores(nmds_fish_occup, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Stream = data_sites$stream_name) #add grouping variable "Management" to dataframe
site.scrs <- cbind(site.scrs, Group = data_sites$site_group) #add grouping variable of cluster grouping to dataframe
spp.scrs <- as.data.frame(scores(fish.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = fish.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

# PLot
nmds.plot.occup <- site.scrs %>% 
  ggplot(aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, 
                 fill = factor(site.scrs$Stream)), size = 2, shape = 21)+
  #coord_fixed() +
  labs(fill = "Stream") + 
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)  +
  scale_color_brewer(palette = "Dark2") + 
  scale_x_continuous(limits = c(-1,1))+ 
  scale_y_continuous(limits = c(-1,0.75))+ 
  scale_fill_brewer(palette = "Dark2") + 
  annotate(geom = "text", x = -1, y = 0.75, size = 4, hjust = 0, label = "2-D stress = 0.11") +
  guides(fill = guide_legend(override.aes=list(shape=21)))
nmds.plot.occup

ggsave(filename = here("figs1","invert_nmds.pdf"), 
       plot = nmds.plot.occup, device = cairo_pdf,  
       units = "in", width = 7, height = 4)

# NMDS data extracition
# Extract nmds site coordinates for axis 1 and 2

min(nmds_fish_occup$points)
mds_fit_occup <- 
  as_tibble(rownames_to_column(as.data.frame(nmds_fish_occup$points), var = "site_id")) %>% 
  mutate(MDS1 = MDS1 + 2, responce = "occupancy") %>% 
  select(site_id, MDS1, responce) %>% 
  arrange(MDS1)
mds_fit_occup



# plot 
nmds_cor <-  mds_fit_occup %>% 
  left_join(gradient, by = "site_id") %>% 
  ggplot(aes(PC1, MDS1, fill = responce)) + 
  geom_point(size = 3, shape = 21)   +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  labs(title = "", x = "Long. gradient (PC1)", y = "Fish NMDS Axis 1", 
       fill = "Fish response") + 
  theme(legend.position = c(0.1,0.9))
nmds_cor

mds_fit_occup <- 
  mds_fit_occup %>% 
  left_join(gradient, by = "site_id")

cor.test(mds_fit_occup$PC1, mds_fit_occup$MDS1, method = "pearson")

