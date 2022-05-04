
# Fish community analyses

## Prep 

# Libraries
library(tidyverse)
library(here)
library(vegan)
library(funrar)

source(here("code", "fx_theme_pub.R"))
theme_set(theme_bw())

## Data

# Electrofishing data
data_field <- read_csv(here("data", "data_field_fish.csv")) %>% 
  select(-life_stage, -count_euthanized, -length_mm, -weight_g, -field_notes) %>% 
  filter(sample_year != "2015") %>%   # remove 2015 pilot data
  filter(sample_year != "2017a") %>%  # remove April 2017 sample event
  filter(site_id != "LR01") %>%       # remove LR01
  filter(stream_name != "Horse")      # remove Horse Creek
data_field

data_sites <- read_csv(here("data", "data_field_sites.csv")) %>% 
  filter(stream_name != "Horse") %>% 
  mutate(site_area = Site_length * Site_width) %>%  # calculate site area
  select(site_id, site_area)  # clean up tbl
data_sites

gradient <- read_csv("data/data_PCA_results.csv")
fish <- read_csv("data/metadata_fishes.csv")


# Fish community metrics ------------------------------------------------------

# Calculate total fish abundance and density
fish_community_summary <- data_field %>% 
  group_by(sample_year, stream_name, site_id) %>%
  summarise(.groups = "drop",
            abund = sum(count)) %>% 
  left_join(data_sites, by = "site_id") %>%
  #ungroup() %>% 
  mutate(density = abund / site_area) %>% 
  select(-site_area) %>% 
  arrange(sample_year, stream_name, site_id)
fish_community_summary

# Next, calculate site diversity metrics

D.Hill <- function(x,a) {  
  tot <- apply(x,1,sum)
  x <- sweep(x,1,tot,"/")
  x <- x^a
  D <- apply(x,1,sum,na.rm=TRUE)
  D <- D^(1/(1-a))
}
D.mg <- function(x){  
  r <- vegan::specnumber(x)
  a <- rowSums(x)
  return((r-1)/log(a))
}

# Make a sites by species matirx (wide data)
data_field_wide <- data_field %>% 
  # group by sample event and species and sum total number fish by site
  group_by(sample_year, stream_name, site_id, taxon_code) %>%
  summarise(.groups = "drop",
            abund = sum(count)) %>% 
  # spread obs to make wide table
  spread(taxon_code, abund, fill = 0) 
data_field_wide

# Set the abundance matrix to its own df
data_field_wide_matrix <- as.matrix(data_field_wide[4:ncol(data_field_wide)])

# Make calculations
(site_richness <-
    enframe(vegan::specnumber(data_field_wide_matrix)) %>%
    rename(richness = value) %>% 
    select(-name) )
(site_shannon <-
    enframe(vegan::diversity(data_field_wide_matrix, index = "shannon")) %>%
    rename(shannon = value) %>% 
    select(-name) )
(site_simpson <-
    enframe(vegan::diversity(data_field_wide_matrix, index = "simpson")) %>%
    rename(simpson = value) %>% 
    select(-name) )
(site_Dhill <-
    enframe(D.Hill(data_field_wide_matrix,0.4)) %>%
    rename(Dhill = value) %>% 
    select(-name) )
(site_mg <-
    enframe(D.mg(data_field_wide_matrix)) %>%
    rename(mg = value) %>% 
    select(-name) )

# Add to summary table
fish_community_summary <- fish_community_summary %>%
  bind_cols(site_richness) %>% 
  bind_cols(site_shannon) %>% 
  bind_cols(site_simpson) %>% 
  bind_cols(site_Dhill) %>% 
  bind_cols(site_mg) 

fish_community_summary %>%
  write_csv(here("out", "fish_structure_summary.csv"))


# Species counts ---------------------------------------------------------

fish_spp_abundnace <- data_field %>% 
  # group and sum total number of fish by site
  group_by(sample_year, stream_name, site_id, taxon_code) %>%
  summarise(.groups = "drop",
            count = sum(count)) %>% 
  # Now, join site area to tbl
  left_join(data_sites, by = "site_id") %>%
  
  # Calculate total fish density  by species
  mutate(density = count / site_area) %>% 
  select(-site_area) %>%  
  arrange(sample_year, stream_name, site_id) 

fish_spp_abundnace %>% 
  write_csv(here("out", "fish_spp_counts.csv"))

# Species relative abundance -----------------------------------------------------------

# Calculate relative abund matrix 
fish_rel_abund <- data_field_wide_matrix %>% 
  funrar::make_relative() %>% 
  round(digits = 3) %>% 
  as_tibble() %>% 
  mutate(site_id = data_field_wide$site_id, 
         stream_name = data_field_wide$stream_name, 
         sample_year = data_field_wide$sample_year)


fish_rel_abund %>% 
  write_csv(here("out", "fish_spp_relative.csv"))




# Vegan NMDS -------------------------


# Make wide data, but this time without year (so means across all years)
data_field_wide_site <- data_field %>%
  group_by(site_id, taxon_code) %>%
  summarise(fish_abundance = mean(count)) %>%
  spread(taxon_code, fish_abundance, fill = 0)

# Vector of site_ids
sitesLabels <- data_field_wide_site$site_id  

# Make data frame for analysis
BCI <- as.matrix(round(data_field_wide_site[,4:ncol(data_field_wide_site)]))
rownames(BCI) <- data_field_wide_site$site_id 

# Set dataframs
mat_fish_abund <- BCI
mat_fish_occup <- as.matrix((BCI > 0) + 0)

data_sites <- read_csv(here("out", "data_PCA_results.csv")) %>% 
  filter(site_id != "LR01") %>% 
  arrange(site_id) %>% 
  select(-site_id, -PC1, -PC2) %>% 
  as.data.frame()
rownames(data_sites) <- rownames(mat_fish_abund)


# Run NMDS
nmds_fish_abund <- metaMDS(mat_fish_abund,
                           distance = "bray", trace = FALSE, 
                           k = 3,maxit = 999, trymax = 500,wascores = TRUE)

nmds_fish_occup <- metaMDS(mat_fish_occup,
                           distance = "bray", trace = FALSE, 
                           k = 3,maxit = 999, trymax = 500,wascores = TRUE)

# Check stress
nmds_fish_abund$stress
nmds_fish_occup$stress

# plot the results
plot(nmds_fish_abund, type = "t")
plot(nmds_fish_occup, type = "t")

# Shepards test/goodness of fit
stressplot(nmds_fish_abund) # Produces a Shepards diagram
stressplot(nmds_fish_occup) # Produces a Shepards diagram

# Plotting NMDS
# Vectors
fish.envfit <- envfit(nmds_fish_abund, data_sites, permutations = 999) # this fits environmental vectors
fish.spp.fit <- envfit(nmds_fish_abund, mat_fish_abund, permutations = 999) # this fits species vectors

site.scrs <- as.data.frame(scores(nmds_fish_abund, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Stream = data_sites$stream_name) #add grouping variable "Management" to dataframe
site.scrs <- cbind(site.scrs, Group = data_sites$site_group) #add grouping variable of cluster grouping to dataframe
spp.scrs <- as.data.frame(scores(fish.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = fish.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

# PLot
nmds.plot.abund <- site.scrs %>% 
  ggplot(aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, 
                 fill = factor(site.scrs$Stream)), size = 2, shape = 21)+
  #coord_fixed() +
  labs(fill = "Stream", shape = "Group") + 
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)  +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  annotate(geom = "text", x = -1.25, y = 1, size = 4, hjust = 0, label = "2-D stress = 0.05") +
  guides(fill = guide_legend(override.aes=list(shape=21)))
nmds.plot.abund

# Vectors
fish.envfit <- envfit(nmds_fish_occup, data_sites, permutations = 999) # this fits environmental vectors
fish.spp.fit <- envfit(nmds_fish_occup, mat_fish_abund, permutations = 999) # this fits species vectors

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
  labs(fill = "Stream", shape = "Group") + 
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), 
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) +
  ggrepel::geom_text_repel(data = sig.spp.scrs, 
                           aes(x=NMDS1, y=NMDS2, label = Species), 
                           cex = 3, direction = "both", segment.size = 0.25, 
                           max.overlaps = 50)  +
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  annotate(geom = "text", x = -1, y = 1, size = 4, hjust = 0, label = "2-D stress = 0.05") +
  guides(fill = guide_legend(override.aes=list(shape=21)))
nmds.plot.occup

# arrange the three plots in a single row
prow <- cowplot::plot_grid(
  nmds.plot.abund + theme(legend.position="none"),
  nmds.plot.occup + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
prow

# extract the legend from one of the plots
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  nmds.plot.abund + theme(legend.box.margin = margin(0, 0, 0, 0))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
nmds_panel <- cowplot::plot_grid(prow, legend, rel_widths = c(3, 0.5))
nmds_panel

ggsave(filename = here("out", "fish_nmds.pdf"), plot = nmds_panel, 
       device = cairo_pdf,  units = "in", width = 10, height = 3.5)



# NMDS data extraction
# Extract nmds site coordinates for axis 1 and 2
min(nmds_fish_abund$points)
mds_fit_abund <- 
  as_tibble(rownames_to_column(as.data.frame(nmds_fish_abund$points), var = "site_id")) %>% 
  mutate(MDS1 = MDS1 + 2, responce = "abundance") %>% 
  select(site_id, MDS1, responce) %>% 
  arrange(MDS1)
mds_fit_abund

min(nmds_fish_occup$points)
mds_fit_occup <- 
  as_tibble(rownames_to_column(as.data.frame(nmds_fish_occup$points), var = "site_id")) %>% 
  mutate(MDS1 = MDS1 + 2, responce = "occupancy") %>% 
  select(site_id, MDS1, responce) %>% 
  arrange(MDS1)
mds_fit_occup



# plot 
nmds_cor <- 
  mds_fit_abund %>% 
  bind_rows(mds_fit_occup) %>% 
  left_join(gradient, by = "site_id") %>% 
  ggplot(aes(PC1, MDS1, fill = responce)) + 
  geom_point(size = 3, shape = 21)   +
  scale_x_continuous(expand=c(0,0), limits=c(-0.25,9), breaks = seq(0,9,1)) +
  labs(title = "", x = "Long. gradient (PC1)", y = "Fish NMDS Axis 1", 
       fill = "Fish response") + 
  theme(legend.position = c(0.1,0.9))
nmds_cor

# ggsave(filename = "figs/fish_nmds_corr.pdf", plot = nmds_cor, 
#        device = cairo_pdf,  units = "in", width = 7, height = 5)

mds_fit_abund %>% 
  bind_rows(mds_fit_occup) %>% 
  left_join(gradient, by = "site_id") %>% 
  filter(responce=="abundance") %>% 
  select(PC1, MDS1)

data <- 
  mds_fit_abund %>% 
  left_join(mds_fit_occup, by = "site_id") %>% 
  left_join(gradient, by = "site_id")
data 

cor.test(data$PC1, data$MDS1.x, method = "pearson")
cor.test(data$PC1, data$MDS1.y, method = "pearson")
cor.test(data$MDS1.x, data$MDS1.y, method = "pearson")

# Calculating Pearson correlation coefficients for dimensions of NMDS plot
species_cor <-
  as.data.frame(
    cor(mat_fish_abund, nmds_fish_abund$points,
        use = "complete.obs",
        method = "pearson"))

species_cor$taxon <- rownames(species_cor)
species_cor %>% 
  ggplot(aes(x = fct_reorder(taxon, MDS1, .desc = TRUE) , y = MDS1)) +
  geom_point() + 
  ylim(-1,1)
