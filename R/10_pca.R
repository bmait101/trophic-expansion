
# Principal components analysis - Longitudinal gradient

## TO DO 

# * custom theme is not working


## prep ---------------------------------
source(here::here("R", "00_prep.R"))
pacman::p_load(GGally, ggpubr, ggpubr, factoextra, dendextend, pvclust)


## data  -------------------------------------------------------
site_data <- read_csv(here("data", "site_covariates.csv"), 
  col_types = cols(
    # site_id = col_factor(levels = label_sites_ids),
    site_name = col_character(),
    drainage_basin = col_factor(levels = NULL),
    state = col_factor(levels = NULL),
    stream_name = col_factor(levels = label_streams),
    lat = col_double(),
    lon = col_double(),
    Sthlr_order = col_integer(),
    Elevation = col_double(),
    Drainage_area = col_double(),
    Dist_N_Platte = col_double(),
    Gradient = col_double(),
    Temp_warm = col_double(),
    Site_length = col_integer(),
    Site_width = col_double()
  )
)

# clean up
site_data_clean <- site_data %>% 
  filter(stream_name != "Horse") %>%  # scouting
  filter(site_id != "LR01") %>%  # scouting
  mutate(site_id = if_else(site_id == "LR00", "LR01", site_id)) |> 
  select(
    site_id, stream_name, 
    Elevation, Sthlr_order, Drainage_area, Dist_N_Platte, Gradient,
    Temp_warm, Site_width
  )

# check it
site_data_clean

# colinearity ---------------

site_data_clean %>% 
  select(-site_id, -stream_name) %>%
  GGally::ggcorr(geom = "circle", nbreaks = 5)

site_data_clean %>% 
  select(-site_id) %>% 
  GGally::ggpairs(axisLabels = "show")


# normality --------------------

# elevation
ggpubr::ggarrange(
  ggpubr::ggdensity(site_data_clean$Elevation, 
                    main = "Density Plot", xlab = "d15N"),
  ggpubr::ggqqplot(site_data_clean$Elevation, 
                   main = "Q-Q Plot", xlab = "d15N"), ncol = 2)
shapiro.test(site_data_clean$Elevation)
# Data normal

# Upstream drainage area
ggpubr::ggarrange(
  ggpubr::ggdensity(log(site_data_clean$Drainage_area +1), 
                    main = "Density Plot", xlab = "d15N"),
  ggpubr::ggqqplot(log(site_data_clean$Drainage_area +1), 
                   main = "Q-Q Plot", xlab = "d15N"), ncol = 2)
shapiro.test(log(site_data_clean$Drainage_area +1))
# logged data is normal

# Dist to North Platte
ggpubr::ggarrange(
  ggpubr::ggdensity(site_data_clean$Dist_N_Platte, 
                    main = "Density Plot", xlab = "d15N"),
  ggpubr::ggqqplot(site_data_clean$Dist_N_Platte, 
                   main = "Q-Q Plot", xlab = "d15N"), ncol = 2)
shapiro.test(site_data_clean$Dist_N_Platte)
# data normal

# Gradient
ggpubr::ggarrange(
  ggpubr::ggdensity(sqrt(site_data_clean$Gradient), 
                    main = "Density Plot", xlab = "d15N"),
  ggpubr::ggqqplot(sqrt(site_data_clean$Gradient), 
                   main = "Q-Q Plot", xlab = "d15N"), ncol = 2)
shapiro.test(sqrt(site_data_clean$Gradient))
# data non-normal; heavily skewed left
# log doesn't help
# sqrt() helps, but still non-normal by Shapiro (p = 0.009); best we can get

# Temperature
ggpubr::ggarrange(
  ggpubr::ggdensity(site_data_clean$Temp_warm, 
                    main = "Density Plot", xlab = "d15N"),
  ggpubr::ggqqplot(site_data_clean$Temp_warm, 
                   main = "Q-Q Plot", xlab = "d15N"), ncol = 2)
shapiro.test(site_data_clean$Temp_warm)
# normal

# Site Width
ggpubr::ggarrange(
  ggpubr::ggdensity(site_data_clean$Site_width, 
                    main = "Density Plot", xlab = "d15N"),
  ggpubr::ggqqplot(site_data_clean$Site_width, 
                   main = "Q-Q Plot", xlab = "d15N"), ncol = 2)
shapiro.test(site_data_clean$Site_width)
# normal

# log transform drainage area,
# sqrt transform gradient. 


## PCA ---------------

# set up df
pca_df <- site_data_clean %>% 
  select(-site_id, -stream_name) %>% 
  mutate(
    Drainage_area = log(Drainage_area +1),
    Gradient = sqrt(Gradient)
  ) %>% 
  as.data.frame()

# Set row names of df for plotting
rownames(pca_df) <- site_data_clean$site_id 

### Cluster analysis

hc1 <-
  pca_df %>%                               # Data
  scale() %>%                              # Scale the data
  dist(method = "euclidean") %>%           # Compute distance matrix
  hclust(method = "ward.D2")               # Hierarchical clustering

# plot
as.dendrogram(hc1) %>%
  dendextend::hang.dendrogram(hang_height = .5) %>%
  dendextend::set("branches_lwd", 1) %>%  # Branches line width
  dendextend::set("branches_k_color", viridis::viridis(3), k = 3) %>%  
  #set("labels_colors", group.colz, k = 3) %>%   
  dendextend::set("labels_cex", 1)  %>% 
  factoextra::fviz_dend(xlab = "Distance", main = "") + 
  theme_bw() 
# three main groups

# add grouping variable
site_data_clean <- site_data_clean %>% 
  mutate(
    site_group = 
      with(site_data_clean,
           ifelse(site_id %in% c(
             "LR00","LR01","MB01","MB02","SW01"),"Up",
                  ifelse(site_id %in% c(
                    "LR02","LR03","LR04","LR05","LR06","SW02","SW03","MB03"),
                    "Mid","Down"))))

### PCA

# run analysis
pca_01 <- prcomp(pca_df, center = TRUE, scale. = TRUE)
summary(pca_01)

# Extract eigenvalues/variances
get_eig(pca_01)

# Visualize eigenvalues/variances
fviz_screeplot(pca_01, addlabels = TRUE, ylim = c(0, 100)) + 
  # theme_Publication() +
  labs(x = "PCA Dimension")

# variable contributions to axes
cowplot::plot_grid(
  # Contributions of variables to PC
  fviz_contrib(pca_01, choice = "var", axes = 1, top = 10, 
               color = "black", fill = "grey") + labs(x = ""),
  # Contributions of variables to PC
  fviz_contrib(pca_01, choice = "var", axes = 2, top = 10, 
               color = "black", fill = "grey") + labs(x = ""),
  # Contributions of variables to PC
  fviz_contrib(pca_01, choice = "var", axes = 3, top = 10, 
               color = "black", fill = "grey") + labs(x = ""),
  ncol = 1, labels = c("A", "B", "C"
  )
)


# Eigen permutation test ----------------------

# Make function for permutation "test" for PCA axes
pca_eigenperm <- 
  function(data, nperm = 1000) {
    pca_out   <- prcomp(data, scale. = TRUE)
    eigenperm <- data.frame(matrix(NA, nperm, ncol(data)))
    n         <- ncol(data)
    data_i    <- data.frame(matrix(NA, nrow(data), ncol(data)))
    for (j in 1:nperm) {
      for (i in 1:n) {
        data_i[,i] <- sample(data[,i], replace = FALSE)
      }
      pca.perm <- prcomp(data_i, scale. = TRUE)
      eigenperm[j,] <- pca.perm$sdev ^2
    }
    colnames(eigenperm) <- colnames(pca_out$rotation)
    eigenperm
  }

# Run permuation PCA test; 1000 runs on data
fa_pca_perm <- pca_eigenperm(pca_df, nperm = 100)

# 
fa_pca_rand95_long <- 
  # 95th percentile of random eigenvalues
  data.frame(Random_Eigenvalues = sapply(fa_pca_perm, quantile, 0.95)) %>%  
  # add PC IDs as discrete variable
  mutate(PC = colnames(pca_01$rotation)) %>%  
  # combine rand95 with real eigenvals
  cbind(Eigenvalues = pca_01$sdev ^ 2) %>% 
  # slice out the first 5 rows, that is PCs 1-5
  slice(1:5) %>% 
  # Gather to long format
  gather(key = Variable, value = Value, -PC) %>% 
  as_tibble()

# Plot obserbed vs. random eigevalues
p.pca.perm <- 
  fa_pca_rand95_long %>% 
  ggplot(aes(PC, Value, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  labs(y = "Eigenvalue", x = "", fill = "", 
       title = "Observed vs. Random PCA Eigenvalues") +
  scale_fill_grey() + 
  theme_Publication() +
  theme(legend.position = c(0.3, 0.9))
p.pca.perm

# save plot
path <- here::here("out", "r1_pca_eigen_perm")
ggsave(glue::glue("{path}.pdf"), plot = p.pca.perm, 
       width = 6, height = 5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)


# Custom PCA biplot ------------------------------------

# Biplot of individuals and variables
p.pca.biplot <- factoextra::fviz_pca_biplot(
  pca_01, 
  habillage = site_data_clean$stream_name,
  geom.ind = c("point", "text"), 
  geom.var = c("arrow", "text"),
  alpha.var = 0.75,
  invisible = "quali",
  axes = c(1, 2), 
  col.var = "grey25", # Variables color
  repel = TRUE     # Avoid text overlapping
) + 
  labs(title = "", 
       # x = "PC1 (75.7%)", 
       # y = "PC2 (12.0%)", 
       color = "Stream System"
       ) +
  theme_Publication() +
  theme(legend.position = c(0.8, 0.9))  +
  guides(shape = 'none') + 
  scale_x_continuous(breaks = seq(-4, 4, 1), 
                     limits = c(-4,4)) + 
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = c(-2,-1,0,1,2)) + 
  scale_color_brewer(palette = "Dark2", 
                     guide = guide_legend(reverse = TRUE,
                                          override.aes = list(alpha = 1, size = 3)))

# save
path <- here::here("out", "r1_pca_biplot")
ggsave(glue::glue("{path}.pdf"), plot = p.pca.biplot, 
       width = 8, height = 5, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)


# export data ------------------------
gradient <- 
  as_tibble(
    rownames_to_column(
      as.data.frame(pca_01$x), 
      var = "site_id")
  ) %>% 
  left_join(site_data_clean, by = "site_id") %>% 
  select(site_id, stream_name, site_group, PC1, PC2) %>% 
  mutate(PC1 = (PC1 + 4), PC2 = PC2 + 2.1) %>% 
  arrange(PC1)

# check it 
gradient

# write 
gradient %>% write_csv(here("out", "r1_PCA_results.csv"))


# Export table for manuscript ------------------

# clean up
table_for_paper <-
  gradient %>%
  left_join(site_data, by = "site_id") %>%
  select(site_group, site_id, lat, lon, PC1, Elevation, Gradient, Dist_N_Platte,
         Sthlr_order, Site_width, Temp_warm, Drainage_area) %>%
  mutate(across(c(Elevation, Dist_N_Platte, Drainage_area), round, 0)) %>%
  mutate(across(c(Temp_warm), round, 1)) %>%
  mutate(across(c(PC1, Site_width), round, 2)) %>%
  mutate(across(c(Gradient), round, 3))

# check it
table_for_paper

# write
table_for_paper %>% write_csv(here("out", "r1_site_enviro_pca_table.csv"))
