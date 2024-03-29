
# Prepare isotope data for niche space analyses

# Set up ------------------------------------------------------------------

source(here::here("R", "00_prep.R"))


## Data --------------------

sia_raw <- read_csv(here("data", "iso_compiled.csv"))

meta_fish <- read_csv(here("data", "metadata_fishes.csv"))

meta_invert <- 
  read_csv(here("data", "metadata_bugs.csv")) |> 
  select(taxon_code, ffg) 

gradient <- read_csv(here("out", "r1_PCA_results.csv"))


# NAs ------------------

map(sia_raw, ~sum(is.na(.)))
NAs <- sia_raw |> filter(is.na(sample_year))

NAs |> ggplot(aes(d13C, d15N)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = sia_sample_id), size=2, max.overlaps=15)

### THESE ARE LIVER DATA: about 65 records to remove

# Tidy Data -------------------

sia_tidy <- sia_raw |> 
  # fix SW site labels from data entry errors (TNC and GR are same site)
  mutate(site_id = if_else(site_id=="SW03-TNC","SW03",
                           if_else(site_id=="SW03-GR","SW03",site_id))) |>
  # Drop liver SI data - only keep muscle tissue data
  drop_na(taxon_code) |> 
  # Remove Horse Creek and other preliminary survey sites
  filter(! stream_name %in% c("Horse")) |> 
  filter(! site_id %in% c("LR01","LR067","MONO","Mono")) |> 
  # rename LR00 to LR01 for simplicity
  mutate(site_id = if_else(site_id == "LR00", "LR01", site_id)) |> 
  # set grouping vars as factors
  mutate_if(is.character, as.factor) |> 
  # remove columns not needed for niche space analyses
  select(-study, -sex)


# Lipid corrections -------------------

sia_tidy <- sia_tidy |> 
  # Apply lipid correction to fish with high lipid content
  mutate(
    d13Cc = if_else(compartment != "baseline" & CN > 3.5,
                   d13C - 3.32 + (0.99 * CN),
                   d13C)
    ) |>
  # Organize columns
  mutate(yr_site = paste(sample_year, site_id, sep = "_")) |> 
  select(sia_sample_id, yr_site, sample_year, stream_name, site_id, sample_hitch,
         compartment, resource, taxon_code, length_mm, d13C, d13Cc, d15N, CN)

# Check data
sia_tidy

# Make subsets for whole community, bugs, and fish ------------------------------

sia_baselines <-sia_tidy |>
  filter(compartment %in% c("baseline"))

sia_inverts <-sia_tidy |>
  filter(compartment == "invert") 

sia_fish <- sia_tidy |>
  filter(compartment == "fish") 

sia_consumers <- sia_tidy |>
  filter(compartment %in% c("invert","fish"))


