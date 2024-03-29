
# Prepare work space for analyses

# Load libraries
pacman::p_load(
  tidyverse, 
  here, 
  cowplot, 
  mgcv, 
  gratia, 
  vegan, 
  funrar, 
  broom, 
  rKIN, 
  lme4, 
  lmerTest
  )

# Load custom plotting theme
source(here("R", "fx_theme_pub.R"))

# Set custom theme as default
theme_set(theme_Publication())

# Create labels for all stream systems
label_streams <- c("Laramie","MedBow","Sweetwater","Horse")

# Create site id labels
label_sites_ids <- 
  c(
    "LR00","LR01","LR02","LR03","LR04","LR05","LR06","LR07","LR08",
    "MB01","MB02","MB03","MB04",
    "SW01","SW02","SW03","SW04",
    "HC00","HC01","HC02"
    )

