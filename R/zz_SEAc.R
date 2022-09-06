library(SIBER)

siber_df <- rkin_df %>% 
  select(iso1=d13C_scl, iso2=d15N_scl, group=yr_site, community=yr_site) %>% 
  mutate(group=1) %>% 
  data.frame()
siber.example <- createSiberObject(siber_df)

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)
seac <- group.ML %>% as.data.frame() %>% rownames_to_column() %>% 
  rename(metric=rowname) %>% 
  pivot_longer(-metric, names_to="yr_site", values_to = "value") %>% 
  separate(yr_site, into = c("yr_site","B"), sep = "\\.") %>% 
  filter(metric=="SEAc") %>% 
  select(-B, -metric) %>%
  rename(SEAc = value)
