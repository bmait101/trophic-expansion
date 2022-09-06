

# Baseline analysis for trophic niche expannsion

sia_baselines %>% distinct(resource)
sia_baselines %>% distinct(taxon_code)

# Carbon
c_base <- 
  bind_cols(
    # Baselines -------------------
    sia_baselines %>% 
      filter(sample_year==2016) %>% 
      filter(resource == "detritus") %>%
      # filter(resource == "terrC3") %>% 
      group_by(site_id) %>% 
      summarise(
        b.n = n(),
        b.meanC = mean(d15N, na.rm = TRUE)
      ),
    # Inverts -------------------
    sia_inverts %>% 
      filter(sample_year == 2016) %>% 
      group_by(site_id) %>% 
      summarise(
        i.n = n(),
        i.meanC = mean(d15N, na.rm = TRUE), 
        i.maxC = max(d15N, na.rm = TRUE), 
        i.minC = min(d15N, na.rm = TRUE)
      ) %>% 
      select(-site_id),
    # Fish  -------------------
    sia_fish %>% 
      filter(sample_year == 2016) %>% 
      group_by(site_id) %>% 
      summarise(
        f.n = n(),
        f.meanC = mean(d15N, na.rm = TRUE)
      ) %>% 
      select(-site_id)) 

c_base <- c_base %>% 
  mutate(
    del.i.mean  = i.meanC - b.meanC, 
    del.i.max = i.maxC - b.meanC, 
    del.i.min = i.minC - b.meanC, 
    del.f.mean   = f.meanC - b.meanC) %>% 
  select(
    site_id, 
    b.n, b.meanC, 
    i.n, i.meanC, del.i.mean, i.maxC, del.i.max, i.minC, del.i.min, 
    f.n, f.meanC, del.f.mean)
c_base

# c_base %>% 
#   mutate(across(where(is.numeric), round, 2)) %>%
#   write_csv(here("results", "baseline_D13C.csv"))

c_base <- c_base %>% 
  left_join(gradient, by = "site_id") 


c_base %>% 
  ggplot(aes(x = PC1, y = del.i.mean)) + 
  geom_point() + 
  ggpubr::stat_cor()

m1 <- gam(del.i.mean  ~ s(PC1), 
          method = "REML",
          data = c_base %>% filter(site_id!="LR08"))
draw(m1, residuals = TRUE)

c_base %>% 
  ggplot(aes(x = PC1, y = del.i.max)) + 
  geom_point() + 
  ggpubr::stat_cor()

c_base %>% 
  ggplot(aes(x = PC1, y = del.i.min)) + 
  geom_point() + 
  ggpubr::stat_cor()

c_base %>% 
  ggplot(aes(x = PC1, y = del.f.mean)) + 
  geom_point() + 
  ggpubr::stat_cor()



m1 <- lm(del.i.mean ~ PC1, data = n_base)






n_base <- 
bind_cols(
  # Baselines -------------------
  sia_baselines %>% 
    filter(sample_year==2016) %>% 
    filter(resource == "terrC3") %>% 
    group_by(site_id) %>% 
    summarise(
      b.n = n(),
      b.meanN = mean(d15N, na.rm = TRUE)
      ),
  # Inverts -------------------
  sia_inverts %>% 
    filter(sample_year == 2016) %>% 
    group_by(site_id) %>% 
    summarise(
      i.n = n(),
      i.meanN = mean(d15N, na.rm = TRUE), 
      i.maxN = max(d15N, na.rm = TRUE), 
      i.minN = min(d15N, na.rm = TRUE)
      ) %>% 
    select(-site_id),
  # Fish  -------------------
  sia_fish %>% 
    filter(sample_year == 2016) %>% 
    group_by(site_id) %>% 
    summarise(
      f.n = n(),
      f.meanN = mean(d15N, na.rm = TRUE)) %>% 
    select(-site_id)) 

n_base <- n_base %>% 
  mutate(
    del.bl.N  = b.meanN - i.meanN, 
    del.i.max = i.maxN - i.meanN, 
    del.i.min = i.minN - i.meanN, 
    del.f.N   = f.meanN - i.meanN) %>% 
  select(
    site_id, 
    b.n, b.meanN, del.bl.N, 
    i.n, i.meanN, i.maxN, del.i.max, i.minN, del.i.min, 
    f.n, f.meanN, del.f.N)
n_base

# n_base %>% 
#   mutate(across(where(is.numeric), round, 2)) %>% 
#   write_csv(here("results", "baseline_D15N.csv"))


n_base %>% 
  left_join(gradient, by = "site_id") %>% 
  ggplot(aes(x = PC1, y = del.bl.N)) + 
  geom_point() + 
  ggpubr::stat_cor()

n_base %>% 
  left_join(gradient, by = "site_id") %>% 
  ggplot(aes(x = PC1, y = del.i.max)) + 
  geom_point() + 
  ggpubr::stat_cor()

n_base %>% 
  left_join(gradient, by = "site_id") %>% 
  ggplot(aes(x = PC1, y = del.i.min)) + 
  geom_point() + 
  ggpubr::stat_cor()

n_base %>% 
  left_join(gradient, by = "site_id") %>% 
  ggplot(aes(x = PC1, y = del.f.N)) + 
  geom_point() + 
  ggpubr::stat_cor()

m1 <- lm(del.bl.N ~ PC1, data = n_base)

