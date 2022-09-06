
source(here::here("R", "00_prep.R"))

library(patchwork)


# Fishes =======================================================================

# KUD --------------------------------------------------------------------------

fish_comm_metrics %>% 
  ggplot(aes(x=PC1,y=KUD)) + 
  geom_point()


# Fit GAMs 
k <- 6
m0 <- gam(KUD ~ 
            s(PC1, k = k),
          method = "REML", 
          data = fish_comm_metrics)
m1 <- gam(KUD ~ 
            s(PC1, k = k) + s(sample_year, bs = "re"),
          method = "REML", 
          data = fish_comm_metrics)
m2 <- gam(KUD ~ 
            s(PC1, k = 6) + 
            s(sample_year, bs = "re")+ 
            s(site_id, bs = "re"),
          method = "REML", 
          data = fish_comm_metrics)

# Diagnostics
k.check(m2)
summary(m2)
appraise(m2)
draw(m2, residuals = TRUE)

# Compare models
AIC(m0, m1, m2) |> arrange(AIC)

# Predict from fitted model
new_data <- with(fish_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(fish_comm_metrics$site_id)[[1]],
  sample_year = levels(fish_comm_metrics$sample_year)[[1]]
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m2, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m1))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions
model_label <- c("s(PC1, 4.3)","'Deviance expl.' == '96.6%'")
p.kud <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey") + 
  geom_line() +
  geom_point(data=fish_comm_metrics, aes(x=PC1,y=KUD, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,8), breaks = seq(0,8,2)) +
  scale_y_continuous(limits=c(0,40), breaks = seq(0,40,10)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Fish trophic diversity ('~KUD[95]~')'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(40, 40*0.9), 
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3) 
p.kud


# CR ------------------------------------------------------------------


fish_comm_metrics %>% 
  ggplot(aes(x=PC1,y=CR)) + 
  geom_point()


# Fit gam
k <- 6
m0 <- gam(CR ~ 
            s(PC1, k = k),
          method = "REML", 
          data = fish_comm_metrics)
m1 <- gam(CR ~ 
            s(PC1, k = k) + s(sample_year, bs = "re"),
          method = "REML", 
          data = fish_comm_metrics)
m2 <- gam(CR ~ 
            s(PC1, k = 6) + 
            s(sample_year, bs = "re")+ 
            s(site_id, bs = "re"),
          method = "REML", 
          data = fish_comm_metrics)

# Diagnostics
k.check(m2)
summary(m2)
appraise(m2)
draw(m2, residuals = TRUE)

# Compare models
AIC(m0, m1, m2) |> arrange(AIC)


# Predict from fitted model
new_data <- with(fish_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  sample_year = levels(fish_comm_metrics$sample_year)[[1]],
  site_id = levels(fish_comm_metrics$site_id)[[1]]
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m2, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m1))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions

model_label <- c("s(PC1, 1.00)","'Deviance expl.' == '98.8%'")
p.CR <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey") + 
  geom_line() +
  geom_point(data=fish_comm_metrics, aes(x=PC1,y=CR, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,8), breaks = seq(0,8,2)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks = seq(0,10,2)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Fish '~{delta}^13*C~'(\u2030) range'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(10, 10*0.9),
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3)

p.CR

# NR ------------------------------------------------------------------


fish_comm_metrics %>% 
  ggplot(aes(x=PC1,y=NR)) + 
  geom_point()


# Fit GAMs 
k <- 6
m0 <- gam(NR ~ 
            s(PC1, k = k),
          method = "REML", 
          data = fish_comm_metrics)
m1 <- gam(NR ~ 
            s(PC1, k = k) + s(sample_year, bs = "re"),
          method = "REML", 
          data = fish_comm_metrics)
m2 <- gam(NR ~ 
            s(PC1, k = 6) + 
            s(sample_year, bs = "re")+ 
            s(site_id, bs = "re"),
          method = "REML", 
          data = fish_comm_metrics)

# Diagnostics
k.check(m2)
summary(m2)
appraise(m2)
draw(m2, residuals = TRUE)

# Compare models
AIC(m0, m1, m2) |> arrange(AIC)


# Predict from fitted model
new_data <- with(fish_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  sample_year = levels(fish_comm_metrics$sample_year)[[1]],
  site_id = levels(fish_comm_metrics$site_id)[[1]]
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m1, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m1))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions
model_label <- c("s(PC1, 1.00)", "'Deviance expl.' == '99.2%'")
p.NR <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+
  geom_line() +
  geom_point(data=fish_comm_metrics, aes(x=PC1,y=NR, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,8), breaks = seq(0,8,2)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks = seq(0,10,2)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Fish '~{delta}^15*N~'(\u2030) range'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(10, 10*0.9),
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3)

p.NR

# Bugs =========================================================================

# Cannot run with random effetcs, more coefs than data

# KUD ---------------------------------------------------------------------------


invert_comm_metrics %>% 
  ggplot(aes(x=PC1,y=KUD)) + 
  geom_point()


# Fit gam
k <- 3
m0 <- gam(KUD ~ 
            s(PC1, k = 3),
          method = "REML", 
          data = invert_comm_metrics)

# Diognostics
k.check(m0)
summary(m0)
appraise(m0)
draw(m0, residuals = TRUE)

# Predict from fitted model
new_data <- with(invert_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(invert_comm_metrics$site_id)[[1]]
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m0, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m1))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions
model_label <- c("s(PC1, 2.68)","'Deviance expl.' == '60.2%'")
p.kud.i <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+ 
  geom_line() + 
  geom_point(data=invert_comm_metrics, aes(x=PC1,y=KUD, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,8), breaks = seq(0,8,2)) +
  scale_y_continuous(limits=c(0,70), breaks = seq(0,70,10)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Invertebrate trophic diversity ('~KUD[95]~')'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(70, 70*0.9), 
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3) 

p.kud.i



# CR ---------------------------------------------------------------------------

invert_comm_metrics %>% 
  ggplot(aes(x=PC1,y=CR)) + 
  geom_point()


# Fit gam
k <- 3
m0 <- gam(CR ~ 
            s(PC1, k = k),
          method = "REML", 
          data = invert_comm_metrics)

k.check(m0)
summary(m0)
appraise(m0)
draw(m0, residuals = TRUE)


# Predict from fitted model
new_data <- with(invert_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(invert_comm_metrics$site_id)[[1]]
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m0, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m1))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions

model_label <- c("s(PC1, 2.13)", "'Deviance expl.' == '37.1%'")
p.CR.i <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+ 
  geom_line() + 
  geom_point(data=invert_comm_metrics, aes(x=PC1,y=CR, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,8), breaks = seq(0,8,2)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks = seq(0,10,2)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Invertebrate '~{delta}^13*C~'(\u2030) range'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 2, y = c(4, 4*0.8), 
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3) 

p.CR.i

# NR ---------------------------------------------------------------------------

invert_comm_metrics %>% 
  ggplot(aes(x=PC1,y=NR)) + 
  geom_point()


# Fit gam
k <- 3
m0 <- gam(NR ~ 
            s(PC1, k = k),
          method = "REML", 
          data = invert_comm_metrics)

k.check(m0)
summary(m0)
appraise(m0)
draw(m0, residuals = TRUE)


# Predict from fitted model
new_data <- with(invert_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(invert_comm_metrics$site_id)[[1]]
))

new_data <- bind_cols(
  new_data, 
  as_tibble(
    predict(
      m0, 
      newdata = new_data, 
      terms="s(PC1)",
      se.fit=TRUE))
)

crit.t <- qt(0.975, df = df.residual(m1))
new_data <- new_data %>% 
  mutate(upper=fit+(crit.t*se.fit), lower=fit-(crit.t*se.fit))

# PLot predictions
model_label <- c("s(PC1, 1.00)", "'Deviance expl.' == '26.5%'")

p.NR.i <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+
  geom_line() +
  geom_point(data=invert_comm_metrics, aes(x=PC1,y=NR, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,8), breaks = seq(0,8,2)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks = seq(0,10,2)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Invertebrate '~{delta}^15*N~'(\u2030) range'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(10, 10*0.9),
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3)
 
p.NR.i


# Panel Plot ===================================================================

# Build panels
panel <- cowplot::plot_grid(
  p.kud.i + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.CR.i + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.NR.i + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.kud + theme(legend.position=c(.75,.2)) + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.CR + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.NR + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A", "B", "C","D","E","F"),
  ncol = 3, nrow = 2, 
  align = "vh"
)

# Print it
panel


# Save it
path <- here::here("out", "r1_td_panel")
ggsave(glue::glue("{path}.pdf"), plot = last_plot(), 
       width = 14, height = 8, device = cairo_pdf)
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png", dpi = 300)


