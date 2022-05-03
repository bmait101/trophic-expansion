
library(tidyverse)
library(mgcv)
library(gratia)
library(lme4)
library(lmerTest)
library(patchwork)



# Whole community ==============================================================


# KUD --------------------------------------------------------------------------

# Fit gam
k <- 3
m1 <- gam(KUD ~ 
            s(PC1, k = k), 
          method = "REML", 
          data = consumer_comm_metrics %>% filter(site_id!="LR00"))

k.check(m1)
summary(m1)
appraise(m1)
draw(m1, residuals = TRUE, select = "s(PC1)")

# Fit linear model
m2 <- lm(KUD ~ PC1, data = consumer_comm_metrics%>% filter(site_id!="LR00"))
summary(m2)

# Compare models
AIC(m1, m2)
BIC(m1, m2)

consumer_comm_metrics %>% 
  filter(site_id!="LR00") %>% 
  ggplot(aes(x=PC1,y=KUD)) + 
  geom_point(color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,9), breaks = seq(0,9,1)) 


# Predict from fitted model
new_data <- with(consumer_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(consumer_comm_metrics$site_id)[[1]]
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
model_label <- c("s(PC1, 1.95)","'Deviance expl.' == '62.2%'")
p.kud.i <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+ 
  geom_line() + 
  geom_point(data=consumer_comm_metrics, aes(x=PC1,y=KUD, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(limits=c(0,9), breaks = seq(0,9,1)) +
  # scale_y_continuous(limits=c(0,35), breaks = seq(0,35,5)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Invertebrate TD ('~KUD[95]~')'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(70, 70*0.9), 
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3) 

p.kud.i



# CR ---------------------------------------------------------------------------

# Fit gam
k <- 3
m1 <- gam(CR ~ 
            s(PC1, k = k),
          method = "REML", 
          data = consumer_comm_metrics)

k.check(m1)
summary(m1)
appraise(m1)
draw(m1, residuals = TRUE, select = "s(PC1)")


# Fit linear model
m2 <- lm(KUD ~ PC1, data = consumer_comm_metrics)
summary(m2)

# Compare models
AIC(m1, m2)
BIC(m1, m2)


# Predict from fitted model
new_data <- with(consumer_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(consumer_comm_metrics$site_id)[[1]]
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

model_label <- c("s(PC1, 1.86)", "'Deviance expl.' == '37.6%'")
p.CR.i <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+ 
  geom_line() + 
  geom_point(data=consumer_comm_metrics, aes(x=PC1,y=CR, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(expand=c(0,0), limits=c(0,9), breaks = seq(0,9,1)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks = seq(0,10,2)) +
  labs(title = "", x = "Long. gradient (PC1)", 
       y = expression('Invertebrate '~{delta}^13*C~'(\u2030) range'), 
       fill = "Stream system") +
  annotate(geom = "text", x = 0.5, y = c(10, 10*0.9), 
           hjust = 0, vjust = 1,  label = model_label, parse = TRUE, size = 3) 

p.CR.i

# NR ---------------------------------------------------------------------------

# Fit gam
k <- 3
m1 <- gam(NR ~ 
            s(PC1, k = k),
          method = "REML", 
          data = consumer_comm_metrics)

k.check(m1)
summary(m1)
appraise(m1)
draw(m1, residuals = TRUE, select = "s(PC1)")


# Fit linear model
m2 <- lm(NR ~ PC1, data = consumer_comm_metrics)
summary(m2)

# Compare models
AIC(m1, m2)
BIC(m1, m2)


# Predict from fitted model
new_data <- with(consumer_comm_metrics, tibble(
  PC1 = seq(min(PC1), max(PC1), length.out = 200), 
  site_id = levels(consumer_comm_metrics$site_id)[[1]]
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
model_label <- c("s(PC1, 1.00)", "'Deviance expl.' == '25.3%'")

p.NR.i <- new_data %>% 
  ggplot(aes(x=PC1,y=fit)) + 
  geom_ribbon(aes(ymin=lower,ymax=upper, x = PC1), alpha=.5, fill="grey")+ 
  geom_line() + 
  geom_point(data=consumer_comm_metrics, aes(x=PC1,y=NR, fill=stream_name), 
             color = "black", size = 2, shape = 21) + 
  scale_color_brewer(palette = "Dark2") + 
  scale_fill_brewer(palette = "Dark2") + 
  scale_x_continuous(expand=c(0,0), limits=c(0,9), breaks = seq(0,9,1)) +
  # scale_y_continuous(expand=c(0,0), limits=c(0,10), breaks = seq(0,10,2)) +
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
  p.kud + theme(legend.position=c(.75,.25)) + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.CR + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  p.NR + theme(legend.position="none") + theme(plot.margin=unit(c(1,1,1,1),"mm")), 
  labels = c("A", "B", "C","D","E","F"),
  ncol = 3, nrow = 2, 
  align = "vh"
)
panel

# Save it
ggsave(here("figs1", "td_panel.pdf"), plot = panel, device = cairo_pdf, 
       units = "in", width = 10, height = 6)



