# Custom ggplot2 theme for publication quality plots

theme_Publication <- function(base_size=12, base_family="Helvetica") {
  # ggthemes for "theme_foundation" - better than building off a default theme
  #library(ggthemes)
  
  # Begin construction of chart
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(  
      # coord_cartesian(clip = 'off'),
      # Plotting/chart region (panel)
      panel.background   = element_rect(fill = "white", color = NA),
      panel.border       = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing      = unit(0.25, "lines"),
      panel.ontop       = FALSE,

      # Plot
      plot.background  = element_rect(colour = NA),
      plot.title       = element_text(hjust = 0),
      plot.subtitle    = element_text(size = rel(1)),
      plot.caption     = element_text(size = rel(0.75)),
      plot.margin      = unit(c(5, 5, 5, 5),"mm"),

      # Strips for paneling
      strip.background  = element_rect(fill = "white", colour = "black", size = 0),
      strip.text.x      = element_text(colour = "black", size = rel(1),
                                       margin = margin(0.1,0,0.1,0, "cm")),
      strip.text.y      = element_text(colour = "black", size = rel(1), angle = -90),

      # Axes
      axis.title.x  = element_text(colour = "black", size = rel(1), margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.text.x   = element_text(colour = "black", size = rel(1), margin = margin(t = 3, r = 0, b = 0, l = 0)),
      axis.title.y  = element_text(colour = "black", size = rel(1), angle = 90, margin = margin(t = 0, r = 3, b = 0, l = 0)),
      axis.text.y   = element_text(colour = "black", size = rel(1), hjust = 1, margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.line     = element_line(colour = "black", size = 1),
      axis.ticks    = element_line(colour = "black", size = 0.5),
      axis.ticks.length = unit(.2, "cm"),

      # Legend
      legend.background   = element_rect(color = NA),
      legend.key          = element_rect(colour = NA),
      legend.key.size     = unit(0.8, "lines"),
      legend.text         = element_text(size = rel(0.8)),
      legend.title        = element_text(size = rel(0.8), face = "bold", hjust = 0),
      legend.margin       = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")
    )
}
