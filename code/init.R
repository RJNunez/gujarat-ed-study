# -- Libraries
library(lubridate)
library(ggplot2)
library(splines)
library(stringr)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
library(directlabels)

# -- My color palette
my_palette        <- c("#456e9d", "#f08127", "#6baeaa", "#509745", "#eac240", "#a66f97", "#ff929d", "#D22B2B", "#252525")
names(my_palette) <- c("blue", "orange", "turquoise", "green", "yellow", "purple", "pink", "red", "black")

# -- Set up for figures
theme_set(theme_minimal(base_size   = 12, 
                        base_family = "Helvetica"))

# -- Modifying plot elements globally
theme_update(
  axis.ticks        = element_line(color = "grey92"),
  axis.ticks.length = unit(.5, "lines"),
  panel.grid.minor  = element_blank(),
  legend.title      = element_text(size = 12),
  legend.text       = element_text(color = "grey30"),
  legend.background = element_rect(color = "black", fill = "#FBFCFC"),
  legend.key        = element_rect(fill = "#FBFCFC"),
  legend.direction  = "horizontal",
  legend.position   = "top",
  plot.title        = element_text(size = 18, face = "bold"),
  plot.subtitle     = element_text(size = 12, color = "grey30"),
  plot.caption      = element_text(size = 9, margin = margin(t = 15)),
  plot.background   = element_rect(fill="#FBFCFC", color = "#FBFCFC"),
  panel.background  = element_rect(fill="#FBFCFC", color = NA),
  strip.text        = element_text(face = "bold", color = "white"),
  strip.background  = element_rect(fill = "#252525"))

# -- To be used in the function below
noleap_yday <- function(x){
  ifelse(lubridate::leap_year(x) & lubridate::month(x)>2,
         lubridate::yday(x)-1,
         lubridate::yday(x))
}

# -- Function to create Fourier basis
fourier_trend <- function(x, k = 3){
  H <- lapply(1:k, function(k){
    cbind(sin(2*pi*k/365*x), cos(2*pi*k/365*x))
  })
  res <- do.call(cbind, H)
  colnames(res) <- paste(rep(c("sin", "cos"), k), rep(1:k, each = 2), sep="_")
  res
}