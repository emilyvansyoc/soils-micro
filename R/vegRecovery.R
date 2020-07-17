## Vegetation recovery

require(tidyverse)

# function for standard error
se <- function(x) sqrt(var(x)/length(x))

# calculate percent difference
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")

# colors
source("./R/RColorBrewer.R")

# get data
dat <- read.table("./data/2018RPM.txt", sep = "\t", header = TRUE) %>% 
  # because the biomass conversion is unreliable, we use "raw" RPM reading as relative biomass
  select(-biomass_kg_plot) %>% 
  # remove Block 4
  filter(!Block == 4)

# LDG finished grazing at 1WK, HDG finished at 24H, NG had no grazing so calculate from PRE
hdg <- dat %>% 
  filter(Treatment == "HI") %>% 
  pivot_wider(names_from = GrazeTime, values_from = reading_rpm, names_prefix = "T") %>% 
  mutate(recov = ((T4WK - T24H) / T24H) * 100)
