## Vegetation recovery

require(tidyverse)

# function for standard error
se <- function(x) sqrt(var(x)/length(x))

# calculate percent difference
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")

# colors
source("./R/RColorBrewer.R")
theme_set(theme_bw())

# get data
dat <- read.table("./data/2018RPM.txt", sep = "\t", header = TRUE) %>% 
  # because the biomass conversion is unreliable, we use "raw" RPM reading as relative biomass
  select(-biomass_kg_plot) %>% 
  # remove Block 4
  filter(!Block == 4)

# calculate percent difference
datpd <- percentDifferences(df = dat,
                            ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                            timeKey = "GrazeTime",
                            timeLevels = c("PRE", "24H", "1WK", "4WK"),
                            level1 = "PRE") 

datsum <- datpd %>% 
  pivot_wider(names_from = diffTimeSeries, values_from = reading_rpm) %>% 
  mutate(PRE = 0) %>%
  pivot_longer(cols = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"), names_to = "Time", values_to = "reading_rpm") %>% 
  group_by(Treatment, Time) %>% 
  summarize(rpm = mean(reading_rpm),
            serpm = se(reading_rpm)) %>% 
  mutate(
         Treatment = case_when(
           Treatment %in% "HI" ~ "HDG",
           Treatment %in% "LO" ~ "LDG",
           Treatment %in% "NO" ~ "NG"
         ))  %>% 
  mutate(Time1 = case_when(
    Time %in% "PRE" ~ "PRE",
    Time %in% "diff_24H" ~ "post24H",
    Time %in% "diff_1WK" ~ "post1WK",
    Time %in% "diff_4WK" ~ "post4WK"
  ),
  Time1 = factor(Time1, ordered = TRUE, levels = c("PRE", "post24H", "post1WK", "post4WK")))

ggplot(data = datsum, aes(x = Time1, y = rpm, group = Treatment, color = Treatment)) +
  geom_point(size = 3) + 
  geom_line(aes(linetype = Treatment), size = 1) +
  geom_errorbar(aes(ymin = rpm - serpm, ymax = rpm + serpm), width = 0.2, color = "black") +
  scale_color_manual(values = treatmentcols) +
  labs(x = "Sampling Time", y = "Percent change in relative biomass") +
  theme_bw()
  
# save
#ggsave("./data/plots/veg-over-time.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


## ---- recovery per day ----
# LDG finished grazing at 1WK, HDG finished at 24H, NG had no grazing so calculate from PRE
hdg <- dat %>% 
  filter(!Treatment == "LO") %>% 
  pivot_wider(names_from = GrazeTime, values_from = reading_rpm, names_prefix = "T") %>% 
  mutate(recov = T4WK-T24H,
         recovday = (T4WK - T24H)  / 28) %>% 
  mutate(Treatment = case_when(
    Treatment %in% "HI" ~ "HDG",
    Treatment %in% "NO" ~ "NG"
  ))
# are there differences?
t.test(recov ~ Treatment, data = hdg) # no significance


ldg<- dat %>% 
  filter(!Treatment == "HI") %>% 
  pivot_wider(names_from = GrazeTime, values_from = reading_rpm, names_prefix = "T") %>% 
  mutate(recov = T4WK-T1WK,
         recovday = (T4WK - T1WK)  / 22) %>% 
  mutate(Treatment = case_when(
    Treatment %in% "LO" ~ "LDG",
    Treatment %in% "NO" ~ "NG"
  ))
t.test(recov ~ Treatment, data = ldg) # no significance

## PLOT

ggplot(data = hdg, aes(x = Treatment, y = recov, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  labs(x = "Treatment", y = "Relative vegetation recovery from 24H to 4WK")

# save
#ggsave("./data/plots/veg-HDGandNG-recov.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


ggplot(data = ldg, aes(x = Treatment, y = recov, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  labs(x = "Treatment", y = "Relative vegetation recovery from 1WK to 4WK")
# save
#ggsave("./data/plots/veg-LDGandNG-recov.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


  