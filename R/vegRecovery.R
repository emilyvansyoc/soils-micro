## Vegetation recovery

require(tidyverse)
require(emmeans)
require(ggpubr)

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

# add PRE at 0
datf <- datpd %>% 
  # add PRE at 0
  pivot_wider(names_from = diffTimeSeries, values_from = reading_rpm) %>% 
  mutate(PRE = 0) %>% 
  pivot_longer(cols = c(PRE, diff_24H, diff_1WK, diff_4WK), names_to = "Time", values_to = "reading_rpm")
  
  
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


## ---- forage utilization summary ----

uthi <- dat %>% 
  filter(Treatment == "HI") %>% 
  pivot_wider(names_from = GrazeTime, values_from = reading_rpm, names_prefix = "T") %>% 
  mutate(ut = ((T24H - TPRE) / TPRE) * 100) %>% 
  summarize(avg = mean(ut),
            se = se(ut))

utlo <- dat %>% 
  filter(Treatment == "LO") %>% 
  pivot_wider(names_from = GrazeTime, values_from = reading_rpm, names_prefix = "T") %>% 
  mutate(ut = ((T1WK - TPRE) / TPRE) * 100) %>% 
  summarize(avg = mean(ut),
            se = se(ut))

## ---- GLM on Treatment-Time ----

# sufficiently normal so no need to log transform

## GLM
mod <- glm(reading_rpm ~ Time * Treatment,
           data = datf,
           family = gaussian(link = "identity"))
# check residuals for normality
shapiro.test(resid(mod))

# add to modelFit
modelFit <-  data.frame(Param = "RPM",
                              deviance = mod$deviance,
                              null.deviance = mod$null.deviance,
                              diff = mod$null.deviance - mod$deviance,
                              df.null = mod$df.null,
                              df.dev = mod$df.residual)  

# get F statistic and DF
anova(mod, test = "F")

# post-hoc test with emmeans
posthocTrt <- data.frame(Param = "RPM",
                         emmeans(mod, pairwise ~ Treatment | Time, type = "response")$contrasts)

# post-hoc test with emmeans
posthocTime <- data.frame(Param = "RPM",
                          emmeans(mod, pairwise ~ Time | Treatment, type = "response")$contrasts)

### Significance
sigtrt <- posthocTrt %>% filter(p.value < 0.05)
sigtime <- posthocTime %>% filter(p.value < 0.05)

# get mean and se for the plot
datsum <- datf %>% 
  group_by(Treatment, Time) %>% 
  summarize(mean = mean(reading_rpm),
            se = se(reading_rpm))  %>% 
  mutate(Time = fct_recode(Time,
                           PRE = "PRE",
                           `24H` = "diff_24H",
                           `1WK` = "diff_1WK",
                           `4WK` = "diff_4WK")) %>% 
  mutate(Treatment = fct_recode(Treatment,
                                HDG = "HI", LDG = "LO", NG = "NO")) %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("PRE", "24H", "1WK", "4WK")))

### rebuild the manuscript plot with ggpubr
ggplot(data = datsum, 
       aes(x = Time, y = mean, group = Treatment)) +
  geom_point(aes(color = Treatment, shape = Treatment), size = 4) +
  geom_line(aes(linetype = Treatment, color = Treatment), size = 1.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se,
                    color = Treatment),
                width = 0.1, size = 1) +
  labs(x = "Sampling Time", y = "Relative vegetation change from PRE (%)") +
  scale_color_manual(values = treatmentcols) + 
  scale_shape_manual(values = c(19, 15, 4)) +
  
  theme(
    # pivot x text to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1)) +
  # add extra white space at top for significance asterisks
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1), add = c(0, 0)))

# save
#ggsave("./data/plots/veg-recovery-v3.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")
