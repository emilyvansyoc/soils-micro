
## Significance from GLM post-hocs
require(tidyverse)
require(ggpubr)
require(RColorBrewer)
theme_set(theme_bw())

# read data from GLMs
source("https://github.com/EmilyB17/soils-micro/raw/master/R/BiogeochemicalGLMs.R")

# function for standard error
se <- function(x) sqrt(var(x)/length(x))

# set brewer colors
brewer.pal(n = 3, name = "Dark2")
treatmentcols <- c("HDG" = "#1B9E77", "LDG" =  "#D95F02", "NG" = "#7570B3")

# get backtransformed data
dat <- datpd %>% 
  pivot_longer(cols = c(NH4_mgkgdrysoil, NPOC_mgkgdrysoil, DON_mgkgdrysoil,
                        grav_mois),
               names_to = "param", values_to = "value") %>% 
  # add PRE at 0
  pivot_wider(names_from = diffTimeSeries, values_from = value) %>% 
  mutate(PRE = 0) %>% 
  pivot_longer(cols = c(PRE, diff_24H, diff_1WK, diff_4WK), names_to = "Time", values_to = "value") %>% 
  # make horizontal
  pivot_wider(names_from = param, values_from = value)

## ---- differences within Treatment ----

sig <- posthocTime %>% filter(p.value < 0.05) 


### PLOT ALL - backtransformed
datv <- dat %>% 
  pivot_longer(cols = c(NH4_mgkgdrysoil, NPOC_mgkgdrysoil, DON_mgkgdrysoil, grav_mois),
               names_to = "param", values_to = "value") %>% 
  group_by(Treatment, Time, param) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            se = se(value)) %>% 
  ungroup() %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"))) %>% 
  mutate(Treatment = fct_recode(Treatment,
                                HDG = "HI", LDG = "LO", NG = "NO"),
         Time = fct_recode(Time,
                           PRE = "PRE",
                           `24H` = "diff_24H",
                           `1WK` = "diff_1WK",
                           `4WK` = "diff_4WK"),
         param = factor(param),
         param = fct_recode(param,
                            DON = "DON_mgkgdrysoil",
                            Moisture = "grav_mois",
                            Ammonium = "NH4_mgkgdrysoil",
                            DOC = "NPOC_mgkgdrysoil"))


# plot all data
ggplot(data = datv, 
       aes(x = Time, y = mean, group = Treatment)) +
  geom_point(aes(color = Treatment), size = 2) +
  geom_line(aes(linetype = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se,
                    linetype = Treatment, color = Treatment),
                width = 0.1) +
  labs(x = "Sampling Time", y = "% change from PRE") +
  scale_color_manual(values = treatmentcols) + 
  facet_wrap(~param, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save
#ggsave("./data/plots/all-time-lineplot.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


