
## Significance from GLM post-hocs
require(tidyverse)
require(ggpubr)
require(RColorBrewer)
theme_set(theme_bw())

# read data from GLMs
source("https://github.com/EmilyB17/soils-micro/raw/master/R/BiogeochemicalGLMs.R")

# get RColorBrewer colors
source("./R/RColorBrewer.R")

# function for standard error
se <- function(x) sqrt(var(x)/length(x))

# set brewer colors
#brewer.pal(n = 3, name = "Dark2")
#treatmentcols <- c("HDG" = "#1B9E77", "LDG" =  "#D95F02", "NG" = "#7570B3")

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
                            `B. Dissolved organic nitrogen` = "DON_mgkgdrysoil",
                            `A. Soil water content` = "grav_mois",
                            `D. Ammonium` = "NH4_mgkgdrysoil",
                            `C. Dissolved organic carbon` = "NPOC_mgkgdrysoil")) %>% 
  # order the variables based on order in the text
  mutate(param = factor(param, ordered = TRUE, levels = c("A. Soil water content",
                                                          "B. Dissolved organic nitrogen",
                                                          "C. Dissolved organic carbon",
                                                          "D. Ammonium")))


# plot all data
ggplot(data = datv, 
       aes(x = Time, y = mean, group = Treatment)) +
  geom_point(aes(color = Treatment, shape = Treatment), size = 3) +
  geom_line(aes(linetype = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se,color = Treatment),
                width = 0.1) +
  labs(x = "Sampling Time", y = "% change from PRE") +
  scale_color_manual(values = treatmentcols) + 
  scale_shape_manual(values = c(19, 15, 4)) +
  facet_wrap(~param, scales = "free") +
  
  theme(
    # pivot x text to 45 degrees
    axis.text.x = element_text(angle = 45, hjust = 1),
    # remove grey background from facet labels
        strip.background = element_rect(
          fill="white", linetype=0
        ), 
    # change font size of facet labels
        strip.text.x = element_text(size = 11, hjust = 0)) +
  # add extra white space at top for significance asterisks
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1), add = c(0, 0)))

# save
#ggsave("./data/plots/all-time-lineplot-v3.png", plot = last_plot(), dpi = 300, height = 4.67, width = 6.48, units = "in")


