### Nitrate troubleshooting

require(tidyverse)
theme_set(theme_bw())

# read data
nit <- read.table("https://github.com/EmilyB17/soils-micro/raw/master/data/2018CN.txt",
                  sep = "\t", header = TRUE) %>% 
  select(Plot, Block, Treatment, GrazeTime, NO3_mgkgdrysoil) %>% 
  mutate(GrazeTime = factor(GrazeTime, ordered = TRUE, levels = c("PRE", "24H", "1WK", "4WK"))) %>% 
  # REMOVE BLOCK 4
  filter(!Block == 4)

# explore data
ggplot(data = nit, aes(x = NO3_mgkgdrysoil)) +
  geom_histogram(binwidth = 1) +
  labs(x = "Nitrate mg/kg dry soil", y = "Count", title = "Nitrate Histogram")

ggplot(data = nit, aes(x = GrazeTime, y = NO3_mgkgdrysoil)) +
  geom_boxplot() +
  labs(x = "GrazeTime", y = "Nitrate", title = "Outliers at 1 and 4WKs")
                 

## ---- Percent change GLM ----

source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")

# calculate percent diff
nitpd <- percentDifferences(df = nit,
                            ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                            timeKey = "GrazeTime",
                            timeLevels = c("PRE", "24H", "1WK", "4WK"),
                            level1 = "PRE") %>% 
  pivot_wider(names_from = diffTimeSeries, values_from = NO3_mgkgdrysoil) %>% 
  mutate(PRE = 0) %>% 
  pivot_longer(cols = c(PRE, diff_24H, diff_1WK, diff_4WK), names_to = "Time", 
               values_to = "NO3_mgkgdrysoil") %>% 
  # LOG TRANSFORM
  mutate(logvalue = log1p(NO3_mgkgdrysoil + 101))

# explore
ggplot(data = nitpd, aes(x = logvalue)) +
  geom_histogram() +
  labs(x = "Log percent change", y = "Count", title = "Nitrate Histogram PD")

ggplot(data = nitpd, aes(x = Time, y = logvalue)) +
  geom_boxplot() +
  labs(x = "Time", y = "Log percent change", title = "Wonkiness in percent change")

shapiro.test(nitpd$logvalue) # strongly not normal

## GLM
mod <- glm(logvalue ~ Time * Treatment,
           data = nitpd,
           family = gaussian(link = "identity"))
## MODEL FIT PARAMETERS
plot(resid(mod))
qqnorm(resid(mod))
qqline(resid(mod))
shapiro.test(resid(mod))

### ---- LMER (Gordon) ----

# LMER structure on un-transformed data

# need Time 1 (PRE) as a separate vector or column
ts <- nit %>% 
  # make wider by GrazeTime
  pivot_wider(names_from = "GrazeTime", names_prefix = "T_", 
              values_from = "NO3_mgkgdrysoil") %>% 
  # make narrow by time but keep PRE separate
  pivot_longer(cols = c(T_24H, T_1WK, T_4WK), names_to = "GrazeTime", values_to = "value")

## LMER
require(multcomp)
require(lme4)

mod1 <- lmer(formula = "value ~ Treatment + T_PRE + GrazeTime + GrazeTime*Treatment + (1|GrazeTime)", data = ts) # get convergence warning; not a good model fit

## MODEL FIT
plot(resid(mod1))
qqnorm(resid(mod1))
qqline(resid(mod1))
shapiro.test(resid(mod1))

#### ---- NONPARAMETRIC ----

## problem: there is no way to get an idea of interaction with non-parametric
# since these tests depend on sign changes and ranking
# instead of analysis of variance

## SO: for supplementary material, perform a Kruskal-Wallis test for each
# Time (vs Treatment) and for each Treatment (over time) to ensure 
# we are not missing trends that could be otherwise important

times <- unique(nit$GrazeTime)
pvalTime <- data.frame()

for(i in 1:length(times)) {
  # kruskal-wallis
  k <- kruskal.test(NO3_mgkgdrysoil ~ Treatment, data = filter(nit, GrazeTime == times[i]))
  # get p value
  pvalTime <- rbind(pvalTime,
                    data.frame(Time = times[i],
                               p = k$p.value))

}

print(pvalTime) # no significance

# differences within each treatment over time
trt <- unique(nit$Treatment)
pvalTrt <- data.frame()

for(i in 1:length(trt)) {
  # kruskal wallis
  k <- kruskal.test(NO3_mgkgdrysoil ~ GrazeTime,
                    data = filter(nit, Treatment == trt[i]))
  # get p value
  pvalTrt <- rbind(pvalTrt,
                   data.frame(Trt = trt[i],
                              p = k$p.value))
}

print(pvalTrt) # no significance

# function for standard error
se <- function(x) sqrt(var(x)/length(x))
# get data for supplementary materials since this isn't presented in % diff
df <- nit %>% 
  group_by(Treatment, GrazeTime) %>% 
  summarize(meannitr = round(mean(NO3_mgkgdrysoil), 3),
            senitr = round(se(NO3_mgkgdrysoil), 3)) 
write.table(df, "./data/nitrate-summary.txt", sep = "\t", row.names = FALSE)
