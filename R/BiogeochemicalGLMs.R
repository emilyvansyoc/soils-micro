### GLMs on biogeochemical data

require(tidyverse)
require(emmeans)
theme_set(theme_bw())
# function for standard error
se <- function(x) sqrt(var(x)/length(x))

# read data
dat <- read.table("https://github.com/EmilyB17/soils-micro/raw/master/data/2018CN.txt",
                  sep = "\t", header = TRUE) %>% 
  # remove nitrate and mineral N
  dplyr::select(-c(NO3_mgkgdrysoil, mineralN_mgkgdrysoil)) %>% 
  mutate(GrazeTime = factor(GrazeTime, ordered = TRUE, levels = c("PRE", "24H", "1WK", "4WK"))) %>% 
  # REMOVE BLOCK 4
  filter(!Block == 4)

# calculate percent difference
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")

# calculate
datpd <- percentDifferences(df = dat,
                            ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                            timeKey = "GrazeTime",
                            timeLevels = c("PRE", "24H", "1WK", "4WK"),
                            level1 = "PRE") 

# add log transformation for normality
datlog <- datpd %>% 
  pivot_longer(cols = c(NH4_mgkgdrysoil, NPOC_mgkgdrysoil, DON_mgkgdrysoil,
                        grav_mois),
               names_to = "param", values_to = "value") %>% 
  # add PRE at 0
  pivot_wider(names_from = diffTimeSeries, values_from = value) %>% 
  mutate(PRE = 0) %>% 
  pivot_longer(cols = c(PRE, diff_24H, diff_1WK, diff_4WK), names_to = "Time", values_to = "value") %>% 
  
  # LOG TRANSFORMATION
  mutate(logvalue = log1p(value + 101)) %>% 
  dplyr::select(-value) %>% 
  # make horizontal
  pivot_wider(names_from = param, values_from = logvalue, names_prefix = "log_")

# get summarized data for significance tables
sum <- datpd %>% 
  group_by(Treatment, diffTimeSeries) %>% 
  summarize(meanDON = mean(DON_mgkgdrysoil),
            seDON = se(DON_mgkgdrysoil),
            meanGrav = mean(grav_mois),
            seGrav = se(grav_mois),
            meanNH4 = mean(NH4_mgkgdrysoil),
            seNH4 = se(NH4_mgkgdrysoil),
            meanDOC = mean(NPOC_mgkgdrysoil),
            seDOC = se(NPOC_mgkgdrysoil)) %>% 
  ungroup()
#write.table(sum, "./data/biophys-percentchange-sum.txt", sep = "\t", row.names = FALSE)


## get empty dataframes to fill as we go
modelFit <- data.frame()
posthocTrt <- data.frame()
posthocTime <- data.frame()

## ---- NPOC (Dissolved organic carbon) ----

npoc <- datlog %>% dplyr::select(Plot, Block, Treatment, Time, log_NPOC_mgkgdrysoil)

# histogram
hist(npoc$log_NPOC_mgkgdrysoil)
# test for normality
shapiro.test(npoc$log_NPOC_mgkgdrysoil)
# QQ plot
qqnorm(npoc$log_NPOC_mgkgdrysoil)
qqline(npoc$log_NPOC_mgkgdrysoil)

## GLM
mod <- glm(log_NPOC_mgkgdrysoil ~ Time * Treatment,
           data = npoc,
           family = gaussian(link = "identity"))
# check residuals for normality
shapiro.test(resid(mod))

# add to modelFit
modelFit <-  rbind(modelFit,
                   data.frame(Param = "NPOC",
                              deviance = mod$deviance,
                              null.deviance = mod$null.deviance,
                              diff = mod$null.deviance - mod$deviance,
                              df.null = mod$df.null,
                              df.dev = mod$df.residual))   

# post-hoc test with emmeans
posthocTrt <- rbind(posthocTrt,
                    data.frame(Param = "NPOC",
                               emmeans(mod, pairwise ~ Treatment | Time, type = "response")$contrasts))

# post-hoc test with emmeans
posthocTime <- rbind(posthocTime,
                     data.frame(Param = "NPOC",
                                emmeans(mod, pairwise ~ Time | Treatment, type = "response")$contrasts))

### ---- NH4 (Ammonium) ----

amm <- datlog %>% 
  dplyr::select(Plot, Block, Treatment, Time, log_NH4_mgkgdrysoil)

# histogram
hist(amm$log_NH4_mgkgdrysoil)
# test for normality
shapiro.test(amm$log_NH4_mgkgdrysoil)
# QQ plot
qqnorm(amm$log_NH4_mgkgdrysoil)
qqline(amm$log_NH4_mgkgdrysoil)

## GLM
mod <- glm(log_NH4_mgkgdrysoil ~ Time * Treatment,
           data = amm,
           family = gaussian(link = "identity"))
# check residuals for normality
shapiro.test(resid(mod))

# add to modelFit
modelFit <-  rbind(modelFit,
                   data.frame(Param = "NH4",
                              deviance = mod$deviance,
                              null.deviance = mod$null.deviance,
                              diff = mod$null.deviance - mod$deviance,
                              df.null = mod$df.null,
                              df.dev = mod$df.residual))   

# post-hoc test with emmeans
posthocTrt <- rbind(posthocTrt,
                    data.frame(Param = "NH4",
                               emmeans(mod, pairwise ~ Treatment | Time, type = "response")$contrasts))

# post-hoc test with emmeans
posthocTime <- rbind(posthocTime,
                     data.frame(Param = "NH4",
                                emmeans(mod, pairwise ~ Time | Treatment, type = "response")$contrasts))

### ---- DON (dissolved organic nitrogen) ----

don <- datlog %>% 
  dplyr::select(Plot, Block, Treatment, Time, log_DON_mgkgdrysoil)

# histogram
hist(don$log_DON_mgkgdrysoil)
# test for normality
shapiro.test(don$log_DON_mgkgdrysoil) # NOT NORMAL
# QQ plot
qqnorm(don$log_DON_mgkgdrysoil)
qqline(don$log_DON_mgkgdrysoil)

## GLM
mod <- glm(log_DON_mgkgdrysoil ~ Time * Treatment,
           data = don,
           family = gaussian(link = "identity"))
# check residuals for normality
shapiro.test(resid(mod))

# add to modelFit
modelFit <-  rbind(modelFit,
                   data.frame(Param = "DON",
                              deviance = mod$deviance,
                              null.deviance = mod$null.deviance,
                              diff = mod$null.deviance - mod$deviance,
                              df.null = mod$df.null,
                              df.dev = mod$df.residual))   

# post-hoc test with emmeans
posthocTrt <- rbind(posthocTrt,
                    data.frame(Param = "DON",
                               emmeans(mod, pairwise ~ Treatment | Time, type = "response")$contrasts))

# post-hoc test with emmeans
posthocTime <- rbind(posthocTime,
                     data.frame(Param = "DON",
                                emmeans(mod, pairwise ~ Time | Treatment, type = "response")$contrasts))

### ---- GRAVIMETRIC MOISTURE ----

grav <- datlog %>% 
  dplyr::select(Plot, Block, Treatment, Time, log_grav_mois)

# histogram
hist(grav$log_grav_mois)
# test for normality
shapiro.test(grav$log_grav_mois) # NOT NORMAL
# QQ plot
qqnorm(grav$log_grav_mois)
qqline(grav$log_grav_mois)

## GLM
mod <- glm(log_grav_mois ~ Time * Treatment,
           data = grav,
           family = gaussian(link = "identity"))
mod <- glm(log.moisture ~ Time * Treatment,
           data = gravpd,
           family = gaussian(link = "identity"))
# check residuals for normality
shapiro.test(resid(mod)) # marginally normal

# add to modelFit
modelFit <-  rbind(modelFit,
                   data.frame(Param = "GravMois",
                              deviance = mod$deviance,
                              null.deviance = mod$null.deviance,
                              diff = mod$null.deviance - mod$deviance,
                              df.null = mod$df.null,
                              df.dev = mod$df.residual))   

# post-hoc test with emmeans
posthocTrt <- rbind(posthocTrt,
                    data.frame(Param = "GravMois",
                               emmeans(mod, pairwise ~ Treatment | Time, type = "response")$contrasts))

# post-hoc test with emmeans
posthocTime <- rbind(posthocTime,
                     data.frame(Param = "GravMois",
                                emmeans(mod, pairwise ~ Time | Treatment, type = "response")$contrasts))

### ---- write output tables ----


#write.table(posthocTime, "./data/GLM-posthoc-Time.txt", sep = "\t", row.names = FALSE)
#write.table(posthocTrt, "./data/GLM-posthoc-Treatment.txt", sep = "\t", row.names = FALSE)
