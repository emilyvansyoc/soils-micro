### Enzymes

require(tidyverse)
require(vegan)
require(ggordiplots)
require(emmeans)
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")

# read cn for C:DOC and N:DON ratios
cn <- read.table("./data/2018CN.txt", sep = "\t", header = TRUE) %>% 
  filter(!Block == 4)

# wrangle enzymes
enz <- read.table("./data/2018Enzymes.txt", sep = "\t", header = TRUE) %>% 
  # remove block 4
  filter(!Block == 4) %>% 
  merge(cn, by = c("Plot", "Block", "Treatment", "GrazeTime")) %>% 
  # make wide
  pivot_wider(names_from = Substrate, values_from = Enzyme_nm_g_hr) %>% 
  # calculate other ratios
  mutate(Ccycl = AG + BG + BX + CBH,
         Ncycl = NAG + LAP,
         CNrat = (AG + BG + BX + CBH) / (NAG + LAP),
         CP = (AG + BG + BX + CBH + 1) / (PHOS + 1),
         NP = (NAG + LAP + 1) / (PHOS + 1),
         CDOC = Ccycl / NPOC_mgkgdrysoil,
         # some DON is zero; need to account for that 
         NDON = case_when(
           DON_mgkgdrysoil == 0 ~ Ncycl /1,
           DON_mgkgdrysoil !=0 ~ Ncycl / DON_mgkgdrysoil),
         CdocNdon = CDOC / NDON) %>% 
  # remove extra cn variables 
  select(-c(NO3_mgkgdrysoil, NH4_mgkgdrysoil, mineralN_mgkgdrysoil,
            NPOC_mgkgdrysoil, DON_mgkgdrysoil, grav_mois))

# calculate percent differences
enzpd <- percentDifferences(df = enz,
                            ids = c("Plot", "GrazeTime", "Block", "Treatment"),
                            timeKey = "GrazeTime",
                            timeLevels = c("PRE", "24H", "1WK", "4WK"),
                            level1 = "PRE") 

# make vertical and add log transformation
enzlog <- enzpd  %>% 
  # make vertical
  pivot_longer(cols = c(AG, BG, BX, CBH, Ccycl, CDOC, CdocNdon, CNrat,
                        CP, LAP, NAG, Ncycl, NDON, NP, PER1, PEROX, PHENOX, PHOS), names_to = "Type", values_to = "value") %>% 
  # add PRE at 0
  pivot_wider(names_from = diffTimeSeries, values_from = value) %>% 
  mutate(PRE = 0) %>% 
  pivot_longer(cols = c(PRE, diff_24H, diff_1WK, diff_4WK), names_to = "Time", values_to = "value") 

### ---- GLM loop ----

## get empty dataframes to fill as we go
modelFit <- data.frame()
posthocTrt <- data.frame()
posthocTime <- data.frame()

types <- unique(enzlog$Type)

for(i in 1:length(types)) {
  
  # perform model
  mod <- glm(value ~ Time * Treatment,
             data = filter(enzlog, Type == types[i]),
             family = gaussian(link = "identity"))
  # get model fit variables
  mf <- data.frame(Type = as.character(types[i]),
                   deviance = mod$deviance,
                   null.deviance = mod$null.deviance,
                   diff = mod$null.deviance - mod$deviance,
                   df.null = mod$df.null,
                   df.dev = mod$df.residual)
  modelFit <- rbind(mf, modelFit)
  # posthoc for treatment differences
  posthocTrt <- rbind(posthocTrt,
                      data.frame(Type = types[i],
                                 emmeans(mod, pairwise ~ Treatment | Time, type = "response")$contrasts))
  
  # posthoc for time differences
  posthocTime <- rbind(posthocTime,
                       data.frame(Type= types[i],
                                  emmeans(mod, pairwise ~ Time | Treatment, type = "response")$contrasts))
  
}

### ---- Significance ----

# function for standard error
se <- function(x) sqrt(var(x)/length(x))

# set brewer colors
brewer.pal(n = 3, name = "Dark2")
treatmentcols <- c("HDG" = "#1B9E77", "LDG" =  "#D95F02", "NG" = "#7570B3")

sig <- posthocTrt %>% filter(p.value < 0.05) # BG, CBH, CCYCL, CDOC*2, CNratio


### Line plot for each treatment over time
# want to plot ALL enzymes 

sigt <- posthocTime %>% filter(p.value < 0.05)

# get significant data 
dat <- enzlog %>% 
  
  # AG, BG, BX, CBH, CCYCLING, CDOCNDON, CNRAT, NAG
  #semi_join(sigt, by = "Type")  %>% 
  # get means and Se for line plot
  group_by(Type, Treatment, Time) %>% 
  summarize(mean = mean(value),
            se = se(value)) %>% 
  # reorder and rename factors
  ungroup() %>% 
  mutate(Time = factor(Time, ordered = TRUE, levels = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"))) %>% 
  mutate(Treatment = fct_recode(Treatment,
                                HDG = "HI", LDG = "LO", NG = "NO"),
         Time = fct_recode(Time,
                           PRE = "PRE",
                           `24H` = "diff_24H",
                           `1WK` = "diff_1WK",
                           `4WK` = "diff_4WK"))



# plot individual enzymes (so we need to remove the ratios)
enzl <- paste("AG", "BG", "BX", "CBH", "NAG", "LAP", "PHOS", "PEROX", "PHENOX")
plotdat <- dat %>% filter(str_detect(enzl, Type)) 




ggplot(data = plotdat, aes(x = Time, y = mean, group = Treatment)) +
  
  geom_hline(yintercept = 0, color = "black") +
  geom_point(aes(color = Treatment)) +
  geom_line(aes(linetype = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se,
                    linetype = Treatment, color = Treatment),
                width = 0.1) +
  facet_wrap(~Type, scales = "free") +
  scale_color_manual(values = treatmentcols) +
  labs(x = "Sampling Time", y = "% change from PRE")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save
#ggsave("./data/plots/enzymes-indiv-lineplot.png", plot = last_plot(), dpi = "print")


# boxplot showing the spike in HDG in BG, CBH, NAG

dat1 <- enzlog %>% 
  
  # AG, BG, BX, CBH, CCYCLING, CDOCNDON, CNRAT, NAG
  #semi_join(sigt, by = "Type")  %>% 
  
  mutate(Time = factor(Time, ordered = TRUE, levels = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"))) %>% 
  mutate(Treatment = fct_recode(Treatment,
                                HDG = "HI", LDG = "LO", NG = "NO"),
         Time = fct_recode(Time,
                           PRE = "PRE",
                           post24H = "diff_24H",
                           post1WK = "diff_1WK",
                           post4WK = "diff_4WK"))

# get all individual enzymes
enzl <- paste("AG", "BG", "BX", "CBH", "NAG", "LAP",  "PEROX", "PHENOX", sep = "|")
plotdat <- dat1 %>% filter(str_detect(enzl, Type))

ggplot(data = filter(plotdat, Time == "post24H"), aes(x = Type, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  labs(x = "Enzyme", y = "% change from PRE to 24H")
# save
#ggsave("./data/plots/enzymes-24h-HDGspike.png", plot = last_plot(), dpi = 600, width = 6.48, height = 4.67, units = "in")

# get only PHOS b/c it's weird
enzl <- paste("PHOS")
plotdat <- dat1 %>% filter(str_detect(enzl, Type))

ggplot(data = filter(plotdat, Time == "post24H"), aes(x = Type, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  labs(x = "Enzyme", y = "% change from 24H to PRE")
# save
#ggsave("./data/plots/enzymes-24h-PHOS-HDGspike.png", plot = last_plot(), dpi = "print")

## group by C substrate (complex vs simple)

# complex substrates
complexC <- paste("BG", "BX", "CBH")
plotdat <- dat1 %>% filter(str_detect(complexC, Type))

ggplot(data = filter(plotdat, Time == "post24H"), aes(x = Type, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  labs(x = "Enzyme", y = "% change from 24H to PRE")
# save
#ggsave("./data/plots/enzymes-24h-complexC-HDGspike.png", plot = last_plot(), dpi = "print")

## add 1WK as a comparison for LDG grazing
complexC <- paste("BG", "BX", "CBH")
t <- c("post24H", "post1WK")
plotdat <- dat1 %>% filter(str_detect(complexC, Type)) %>% filter(Time %in% t)

ggplot(data = plotdat, aes(x = Type, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  facet_wrap(~Time) +
  labs(x = "Enzyme", y = "% change from PRE")
# save
#ggsave("./data/plots/enzymes-24h-1WK-complexC-HDGspike.png", plot = last_plot(), dpi = "print")

## plot ratios
rats <- paste("Ccycl", "CdocNdon", "CNrat", sep = "|")
plotdat <- dat %>% filter(str_detect(rats, Type))

ggplot(data = plotdat, aes(x = Time, y = mean, group = Treatment)) +
  geom_point(aes(color = Treatment)) +
  geom_line(aes(linetype = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se,
                    linetype = Treatment, color = Treatment),
                width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(~Type) +
  scale_color_manual(values = treatmentcols) +
  labs(x = "Sampling Time", y = "% change from PRE")



## ---- Regression Plots ----

## plot from Sinsabaugh paper

theme_set(theme_bw())

# get data
enzln <- enz %>% 
  # make vertical
  pivot_longer(cols = c(AG, BG, BX, CBH, Ccycl, CDOC, CdocNdon, CNrat,
                        CP, LAP, NAG, Ncycl, NDON, NP, PER1, PEROX, 
                        PHENOX, PHOS), names_to = "Type", values_to = "value") %>% 
  
  # calculate natural log (log() in R, ln in calculators) 
  mutate(lnvalue = log(value)) %>% 
  select(-value) %>% 
  # make horizontal
  pivot_wider(names_from = Type, values_from = lnvalue, names_prefix = "ln")

# get the linear regression
modNcycl <- summary(lm(lnNcycl ~ lnBG, data = enzln))


# plot BG vs Ncycl

ggplot(data = enzln, aes(x = lnNcycl, y = lnBG)) + 
  geom_point(aes(color = Treatment, shape = GrazeTime)) +
  # plot 1:1 line (perfect diagonal)
  ggtitle("BG vs NAP + LAP")

# get the linear regression - first remove NAs
enzna <- enzln %>% drop_na() %>% filter(!enzln$lnPHOS == "-Inf")
modP <- summary(lm(lnPHOS ~ lnBG, data = enzna))

# plot BG vs P
ggplot(data = enzna, aes(x = lnPHOS, y = lnBG)) + 
  geom_point() +
  geom_point(aes(color = Treatment, shape = GrazeTime)) +
  geom_smooth(method = "lm") +
  geom_text(label = "r2 = 0.004", x = 2, y = 6.3) +
  ggtitle("BG vs PHOS")

