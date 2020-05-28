### Enzymes

require(tidyverse)
require(vegan)
require(ggordiplots)
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


### TIME

sigt <- posthocTime %>% filter(p.value < 0.05)

# get significant data 
dat <- enzlog %>% 
  
  # AG, BG, BX, CBH, CCYCLING, CDOCNDON, CNRAT, NAG
  semi_join(sigt, by = "Type")  %>% 
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



# plot individual enzymes
enzl <- paste("AG", "BG", "BX", "CBH", "NAG")
plotdat <- dat %>% filter(str_detect(enzl, Type))

ggplot(data = plotdat, aes(x = Time, y = mean, group = Treatment)) +
  geom_point(aes(color = Treatment)) +
  geom_line(aes(linetype = Treatment, color = Treatment)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se,
                    linetype = Treatment, color = Treatment),
                width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  facet_wrap(~Type) +
  scale_color_manual(values = treatmentcols) +
  labs(x = "Sampling Time", y = "% change from PRE")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save
#ggsave("./data/plots/enzymes-indiv-lineplot.png", plot = last_plot(), dpi = "print")


# boxplot showing the spike in HDG in BG, CBH, NAG

dat1 <- enzlog %>% 
  
  # AG, BG, BX, CBH, CCYCLING, CDOCNDON, CNRAT, NAG
  semi_join(sigt, by = "Type")  %>% 
  
  mutate(Time = factor(Time, ordered = TRUE, levels = c("PRE", "diff_24H", "diff_1WK", "diff_4WK"))) %>% 
  mutate(Treatment = fct_recode(Treatment,
                                HDG = "HI", LDG = "LO", NG = "NO"),
         Time = fct_recode(Time,
                           PRE = "PRE",
                           post24H = "diff_24H",
                           post1WK = "diff_1WK",
                           post4WK = "diff_4WK"))

# get only BG, CBH, NAG
enzl <- paste("Ccycl", "CdocNdon", "CNrat", sep = "|")
plotdat <- dat1 %>% filter(!str_detect(enzl, Type))

ggplot(data = filter(plotdat, Time == "post24H"), aes(x = Type, y = value, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = treatmentcols) +
  labs(x = "Enzyme", y = "% change from 24H to PRE")
# save
#ggsave("./data/plots/enzymes-24h-HDGspike.png", plot = last_plot(), dpi = "print")


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
