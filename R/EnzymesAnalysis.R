### Enzymes

require(tidyverse)
require(vegan)
if(!require(ggordiplots)) {
  install.packages("remotes")
  remotes::install_github("jfq3/ggordiplots")
  require(ggordiplots)
}
require(emmeans)
source("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/R/CustomFunctions/FUNCTION_percentDifferences.R")
source("./R/RColorBrewer.R")

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

## get info for summary table
sum <- enzpd %>% 
  select(Plot, Block, Treatment, diffTimeSeries, AG, BG, BX, CBH, 
         LAP, NAG, PEROX, PHENOX, PHOS) %>% 
  pivot_longer(cols = c(AG, BG, BX, CBH, LAP, NAG, PEROX, PHENOX, PHOS),
               names_to = "Type", values_to = "percentdiff") %>% 
  group_by(Treatment, diffTimeSeries, Type) %>% 
  summarize(mean = round(mean(percentdiff), 2),
            se = round(se(percentdiff), 2)) %>% 
  pivot_wider(names_from = Type, values_from = c(mean, se))

#write.table(sum, file = "./data/enzymes-summary.txt", sep = "\t", row.names = FALSE)
  

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
enzl <- paste("AG", "BG", "BX", "CBH", "NAG", "LAP", "PEROX", "PHENOX")
plotdat <- dat %>% filter(str_detect(enzl, Type)) 


# PLOT

ggplot(data = plotdat, aes(x = Time, y = mean, group = Treatment)) +
  
  #geom_hline(yintercept = 0, color = "black") +
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

# for simplicity's sake, divide into C cycling and other
ccy <- c("AG", "BG", "BX", "CBH")
plotdat <- dat %>% filter(Type %in% ccy) 

# PLOT
ggplot(data = plotdat, aes(x = Time, y = mean, group = Treatment)) +
  
  #geom_hline(yintercept = 0, color = "black") +
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
#ggsave("./data/plots/enzymes-Ccycl-lines.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")

other <- c("LAP", "NAG", "PEROX", "PHENOX")
plotdat <- dat %>% filter(Type %in% other) 

# PLOT
ggplot(data = plotdat, aes(x = Time, y = mean, group = Treatment)) +
  
  #geom_hline(yintercept = 0, color = "black") +
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
#ggsave("./data/plots/enzymes-other-lines.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")



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

## ---- EEA PERMANOVAs ----

# get bacterial and fungal community data
require(phyloseq)
load("./data/phyloseq16S_ITS.RData")

# get enzyme data and calculate relative abundnace with vegan's decostand
enzr <- enz %>% 
  select(AG, BG, BX, CBH, LAP, NAG, PER1, PEROX, PHENOX, PHOS)

enzd <- decostand(enzr, method = "total", MARGIN = 2, na.rm = TRUE)

## get bacterial and fungal Phyla for NMDS metadata (only bulk samples)
tgf <- tax_glom(subset_samples(psITS18, psITS18@sam_data$sample_type != "rhizospheric"), 
                taxrank = rank_names(psITS18)[2], NArm = TRUE)
taxa_names(tgf) <- tgf@tax_table@.Data[,2]
# make into dataframe
df <- as.data.frame(otu_table(tgf))
# get relative abundance
dc <- decostand(df, method = "total", MARGIN = 2)
colnames(dc) <- gsub("p__", "fungi_", colnames(dc))
# same for 16S
tgs <- tax_glom(subset_samples(ps16S18, ps16S18@sam_data$sample_type != "rhizospheric"),
                taxrank = rank_names(ps16S18)[2], NArm = TRUE)
taxa_names(tgs) <- tgs@tax_table@.Data[,2]
# make into dataframe
dfs <- as.data.frame(otu_table(tgs))
dcs <- decostand(dfs, method = "total", MARGIN = 2)
colnames(dcs) <- paste0("bac_", colnames(dcs))
dcs <- dcs %>% rownames_to_column(var = "name")

# ITS has 71 samples while 16S has 72 - make the missing sample of ITS into "NAs"
which(!row.names(dcs) %in% row.names(dc))
row.names(dcs[35,]) # PRE_8_B
dc["PRE_8_B",] <- 0
dc <- dc[order(row.names(dc)),]
dc <- dc %>% rownames_to_column(var = "name")

# combine with other sample data
sd <- sd16S18 %>% 
  filter(sample_type == "bulk") %>% 
  select(-c(AG, BG, BX, CBH, LAP, NAG, PER1, PEROX, PHENOX, PHOS, biomass_kg_plot)) %>% 
  rownames_to_column(var = "name") %>% 
  left_join(dc, by = "name") %>% 
  left_join(dcs, by = "name") %>% 
  column_to_rownames(var = "name")

# PERMANOVA testing
dis <- vegdist(enzd)
adonis(dis ~ Treatment * GrazeTime, data = sd, permutations = 999) # GrazeTime is significant
# pairwise for GrazeTime
pairwise.adonis(x= dis, factors = sd$GrazeTime, p.adjust.m = "bonferroni") # 1WK vs 4WK

# ordinate
ord <- metaMDS(enzd, distance = "bray", k = 2, trymax = 1000, previousBest = TRUE)
stressplot(ord, dis) # OK with 2 dimensions
eft <- envfit(ord, env = sd, perm = 1000, na.rm = TRUE)

## PLOT NMDS
## 1.  get NMDS coordinates from ordination object
coords <- as.data.frame(scores(ord, display = "sites")) %>% 
  # and add grouping variable
  mutate(Treatment = sd$Treatment,
         Time = sd$GrazeTime)

## 2.  get significant envfit objects to plot as arrows
spp <- as.data.frame(scores(eft, display = "vectors")) 
spp$species <- rownames(spp)
# subset to only p < 0.05 vectors to de-clutter 
sigs <- data.frame(eft$vectors$arrows) %>% 
  filter(eft$vectors$pvals < 0.05,)
sigspecies <- rownames(sigs)
# subset
sigspecies <- spp %>% filter(species %in% sigspecies)

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sd$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

## Plot by TIME - WITH LABELS
## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, color = Time), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  # ENVFIT ARROWS
  geom_segment(data = sigspecies,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # ARROW TEXT
  geom_text(data = sigspecies, aes(x = NMDS1, y = NMDS2, label = species),
            size = 3) +
  # SCALE CORRECTLY
  coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Sampling Time") +
  # add constant colors for each Time
  scale_color_manual(values = timecols)

## ---- Edaphic PERMANOVA ----

# Euclidean distance PERMANOVA with soils variables

# normalize to relative abundance
sdf <- cn %>% 
  select(-c(NO3_mgkgdrysoil, mineralN_mgkgdrysoil, Plot, Block, Treatment, GrazeTime))
sdf <- decostand(sdf, method = "total", MARGIN = 2)
# make distance matrix
dis <- vegdist(sdf, method = "euclidean")

# get sample data
sd <- sd16S18 %>% 
  filter(sample_type == "bulk") %>% 
  select(-c(biomass_kg_plot, sample_type, NO3_mgkgdrysoil, mineralN_mgkgdrysoil)) %>% 
  rownames_to_column(var = "name") %>% 
  left_join(dc, by = "name") %>% 
  left_join(dcs, by = "name") %>% 
  column_to_rownames(var = "name")

# PERMANOVA
adonis(dis ~ Treatment * GrazeTime, data = sd) # GrazeTime is significant
# pairwise for GrazeTime
pairwise.adonis(x= dis, factors = sd$GrazeTime, p.adjust.m = "bonferroni") # 1WK vs 4WK, 24H vs 4WK, PRE vs 24H

# Ordinate
ord <- metaMDS(sdf,  k = 2, trymax = 1000, previousBest = TRUE)
stressplot(ord, dis)
eft <- envfit(ord, env = sd, perm = 1000, na.rm  =TRUE)

# Plot by GrazeTime
## PLOT NMDS
## 1.  get NMDS coordinates from ordination object
coords <- as.data.frame(scores(ord, display = "sites")) %>% 
  # and add grouping variable
  mutate(Treatment = sd$Treatment,
         Time = sd$GrazeTime)

## 2.  get significant envfit objects to plot as arrows
spp <- as.data.frame(scores(eft, display = "vectors")) 
spp$species <- rownames(spp)
# subset to only p < 0.05 vectors to de-clutter 
sigs <- data.frame(eft$vectors$arrows) %>% 
  filter(eft$vectors$pvals < 0.05,)
sigspecies <- rownames(sigs)
# subset
sigspecies <- spp %>% filter(species %in% sigspecies)

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sd$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

## Plot by TIME - WITH LABELS
## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, color = Time), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  # ENVFIT ARROWS
  geom_segment(data = sigspecies,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # ARROW TEXT
  geom_text(data = sigspecies, aes(x = NMDS1, y = NMDS2, label = species),
            size = 3) +
  # SCALE CORRECTLY
  coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Sampling Time") +
  # add constant colors for each Time
  scale_color_manual(values = timecols)

### ---- Mantel correlations ----

## EEA vs soil abiotic profiles
enzdis <- vegdist(enzd, method = "bray")
sdis <- vegdist(sdf, method = "euclidean")

mantel(xdis = sdis, ydis = enzdis, method = "spearman", na.rm = TRUE) # not significant
plot(sdis ~ enzdis)

# bacterial vs soil aboitic
bacdis <- distance(subset_samples(ps16S18T, sample_type == "bulk"), method = "bray")
mantel(xdis = sdis, ydis = bacdis, method = "spearman") # not significant

# bacterial vs EEA
mantel(xdis = bacdis, ydis = enzdis, method = "spearman") # not significant

# fungal vs soil abiotic
# fungi is missing one sample ("PRE_8_B") 
fundis <- distance(subset_samples(psITS18T, sample_type == "bulk"), method = "bray")
mantel(xdis = vegdist(sdf[-35,], method = "euclidean"), ydis = fundis, method = "spearman") # not significant

# fungal vs EEA
mantel(xdis = vegdist(enzd[-35,], method = "bray"), ydis = fundis, method = "spearman") # not significant

# bact vs fungal
bacs <- subset_samples(ps16S18T, rownames(ps16S18T@sam_data) != "PRE_8_B")
bacs1 <- subset_samples(bacs, bacs@sam_data$sample_type == "bulk")
bacd <- as.data.frame(otu_table(bacs1))
mantel(xdis = fundis, ydis = vegdist(bacd, method = "bray"), method = "spearman") ## significant
plot(fundis ~ vegdist(bacd, method = "bray"))
