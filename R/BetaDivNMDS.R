## BETA DIVERSITY (NMDS AND PERMANOVA)


## READ IN AND PREPARE DATA
if(!require(ggordiplots)) {
  require(devtools)
  devtools::install_github("jfq3/ggordiplots")
  require(ggordiplots)
}
require(vegan)
require(tidyverse)
require(phyloseq)
require(pairwiseAdonis)
require(RColorBrewer)
# load phyloseq objects and metadata
if(!require(repmis)) {
  install.packages("repmis")
  require(repmis)
}
source_data("https://github.com/EmilyB17/soils-micro/blob/master/data/phyloseq16S_ITS.RData?raw=true")

set.seed(123)

# set brewer colors
#brewer.pal(n = 3, name = "Dark2")
#treatmentcols <- c("HDG" = "#1B9E77", "LDG" =  "#D95F02", "NG" = "#7570B3")
source("./R/RColorBrewer.R")
# quick shortcut to avoid re-naming everything
treatmentcols <- c("HI" = "#1B9E77", "LO" =  "#D95F02", "NO" = "#7570B3")

#### ---- ITS -----

# Add EEA C cycling
# also make sure ADONIS variables are factors
sdITS18 <- sdITS18 %>% 
  mutate(Ccycl = AG + BG + BX + CBH,
         Ncycl = NAG + LAP,
         CNrat = Ccycl / Ncycl,
         Plot = as.factor(Plot),
         Block = as.factor(Block))

## conglomerate OTUs by Phyla and add to metadata
# agglomerate data by Phyla
tg <- tax_glom(psITS18, taxrank = rank_names(psITS18)[2], NArm = TRUE)
# name the phyla
taxa_names(tg) <- tg@tax_table@.Data[,2]
# make a dataframe
df <- as.data.frame(otu_table(tg))
# create relative abundance of Phyla for each sample
dc <- decostand(df, method = "total", MARGIN = 1)
# combine with sample data
sdITS18 <- cbind(dc, sdITS18) 
# remove Plot and Block (not interested in differences)
sdITS18$Plot <- NULL
sdITS18$Block <- NULL
# make sure GrazeTime is ordered
sdITS18$GrazeTime <- factor(sdITS18$GrazeTime, ordered = TRUE,
                            levels = c("PRE", "24H", "1WK", "4WK"))

# statistical testing
dis <- distance(psITS18T, method = "bray")
adonis(dis ~ Treatment * GrazeTime + sample_type,
       data = sdITS18, permutations = 999) # all are significant except Interactions
pairwise.adonis2(dis ~  Treatment * GrazeTime + sample_type,
                             data = sdITS18,
                             p.adjust.m = "bon", perm = 1000)
 
 
## BULK SOIL ONLY
disb <- distance(subset_samples(psITS18T, sample_type == "bulk"), method = "bray")
adonis(disb ~ Treatment * GrazeTime,
       data = filter(sdITS18, sample_type == "bulk"), permutations = 999) #Treatment is significant
# because treatment is significant, look at pairwise treatment comparisons
pairwise.adonis2(disb ~  Treatment * GrazeTime,
                 filter(sdITS18, sample_type == "bulk"),
                 p.adjust.m = "bon", perm = 1000)

## RHIZOSPHERIC SOIL ONLY
disr <- distance(subset_samples(psITS18T, sample_type == "rhizospheric"), method = "bray")
adonis(disr ~ Treatment * GrazeTime,
       data = filter(sdITS18, sample_type == "rhizospheric"), permutations = 999)  # Treatment is also significant
# because treatment is significant, look at pairwise treatment comparisons
pairwise.adonis2(disr ~  Treatment * GrazeTime,
                 filter(sdITS18, sample_type == "rhizospheric"),
                 p.adjust.m = "bon", perm = 1000) # no differences 

## since bulk and rhizospheric had the same differences, plot NMDS together

# ordinate
ord <- ordinate(psITS18T, method = "NMDS", distance = "bray", k = 2, 
                trymax = 1000, previousBest = TRUE) 
stressplot(ord, dis) #OK with 2 dimensions
eft <- envfit(ord, env = sdITS18, perm = 1000, na.rm = TRUE)
eft #BG, PEROX, Ascomycota, Basidiomycota, factors: Plot, Block, Treatment


## 1.  get NMDS coordinates from ordination object
coords <- as.data.frame(scores(ord, display = "sites")) %>% 
  # and add grouping variable
  mutate(Treatment = sdITS18$Treatment,
         Time = sdITS18$GrazeTime,
         Type = sdITS18$sample_type)

## 2.  get significant envfit objects to plot as arrows
spp <- as.data.frame(scores(eft, display = "vectors")) 
spp$species <- rownames(spp)
# subset to only p < 0.05 vectors to de-clutter 
sigs <- data.frame(eft$vectors$arrows[eft$vectors$pvals < 0.05,])
sigspecies <- rownames(sigs)
# subset
sigspecies <- spp %>% filter(species %in% sigspecies)

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sdITS18$Treatment, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

## 4. PLOT BY TREATMENT
## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, shape = Time, color = Treatment), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  # ENVFIT ARROWS
  geom_segment(data = sigspecies,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # ARROW TEXT
  #geom_text(data = sigspecies, aes(x = NMDS1, y = NMDS2, label = species),
   #         size = 3) +
  # SCALE CORRECTLY
  coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Treatment", shape = "Time") +
  scale_color_manual(values = treatmentcols)
# save
#ggsave("./data/plots/NMDS-fungi-Treatment-nolabels.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


## 5. PLOT BY GRAZETIME
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sdITS18$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

## 6. PLOT BY TIME
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sdITS18$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, color = Time, shape = Treatment), size = 1) +
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
  labs(color = "Sampling Time")
# save
#ggsave("./data/plots/NMDS-fungi-Time-labels.png", plot = last_plot(), dpi = "print")


## SOIL TYPE

# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sdITS18$sample_type, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 
## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, color = Type), size = 1) +
  # ENVFIT ARROWS
  geom_segment(data = sigspecies,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # GROUP NAMES AT ELLIPSE CENTER 
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  
  # ARROW TEXT
  geom_text(data = sigspecies, aes(x = NMDS1, y = NMDS2, label = species),
           size = 3) +
  # SCALE CORRECTLY
  coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Soil Type")
# save
#ggsave("./data/plots/NMDS-fungi-soil-type.png", plot = last_plot(), dpi = 600, width = 6.52, height = 4.83, units = "in")


#### ---- 16S ----

# Add EEA C cycling
# also make sure ADONIS variables are factors
sd16S18 <- sd16S18 %>% 
  mutate(Ccycl = AG + BG + BX + CBH,
         Ncycl = NAG + LAP,
         CNrat = Ccycl / Ncycl,
         Plot = as.factor(Plot),
         Block = as.factor(Block))

## conglomerate OTUs by Phyla and add to metadata
# agglomerate data by Phyla
tg <- tax_glom(ps16S18, taxrank = rank_names(ps16S18)[2], NArm = TRUE)
# name the phyla
taxa_names(tg) <- tg@tax_table@.Data[,2]
# make a dataframe
df <- as.data.frame(otu_table(tg))
# create relative abundance of Phyla for each sample
dc <- decostand(df, method = "total", MARGIN = 1)
# combine with sample data
sd16S18 <- cbind(dc, sd16S18) 
# remove Plot and Block (not interested in differences)
sd16S18$Plot <- NULL
sd16S18$Block <- NULL
# make sure GrazeTime is ordered for plots
sd16S18$GrazeTime <- factor(sd16S18$GrazeTime, ordered = TRUE,
                            levels = c("PRE", "24H", "1WK", "4WK"))

# statistical testing
dis <- distance(ps16S18T, method = "bray")

adonis(dis ~ Treatment * GrazeTime + sample_type,
       data = sd16S18, permutations = 999) # only GrazeTime
pairwise.adonis2(dis ~  GrazeTime * Treatment,
                 data = sd16S18,
                 p.adjust.m = "bon", perm = 1000) # all GrazeTime comparisons are significant
pairwise.adonis2(dis ~  Treatment * GrazeTime,
                 data = sd16S18,
                 p.adjust.m = "bon", perm = 1000)  # NG vs HDG difference; likely in the 1WK-4WK time above 

## BULK SOIL ONLY
disb <- distance(subset_samples(ps16S18T, sample_type == "bulk"), method = "bray")
adonis(disb ~ Treatment * GrazeTime,
       data = filter(sd16S18, sample_type == "bulk"), permutations = 999) #GrazeTime is significant
# because GrazeTime is significant, look at pairwise treatment comparisons
pairwise.adonis2(disb ~  GrazeTime * Treatment,
                 filter(sd16S18, sample_type == "bulk"),
                 p.adjust.m = "bon", perm = 1000) # PRE vs 24H, PRE vs 1WK

## RHIZOSPHERIC SOIL ONLY
disr <- distance(subset_samples(ps16S18T, sample_type == "rhizospheric"), method = "bray")
adonis(disr ~ Treatment * GrazeTime,
       data = filter(sd16S18, sample_type == "rhizospheric"), permutations = 999)  # also GrazeTime

# because GrazeTime is significant, look at pairwise treatment comparisons
pairwise.adonis2(disr ~  GrazeTime * Treatment,
                 filter(sd16S18, sample_type == "rhizospheric"),
                 p.adjust.m = "bon", perm = 1000) # PRE vs 24H, 24H vs 4WK, PRE vs 1WK, 24H vs 1WK

# ordinate
ord <- ordinate(ps16S18T, method = "NMDS", distance = "bray", k = 2, 
                trymax = 1000, previousBest = TRUE) 
stressplot(ord, dis) #OK with 2 dimensions
eft <- envfit(ord, env = sd16S18, perm = 1000, na.rm = TRUE)
eft #BG, PEROX, Ascomycota, Basidiomycota, factors: Plot, Block, Treatment


## 1.  get NMDS coordinates from ordination object
coords <- as.data.frame(scores(ord, display = "sites")) %>% 
  # and add grouping variable
  mutate(Treatment = sd16S18$Treatment,
         Time = sd16S18$GrazeTime,
         Type = sd16S18$sample_type)

## 2.  get significant envfit objects to plot as arrows
spp <- as.data.frame(scores(eft, display = "vectors")) 
spp$species <- rownames(spp)
# subset to only p < 0.05 vectors to de-clutter 
sigs <- data.frame(eft$vectors$arrows[eft$vectors$pvals < 0.05,])
sigspecies <- rownames(sigs)
# subset
sigspecies <- spp %>% filter(species %in% sigspecies)

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sd16S18$Treatment, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

## 4. PLOT BY TREATMENT
## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, shape = Time, color = Treatment), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  # ENVFIT ARROWS
  geom_segment(data = sigspecies,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # ARROW TEXT
  #geom_text(data = sigspecies, aes(x = NMDS1, y = NMDS2, label = species),
   #         size = 3) +
  # SCALE CORRECTLY
  coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Treatment", shape = "Time")



## 5. PLOT BY GRAZETIME
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sd16S18$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 

## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, shape = Treatment, color = Time), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
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
  labs(color = "Time", shape = "Treatment") +
  scale_color_manual(values = timecols)

# save
#ggsave("./data/plots/nmds-bacteria-time.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


## 6. PLOT BY SOIL TYPE
# save gg_ordiplot object to get ellipse values
plot <-  gg_ordiplot(ord, groups = sd16S18$sample_type, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean <- plot$df_mean.ord
# pull NMDS coordinates
ord.data <- plot$df_ord 


## create in ggplot2
ggplot(data = coords, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords, aes(x = NMDS1, y = NMDS2, color = Type), size = 1) +
  # ENVFIT ARROWS
  geom_segment(data = sigspecies,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  # GROUP NAMES AT ELLIPSE CENTER 
  annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  
  # ARROW TEXT
  geom_text(data = sigspecies, aes(x = NMDS1, y = NMDS2, label = species),
           size = 3) +
  # SCALE CORRECTLY
  coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Soil Type")
# save
#ggsave("./data/plots/NMDS-bacteria-soil-type.png", plot = last_plot(), dpi = 600, width = 6.52, height = 4.83, units = "in")



