## revisions: figures 5 and 6
# Emily Van Syoc
# 11/2021

# load packages
require(vegan)
require(tidyverse)
require(phyloseq)
require(ggordiplots)
require(RColorBrewer)
require(ggpubr)

# load data
load("./data/phyloseq16S_ITS.RData")


# set brewer colors
brewer.pal(n = 8, name = "Dark2")
treatmentcols <- c("HDG" = "#1B9E77", "LDG" =  "#D95F02", "NG" = "#7570B3")
treatmentcols2 <- c("HI" = "#1B9E77", "LO" =  "#D95F02", "NO" = "#7570B3")
timecols <- c("PRE" = "#E7298A", "24H" = "#66A61E", "1WK" = "#E6AB02", "4WK" = "#A6761D")

## ---- ITS wrangle ----

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
disITS <- distance(psITS18T, method = "bray")

# ordinate
ordITS <- ordinate(psITS18T, method = "NMDS", distance = "bray", k = 2, 
                trymax = 1000, previousBest = TRUE) 

## ---- 16S Wrangle ----

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
dis16S <- distance(ps16S18T, method = "bray")

# ordinate
ord16S <- ordinate(ps16S18T, method = "NMDS", distance = "bray", k = 2, 
                trymax = 1000, previousBest = TRUE) 

## ---- Figure 5: bulk vs rhizo ----

## ITS


## 1.  get NMDS coordinates from ordination object
coordsITS <- as.data.frame(scores(ordITS, display = "sites")) %>% 
  # and add grouping variable
  mutate(Treatment = sdITS18$Treatment,
         Time = sdITS18$GrazeTime,
         Type = sdITS18$sample_type)

# plot by Soiltype
# save gg_ordiplot object to get ellipse values
plotITS <-  gg_ordiplot(ordITS, groups = sdITS18$sample_type, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ellITS <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.meanITS <- plot$df_mean.ord
# pull NMDS coordinates
ord.dataITS <- plot$df_ord 
## create in ggplot2
p1 <- ggplot(data = coordsITS, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ellITS, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coordsITS, aes(x = NMDS1, y = NMDS2, color = Type), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +

  # SCALE CORRECTLY
 # coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Soil Type")
# save
#ggsave("./data/plots/NMDS-fungi-soil-type.png", plot = last_plot(), dpi = 600, width = 6.52, height = 4.83, units = "in")


### 16S

## 1.  get NMDS coordinates from ordination object
coords16S <- as.data.frame(scores(ord16S, display = "sites")) %>% 
  # and add grouping variable
  mutate(Treatment = sd16S18$Treatment,
         Time = sd16S18$GrazeTime,
         Type = sd16S18$sample_type)

# plot by Soiltype
# save gg_ordiplot object to get ellipse values
plot16S <-  gg_ordiplot(ord16S, groups = sd16S18$sample_type, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell16S <- plot$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean16S <- plot$df_mean.ord
# pull NMDS coordinates
ord.data16S <- plot$df_ord 
## create in ggplot2
p2 <- ggplot(data = coords16S, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell16S, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords16S, aes(x = NMDS1, y = NMDS2, color = Type), size = 1) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  
  # SCALE CORRECTLY
  #coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Soil Type")
# save
#ggsave("./data/plots/NMDS-fungi-soil-type.png", plot = last_plot(), dpi = 600, width = 6.52, height = 4.83, units = "in")

### arrange
ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "right")
# save
ggsave("./data/plots/NMDS-new-fig5.png", plot = last_plot(), dpi = 600, width = 8, height = 4, units = "in")

### ---- Fig 6: Fungi ----

# plot by TIME

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plotITS <-  gg_ordiplot(ordITS, groups = sdITS18$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ellITS <- plotITS$df_ellipse
# get label coordinates for ellipse centers
NMDS.meanITS <- plot$df_mean.ord
# pull NMDS coordinates
ord.dataITS <- plot$df_ord 

## also need to rename factor levels
#coordsITS$Treatment <- fct_recode(coordsITS$Treatment, HDG = "HI", LDG = "LO", NG = "NO")

## create in ggplot2
timeITS <- ggplot(data = coordsITS, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ellITS, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coordsITS, aes(x = NMDS1, y = NMDS2, shape = Treatment, color = Time), size = 1.3) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +

  # SCALE CORRECTLY
  #coord_fixed() +
  # THEME
  theme_bw() +
  scale_color_manual(values = timecols) +
  # CHANGE LEGEND TITLE
  labs(color = "Time", shape = "Treatment")
# save
#ggsave("./data/plots/NMDS-fungi-Treatment-labels.png", plot = last_plot(), dpi = "print")

## Treatment

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plotITS1 <-  gg_ordiplot(ordITS, groups = sdITS18$Treatment, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ellITS1 <- plotITS1$df_ellipse
# get label coordinates for ellipse centers
NMDS.meanITS1 <- plotITS1$df_mean.ord
# pull NMDS coordinates
ord.dataITS1 <- plotITS1$df_ord 

## also need to rename factor levels
coords$Treatment <- fct_recode(coords$Treatment, HDG = "HI", LDG = "LO", NG = "NO")
df_ellITS1$Group <-fct_recode(df_ellITS1$Group, HDG = "HI", LDG = "LO", NG = "NO")
NMDS.meanITS1$Group <- fct_recode(NMDS.meanITS1$Group, HDG = "HI", LDG = "LO", NG = "NO")
## 4. PLOT BY TREATMENT
## create in ggplot2
treatITS <- ggplot(data = coordsITS, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ellITS1, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coordsITS, aes(x = NMDS1, y = NMDS2, shape = Time, color = Treatment), size = 1.3) +
  # GROUP NAMES AT ELLIPSE CENTER 
 # annotate("text",x = NMDS.meanITS1$x, y = NMDS.meanITS1$y,label=NMDS.meanITS1$Group, size = 5) +
  # SCALE CORRECTLY
  #coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Treatment", shape = "Time") +
  # colors for each treatment
  scale_color_manual(values = treatmentcols)

## ---- Fig 6: Bacteria ----

# plot by TIME

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plot16S <-  gg_ordiplot(ord16S, groups = sd16S18$GrazeTime, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell16S <- plot16S$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean16S <- plot$df_mean.ord
# pull NMDS coordinates
ord.data16S <- plot$df_ord 

## also need to rename factor levels
coords16S$Treatment <- fct_recode(coords16S$Treatment, HDG = "HI", LDG = "LO", NG = "NO")

## create in ggplot2
time16S <- ggplot(data = coords16S, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell16S, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords16S, aes(x = NMDS1, y = NMDS2, shape = Treatment, color = Time), size = 1.3) +
  # GROUP NAMES AT ELLIPSE CENTER 
  #annotate("text",x = NMDS.mean$x, y = NMDS.mean$y,label=NMDS.mean$Group, size = 5) +
  
  # SCALE CORRECTLY
  #coord_fixed() +
  # THEME
  theme_bw() +
  scale_color_manual(values = timecols) +
  # CHANGE LEGEND TITLE
  labs(color = "Time", shape = "Treatment")
# save
#ggsave("./data/plots/NMDS-fungi-Treatment-labels.png", plot = last_plot(), dpi = "print")

## Treatment

## 3. Pull ellipse information using gg_ordiplot package
# save gg_ordiplot object to get ellipse values
plot16S1 <-  gg_ordiplot(ord16S, groups = sd16S18$Treatment, label = FALSE, plot = FALSE)
# get ellipse coordinates
df_ell16S1 <- plot16S1$df_ellipse
# get label coordinates for ellipse centers
NMDS.mean16S1 <- plot16S1$df_mean.ord
# pull NMDS coordinates
ord.data16S1 <- plot16S1$df_ord 

## also need to rename factor levels
#coords$Treatment <- fct_recode(coords$Treatment, HDG = "HI", LDG = "LO", NG = "NO")
df_ell16S1$Group <-fct_recode(df_ell16S1$Group, HDG = "HI", LDG = "LO", NG = "NO")
#NMDS.mean16S1$Group <- fct_recode(NMDS.mean16S1$Group, HDG = "HI", LDG = "LO", NG = "NO")
## 4. PLOT BY TREATMENT
## create in ggplot2
treat16S <- ggplot(data = coords16S, aes(x = NMDS1, y = NMDS2)) + # label axises automatically
  # ELLIPSES
  geom_path(data = df_ell16S1, aes(x = x, y = y, color = Group), show.legend = FALSE) +
  # ORDINATION POINTS
  geom_point(data = coords16S, aes(x = NMDS1, y = NMDS2, shape = Time, color = Treatment), size = 1.3) +
  # GROUP NAMES AT ELLIPSE CENTER 
  # annotate("text",x = NMDS.mean16S1$x, y = NMDS.mean16S1$y,label=NMDS.mean16S1$Group, size = 5) +
  # SCALE CORRECTLY
  #coord_fixed() +
  # THEME
  theme_bw() +
  # CHANGE LEGEND TITLE
  labs(color = "Treatment", shape = "Time") +
  # colors for each treatment
  scale_color_manual(values = treatmentcols)

## ---- Fig 6 all ----

g1 <- ggarrange(timeITS, time16S, nrow = 1,
          common.legend = TRUE, legend = "right")

g2 <- ggarrange(treatITS, treat16S, nrow = 1,
                common.legend = TRUE, legend = "right")

ggarrange(g2, g1, nrow = 2)
# save
ggsave("./data/plots/NMDS-newfig6.png", plot = last_plot(), dpi = 600, width = 8, height = 6, units = "in")
