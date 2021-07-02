## Alpha Diversity

require(tidyverse)
source("./R/RColorBrewer.R")
theme_set(theme_bw())

## note: block 4 has already been removed
# this data is rarefied
dat <- read.table("https://raw.githubusercontent.com/EmilyB17/grazing_soil_microbes/master/data/Metabarcoding_AlphaDiversity_All2.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = TRUE) %>% 
  # a priori: Shannon and Observed
  select(seq, Year, Treatment, GrazeTime, ID, Shannon, Observed) %>% 
  # get Plot and soil type info
  mutate(
    Plot = sapply(str_split(ID, "_"), `[`, 2),
    soil_type = sapply(str_split(ID, "_"), `[`, 3)
  )  %>% 
  # only year 2
  filter(Year == "2018")

# function for standard error
se <- function(x) sqrt(var(x)/length(x))


### ---- ITS: Observed ----

# subset data
its <- dat %>% 
  filter(seq == "ITS")

# test for normality
hist(its$Observed)
shapiro.test(its$Observed) # good to go

# ANOVA
mod <- aov(Observed ~ Treatment * GrazeTime, data = its) 

# t test for soil type
t <- t.test(Observed ~ soil_type, data = its) ## BULK IS HIGHER

# plot
ggplot(data = its, aes(x = soil_type, y = Observed, fill = soil_type)) +
  geom_boxplot() +
  labs(fill = "Soil Type", x = "Soil Type", y = "Observed ASVs")
# save
#ggsave("./data/plots/obeservedotus-fungi.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")


### subset bulk & rhizospheric
summary(aov(Observed ~ Treatment * GrazeTime, 
            data = filter(its, soil_type == "B")))
summary(aov(Observed ~ Treatment * GrazeTime, 
            data = filter(its, soil_type == "R")))

## get indices for supplementary data
sum <- its %>% 
  group_by(Treatment, GrazeTime, soil_type) %>% 
  summarize(meanOb = round(mean(Observed), 3),
            seOb = round(se(Observed), 3),
            meanShan = round(mean(Shannon), 3),
            seShan = round(se(Shannon), 3)
            
            )

#write.table(sum, "./data/alpha-indices-its.txt", sep = "\t", row.names = FALSE)

# get just bulk & rhizo
br <- its %>% 
  group_by(soil_type) %>% 
  summarize(meanOb = round(mean(Observed), 3),
            seOb = round(se(Observed), 3),
            meanShan = round(mean(Shannon), 3),
            seShan = round(se(Shannon), 3)
            
  )


### ---- ITS: Shannon ----

# test for normality
hist(its$Shannon)
shapiro.test(its$Shannon) # not normal

# GLM
summary(glm(Shannon ~ Treatment * GrazeTime, data = its))

# Wilcox rank-sum test
wilcox.test(Shannon ~ soil_type, data = its)

### subset bulk & rhizospheric
summary(glm(Shannon ~ Treatment * GrazeTime,
            data = filter(its, soil_type == "B")))
summary(glm(Shannon ~ Treatment * GrazeTime,
            data = filter(its, soil_type == "R")))

# plot
ggplot(data = its, aes(x = soil_type, y = Shannon, fill = soil_type)) +
  geom_boxplot() +
  labs(fill = "Soil Type", x = "Soil Type", y = "Shannon H'")
# save
#ggsave("./data/plots/shannon-fungi.png", plot = last_plot(), dpi = 600, height = 4.67, width = 6.48, units = "in")



### --- 16S: Observed ----

# subset data
bac <- dat %>% filter(seq == "16S")

# test for normality
hist(bac$Observed)
shapiro.test(bac$Observed) # not normal

# GLM
summary(glm(Observed ~ Treatment * GrazeTime, data = bac))

# Wilcox
wilcox.test(Observed ~ soil_type, data = bac, exact = FALSE)

### subset bulk & rhizospheric
summary(glm(Observed ~ Treatment * GrazeTime,
            data = filter(bac, soil_type == "B")))
summary(glm(Observed ~ Treatment * GrazeTime,
            data = filter(bac, soil_type == "R")))
# get mean & se for summary tables
sum <- bac %>% 
  group_by(Treatment, GrazeTime, soil_type) %>% 
  summarize(meanOb = round(mean(Observed), 3),
            seOb = round(se(Observed), 3),
            meanShan = round(mean(Shannon), 3),
            seShan = round(se(Shannon), 3)
            
  )

#write.table(sum, "./data/alpha-indices-16S.txt", sep = "\t", row.names = FALSE)


### ---- 16S: Shannon ----

# test for normality
hist(bac$Shannon)
shapiro.test(bac$Shannon) # not normal

# GLM
summary(glm(Shannon ~ Treatment * GrazeTime, data = bac))

# Wilcox
wilcox.test(Shannon ~ soil_type, data = bac, exact = FALSE)

### subset bulk & rhizospheric
summary(glm(Shannon ~ Treatment * GrazeTime,
            data = filter(bac, soil_type == "B")))
summary(glm(Shannon ~ Treatment * GrazeTime,
            data = filter(bac, soil_type == "R")))

### ---- write means/medians ----

# get mean for ITS Observed, everything else median/IQR

sumITS <- dat %>% 
  filter(seq == "ITS") %>% 
  group_by(Treatment, GrazeTime) %>% 
  summarize(meanOb = mean(Observed, na.rm = TRUE),
            sdOb = sd(Observed))

write.table(sumITS, file = "./data/results/alpha-div-its-observed.txt", sep = "\t", row.names = FALSE)

sumAll <- dat %>% 
  group_by(seq, Treatment, GrazeTime) %>% 
  summarize(medShannon = median(Shannon),
            iqrShannon = IQR(Shannon),
            medOb = median(Observed),
            iqrOb = IQR(Observed))

write.table(sumAll, file = "./data/results/alpha-div-its-16s-obs-shannon.txt",
            sep = "\t", row.names = FALSE)

# soil type
sumtypei <- dat %>% 
  filter(seq == "ITS") %>% 
  group_by(soil_type) %>% 
  summarize(meanOb = mean(Observed),
            sdOb = sd(Observed))
write.table(sumtypei, file = "./data/results/alpha-div-its-obs-soiltype.txt",
            sep = "\t", row.names = FALSE)

sumtype <- dat %>% 
  group_by(seq, soil_type) %>% 
  summarize(medShannon = median(Shannon),
            iqrShannon = IQR(Shannon),
            medOb = median(Observed),
            iqrOb = IQR(Observed))
write.table(sumtype, file = "./data/results/alpha-div-its-16s-ob-sha-soiltype.txt",
            sep = "\t", row.names = FALSE)
