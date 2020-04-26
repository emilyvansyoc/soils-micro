## read in all univariate data from SQL for 2018 grazing year

require(tidyverse) # gonna need it

# lu tables
luID <- read.table("./data/Bean_SoilMicrobe_Database/2018luLab_ID.txt", 
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
luPlot <- read.table("./data/Bean_SoilMicrobe_Database/2018luPlots.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)[-13,] %>% 
  mutate(Plot = as.integer(Plot))

luSampPer <- read.table("./data/Bean_SoilMicrobe_Database/2018luSamplePeriods.txt",
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# data
cn <- read.table("./data/Bean_SoilMicrobe_Database/2018tblminN_NPOC_DON.txt",
                 sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  left_join(luID, by = c("lab_ID" = "lab_id")) %>% 
  left_join(luPlot, by = "Plot") %>% 
  left_join(luSampPer, by = "SampDate") %>% 
  # clean up
  select(-c(lab_ID, SampDate))

write.table(cn, "./data/2018CN.txt", sep = "\t", row.names = FALSE)

# data - enzymes
enz <- read.table("./data/Bean_SoilMicrobe_Database/2018tblEnzymes_Vertical.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  select(-RunDay) %>% 
  left_join(luID, by = c("lab_ID" = "lab_id")) %>% 
  left_join(luPlot, by = "Plot") %>% 
  left_join(luSampPer, by = "SampDate") %>% 
  # clean up
  select(-c(lab_ID, SampDate))

write.table(enz, "./data/2018Enzymes.txt", sep = "\t", row.names = FALSE)

# data - grav moisture
grav <- read.table("./data/Bean_SoilMicrobe_Database/2018tblGravMoisture.txt",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  select(-grav_mois) %>% 
  rename(gravmois = decimalgravmois) %>% 
  left_join(luID, by = c("lab_ID" = "lab_id")) %>% 
  left_join(luPlot, by = "Plot") %>% 
  left_join(luSampPer, by = "SampDate") %>% 
  # clean up
  select(-c(lab_ID, SampDate))

write.table(grav, file = "./data/2018gravmois.txt", sep = "\t", row.names = FALSE)

# data - RPM/biomass
rpm <- read.table("./data/Bean_SoilMicrobe_Database/2018tblRPM.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  left_join(luID, by = c("Plot", "SampDate")) %>% 
  left_join(luPlot, by = "Plot") %>% 
  left_join(luSampPer, by = "SampDate") %>% 
  # clean up
  select(-c(lab_id, SampDate))

write.table(rpm, "./data/2018RPM.txt", sep = "\t", row.names =  FALSE)
