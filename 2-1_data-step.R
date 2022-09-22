# Set up Aim 2 data for tabulation and analysis

library(tidyverse)
library(mice)

adat0 <- readRDS("/data/analytic/PRACTICE-analytic-data_2021-08-23.rds")

adat <- readRDS("/data/analytic/PRACTICE-analytic-data_2021-08-23.rds") %>%
  filter(studyperiod != "Washout 2020", studyperiod != "Washout 2019")

# Select table 1 data
adat2 <- adat %>% filter(cancer %in% c("bre", "pro", "ucc", "nsc"))

adat3 <- adat2 %>% filter(surv.days.a2 <= 60 & surv.ind == 1)

adat4 <- adat3 %>% filter(targeted.yes == 0)

saveRDS(adat4, paste("/data/analytic/A2_PRACTICE-analytic-data_", Sys.Date(), ".rds", sep=""))

# Select imputed data for analysis and merge with myelosuppressive indicator
adat.merge <- adat0 %>% select(patientid, myelo.yes, targeted.yes, surv.days.a2)

adat.mi <- readRDS("../../data/analytic/PRACTICE-imputed-data-derv_2021-08-23.rds") %>%
  mice::complete("long", include = T) %>%
  # Add variables left out of imputation data set
  left_join(adat.merge, by="patientid") %>%
  as.mids() %>%
  # Apply aim2 exclusion criteria.
  mice::filter(cancer %in% c("bre", "pro", "ucc", "nsc"),
         surv.days.a2 <= 60 & surv.ind90 == 1,
         targeted.yes == 0, 
         studyperiod != "Washout 2020", studyperiod != "Washout 2019")

# Save
saveRDS(adat.mi, paste("/data/analytic/A2_PRACTICE-imputed-data-derv_", Sys.Date(), ".rds", sep=""))
