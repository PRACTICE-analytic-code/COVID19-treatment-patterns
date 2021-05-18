# Aim 2 sensitivity analysis code
# Author: Daniel Vader

library(dplyr)
library(mice)
library(miceadds)
library(stdReg)
library(lmtest)
library(survival)

################################################################################
# Prep data
################################################################################
adat <- readRDS("/data/analytic/A1SA_PRACTICE-analytic-data_2021-04-13.rds")

# Select imputed data for analysis and merge with myelosuppressive indicator
adat.merge <- readRDS("/data/analytic/PRACTICE-analytic-data_2021-03-30.rds") %>% 
  select(patientid, myelo.yes, targeted.yes, surv.days.a2)

adat.mi <- adat %>%
  mice::complete("long", include = T) %>%
  # Add variables left out of imputation data set
  left_join(adat.merge, by="patientid") %>%
  as.mids() %>%
  # Apply aim2 exclusion criteria.
  mice::filter(cancer %in% c("bre", "pro", "ucc", "nsc"),
               surv.days.a2 <= 60 & surv.ind90 == 1,
               targeted.yes == 0, 
               studyperiod2 != "washout2020", studyperiod != "washout2019")

# Save
saveRDS(adat.mi, paste("/data/analytic/A2SA1_PRACTICE-imputed-data-derv_", Sys.Date(), ".rds", sep=""))

################################################################################
source("2-5_functions.R")

adat <- readRDS("/data/analytic/A2SA1_PRACTICE-imputed-data-derv_2021-04-13.rds")

# Convert to list
adat.l <- miceadds::mids2datlist(adat)

# Fit univariable model ########################################################
m.uni <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

# Fit model w/ no cancer*period interaction ####################################
m.nint <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                age + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                diagnosisday + studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

# Fit cancer interaction model #################################################
m.wint <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                age + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                diagnosisday + studypf +
                
                cnsc*studypf +
                cpro*studypf +
                cucc*studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

# Fit race/eth interaction model ###############################################
m.wint.reth <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                age + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                diagnosisday + studypf +
                
                rblack*studypf +
                rhisp*studypf +
                rother*studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

# Fit age interaction model ####################################################
m.wint.age <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                age + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                diagnosisday + studypf +
                
                age*studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

################################################################################
# Calculate predicted probs.
################################################################################

# Univariable model estimates
pp.uni <- glmpp(m.uni, adat.l) # Marginal
write.csv(pp.uni, paste0("tables/intermediary/a2sa1_pp_univariable_", Sys.Date(), ".csv"))

# Primary model estimates
pp.wint <- glmpp(m.wint, adat.l) # Marginal
pp.wint.c <- glmpp(m.wint, adat.l, cp = T) # Conditional on cancer

pp.wint2 <- bind_rows(pp.wint, pp.wint.c) # Join marginal and conditional
write.csv(pp.wint2, paste0("tables/intermediary/a2sa1_pp_cancer-interaction-full_", Sys.Date(), ".csv"))

# Race/ethnicity model estimates
pp.wint.reth <- glmpp(m.wint.reth, adat.l) # Marginal
pp.wint.reth.c <- glmpp(m.wint.reth, adat.l, cp = T, cvar="reth") # Conditional on race/eth

pp.wint.reth2 <- bind_rows(pp.wint.reth, pp.wint.reth.c) # Join marginal and conditional
write.csv(pp.wint.reth2, paste0("tables/intermediary/a2sa1_pp_reth-interaction-full_", Sys.Date(), ".csv"))

# Age model estimates
pp.wint.age <- glmpp(m.wint.age, adat.l) # Marginal
pp.wint.age.c <- glmpp(m.wint.age, adat.l, cp = T, cvar="age") # Conditional on age

pp.wint.age2 <- bind_rows(pp.wint.age, pp.wint.age.c) # Join marginal and conditional
write.csv(pp.wint.age2, paste0("tables/intermediary/a2sa1_pp_age-interaction-full_", Sys.Date(), ".csv"))