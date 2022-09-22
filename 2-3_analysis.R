# Aim 2 Analysis
# Author: Daniel Vader
# Date created: 2021-02-22

library(dplyr)
library(mice)
library(miceadds)
library(stdReg)
library(splines)
library(lmtest)
library(sandwich)
library(mitml)

# Functions for predicted probability generation and effect/CI estimation for
  # glm with robust standard errors.
source("2-5_functions.R")

# Aim 2 Analysis
adat <- readRDS("../../data/analytic/A2_PRACTICE-imputed-data-derv_2021-08-23.rds")

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

# Pool and save table of coefficients 
m.uni.p <- sglm(m.uni) %>% select(results, se, '(lower', 'upper)')
write.csv(m.uni.p, paste0("tables/intermediary/aim2_univariable_fm_", Sys.Date(), ".csv"))

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

# Pool and save table of coefficients 
m.nint.p <- sglm(m.nint) %>% select(results, se, '(lower', 'upper)')
write.csv(m.nint.p, paste0("tables/intermediary/aim2_no-cancer-interaction_fm_", Sys.Date(), ".csv"))


# Fit model w/ no cancer*period interaction: test splines - diagnosisday
m.nint.sd <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                age + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                ns(diagnosisday, df=4) + studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})


# Fit model w/ no cancer*period interaction: test splines - age

m.nint.sa3 <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                ns(age, df=3) + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                diagnosisday + studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

# Fit model w/ no cancer*period interaction: test splines - age
m.nint.sa4 <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                ns(age, df=4) + gendm + 
                igov + iother + # commercial ref 
                ecog.c + opioid_adv + 
                cnsc + cpro + cucc + # cbre ref
                diagnosisday + studypf, 
              family=binomial(link="logit"), 
              cluster=adat$data$practiceid
  )
})

testspline(m.nint, m.nint.sd) # Not sig
testspline(m.nint, m.nint.sa3) # Sig 
testspline(m.nint.sa3, m.nint.sa4) # Not sig



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

# Pool and save table of coefficients
m.wint.p <- sglm(m.wint) %>% select(results, se, '(lower', 'upper)')
write.csv(m.wint.p, paste0("tables/intermediary/aim2_cancer-interaction_fm_", Sys.Date(), ".csv"))

# compare (need to use mitlm since mice can't handle glm.cluster)
mitml::testModels(m.wint, m.nint)

# Fit cancer interaction model with splines for age ##################
m.wint.sa3 <- lapply(adat.l, FUN=function(data){
  glm.cluster(data=data, 
              formula=myelo.yes ~ 
                rblack + rhisp + rother + # rwhite ref
                ns(age, df=3) + gendm + 
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
sglm(m.wint.sa3)

testspline(m.wint, m.wint.sa3) # Sig

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

# Pool and save table of coefficients
m.wint.reth.p <- sglm(m.wint.reth) %>% select(results, se, '(lower', 'upper)')
write.csv(m.wint.reth.p, paste0("tables/intermediary/aim2_reth-interaction_fm_", Sys.Date(), ".csv"))

# compare (need to use mitlm since mice can't handle glm.cluster)
mitml::testModels(m.wint.reth, m.nint)

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

# Pool and save table of coefficients
m.wint.age.p <- sglm(m.wint.age) %>% select(results, se, '(lower', 'upper)')
write.csv(m.wint.age.p, paste0("tables/intermediary/aim2_age-interaction_fm_", Sys.Date(), ".csv"))

# compare (need to use mitlm since mice can't handle glm.cluster)
mitml::testModels(m.wint.age, m.nint)

################################################################################
# Calculate predicted probs.
################################################################################

# Univariable model estimates
pp.uni <- glmpp(m.uni, adat.l) # Marginal
write.csv(pp.uni, paste0("tables/intermediary/aim2_pp_univariable_", Sys.Date(), ".csv"))

# Primary model estimates
pp.wint <- glmpp(m.wint, adat.l) # Marginal
pp.wint.c <- glmpp(m.wint, adat.l, cp = T) # Conditional on cancer

pp.wint2 <- bind_rows(pp.wint, pp.wint.c) # Join marginal and conditional
#saveRDS(pp.wint2, paste0("models/aim2_pp_cancer-interaction-full_", Sys.Date(), ".rds"))
write.csv(pp.wint2, paste0("tables/intermediary/aim2_pp_cancer-interaction-full_", Sys.Date(), ".csv"))

# Race/ethnicity model estimates
pp.wint.reth <- glmpp(m.wint.reth, adat.l) # Marginal
pp.wint.reth.c <- glmpp(m.wint.reth, adat.l, cp = T, cvar="reth") # Conditional on cancer

pp.wint.reth2 <- bind_rows(pp.wint.reth, pp.wint.reth.c) # Join marginal and conditional
write.csv(pp.wint.reth2, paste0("tables/intermediary/aim2_pp_reth-interaction-full_", Sys.Date(), ".csv"))

# Age model estimates
pp.wint.age <- glmpp(m.wint.age, adat.l) # Marginal
pp.wint.age.c <- glmpp(m.wint.age, adat.l, cp = T, cvar="age") # Conditional on cancer

pp.wint.age2 <- bind_rows(pp.wint.age, pp.wint.age.c) # Join marginal and conditional
write.csv(pp.wint.age2, paste0("tables/intermediary/aim2_pp_age-interaction-full_", Sys.Date(), ".csv"))


