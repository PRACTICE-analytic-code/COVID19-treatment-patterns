# Aim 1 analytic code
# Author: Daniel Vader

#.libPaths(c(.libPaths(), "rlibs"))

library(tidyverse)
library(mice)
library(survival)
library(splines)

# Read in analytic data set
adat <- readRDS("../../data/analytic/A1_PRACTICE-imputed-data-derv_2021-04-05.rds")
  

################################################################################
# Fit models
################################################################################

# Fit model w/ no cancer*period interaction
m.nint <- with(adat, 
          survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                              rblack + rhisp + rother + # rwhite ref
                              age + gendm + 
                              igov + iother + # commercial ref 
                              ecog.c + opioid_adv + 
                              ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                              diagnosisday + studypf, 
                            robust=T, cluster=practiceid,
                            ties = "breslow",
                            x = T
                          )
          )

saveRDS(m.nint, "models/aim1_no-cancer-interaction.rds")

# Fit model w/ cancer*period interaction
m.wint <- with(adat, 
                 survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                     rblack + rhisp + rother + # rwhite ref
                                     age + gendm + 
                                     igov + iother + # commercial ref 
                                     ecog.c + opioid_adv + 
                                     ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                     diagnosisday + studypf +
                                   
                                     ccrc*studypf + 
                                     cnsc*studypf + 
                                     cpan*studypf +
                                     cpro*studypf + 
                                     crcc*studypf + 
                                     cucc*studypf
                                 , 
                                 robust=T, cluster=practiceid,
                                 ties = "breslow",
                                 x = T
                 )
)
saveRDS(m.wint, "models/aim1_cancer-interaction.rds")

# Fit model w/ raceeth*period interaction
m.wint.reth <- with(adat, 
               survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                 rblack + rhisp + rother + # rwhite ref
                                 age + gendm + 
                                 igov + iother + # commercial ref 
                                 ecog.c + opioid_adv + 
                                 ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                 diagnosisday + studypf +

                                 rblack*studypf +
                                 rhisp*studypf +
                                 rother*studypf
                               , 
                               robust=T, cluster=practiceid,
                               ties = "breslow",
                               x = T
               )
)
saveRDS(m.wint.reth, "models/aim1_cancer-interaction-reth.rds")

# Fit model age*period interaction
m.wint.age <- with(adat, 
                     survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                       rblack + rhisp + rother + # rwhite ref
                                       age + gendm + 
                                       igov + iother + # commercial ref 
                                       ecog.c + opioid_adv + 
                                       ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                       diagnosisday + studypf +
                                       
                                       age*studypf
                                     , 
                                     robust=T, cluster=practiceid,
                                     ties = "breslow",
                                     x = T
                     )
)
saveRDS(m.wint.age, "models/aim1_age-interaction.rds")


################################################################################
# Export model fit statistics
################################################################################
m.nint <- readRDS("models/aim1_no-cancer-interaction.rds")
m.wint <- readRDS("models/aim1_cancer-interaction.rds")

# Get and store fit statistics function
g.aicbic <- function(m){
  a <- m$analyses
  t <- tibble(aic=numeric(), bic=numeric())
  
  for(i in 1:length(a)){
    aic <- AIC(a[[i]])
    bic <- BIC(a[[i]])
    t <- t %>% add_row(aic=aic, bic=bic)
  }
  return(t)
}

f.wint <- g.aicbic(m.wint)
f.nint <- g.aicbic(m.nint)

ftable <- bind_cols(f.nint, f.wint, f.nint-f.wint)
colnames(ftable) <- c("AIC:no cancer interaction", "BIC:no cancer interaction",
                        "AIC:cancer interaction", "BIC:cancer interaction",
                        "AIC:noi - AIC:yesi", "BIC:noi - BIC:yesi")

write.csv(ftable, "tables/aim1_fitstats.csv")



################################################################################
# Tabulate estimates
################################################################################

# Table that reports the hazard ratios for time periods compared to ref
# No cancer interaction model
p.nint <- pool(m.nint) # pool results
est.nint.fm <- summary(p.nint, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.nint.fm, "tables/intermediary/aim1_no-cancer-interaction_fm_2021-04-01.csv")

# Cancer interaction model
p.wint <- pool(m.wint) # pool results
est.wint.fm <- summary(p.wint, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.fm, "tables/intermediary/aim1_cancer-interaction_fm_2021-04-01.csv")

# Race/ethnicity interaction model
p.wint.reth <- pool(m.wint.reth) # pool results
est.wint.reth.fm <- summary(p.wint.reth, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.reth.fm, "tables/intermediary/aim1_reth-interaction_fm_2021-04-01.csv")

# Age interaction model
p.wint.age <- pool(m.wint.age) # pool results
est.wint.age.fm <- summary(p.wint.age, conf.int=T, exponentiate = T) %>%
  select(term, estimate, "2.5 %", "97.5 %") %>% 
  rename(Term = term, HR = estimate)
write.csv(est.wint.age.fm, "tables/intermediary/aim1_age-interaction_fm_2021-04-01.csv")



################################################################################
# Grab predicted probabilities of survival at 30 days by: 
   # year, period, and cancer (only in cancer interaction model)
################################################################################
m.nint <- readRDS("models/aim1_no-cancer-interaction.rds")
m.wint <- readRDS("models/aim1_cancer-interaction.rds")
m.wint.reth <- readRDS("models/aim1_cancer-interaction-reth.rds")
m.wint.age <- readRDS("models/aim1_age-interaction.rds")

adat <- readRDS("../../data/analytic/A1_PRACTICE-imputed-data-derv_2021-04-05.rds") 

source("1-5_functions.R")

# No interaction model estimates ###############################################
pp.nint <- survpp(m.nint, adat, cp=F)
write.csv(pp.nint, "tables/intermediary/aim1_pp_no-cancer-interaction-marginal.csv")

# Primary model estimates (cancer interaction) #################################
pp.wint <- survpp(m.wint, adat) # Marginal
pp.wint.c <- survpp(m.wint, adat, cp = T) # Conditional on cancer

pp.wint2 <- bind_rows(pp.wint, pp.wint.c) # Join marginal and conditional
write.csv(pp.wint2, "tables/intermediary/aim1_pp_cancer-interaction-full-p.csv")


# Race / ethnicity estimates ###################################################
pp.wint.reth <- survpp(m.wint.reth, adat) # Marginal
pp.wint.reth.c <- survpp(m.wint.reth, adat, cp = T, cvar="reth") # Conditional on race/eth

pp.wint.reth2 <- bind_rows(pp.wint.reth, pp.wint.reth.c) # Join marginal and conditional
write.csv(pp.wint.reth2, "tables/intermediary/aim1_pp_cancer-interaction-reth-full.csv")


# Age estimates ################################################################
pp.wint.age <- survpp(m.wint.age, adat) # Marginal
pp.wint.age.c <- survpp(m.wint.age, adat, cp = T, cvar="age") # Conditional on age

pp.wint.age2 <- bind_rows(pp.wint.age, pp.wint.age.c) # Join marginal and conditional
write.csv(pp.wint.age2, "tables/intermediary/aim1_pp_age-interaction-full.csv")