# Aim 1 sensitivity analysis code
# Author: Daniel Vader

library(dplyr)
library(readr)
library(lubridate)
library(mice)
library(survival)

################################################################################
# SA 1 Exposure defined using state stay at home order.
################################################################################

# Data Prep ####################################################################
# N = 16259
tdat <- readRDS("/data/analytic/PRACTICE-analytic-data_2021-08-23.rds") %>%
   select(patientid, state, diagnosisdate)
adat <- readRDS("/data/analytic/PRACTICE-imputed-data-derv_2021-08-23.rds") %>%
   mice::complete("long", include = T) %>%
   left_join(tdat, by="patientid") %>%
   # Clasify regions
   mutate(region = ifelse(state %in% c("CT", "ME", "MA", "NH", "RI", "Vt",
                                       "NJ", "NY", "PA"), 1,
                   ifelse(state %in% c("IL", "IN", "MI", "OH", "WI", "IA", "KS",
                                       "MN", "MO", "NE", "ND", "SD"), 2,
                   ifelse(state %in% c("DE", "FL", "GA", "MD", "NC", "SC", "VA", "DC",
                                       "AL", "KY", "MS", "TN", "AR", "LA", "OK",
                                       "TX"), 3,
                   ifelse(state %in% c("AZ", "CO", "ID", "MT", "NV", "NM", "UT",
                                       "WY", "AK", "CA", "HI", "OR", "WA"), 4, 9
                          )))),
          regionf = factor(region, levels=c(1,2,3,4,9), 
                           labels=c("Northeast", "Midwest", "South", "West", 
                                    "Other"))
          )
saveRDS(filter(adat, practicetype == "COMMUNITY", 
               !(studyperiod %in% c("Washout 2019", "Washout 2020"))),
        paste0("/data/analytic/A1SA2_PRACTICE-analytic-data_", Sys.Date(), ".rds"))

# Load state stay at home data
   # Data from: "See Which States and Cities Have Told Residents to Stay at 
   # Home" by Sarah Mervosh, Denise Lu, and Vanessa Swales. New York Times. 
   # April 2020. Accessed 8-29-2020
sah <- read_csv("/data/stay-at-home-orders/sahd.csv") %>%
  mutate(sah.date = mdy(date),
         # Switch state names to abbreviations
         state = ifelse(state == "Puerto Rico", "PR",
                        state.abb[match(state, state.name)])
         ) %>%
  select(state, sah.date)

# Drop academic practices (no state recorded). N = 14747
sa1.dat <- adat %>% filter(practicetype == "COMMUNITY")
table(sa1.dat$.imp)

# Drop subjects with missing state. N= 14669
sa1.dat <- sa1.dat %>% filter(!is.na(state))
table(sa1.dat$.imp)

# Join with stay at home order info.
sa1.dat <- left_join(sa1.dat, sah, by="state") %>% 
  mutate(sah.date = if_else(is.na(sah.date), mdy("04-07-2020"), sah.date),
         # Set study period relative to stay at home order date. Patients in 
            # washout if diagnosis occurred within 30 days of state stay at home
            # order.
         period = ifelse(yday(diagnosisdate) < yday(sah.date),
                         ifelse(yday(sah.date) - yday(diagnosisdate) <= 30, 
                                "washout", "pre"),
                         "post"
                         ),
         studyperiod2 = ifelse(period == "pre",
                              ifelse(year(diagnosisdate) == "2019", 
                                     "pre2019", "pre2020"),
                              ifelse(period == "post",
                                ifelse(year(diagnosisdate) == "2019",
                                       "post2019", "post2020"),
                                ifelse(year(diagnosisdate) == "2019",
                                       "washout2019", "washout2020")
                                )
                              ),
         studypf = ifelse(studyperiod2 == "pre2019", 1,
                          ifelse(studyperiod2 == "post2019", 2,
                                 ifelse(studyperiod2 == "pre2020", 3, 
                                        ifelse(studyperiod2 == "post2020", 4, NA)))),
         studypf = factor(studypf, 
                          levels=c(1,2,3,4), 
                          labels=c("pre2019", "post2019", "pre2020", "post2020"))
         ) %>% as.mids()

saveRDS(sa1.dat, paste0("/data/analytic/A1SA_PRACTICE-analytic-data_", Sys.Date(), ".rds"))

# Fit models ###################################################################
# N = 12589 after excluding washout
adat <- readRDS("/data/analytic/A1SA_PRACTICE-analytic-data_2021-10-18.rds") %>%
   filter(!(studyperiod2 %in% c("washout2019", "washout2020")))
source("1-5_functions.R")

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

saveRDS(m.nint, "models/a1sa1_no-cancer-interaction.rds")

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
saveRDS(m.wint, "models/a1sa1_cancer-interaction.rds")

# Compare
pool.compare(m.wint, m.nint)

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
saveRDS(m.wint.reth, "models/a1sa1_cancer-interaction-reth.rds")

# Compare
pool.compare(m.wint.reth, m.nint)

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
saveRDS(m.wint.age, "models/a1sa1_age-interaction.rds")

# Compare
pool.compare(m.wint.age, m.nint)

# Calculate 30 day survival probs ##############################################

# Unadjusted estimates
adat.u <- adat$data
tti90.fit <- survfit(Surv(adat.u$surv.days90, adat.u$surv.ind90) ~ adat.u$studypf)

# Extract unadjusted PP at 30 days by time period and overall
tti.fit.s <- summary(tti90.fit, time=30) 
tti.fit.sp <- cbind(tti.fit.s$strata,
                    1 -tti.fit.s$surv, 
                    1-tti.fit.s$upper, 
                    1-tti.fit.s$lower) %>% as.data.frame()
colnames(tti.fit.sp) <- c("period", "est", "lci", "uci")
write.csv(tti.fit.sp, "tables/intermediary/a1sa1_unadjusted-30daysurv-probs.csv", na = "")  

# Primary model estimates (cancer interaction)
pp.wint <- survpp(m.wint, adat) # Marginal
pp.wint.c <- survpp(m.wint, adat, cp = T) # Conditional on cancer

pp.wint2 <- bind_rows(pp.wint, pp.wint.c) # Join marginal and conditional
write.csv(pp.wint2, "tables/intermediary/a1sa1_pp_cancer-interaction-full.csv")

# Race / ethnicity estimates 
pp.wint.reth <- survpp(m.wint.reth, adat) # Marginal
pp.wint.reth.c <- survpp(m.wint.reth, adat, cp = T, cvar="reth") # Conditional on race/eth

pp.wint.reth2 <- bind_rows(pp.wint.reth, pp.wint.reth.c) # Join marginal and conditional
write.csv(pp.wint.reth2, "tables/intermediary/a1sa1_pp_reth-interaction-full.csv")

# Age estimates 
pp.wint.age <- survpp(m.wint.age, adat) # Marginal
pp.wint.age.c <- survpp(m.wint.age, adat, cp = T, cvar="age") # Conditional on age

pp.wint.age2 <- bind_rows(pp.wint.age, pp.wint.age.c) # Join marginal and conditional
write.csv(pp.wint.age2, "tables/intermediary/a1sa1_pp_age-interaction-full.csv")



################################################################################
# SA 1a Region interaction
################################################################################
adat <- readRDS("/data/analytic/A1SA2_PRACTICE-analytic-data_2021-10-18.rds") %>%
   as.mids()
# N = 12810
source("1-5_functions.R")

# Fit model w/ no interaction
m.nint.r <- with(adat, 
                 survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                    rblack + rhisp + rother + # rwhite ref
                                    age + gendm + 
                                    igov + iother + # commercial ref 
                                    ecog.c + opioid_adv + 
                                    ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                    diagnosisday + regionf + studypf 
                                 , 
                                 robust=T, cluster=practiceid,
                                 ties = "breslow",
                                 x = T
                 )
)

# Fit model w/ region*period interaction
m.wint.r <- with(adat, 
               survival::coxph(Surv(surv.days90, surv.ind90) ~ 
                                  rblack + rhisp + rother + # rwhite ref
                                  age + gendm + 
                                  igov + iother + # commercial ref 
                                  ecog.c + opioid_adv + 
                                  ccrc + cnsc + cpan + cpro + crcc + cucc + # cbre ref
                                  diagnosisday + regionf + studypf +
                                  
                                  regionf*studypf
                               , 
                               robust=T, cluster=practiceid,
                               ties = "breslow",
                               x = T
               )
)
saveRDS(m.wint.r, "models/a1sa1a_region-interaction.rds")

pool.compare(m.wint.r, m.nint.r)


pp.wint.r <- survpp(m.wint.r, adat) # Marginal
pp.wint.r.c <- survpp(m.wint.r, adat, cp = T, cvar="region") # Conditional on race/eth

pp.wint.r2 <- bind_rows(pp.wint.r, pp.wint.r.c) # Join marginal and conditional
write.csv(pp.wint.r2, "tables/intermediary/a1sa1a_pp_region-interaction-full.csv")


